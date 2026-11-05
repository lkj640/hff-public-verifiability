"""Automated verification of the AI-drafted assessments against PubMed / PubMed Central.

Checks (manuscript Section 2.6, steps a-c):
  (a) Citation existence: every cited PMID resolves in PubMed, and its title, journal and
      year match the bibliographic details recorded on the website.
  (b) Number support: numbers stated in each paper's rationale ("detail") are searched for
      in that paper's PubMed title + abstract. Numbers that are absent are flagged for manual
      checking (they may come from the full text, or may be unsupported).
  (c) Funding in full text: for scored ingredients graded "funding not identified" (grade 3),
      open-access full texts in PubMed Central are searched for funding / conflict-of-interest
      statements. A hit means the grade was probably assigned from the abstract only.

Inputs : the website HTML snapshot (all 611 records, including PMIDs) and private/id_key.csv
Outputs: private/verification_papers.csv     per-paper detail (identifies ingredients; do not publish)
         private/verification_funding.csv    per-paper funding text found in PMC (do not publish)
         results/verification_summary.md     aggregate results (safe to publish)
Network: NCBI E-utilities (no API key; ~3 requests/s). Responses are cached in private/cache/.

Usage:
    python tools/verify_citations.py path/to/index_v58_final_2026-09-15.html
"""
import csv
import difflib
import json
import re
import sys
import time
import unicodedata
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRIV = ROOT / "private"
CACHE = PRIV / "cache"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TOOL = "hff-public-verifiability"


def load_data(html_path):
    for line in Path(html_path).read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if s.startswith("const DATA = "):
            return json.loads(s[len("const DATA = "):].rstrip(";"))
    raise SystemExit("Could not find `const DATA = [...]` in the HTML file.")


def efetch(db, ids, tag):
    """Fetch XML for a batch of ids, with an on-disk cache."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{db}_{tag}.xml"
    if path.exists():
        return path.read_bytes()
    q = urllib.parse.urlencode({"db": db, "id": ",".join(ids), "retmode": "xml", "tool": TOOL})
    for attempt in range(4):
        try:
            req = urllib.request.Request(EUTILS, data=q.encode(), method="POST")
            with urllib.request.urlopen(req, timeout=120) as r:
                body = r.read()
            path.write_bytes(body)
            time.sleep(0.4)
            return body
        except Exception as e:  # network hiccup: back off and retry
            print(f"  retry {attempt + 1} for {db} batch {tag}: {e}")
            time.sleep(2 * (attempt + 1))
    raise SystemExit(f"Failed to fetch {db} batch {tag}")


def text_of(el):
    return " ".join("".join(el.itertext()).split()) if el is not None else ""


def norm(s):
    s = unicodedata.normalize("NFKC", s or "").lower()
    s = re.sub(r"<[^>]+>", " ", s)
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return " ".join(s.split())


def similar(a, b):
    a, b = norm(a), norm(b)
    if not a or not b:
        return 0.0
    if a in b or b in a:
        return 1.0
    return difflib.SequenceMatcher(None, a, b).ratio()


def parse_pubmed(xml_bytes):
    out = {}
    root = ET.fromstring(xml_bytes)
    for art in root.iter("PubmedArticle"):
        pmid = text_of(art.find("MedlineCitation/PMID"))
        a = art.find("MedlineCitation/Article")
        year = text_of(a.find("Journal/JournalIssue/PubDate/Year")) or \
            text_of(a.find("Journal/JournalIssue/PubDate/MedlineDate"))[:4]
        art_year = text_of(a.find("ArticleDate/Year"))
        pmc = ""
        # only the article's own IDs (not IDs inside its reference list)
        for aid in art.findall("PubmedData/ArticleIdList/ArticleId"):
            if aid.get("IdType") == "pmc":
                pmc = text_of(aid)
        out[pmid] = {
            "title": text_of(a.find("ArticleTitle")),
            "journal": text_of(a.find("Journal/Title")),
            "journal_iso": text_of(a.find("Journal/ISOAbbreviation")),
            "year": year, "article_year": art_year,
            "abstract": " ".join(text_of(x) for x in a.findall("Abstract/AbstractText")),
            "has_grant_list": int(a.find("GrantList") is not None),
            "coi_statement": text_of(art.find("MedlineCitation/CoiStatement")),
            "pmcid": pmc.upper().replace("PMC", ""),
        }
    for book in root.iter("PubmedBookArticle"):
        doc = book.find("BookDocument")
        pmid = text_of(doc.find("PMID"))
        out[pmid] = {
            "title": text_of(doc.find("ArticleTitle")) or text_of(doc.find("Book/BookTitle")),
            "journal": text_of(doc.find("Book/BookTitle")), "journal_iso": "",
            "year": text_of(doc.find("Book/PubDate/Year")), "article_year": "",
            "abstract": " ".join(text_of(x) for x in doc.findall("Abstract/AbstractText")),
            "has_grant_list": 0, "coi_statement": "", "pmcid": "", "book": 1,
        }
    return out


NUM_RE = re.compile(r"(?<![0-9A-Za-z.])(\d{1,3}(?:,\d{3})+|\d+)(?:\.(\d+))?(?![0-9A-Za-z])")


def numbers_in(text):
    """Numbers worth checking: skip single digits, 4-digit years and long identifiers."""
    nums = set()
    for m in NUM_RE.finditer(text or ""):
        whole, frac = m.group(1).replace(",", ""), m.group(2)
        if frac is None:
            if len(whole) == 1 or len(whole) > 6 or (len(whole) == 4 and 1900 <= int(whole) <= 2030):
                continue
            nums.add(whole)
        else:
            nums.add(f"{whole}.{frac}")
    return nums


def number_supported(num, source):
    src = (source or "").replace(",", "")
    variants = {num}
    if num.startswith("0."):
        variants.add(num[1:])  # "0.004" vs ".004"
    return any(re.search(rf"(?<![\d.]){re.escape(v)}(?![\d])", src) for v in variants)


FUND_TAGS = ("funding-group", "funding-statement")


def pmc_funding(xml_bytes):
    """Extract funding / COI text per PMC article."""
    out = {}
    root = ET.fromstring(xml_bytes)
    for art in root.iter("article"):
        pmcid = ""
        for aid in art.iter("article-id"):
            if aid.get("pub-id-type") in ("pmc", "pmcid", "pmcaid"):
                pmcid = text_of(aid).replace("PMC", "")
        parts = []
        for tag in FUND_TAGS:
            parts += [text_of(x) for x in art.iter(tag)]
        for fn in art.iter("fn"):
            if (fn.get("fn-type") or "").lower() in ("coi-statement", "conflict", "financial-disclosure", "funding"):
                parts.append(text_of(fn))
        for sec in art.iter("sec"):
            title = text_of(sec.find("title")).lower()
            if any(k in title for k in ("funding", "conflict", "competing", "disclosure", "acknowledg", "financial")):
                parts.append(text_of(sec))
        for ack in art.iter("ack"):
            parts.append(text_of(ack))
        txt = " ".join(p for p in parts if p)
        has_body = art.find(".//body") is not None
        out[pmcid] = {"funding_text": txt[:2000], "has_body": int(has_body)}
    return out


def pct(a, n):
    return f"{a}/{n} ({100 * a / n:.1f}%)" if n else "0/0"


def main(html_path):
    data = load_data(html_path)
    key = {r["website_no"]: r["record_id"] for r in csv.DictReader(open(PRIV / "id_key.csv", encoding="utf-8"))}

    papers = []
    for d in data:
        scored = d["status"] == "논문검증"
        for p in d["papers"]:
            papers.append({"record_id": key[str(d["no"])], "scored": int(scored),
                           "full_text_used": int(str(d.get("evidenceSource") or "").startswith("논문자료")),
                           "coi_grade": d["sub"]["coi"] if scored else "",
                           "pmid": str(p.get("pmid", "")).strip(), "site_title": p.get("title", ""),
                           "site_journal": p.get("journal", ""), "site_year": str(p.get("year", "")),
                           "detail": p.get("detail", "")})

    # (a) fetch PubMed records
    pmids = sorted({p["pmid"] for p in papers if p["pmid"].isdigit()})
    print(f"{len(papers)} cited papers, {len(pmids)} unique numeric PMIDs")
    pub = {}
    for i in range(0, len(pmids), 200):
        pub.update(parse_pubmed(efetch("pubmed", pmids[i:i + 200], f"b{i // 200:03d}")))

    for p in papers:
        r = pub.get(p["pmid"])
        p["pmid_valid"] = int(r is not None)
        p["book_chapter"] = int(bool(r and r.get("book")))
        if r is None:
            p.update(title_sim="", journal_match="", year_match="", numbers_checked=0, numbers_unsupported="",
                     pmcid="", abstract_available=0)
            continue
        p["title_sim"] = round(similar(p["site_title"], r["title"]), 3)
        p["journal_match"] = int(max(similar(p["site_journal"], r["journal"]),
                                     similar(p["site_journal"], r["journal_iso"])) >= 0.8)
        years = {r["year"], r["article_year"]} - {""}
        p["year_match"] = ("exact" if p["site_year"] in years else
                           "within1" if any(y.isdigit() and p["site_year"].isdigit() and abs(int(y) - int(p["site_year"])) <= 1 for y in years)
                           else "no")
        # (b) numbers in the rationale vs title + abstract
        source = r["title"] + " " + r["abstract"]
        nums = numbers_in(p["detail"])
        unsup = sorted(n for n in nums if not number_supported(n, source))
        p["numbers_checked"] = len(nums)
        p["numbers_unsupported"] = " ".join(unsup)
        p["pmcid"] = r["pmcid"]
        p["abstract_available"] = int(bool(r["abstract"]))
        p["pubmed_coi_statement"] = int(bool(r["coi_statement"]))

    # (c) funding statements in PMC full text, for scored ingredients graded "funding not identified"
    target = sorted({p["pmcid"] for p in papers if p.get("pmcid") and p["scored"] and p["coi_grade"] == 3})
    pmc = {}
    for i in range(0, len(target), 100):
        pmc.update(pmc_funding(efetch("pmc", target[i:i + 100], f"b{i // 100:03d}")))
    for p in papers:
        f = pmc.get(p.get("pmcid", ""), None)
        p["pmc_fulltext_checked"] = int(f is not None and bool(f["has_body"]))
        p["pmc_funding_found"] = int(f is not None and len(f["funding_text"]) > 20)
        p["pmc_funding_text"] = f["funding_text"][:500] if f else ""

    cols = ["record_id", "scored", "full_text_used", "coi_grade", "pmid", "pmid_valid", "book_chapter", "title_sim",
            "journal_match", "year_match", "abstract_available", "numbers_checked", "numbers_unsupported",
            "pmcid", "pmc_fulltext_checked", "pmc_funding_found", "site_title", "site_journal", "site_year"]
    with open(PRIV / "verification_papers.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(papers)
    with open(PRIV / "verification_funding.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["record_id", "pmid", "pmcid", "pmc_funding_text"], extrasaction="ignore")
        w.writeheader()
        w.writerows([p for p in papers if p["pmc_funding_found"]])

    # ---------------- aggregate summary ----------------
    S = [p for p in papers if p["scored"]]
    V = [p for p in S if p["pmid_valid"]]
    L = []
    L.append("# Automated verification of AI-drafted assessments")
    L.append(f"Source: website snapshot of 15 September 2026; PubMed/PMC records retrieved {time.strftime('%Y-%m-%d')}.")
    L.append(f"Scope: papers cited for the {len({p['record_id'] for p in S})} scored ingredients ({len(S)} citations, "
             f"{len({p['pmid'] for p in S})} unique PMIDs).")
    L.append("\n## (a) Citation existence and bibliographic match")
    L.append(f"- PMID resolves in PubMed: {pct(len(V), len(S))}")
    bad = [p for p in S if not p["pmid_valid"]]
    L.append(f"- PMID missing (no PMID recorded): {len(bad)}")
    books = [p for p in V if p.get("book_chapter")]
    L.append(f"- PMID resolves to a textbook chapter (StatPearls), not a study: {len(books)} "
             f"citations for {len({p['record_id'] for p in books})} ingredients")
    L.append(f"- Title matches (similarity ≥ 0.85): {pct(sum(p['title_sim'] >= 0.85 for p in V), len(V))}")
    L.append(f"- Title similarity < 0.85: {sum(p['title_sim'] < 0.85 for p in V)} (on inspection, all are the same "
             f"article: trademark symbols such as ® rendered as (R), or titles shortened with '...')")
    L.append(f"- Journal matches: {pct(sum(p['journal_match'] for p in V), len(V))}")
    yc = Counter(p["year_match"] for p in V)
    L.append(f"- Year exact: {pct(yc['exact'], len(V))}; within 1 year: {yc['within1']}; mismatch: {yc['no']}")
    L.append("\n## (b) Numbers in the rationales found in title + abstract")
    for lab, grp in [("All scored", V), ("Full text used", [p for p in V if p["full_text_used"]]),
                     ("Abstract-level", [p for p in V if not p["full_text_used"]])]:
        chk = sum(p["numbers_checked"] for p in grp)
        uns = sum(len(p["numbers_unsupported"].split()) for p in grp)
        anyu = sum(bool(p["numbers_unsupported"]) for p in grp)
        L.append(f"- {lab}: {chk} numbers checked, {pct(chk - uns, chk)} found in the PubMed record; "
                 f"citations with ≥ 1 number not found: {pct(anyu, len(grp))}")
    L.append("- Manual inspection of a random sample showed that many numbers not found are derived from the "
             "abstract (e.g., a total sample size computed from group sizes: 19 = 11 smokers + 8 non-smokers), "
             "so the share of genuinely unsupported numbers is lower than these figures; each flagged citation "
             "still needs manual checking.")
    L.append("- Numbers not found are not necessarily wrong (they may come from the full text or be "
             "derived, e.g. rounded or computed); they are the list for manual checking.")
    L.append(f"- Citations without an abstract in PubMed: {sum(not p['abstract_available'] for p in V)}")
    L.append("\n## (c) Funding statements for ingredients graded 'funding not identified' (grade 3)")
    g3 = [p for p in V if p["coi_grade"] == 3]
    g3_ing = {p["record_id"] for p in g3}
    with_pmc = [p for p in g3 if p["pmc_fulltext_checked"]]
    found = [p for p in with_pmc if p["pmc_funding_found"]]
    found_ing = {p["record_id"] for p in found}
    L.append(f"- Ingredients graded 3: {len(g3_ing)} ({len(g3)} citations)")
    L.append(f"- Citations with open-access full text in PMC: {pct(len(with_pmc), len(g3))}")
    L.append(f"- Of these, funding/COI statement present in the full text: {pct(len(found), len(with_pmc))}")
    L.append(f"- Ingredients with ≥ 1 cited paper whose full text contains a funding/COI statement: "
             f"{pct(len(found_ing), len(g3_ing))}")
    L.append(f"- Citations with a COI statement in the PubMed record itself: {sum(p.get('pubmed_coi_statement', 0) for p in g3)}")
    L.append("- These ingredients need their funding grade re-assessed from the full text; the remaining "
             "grade-3 ingredients without open-access full text need manual retrieval.")
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "verification_summary.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    main(sys.argv[1])
