"""Build the de-identified analysis dataset from the evidence website's HTML snapshot.

The website (index_v58_final_2026-09-15.html) embeds all ingredient-level assessments
as a JSON array (`const DATA = [...]`). This script reads that array and writes:

  data/ingredients.csv        one row per recognition record (de-identified for the
                              individually recognized pathway)
  data/retained_records.csv   retained publications (PMIDs) for notified ingredients
  private/id_key.csv          mapping from pseudonymous IDs to the original records
                              (NOT for publication; listed in .gitignore)

Usage:
    python tools/extract_from_website.py path/to/index_v58_final_2026-09-15.html
"""
import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

CATEGORY_EN = {
    "체지방": "Body fat", "피부": "Skin", "혈행": "Blood circulation", "면역기능": "Immune function",
    "관절": "Joint", "배변활동": "Bowel movement", "항산화": "Antioxidant", "인지기능": "Cognitive function",
    "눈건강": "Eye health", "갱년기": "Menopause", "혈당": "Blood glucose", "간건강": "Liver health",
    "소화/위건강": "Digestion/gastric", "근력": "Muscle strength", "피로개선": "Fatigue", "수면": "Sleep",
    "혈압": "Blood pressure", "스트레스": "Stress", "전립선": "Prostate", "기관지": "Bronchial health",
    "에너지·대사": "Energy metabolism", "뼈건강": "Bone health", "여성건강": "Women's health",
    "운동수행능력": "Exercise performance", "키성장": "Height growth", "잇몸건강": "Gum health", "모발": "Hair",
}
DOMAINS = ["design", "sample", "journal", "repro", "effect", "coi", "safety"]


def load_data(html_path):
    for line in Path(html_path).read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if s.startswith("const DATA = "):
            return json.loads(s[len("const DATA = "):].rstrip(";"))
    raise SystemExit("Could not find `const DATA = [...]` in the HTML file.")


def funding_category(grade):
    if grade <= 2:
        return "manufacturer_sponsored_or_led"
    if grade == 3:
        return "not_reported"
    return "independent_documented"


def main(html_path):
    data = load_data(html_path)
    (ROOT / "data").mkdir(exist_ok=True)
    (ROOT / "private").mkdir(exist_ok=True)

    rows, key, records = [], [], []
    ir_n = 0
    for d in sorted(data, key=lambda x: int(x["no"])):
        appno = str(d["appno"])
        notified = appno.startswith("고시형")
        scored = d["status"] == "논문검증"
        year_m = re.match(r"^(20\d\d)-\d+$", appno)
        if notified:
            code = appno.replace("고시형", "").strip()
            code = "" if "없음" in code else code
            rid = f"N-{code}" if code else f"N-X{int(d['no'])}"
        else:
            ir_n += 1
            rid = f"IR-{ir_n:03d}"
        row = {
            "record_id": rid,
            "pathway": "notified" if notified else "individually_recognized",
            "hff_code": code if notified else "",
            "ingredient_en": (d.get("nameEn") or "") if notified else "",
            "category_en": CATEGORY_EN.get(d["category"], d["category"]),
            "recognition_year": year_m.group(1) if (year_m and not notified) else "",
            "code_style_number_without_notified_label": int(
                (not notified) and bool(re.match(r"^\d-\d+$", appno))),
            "scored": int(scored),
            "n_retained_articles": len(d["papers"]),
            "verified_date": d["verifiedAt"],
        }
        for k in DOMAINS:
            row[f"grade_{k}"] = d["sub"][k] if scored else ""
        row["composite_score"] = int(d["score"]) if scored else ""
        row["funding_category"] = funding_category(d["sub"]["coi"]) if scored else ""
        rows.append(row)
        key.append({"record_id": rid, "website_no": d["no"], "recognition_no": appno, "name_ko": d["name"]})
        if notified:
            for p in d["papers"]:
                records.append({"record_id": rid, "pmid": p.get("pmid", ""), "year": p.get("year", ""),
                                "journal": p.get("journal", ""), "title": p.get("title", ""),
                                "url": p.get("url", "")})

    def write(path, items):
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(items[0].keys()))
            w.writeheader()
            w.writerows(items)

    write(ROOT / "data" / "ingredients.csv", rows)
    write(ROOT / "data" / "retained_records_notified.csv", records)
    write(ROOT / "private" / "id_key.csv", key)
    print(f"{len(rows)} records -> data/ingredients.csv; {len(records)} notified retained records")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    main(sys.argv[1])
