"""Reproduce every statistic, Table 2, and Figures 1-3 of the manuscript from data/ingredients.csv.

Usage:
    python analysis/reproduce.py

Outputs:
    results/results.md         all numbers reported in the manuscript, by section
    results/table2.csv         Table 2
    results/table3.csv         Table 3 (top/bottom ten notified ingredients)
    results/categories.csv     identification rate by functional category (Figure 3)
    figures/figure1.png ... figure3.png
"""
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.stats import chi2_contingency, mannwhitneyu
from statsmodels.stats.proportion import proportion_confint

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
RES.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

WEIGHTS = {"effect": 20, "design": 15, "sample": 15, "repro": 15, "coi": 15, "journal": 10, "safety": 10}
LABELS = {"design": "Study design", "journal": "Journal standing", "sample": "Participants & attrition",
          "safety": "Safety reporting", "effect": "Effect size & statistical rigor",
          "coi": "Funding independence", "repro": "Replication"}

df = pd.read_csv(ROOT / "data" / "ingredients.csv", dtype={"hff_code": str, "recognition_year": "Int64"})
IR = df[df.pathway == "individually_recognized"]
NT = df[df.pathway == "notified"]
S, sIR, sNT = df[df.scored == 1], IR[IR.scored == 1], NT[NT.scored == 1]

# Sanity check: composite score = sum(grade / 5 * weight)
calc = sum(S[f"grade_{k}"] / 5 * w for k, w in WEIGHTS.items())
assert (calc.round().astype(int) == S.composite_score.astype(int)).all(), "composite score mismatch"

out = []
def say(s=""): out.append(s)
def pct(a, n): return f"{100 * a / n:.1f}%"
def ci(a, n):
    lo, hi = proportion_confint(a, n, method="wilson")
    return f"{100 * a / n:.1f}% (95% CI {100 * lo:.1f}–{100 * hi:.1f})"
def rd(a, n, b, m):
    p1, p2 = a / n, b / m
    se = math.sqrt(p1 * (1 - p1) / n + p2 * (1 - p2) / m)
    d = p1 - p2
    return f"{100 * d:.1f} percentage points (95% CI {100 * (d - 1.96 * se):.1f}–{100 * (d + 1.96 * se):.1f})"
def chi(a, n, b, m):  # Yates-corrected chi-square (scipy default for 2x2)
    stat, p, *_ = chi2_contingency([[a, n - a], [b, m - b]])
    return f"χ² = {stat:.1f}, p = {p:.2g}"
def mwu(x, y):
    r = mannwhitneyu(x, y, alternative="two-sided")
    return f"U = {r.statistic:,.0f}, p = {r.pvalue:.3g}"

def criteria(frame):
    c = {"rct": frame.grade_design >= 3, "dbrct": frame.grade_design >= 4, "n50": frame.grade_sample >= 3,
         "peer": frame.grade_journal >= 3, "repl": frame.grade_repro >= 4, "indep": frame.grade_coi >= 4}
    return c

nIR, nNT, nAll = len(IR), len(NT), len(df)
say("# Reproduced results")
say(f"Dataset: {nAll} recognition records ({nIR} individually recognized, {nNT} notified); "
    f"assessments verified {df.verified_date.min()} to {df.verified_date.max()}.")

say("\n## 2.2 Study population")
yr = IR.recognition_year
say(f"Individually recognized: {nIR} ({pct(nIR, nAll)}); notified: {nNT} ({pct(nNT, nAll)}).")
say(f"Year-coded recognition numbers: {yr.notna().sum()} records ({int(yr.min())}–{int(yr.max())}); "
    f"recognized 2022 onwards: {(yr >= 2022).sum()} ({pct((yr >= 2022).sum(), nIR)} of pathway); "
    f"no year-coded number: {yr.isna().sum()}.")
say(f"Code-style numbers without notified label: {int(IR.code_style_number_without_notified_label.sum())}.")

say("\n## 3.1 Identification of human-trial literature")
say(f"All: {len(S)}/{nAll} = {ci(len(S), nAll)}; not scored: {nAll - len(S)} ({pct(nAll - len(S), nAll)}).")
say(f"Notified: {len(sNT)}/{nNT} = {ci(len(sNT), nNT)}; individually recognized: {len(sIR)}/{nIR} = {ci(len(sIR), nIR)}.")
say(f"Risk difference: {rd(len(sNT), nNT, len(sIR), nIR)}; {chi(len(sIR), nIR, len(sNT), nNT)}.")
post = IR[yr >= 2022]; pre = IR[yr < 2022]; noyr = IR[yr.isna()]
say(f"By recognition period (individually recognized): 2022+ {post.scored.sum()}/{len(post)} ({pct(post.scored.sum(), len(post))}); "
    f"before 2022 {pre.scored.sum()}/{len(pre)} ({pct(pre.scored.sum(), len(pre))}); "
    f"{chi(post.scored.sum(), len(post), pre.scored.sum(), len(pre))}; "
    f"no year-coded number {noyr.scored.sum()}/{len(noyr)} ({pct(noyr.scored.sum(), len(noyr))}).")

say("\n## 3.2 Score distribution")
sc = S.composite_score.astype(int)
bands = (sc // 10).clip(upper=9).value_counts().sort_index()
say(f"Mean {sc.mean():.1f} (SD {sc.std():.1f}; median {sc.median():.1f}; range {sc.min()}–{sc.max()}).")
say(f"60–69 band: n = {bands.get(6, 0)} ({pct(bands.get(6, 0), len(S))}); ≥ 90: n = {(sc >= 90).sum()} ({pct((sc >= 90).sum(), len(S))}).")
say(f"Individually recognized {sIR.composite_score.mean():.1f} ± {sIR.composite_score.std():.1f}; "
    f"notified {sNT.composite_score.mean():.1f} ± {sNT.composite_score.std():.1f}; "
    f"{mwu(sIR.composite_score, sNT.composite_score)}.")

say("\n## 3.3 Domain-level results (Figure 1)")
dom_rows = []
for k in WEIGHTS:
    g = S[f"grade_{k}"]
    dom_rows.append((k, 100 * (g <= 2).mean(), 100 * (g == 3).mean(), 100 * (g >= 4).mean(), 100 * g.mean() / 5))
dom_rows.sort(key=lambda r: -r[4])
for k, lo, mid, hi, att in dom_rows:
    say(f"{LABELS[k]}: attainment {att:.1f}%; low (1–2) {lo:.1f}%, grade 3 {mid:.1f}%, high (4–5) {hi:.1f}%.")
c = criteria(S)
say(f"Double-blind RCT: {c['dbrct'].sum()} ({pct(c['dbrct'].sum(), len(S))}).")
r1 = (S.grade_repro == 1).sum()
say(f"No independent replication (grade 1): {r1} = {ci(r1, len(S))}.")
for cat in ["not_reported", "manufacturer_sponsored_or_led", "independent_documented"]:
    a = (S.funding_category == cat).sum()
    say(f"Funding {cat}: {a} = {ci(a, len(S))}.")
say(f"Replication grade: IR {sIR.grade_repro.mean():.2f} vs notified {sNT.grade_repro.mean():.2f}; {mwu(sIR.grade_repro, sNT.grade_repro)}.")
say(f"Funding grade: IR {sIR.grade_coi.mean():.2f} vs notified {sNT.grade_coi.mean():.2f}; {mwu(sIR.grade_coi, sNT.grade_coi)}.")
mI = (sIR.funding_category == "manufacturer_sponsored_or_led").sum()
mN = (sNT.funding_category == "manufacturer_sponsored_or_led").sum()
say(f"Manufacturer-sponsored: IR {pct(mI, len(sIR))} vs notified {pct(mN, len(sNT))}; {chi(mI, len(sIR), mN, len(sNT))}.")

say("\n## 3.4 Cumulative criteria (Figure 2)")
steps = [("Individually recognized ingredients (all)", IR), ("Human trial identified in public literature", sIR)]
cur = sIR
for lab, key in [("Randomized controlled trial or higher", "rct"), ("Double-blind RCT or higher", "dbrct"),
                 ("+ ≥50 participants", "n50"), ("+ Peer-reviewed journal", "peer"),
                 ("+ Independent replication†", "repl"), ("+ Funding independence documented", "indep")]:
    cur = cur[criteria(cur)[key]]
    steps.append((lab, cur))
for lab, fr in steps:
    say(f"{lab}: {len(fr)} = {ci(len(fr), nIR)}.")
cN = criteria(sNT)
j4N = (cN["dbrct"] & cN["n50"] & cN["peer"] & cN["repl"]).sum()
j4I = len(steps[6][1])
say(f"Four criteria jointly: IR {j4I}/{nIR} ({pct(j4I, nIR)}) vs notified {j4N}/{nNT} ({pct(j4N, nNT)}); "
    f"risk difference {rd(j4N, nNT, j4I, nIR)}; {chi(j4I, nIR, j4N, nNT)}.")

say("\n## 3.5 Functional categories (Figure 3)")
cat = df.groupby("category_en").agg(n=("scored", "size"), identified=("scored", "sum")).reset_index()
cat["rate"] = 100 * cat.identified / cat.n
cat = cat[cat.n >= 11].sort_values(["rate", "n"], ascending=[False, False])
cat.to_csv(RES / "categories.csv", index=False)
for _, r in cat.iterrows():
    say(f"{r.category_en} (n = {r.n}): {r.identified}/{r.n} = {r.rate:.1f}%.")

say("\n## 3.6 Notified ingredients (Table 3)")
say(f"Scored notified: n = {len(sNT)}, {sNT.composite_score.mean():.1f} ± {sNT.composite_score.std():.1f}.")
one = (sNT.n_retained_articles == 1).sum()
noind = (sNT.funding_category != "independent_documented").sum()
say(f"Single-article evidence base: {pct(one, len(sNT))}; no documented independent funding: {noind}/{len(sNT)} ({pct(noind, len(sNT))}).")
t3 = sNT.sort_values("composite_score", ascending=False)
t3 = pd.concat([t3.head(10), t3.tail(10)])[["hff_code", "ingredient_en", "category_en", "n_retained_articles", "composite_score"]]
t3.to_csv(RES / "table3.csv", index=False)

# Table 2
def col(fr_all, fr_s):
    n, m = len(fr_all), len(fr_s)
    return {
        "Human-trial literature identified, n (%)": f"{m} ({pct(m, n)})",
        "Composite score, mean ± SD": f"{fr_s.composite_score.mean():.1f} ± {fr_s.composite_score.std():.1f}",
        "Replication grade, mean (1–5 scale)": f"{fr_s.grade_repro.mean():.2f}",
        "Funding-independence grade, mean (1–5 scale)": f"{fr_s.grade_coi.mean():.2f}",
        "Manufacturer-sponsored or -led, % of scored": pct((fr_s.funding_category == "manufacturer_sponsored_or_led").sum(), m),
        "Funding not reported, % of scored": pct((fr_s.funding_category == "not_reported").sum(), m),
        "Single-article evidence base, % of scored": pct((fr_s.n_retained_articles == 1).sum(), m),
        "Independent funding documented, % of scored": pct((fr_s.funding_category == "independent_documented").sum(), m),
    }
t2 = pd.DataFrame({f"Individually recognized (n = {nIR})": col(IR, sIR), f"Notified (n = {nNT})": col(NT, sNT),
                   f"All ingredients (N = {nAll})": col(df, S)})
t2.loc["Four evidence criteria met jointly, % of total"] = [pct(j4I, nIR), pct(j4N, nNT), "—"]
t2.to_csv(RES / "table2.csv")
say("\n## Table 2")
say(t2.to_markdown() if hasattr(t2, "to_markdown") else t2.to_string())

(RES / "results.md").write_text("\n".join(out) + "\n", encoding="utf-8")

# ---------------- Figures ----------------
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
                     "axes.spines.top": False, "axes.spines.right": False})
RED, GREY, BLUE, NAVY, LIGHT = "#b03a2e", "#c8c8c2", "#5b8bd6", "#1f3a6e", "#bcd0ef"

# Figure 1
fig, ax = plt.subplots(figsize=(7.2, 3.5), dpi=300)
for i, (k, lo, mid, hi, att) in enumerate(dom_rows):
    y = len(dom_rows) - 1 - i
    ax.barh(y, -lo, left=-mid / 2, color=RED, hatch="///", edgecolor="white", height=0.6, lw=0)
    ax.barh(y, mid, left=-mid / 2, color=GREY, height=0.6)
    ax.barh(y, hi, left=mid / 2, color=BLUE, height=0.6)
    ax.text(-mid / 2 - lo - 1, y, f"{lo:.0f}" if lo >= 1 else f"{lo:.1f}", ha="right", va="center",
            color=RED, fontsize=8, fontweight="bold")
    ax.text(mid / 2 + hi + 1, y, f"{hi:.0f}", ha="left", va="center", color=NAVY, fontsize=8, fontweight="bold")
    ax.text(1.0, y, f"{att:.1f}%", transform=ax.get_yaxis_transform(), ha="right", va="center",
            fontsize=9, fontweight="bold")
ax.text(1.0, len(dom_rows) - 0.45, "Attainment", transform=ax.get_yaxis_transform(), ha="right",
        fontsize=8, fontweight="bold")
ax.set_yticks(range(len(dom_rows)))
ax.set_yticklabels([LABELS[r[0]] for r in dom_rows][::-1], fontsize=9)
ax.axvline(0, color="black", lw=0.8)
ax.set_xlim(-65, 115)
ax.set_xticks([-60, -40, -20, 0, 20, 40, 60, 80])
ax.set_xticklabels(["60", "40", "20", "0", "20", "40", "60", "80"], fontsize=9)
ax.set_xlabel(f"% of ingredients (n = {len(S)}); left: low grades (1–2), centre: grade 3, right: high grades (4–5)", fontsize=9)
ax.legend(handles=[Patch(facecolor=RED, hatch="///", edgecolor="white", label="Low (1–2)"),
                   Patch(facecolor=GREY, label="Intermediate (3)"), Patch(facecolor=BLUE, label="High (4–5)")],
          loc="upper center", bbox_to_anchor=(0.6, -0.2), ncol=3, frameon=False, fontsize=8)
plt.tight_layout()
fig.savefig(FIG / "figure1.png", dpi=300)
plt.close(fig)

# Figure 2
fig = plt.figure(figsize=(7.2, 3.6), dpi=300)
gs = fig.add_gridspec(1, 2, width_ratios=[3.6, 1], wspace=0.5)
ax = fig.add_subplot(gs[0])
colors = [LIGHT, LIGHT, BLUE, BLUE, BLUE, BLUE, NAVY, NAVY]
for i, ((lab, fr), colr) in enumerate(zip(steps, colors)):
    y = len(steps) - 1 - i
    p = 100 * len(fr) / nIR
    ax.barh(y, p, color=colr, height=0.62)
    txt = f"{len(fr)} ({p:.1f}%)"
    if p > 30:
        ax.text(p - 1.5, y, txt, ha="right", va="center", fontsize=8, fontweight="bold",
                color="white" if colr == BLUE else "black")
    else:
        ax.text(p + 1.5, y, txt, ha="left", va="center", fontsize=8, fontweight="bold")
ax.set_yticks(range(len(steps)))
ax.set_yticklabels([s[0] for s in steps][::-1], fontsize=8)
ax.set_xlim(0, 100)
ax.set_xlabel(f"% of {nIR} individually recognized ingredients", fontsize=9)
ax.tick_params(axis="x", labelsize=8)
ax.set_title("(a) Cumulative attainment of evidence criteria", fontsize=9, fontweight="bold", loc="left")
ax2 = fig.add_subplot(gs[1])
vals = [100 * j4I / nIR, 100 * j4N / nNT]
ax2.bar([0, 1], vals, color=[BLUE, NAVY], width=0.55)
for i, v in enumerate(vals):
    ax2.text(i, v + 1.5, f"{v:.1f}%", ha="center", fontsize=8, fontweight="bold")
ax2.set_xticks([0, 1])
ax2.set_xticklabels([f"Individually\nrecognized\n(n={nIR})", f"Notified\n(n={nNT})"], fontsize=7)
ax2.set_ylim(0, 70)
ax2.set_ylabel("%", fontsize=8)
ax2.tick_params(axis="y", labelsize=8)
ax2.set_title("(b) Four criteria\nmet jointly‡", fontsize=9, fontweight="bold", loc="left")
fig.subplots_adjust(left=0.30, right=0.98, bottom=0.17, top=0.86)
fig.savefig(FIG / "figure2.png", dpi=300)
plt.close(fig)

# Figure 3
overall = 100 * len(S) / nAll
fig, ax = plt.subplots(figsize=(6.2, 5.2), dpi=300)
rows = list(cat.itertuples())
for i, r in enumerate(rows):
    y = len(rows) - 1 - i
    ax.barh(y, r.rate, color=BLUE if r.rate >= overall else RED, height=0.65)
    ax.text(r.rate + 1.2, y, f"{r.rate:.1f}%", va="center", fontsize=8)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([f"{r.category_en} (n={r.n})" for r in rows][::-1], fontsize=8)
ax.axvline(overall, ls="--", color="black", lw=0.9)
ax.text(overall + 0.8, len(rows) - 0.3, f"Overall {overall:.1f}%", fontsize=8)
ax.set_xlim(0, 118)
ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.tick_params(axis="x", labelsize=8)
ax.set_xlabel("Ingredients with human-trial literature identified (%)", fontsize=9)
plt.tight_layout()
fig.savefig(FIG / "figure3.png", dpi=300)
plt.close(fig)

print("\n".join(out))
print("\nWrote results/ and figures/.")
