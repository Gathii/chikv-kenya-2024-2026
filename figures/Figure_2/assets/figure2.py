# -*- coding: utf-8 -*-
"""
Figure 2

"""
import os
import csv
import datetime as dt
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

#conf
MIN_DENOM = 5          # hide weekly positivity where the denominator is tiny
BAR_WIDTH = 0.85

try:
    from scipy.signal import savgol_filter
    def smooth(y):
        n = len(y)
        w = min(9, n if n % 2 else n - 1)
        return np.asarray(y, float) if w < 3 else np.clip(savgol_filter(y, w, 2), 0, None)
except Exception:
    def smooth(y, w=5):
        y = np.asarray(y, float)
        return np.clip(np.convolve(y, np.ones(w) / w, mode="same"), 0, None)


def tint(hexcol, f=0.28):
    """Blend toward white -> a solid light colour. Avoids alpha."""
    h = hexcol.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    mix = lambda c: (c * f + 255 * (1 - f)) / 255
    return (mix(r), mix(g), mix(b))


HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "CHIKV_pos-neg_ajtmh_alltime_latest.txt")
OUT = HERE

#data
#Southcoast = Mombasa (sites LKN, MNB) | Northcoast = Lamu (LAM, KNM) | Hagadera = Garissa (GSA)
SERIES = [("Southcoast", "South Coast", "#058D59", "-"),
          ("Northcoast", "North Coast", "#CFCC1B", "--"),
          ("Hagadera",   "Hagadera",    "#C04627", ":")]

rows = []
with open(SRC, encoding="utf-8-sig", newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        r = {k: (v or "").strip() for k, v in r.items()}
        if r.get("Date_of_collection"):
            rows.append(r)


def week_start(iso):
    d = dt.date.fromisoformat(iso)
    return d - dt.timedelta(days=d.weekday())          # Monday of that week


agg = defaultdict(lambda: defaultdict(lambda: [0, 0]))  # region -> week -> [tested, positive]
epiweek = {}
for r in rows:
    w = week_start(r["Date_of_collection"])
    a = agg[r["Region"]][w]
    a[0] += 1
    if r["Results"] == "Positive":
        a[1] += 1
    if r.get("Epi_week", "").isdigit():
        epiweek.setdefault(w, int(r["Epi_week"]))

#continuous weekly grid - weeks with no surveillance samples stay visible as gaps
first = min(w for reg in agg for w in agg[reg])
last = max(w for reg in agg for w in agg[reg])
weeks = [first + dt.timedelta(weeks=i) for i in range((last - first).days // 7 + 1)]
x = np.arange(len(weeks))

tested, positive, negative = {}, {}, {}
for key, _lab, _c, _ls in SERIES:
    tested[key] = np.array([agg[key].get(w, [0, 0])[0] for w in weeks], float)
    positive[key] = np.array([agg[key].get(w, [0, 0])[1] for w in weeks], float)
    negative[key] = tested[key] - positive[key]

total_tested = sum(tested[k] for k, *_ in SERIES)

#tick labels every 4 weeks
ticks, ticklabels = [], []
for i, w in enumerate(weeks):
    if i % 4 == 0:
        ew = epiweek.get(w, w.isocalendar()[1])
        ticks.append(i)
        ticklabels.append(f"{w.year}\n{w.strftime('%b')}\nW{ew:02d}")

#styles
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Nimbus Sans", "DejaVu Sans"],
    "font.size": 8,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})


def panel_tag(ax, letter):
    ax.text(-0.062, 1.02, letter, transform=ax.transAxes,
            fontsize=10, fontweight="bold", va="bottom", ha="left")


def tidy(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3, width=0.8)
    ax.margins(x=0.01)


def save_all(fig, base):
    fig.savefig(base + ".png", dpi=600)
    fig.savefig(base + ".pdf")
    fig.savefig(base + ".eps")
    try:
        fig.savefig(base + ".tif", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    except Exception:
        fig.savefig(base + ".tif", dpi=600)
    print("wrote", os.path.basename(base) + ".{png,pdf,eps,tif}")


#figure
fig, (axA, axB, axC) = plt.subplots(
    3, 1, figsize=(7.2, 7.5), sharex=True,
    gridspec_kw={"height_ratios": [1.05, 1.0, 0.72], "hspace": 0.16})

fig.subplots_adjust(left=0.095, right=0.985, top=0.962, bottom=0.104)

#A: weekly confirmed cases
for key, lab, col, ls in SERIES:
    axA.plot(x, positive[key], lw=0, marker="o", ms=2.2,
             mfc=tint(col), mec="none", zorder=2)
    axA.plot(x, smooth(positive[key]), color=col, lw=1.8, ls=ls, label=lab,
             solid_capstyle="round", zorder=3)
axA.set_ylabel("Confirmed chikungunya\ncases per week (n)", fontsize=9, linespacing=1.3)
axA.set_ylim(bottom=0)
axA.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=7))
tidy(axA)
leg = axA.legend(frameon=False, fontsize=8, loc="upper right", handlelength=2.4,
                 title="Surveillance site")
leg.get_title().set_fontsize(8)
panel_tag(axA, "A")

#B: weekly samples tested, stacked and split by result
bottom = np.zeros(len(x))
handles = []
for key, lab, col, _ls in SERIES:
    for vals, fc, tag in ((positive[key], col, "positive"),
                          (negative[key], tint(col, 0.30), "negative")):
        if vals.sum() == 0:
            continue
        axB.bar(x, vals, BAR_WIDTH, bottom=bottom, color=fc,
                edgecolor="white", linewidth=0.15, zorder=2)
        bottom += vals
        handles.append(Patch(facecolor=fc, edgecolor="none", label=f"{lab}, {tag}"))
axB.set_ylabel("Samples tested\nper week (n)", fontsize=9, linespacing=1.3)
axB.set_ylim(0, float(total_tested.max()) * 1.30)
axB.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=6))
tidy(axB)
axB.legend(handles=handles, frameon=False, fontsize=6.6, loc="upper right",
           ncol=3, handlelength=1.5, columnspacing=1.0, labelspacing=0.3,
           handletextpad=0.5)
panel_tag(axB, "B")

#C: weekly test positivity
for key, lab, col, ls in SERIES:
    pct = np.where(tested[key] >= MIN_DENOM,
                   100 * positive[key] / np.maximum(tested[key], 1), np.nan)
    axC.plot(x, pct, color=col, lw=1.4, ls=ls, label=lab,
             solid_capstyle="round", zorder=3)
for key, _lab, col, _ls in SERIES:
    sampled = np.nonzero(tested[key] > 0)[0]
    if sampled.size == 0:
        continue
    overall = 100 * positive[key].sum() / max(tested[key].sum(), 1)
    axC.plot([sampled.min() - 0.5, sampled.max() + 0.5], [overall, overall],
             color=col, lw=0.6, ls=(0, (1, 2.5)), zorder=1)
axC.set_ylabel("Test positivity (%)", fontsize=9, linespacing=1.3)
axC.set_ylim(0, 104)
axC.set_yticks([0, 25, 50, 75, 100])
tidy(axC)
axC.text(0.615, 0.93,
         f"weeks with ≥{MIN_DENOM} tests at that site\ndotted line = site average over its sampled weeks",
         transform=axC.transAxes, fontsize=6.4, color="0.35",
         va="top", ha="left", linespacing=1.35)
panel_tag(axC, "C")

axC.set_xlabel("Epidemiological week", fontsize=9)
axC.set_xticks(ticks)
axC.set_xticklabels(ticklabels, fontsize=6)

save_all(fig, os.path.join(OUT, "Figure_2"))
plt.close(fig)

#numbers for text
print()
print("%-12s %7s %7s %7s   %s" % ("site", "tested", "cases", "pos.%", "peak-case week (cases/tested, positivity)"))
for key, lab, _c, _ls in SERIES:
    t, p = tested[key].sum(), positive[key].sum()
    i = int(np.argmax(positive[key]))
    print("%-12s %7d %7d %6.1f%%   %s: %d/%d = %.0f%%"
          % (lab, t, p, 100 * p / max(t, 1), weeks[i], positive[key][i],
             tested[key][i], 100 * positive[key][i] / max(tested[key][i], 1)))
print()
print("weeks on axis: %d (%s to %s); weeks with any sample: %d"
      % (len(weeks), weeks[0], weeks[-1], int((total_tested > 0).sum())))
print("TOTAL tested %d, cases %d" % (total_tested.sum(),
                                     sum(positive[k].sum() for k, *_ in SERIES)))
