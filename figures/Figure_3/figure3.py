"""
Figure 3 — AJTMH submission.

"""

import sys, os, json, csv, re
from datetime import datetime
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
sys.setrecursionlimit(300000)

BASE = os.path.dirname(os.path.abspath(__file__))
RES  = os.path.join(BASE, "results_20260819")
NWK  = os.path.join(RES, "tree", "ecsa_tree_01.nwk")
BLJ  = os.path.join(RES, "tree", "branch_lengths.json")
META = os.path.join(RES, "metadata.tsv")
TRAJ = os.path.join(RES, "traits.json")
AAJ  = os.path.join(RES, "aa_muts.json")

#colours
#Panel A
REGION_COL = {
    "Africa"        : "#6C3483",   # dark purple
    "Asia"          : "#17A589",   # teal
    "Europe"        : "#A569BD",   # medium purple
    "Indian Ocean"  : "#2E86C1",   # steel blue
    "South America" : "#E91E8C",   # pink/magenta
    "North America" : "#CA6F1E",   # burnt orange
    "?"             : "#AAAAAA",
}
KENYA_COL  = "#111111"   # all Kenya — black
INDIA_COL  = "#1A3A8F"   # OR715104 dark navy

#Panel B tip colours by division
DIV_COL = {
    "Garissa" : "#E87722",
    "Mombasa" : "#2E8B57",
    "Lamu"    : "#F5C518",
}
CLADE_BRANCH = "#5D4037"   # dark brown for Panel B branches

#load data
bl     = json.load(open(BLJ))["nodes"]
md     = {r["strain"].strip(): r
          for r in csv.DictReader(open(META, encoding="utf-8"), delimiter="\t")}
traits = json.load(open(TRAJ))["nodes"]
aa_muts= json.load(open(AAJ))["nodes"]

def gmeta(name, field):
    return md.get(name or "", {}).get(field, "").strip()

def to_numdate(name):
    """Date string → decimal year"""
    d = bl.get(name or "", {})
    s = d.get("date", "")
    try:
        fmt = "%Y-%m-%d" if len(s) >= 10 else "%Y-%m" if len(s) == 7 else "%Y"
        dt  = datetime.strptime(s[:len(fmt)], fmt)
        return dt.year + (dt.timetuple().tm_yday - 1) / 365.25
    except Exception:
        pass
    ci = d.get("num_date_confidence")
    if ci and len(ci) == 2:
        return (ci[0] + ci[1]) / 2
    return d.get("numdate") or d.get("num_date")

def node_region(name):
    r = gmeta(name, "region")
    if r: return r
    return traits.get(name or "", {}).get("region", "?") or "?"

def node_country(name):
    c = gmeta(name, "country")
    if c: return c
    return traits.get(name or "", {}).get("country", "?") or "?"

#newick parser
class Node:
    __slots__ = ("name","children","parent","x","y")
    def __init__(self):
        self.name=None; self.children=[]; self.parent=None; self.x=0.0; self.y=0.0

def parse_nwk(s):
    s = s.strip().rstrip(";"); p = [0]
    def lab():
        m = re.match(r"[^,():;]*", s[p[0]:]); l = m.group(0); p[0]+=len(l); return l
    def skip_bl():
        if p[0]<len(s) and s[p[0]]==":":
            p[0]+=1; m=re.match(r"[-0-9.eE+]*", s[p[0]:]); p[0]+=len(m.group(0))
    def pc():
        n=Node()
        if s[p[0]]=="(":
            p[0]+=1
            while True:
                c=pc(); c.parent=n; n.children.append(c)
                if s[p[0]]==",": p[0]+=1; continue
                if s[p[0]]==")": p[0]+=1; break
            n.name=lab() or None
        else:
            n.name=lab() or None
        skip_bl(); return n
    return pc()

root = parse_nwk(open(NWK).read())

#assign x = decimal year
def assign_x(n):
    nd = to_numdate(n.name)
    n.x = nd if nd else (n.parent.x if n.parent else 2004.0)
    for c in n.children: assign_x(c)
assign_x(root)

#assign y = leaf-counter order
lc = [0]
def assign_y(n):
    if not n.children:
        n.y = lc[0]; lc[0] += 1
    else:
        for c in n.children: assign_y(c)
        n.y = np.mean([c.y for c in n.children])
assign_y(root)
N_tips = lc[0]

all_nodes = []
def collect(n):
    all_nodes.append(n)
    for c in n.children: collect(c)
collect(root)
tips  = [n for n in all_nodes if not n.children]
itern = [n for n in all_nodes if n.children]

#Kenyan 2024-2026 tips for clade shading
ke_all   = {t.name for t in tips if node_country(t.name) == "Kenya"}
ke2425   = {t.name for t in tips
            if node_country(t.name) == "Kenya"
            and gmeta(t.name, "date")[:4] in ("2024","2025","2026")}
#functional node lookup
def _path(n):
    p = []
    while n: p.append(n); n = n.parent
    return p

def mrca(nodes):
    """Most recent common ancestor of a list of nodes."""
    nodes = [x for x in nodes if x is not None]
    if not nodes: return None
    common = set(id(x) for x in _path(nodes[0]))
    for x in nodes[1:]:
        common &= set(id(y) for y in _path(x))
    best, bd = None, -1
    for x in all_nodes:
        if id(x) in common and len(_path(x)) > bd:
            best, bd = x, len(_path(x))
    return best

def subtree_tips(n):
    out = []
    def w(x):
        if x.children:
            for c in x.children: w(c)
        else:
            out.append(x)
    if n: w(n)
    return out

STUDY_PREFIX = ("AFI-", "nvrl-")
study_tips   = [t for t in tips if t.name and t.name.startswith(STUDY_PREFIX)]
garissa_tips = [t for t in study_tips if t.name.startswith("nvrl-")]
coastal_tips = [t for t in study_tips if t.name.startswith("AFI-")]

clade_root = mrca(study_tips)
clade_set  = set()
if clade_root:
    def cc(n): clade_set.add(n.name); [cc(c) for c in n.children]
    cc(clade_root)

#fig setup
mm = 1/25.4
FIG_W, FIG_H = 174*mm, 220*mm
DPI = 300
from matplotlib import font_manager as _fm
_avail = {f.name for f in _fm.fontManager.ttflist}
FONT = next((f for f in ("Arial", "Helvetica", "Nimbus Sans", "DejaVu Sans") if f in _avail),
            "sans-serif")
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})
FS   = {"axis":6, "panel":9, "legend":5.5, "annot":4.5, "mut":4.0}

fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=DPI)
ax_A = fig.add_axes([0.02, 0.06, 0.46, 0.91])   # Panel A: left 48%
ax_B = fig.add_axes([0.54, 0.06, 0.44, 0.91])   # Panel B: right 44%

# PANEL A
# =======
def branch_col_A(parent_name):
    """Branch colour = parent node's inferred region."""
    ctry = node_country(parent_name)
    if ctry == "Kenya": return KENYA_COL, 0.55
    reg = node_region(parent_name)
    return REGION_COL.get(reg, REGION_COL["?"]), 0.35

#draw branches
for n in all_nodes:
    if n.parent is None: continue
    c, lw = branch_col_A(n.parent.name)
    ax_A.plot([n.parent.x, n.x], [n.y, n.y],
              lw=lw, color=c, zorder=1, rasterized=True,
              solid_capstyle="butt")
for n in itern:
    c, lw = branch_col_A(n.name)
    ys = [ch.y for ch in n.children]
    ax_A.plot([n.x, n.x], [min(ys), max(ys)],
              lw=lw, color=c, zorder=1, rasterized=True,
              solid_capstyle="butt")

#draw tips
for t in tips:
    ctry = node_country(t.name)
    if ctry == "Kenya":
        sz = 5.0 if t.name in ke2425 else 3.5
        ax_A.scatter(t.x, t.y, s=sz, c=KENYA_COL,
                     edgecolors="none", zorder=3, rasterized=True)
    else:
        reg = node_region(t.name)
        c   = REGION_COL.get(reg, REGION_COL["?"])
        ax_A.scatter(t.x, t.y, s=2.5, c=c,
                     edgecolors="none", zorder=3, rasterized=True)

# # OR715104 star marker
or715 = next((n for n in tips if n.name == "OR715104"), None)
if or715:
    ax_A.scatter(or715.x, or715.y, s=22, marker="*",
                 c=INDIA_COL, edgecolors="#0A1F5C", linewidths=0.3,
                 zorder=5, rasterized=False)

#Kenya 2024-26 clade highlight
if clade_root:
    ct = [t for t in tips if t.name in clade_set]
    ylo, yhi = min(t.y for t in ct)-1, max(t.y for t in ct)+1
    all_xs = [n.x for n in all_nodes if n.x]
    xmax_d = max(all_xs)
    ax_A.fill_betweenx([ylo, yhi], clade_root.x-0.2, xmax_d+0.1,
                       alpha=0.08, color="#E87722", zorder=0, linewidth=0)
    ax_A.annotate("Kenya\n2024–26",
                  xy=(xmax_d+0.1, (ylo+yhi)/2),
                  xytext=(xmax_d+0.8, (ylo+yhi)/2),
                  fontsize=FS["annot"], fontfamily=FONT,
                  va="center", color="#B03A2E",
                  arrowprops=dict(arrowstyle="-", color="#B03A2E", lw=0.4))

#x-axis: full date range, ticks every 10 years
all_xs    = [n.x for n in all_nodes if n.x]
xmin_data = min(all_xs)
xmax_data = max(all_xs)
ax_A.set_xlim(xmin_data - 1, xmax_data + 2.5)
ax_A.set_ylim(-2, N_tips + 2)
x0 = int(xmin_data / 10) * 10
xticks = list(range(x0, int(xmax_data) + 12, 10))
ax_A.set_xticks(xticks)
ax_A.set_xticklabels([str(y) for y in xticks],
                     fontsize=FS["axis"], fontfamily=FONT)
ax_A.set_yticks([])
ax_A.set_xlabel("Year", fontsize=FS["axis"], fontfamily=FONT, labelpad=2)
ax_A.spines[["top","right","left"]].set_visible(False)
ax_A.spines["bottom"].set_linewidth(0.4)
ax_A.tick_params(axis="x", length=2, width=0.4, pad=2)

#IOL emergence annotation
ax_A.annotate("IOL emergence\n~2004",
              xy=(2004.5, N_tips*0.92),
              xytext=(1996, N_tips*0.96),
              fontsize=FS["annot"], fontfamily=FONT, color="#555555",
              arrowprops=dict(arrowstyle="-|>", color="#555555",
                              lw=0.4, mutation_scale=5))

#Panel A legend
handles_A = [
    mpatches.Patch(fc=REGION_COL["Asia"],          label="Asia",              lw=0),
    mpatches.Patch(fc=REGION_COL["Africa"],        label="Africa (non-Kenya)",lw=0),
    mpatches.Patch(fc=REGION_COL["Europe"],        label="Europe",            lw=0),
    mpatches.Patch(fc=REGION_COL["South America"], label="South America",     lw=0),
    mpatches.Patch(fc=REGION_COL["Indian Ocean"],  label="Indian Ocean",      lw=0),
    Line2D([0],[0], marker="o", lw=0, mfc=KENYA_COL,
           mec="none", ms=4, label="Kenya (all outbreaks)"),
    Line2D([0],[0], marker="*", lw=0, mfc=INDIA_COL,
           mec="#0A1F5C", mew=0.3, ms=6, label="OR715104 (India, Nov 2023)"),
]
ax_A.legend(handles=handles_A, fontsize=FS["legend"], frameon=False,
            loc="lower right", handlelength=0.9,
            prop={"family": FONT, "size": FS["legend"]})

ax_A.text(-0.04, 1.005, "A", transform=ax_A.transAxes,
          fontsize=FS["panel"], fontweight="bold", fontfamily=FONT, va="top")

#PANEL B — zoomed Kenya clade plus sister tip
#============================================
#Panel B = the clade plus its sister lineage
b_root = clade_root.parent if (clade_root and clade_root.parent) else clade_root
b_set  = set()
def cb(n): b_set.add(n.name); [cb(c) for c in n.children]
cb(b_root)

b_all   = [n for n in all_nodes if n.name in b_set]
b_tips  = [n for n in b_all if not n.children]
b_itern = [n for n in b_all if n.children]

#local y positions
b_sorted = sorted(b_tips, key=lambda n: n.y)
loc_y    = {n.name: i for i, n in enumerate(b_sorted)}
N_B      = len(b_sorted)

def bly(n):
    if not n.children:
        return loc_y.get(n.name, n.y)
    ch = [c for c in n.children if c.name in b_set]
    return np.mean([bly(c) for c in ch]) if ch else n.y

#year gridlines
b_xs  = [n.x for n in b_all if n.x]
xlo_B = min(b_xs) - 0.5
xhi_B = max(b_xs) + 0.4
for yr in range(int(xlo_B), int(xhi_B) + 2):
    ax_B.axvline(yr, color="#E8E8E8", lw=0.5, zorder=0)

#branches
BLW = 1.3
for n in b_all:
    if n.parent is None or n.parent.name not in b_set: continue
    ax_B.plot([n.parent.x, n.x], [bly(n), bly(n)],
              lw=BLW, color=CLADE_BRANCH, zorder=2,
              solid_capstyle="butt", rasterized=False)
for n in b_itern:
    ch = [c for c in n.children if c.name in b_set]
    if len(ch) > 1:
        ys = [bly(c) for c in ch]
        ax_B.plot([n.x, n.x], [min(ys), max(ys)],
                  lw=BLW, color=CLADE_BRANCH, zorder=2,
                  solid_capstyle="butt", rasterized=False)

#tips coloured circles
TIP_S = 16
for t in b_tips:
    if t.name == "OR715104":
        fc, ec = INDIA_COL, "#0A1F5C"
    else:
        div = gmeta(t.name, "division")
        fc  = DIV_COL.get(div, "#888888")
        ec  = "#333333"
    ax_B.scatter(t.x, bly(t), s=TIP_S, c=fc,
                 edgecolors=ec, linewidths=0.25,
                 zorder=4, rasterized=False)

# OR715104 label (right of tip)
if or715:
    ort = next((n for n in b_tips if n.name == "OR715104"), None)
    if ort:
        ax_B.text(ort.x + 0.07, bly(ort), "OR715104 (India, Nov 2023)",
                  fontsize=3.8, fontfamily=FONT, va="center",
                  ha="left", color=INDIA_COL, zorder=6)

#branch mutation labels
#Label the branches inside panel B carrying an amino-acid change, ranked by how much
#of the clade they subtend, so the figure follows the tree rather than a fixed list.
MIN_TIPS_TO_LABEL = 4
MAX_LABELS        = 3
_cand = []
for _bn in b_set:
    _m = aa_muts.get(_bn, {}).get("aa_muts", {})
    _txt = "; ".join("%s: %s" % (g, ", ".join(v)) for g, v in _m.items() if v)
    if not _txt or _bn == b_root.name:
        continue
    _nd = next((x for x in b_all if x.name == _bn), None)
    _nt = len(subtree_tips(_nd))
    if _nt >= MIN_TIPS_TO_LABEL:
        _cand.append((_nt, _bn, _txt))
BRANCH_LABELS = {bn: txt for _nt, bn, txt in
                 sorted(_cand, key=lambda r: -r[0])[:MAX_LABELS]}
for nname, label in BRANCH_LABELS.items():
    mn = next((n for n in b_all if n.name == nname), None)
    if mn is None or mn.parent is None: continue
    mid_x = (mn.parent.x + mn.x) / 2
    ax_B.text(mid_x, bly(mn) + 1.0, label,
              fontsize=4.5, fontfamily=FONT,
              ha="center", va="bottom", color="#222222", zorder=5,
              bbox=dict(boxstyle="round,pad=0.10", fc="white",
                        ec="#CCCCCC", alpha=0.85, lw=0.3))

#sub-lineage annotations
def clade_annot(node_name, label, x_off=-1.2, y_off=0, ha="right"):
    sn = next((n for n in b_all if n.name == node_name), None)
    if sn is None: return
    sub = set()
    def sc(x): sub.add(x.name); [sc(c) for c in x.children]
    sc(sn)
    ys = [bly(n) for n in b_tips if n.name in sub]
    if not ys: return
    ymid = np.mean(ys)
    ty = min(max(ymid + y_off, 2.0), N_B - 2.0)   # keep the label inside the axes
    ax_B.annotate(label,
                  xy=(sn.x, ymid),
                  xytext=(sn.x + x_off, ty),
                  fontsize=FS["annot"], fontfamily=FONT, va="center", ha=ha,
                  color="#333333",
                  arrowprops=dict(arrowstyle="-|>", color="#777777",
                                  lw=0.45, mutation_scale=5))

#Sub-lineage annotations from the tree.
#
_coastal_mrca = mrca(coastal_tips)
_gar_basal    = mrca([t for t in garissa_tips
                      if t.name in ("nvrl-605", "nvrl-606", "nvrl-608", "nvrl-614")])
_gar_nested   = mrca([t for t in garissa_tips if t.name in ("nvrl-603", "nvrl-616")])

if _coastal_mrca is not None:
    clade_annot(_coastal_mrca.name, "Coastal sub-lineage", x_off=-0.35, y_off=45)
if _gar_basal is not None:
    clade_annot(_gar_basal.name, "Garissa cluster 1 (n=4)", x_off=-0.35, y_off=13)
if _gar_nested is not None:
    clade_annot(_gar_nested.name,
                "Garissa cluster 2 (n=2),\nnested within coastal", x_off=-0.35, y_off=25)

#single-introduction annotation
intro = clade_root
if intro:
    ax_B.annotate("Single introduction\ntMRCA Nov 2022\n(95% CI Apr 2022–Oct 2024)",
                  xy=(intro.x, bly(intro)),
                  xytext=(intro.x - 1.5, bly(intro) + 60),
                  fontsize=4.0, fontfamily=FONT, va="bottom", color="#444444",
                  arrowprops=dict(arrowstyle="-|>", color="#888888",
                                  lw=0.4, mutation_scale=4))

#axes
ax_B.set_xlim(xlo_B, xhi_B + 1.5)
ax_B.set_ylim(-4, N_B + 3)
yr_B = list(range(int(xlo_B)+1, int(xhi_B)+2))
ax_B.set_xticks(yr_B)
ax_B.set_xticklabels([str(y) for y in yr_B],
                     fontsize=FS["axis"], fontfamily=FONT)
ax_B.set_yticks([])
ax_B.set_xlabel("Date", fontsize=FS["axis"], fontfamily=FONT, labelpad=2)
ax_B.spines[["top","right","left"]].set_visible(False)
ax_B.spines["bottom"].set_linewidth(0.4)
ax_B.tick_params(axis="x", length=2, width=0.4, pad=2)

#Panel B legend
handles_B = [
    mpatches.Patch(fc=DIV_COL["Garissa"], ec="#333333", lw=0.3, label="Garissa"),
    mpatches.Patch(fc=DIV_COL["Lamu"],    ec="#333333", lw=0.3, label="Lamu"),
    mpatches.Patch(fc=DIV_COL["Mombasa"], ec="#333333", lw=0.3, label="Mombasa"),
    Line2D([0],[0], marker="o", lw=0, mfc=INDIA_COL,
           mec="#0A1F5C", mew=0.3, ms=4, label="India (OR715104)"),
]
ax_B.legend(handles=handles_B, fontsize=FS["legend"], frameon=False,
            loc="upper left", handlelength=0.8,
            prop={"family": FONT, "size": FS["legend"]})

ax_B.text(-0.05, 1.005, "B", transform=ax_B.transAxes,
          fontsize=FS["panel"], fontweight="bold", fontfamily=FONT, va="top")

#save fig
OUT_TIFF = os.path.join(BASE, "Figure_3_revised.tiff")
OUT_PDF  = os.path.join(BASE, "Figure_3_revised.pdf")
OUT_PNG  = os.path.join(BASE, "Figure_3_revised.png")
fig.savefig(OUT_TIFF, dpi=DPI, format="tiff",
            bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
fig.savefig(OUT_PDF,  dpi=DPI, format="pdf", bbox_inches="tight")
fig.savefig(OUT_PNG,  dpi=DPI, format="png", bbox_inches="tight")
plt.close(fig)
print(f"Saved:\n  {OUT_TIFF}\n  {OUT_PDF}\n  {OUT_PNG}")
