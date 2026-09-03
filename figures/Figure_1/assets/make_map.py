# -*- coding: utf-8 -*-
"""Figure 2 - Map of Kenya showing CHIKV outbreak sample sites + Somalia snapshot.
AJTMH, 600 dpi. Traditional (non-contested) Kenyan boundaries (GADM, de-facto admin)."""
import os, sys
import numpy as np
import geopandas as gpd
import shapely.geometry as sg
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle
from matplotlib.lines import Line2D

DATA = sys.argv[1] if len(sys.argv) > 1 else "."
OUT  = sys.argv[2] if len(sys.argv) > 2 else "."
DPI_PREVIEW = int(sys.argv[3]) if len(sys.argv) > 3 else 150

CRS = "EPSG:32737"  # UTM zone 37N

#colors (from geoloc.txt)
LAND      = "#e5f5f9"   # default county / Somalia fill
BORDER    = "#bdbdbd"   # county boundary
WATER     = "#2b8cbe"   # ocean + lakes
NAT_LINE  = "#636363"   # Kenya national outline
SOM_LINE  = "#969696"   # Somalia outline
CTY = {"Mombasa": "#756bb1", "Garissa": "#99d8c9", "Lamu": "#fc9272"}
MARK = {  # site: (lon, lat, color, label)
    "Hagadera": (40.37555876966763,  0.0103044377418767, "#e41a1c", "Hagadera Refugee Camp, Garissa"),
    "Lamu":     (40.900849860702365, -2.2710688221588584, "#7fc97f", "Lamu Hospital"),
    "Likoni":   (39.66160625482952,  -4.089228060207857, "#1230b6", "Likoni Hospital"),
}

#load data
ken = gpd.read_file(os.path.join(DATA, "gadm41_KEN_1.json")).to_crs(CRS)
som = gpd.read_file(os.path.join(DATA, "gadm41_SOM_0.json")).to_crs(CRS)
lk  = gpd.read_file(os.path.join(DATA, "lakes", "ne_10m_lakes.shp"))
lk  = lk[lk.intersects(sg.box(33.5, -5, 42.5, 5.5))].to_crs(CRS)
oc  = gpd.read_file(os.path.join(DATA, "ocean", "ne_10m_ocean.shp"))
oc  = oc.clip(sg.box(32, -7, 44, 7)).to_crs(CRS)

#styles
plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.linewidth": 0.8,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})

LAND_BG = "#f7f7f7"
fig, ax = plt.subplots(figsize=(6.5, 6.5))

#neutral land backdrop
kb = ken.total_bounds
W = kb[2] - kb[0]; H = kb[3] - kb[1]
ax.add_patch(Rectangle((kb[0] - W, kb[1] - H), 3 * W, 3 * H,
             facecolor=LAND_BG, edgecolor="none", zorder=0))
#Indian Ocean
oc.buffer(800).plot(ax=ax, facecolor=WATER, edgecolor="none", zorder=1)

#Somalia snapshot
som.plot(ax=ax, facecolor=LAND, edgecolor=SOM_LINE, linewidth=0.6, zorder=1)

#Kenya counties
def cty_color(n): return CTY.get(n, LAND)
ken["c"] = ken["NAME_1"].map(cty_color)
ken.plot(ax=ax, color=ken["c"], edgecolor=BORDER, linewidth=0.4, zorder=2)
#Kenya national outline
ken.dissolve().boundary.plot(ax=ax, color=NAT_LINE, linewidth=1.0, zorder=3)

#Lakes (water bodies) over land
lk.plot(ax=ax, facecolor=WATER, edgecolor=WATER, linewidth=0, zorder=4)

#markers
pts = gpd.GeoSeries([sg.Point(MARK[k][0], MARK[k][1]) for k in MARK],
                    crs="EPSG:4326").to_crs(CRS)
for (k, p) in zip(MARK, pts):
    ax.plot(p.x, p.y, marker="o", ms=9, mfc=MARK[k][2], mec="white",
            mew=1.2, zorder=6)

xmin, ymin, xmax, ymax = ken.total_bounds
padx = (xmax - xmin) * 0.06
pady = (ymax - ymin) * 0.06
ax.set_xlim(xmin - padx, xmax + padx + (xmax - xmin) * 0.12)
ax.set_ylim(ymin - pady * 1.6, ymax + pady)
ax.set_aspect("equal")

ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values():
    s.set_visible(True); s.set_linewidth(0.8); s.set_edgecolor("black")

#Somalia label
som_lbl = gpd.GeoSeries([sg.Point(41.2, 2.2)], crs="EPSG:4326").to_crs(CRS).iloc[0]
ax.text(som_lbl.x, som_lbl.y, "SOMALIA", fontsize=11, style="italic",
        color="#525252", ha="center", va="center", rotation=-38, zorder=7)
ken_lbl = gpd.GeoSeries([sg.Point(37.4, 0.8)], crs="EPSG:4326").to_crs(CRS).iloc[0]
ax.text(ken_lbl.x, ken_lbl.y, "KENYA", fontsize=13, weight="bold",
        color="#9e9e9e", ha="center", va="center", zorder=7)

#scale bar (km), lower-right over ocean
sb_len = 200_000  # 200 km in metres
sb_h = (ymax - ymin) * 0.009
sb_x0 = (xmax + (xmax - xmin) * 0.16) - sb_len
sb_y0 = (ymin - pady * 1.6) + (ymax - ymin) * 0.03
ax.add_patch(Rectangle((sb_x0 - sb_len * 0.10, sb_y0 - sb_h * 1.2),
             sb_len * 1.20, sb_h * 4.6, facecolor="white",
             edgecolor="#bdbdbd", lw=0.5, zorder=8))
for i in range(2):
    ax.add_patch(Rectangle((sb_x0 + i * sb_len / 2, sb_y0), sb_len / 2, sb_h,
                 facecolor=("black" if i == 0 else "white"),
                 edgecolor="black", lw=0.6, zorder=9))
for v, lab in [(0, "0"), (sb_len / 2, "100"), (sb_len, "200 km")]:
    ax.text(sb_x0 + v, sb_y0 + sb_h * 1.5, lab, fontsize=6, ha="center",
            va="bottom", zorder=9)

#north arrow
na_x = xmax + (xmax - xmin) * 0.07
na_y = ymax - (ymax - ymin) * 0.05
ax.annotate("N", xy=(na_x, na_y), xytext=(na_x, na_y - (ymax - ymin) * 0.07),
            arrowprops=dict(facecolor="black", width=3, headwidth=9),
            ha="center", va="center", fontsize=10, fontweight="bold", zorder=8)

#legend
site_handles = [
    Line2D([0], [0], marker="o", ls="", mfc=MARK[k][2], mec="white", mew=1.0,
           ms=9, label=MARK[k][3]) for k in ["Hagadera", "Lamu", "Likoni"]
]
cty_handles = [
    Line2D([0], [0], marker="s", ls="", mfc=CTY["Garissa"], mec=BORDER, ms=10, label="Garissa County"),
    Line2D([0], [0], marker="s", ls="", mfc=CTY["Lamu"], mec=BORDER, ms=10, label="Lamu County"),
    Line2D([0], [0], marker="s", ls="", mfc=CTY["Mombasa"], mec=BORDER, ms=10, label="Mombasa County"),
    Line2D([0], [0], marker="s", ls="", mfc=WATER, mec=WATER, ms=10, label="Water bodies"),
]
leg1 = ax.legend(handles=site_handles, title="Outbreak sample sites",
                 loc="lower left", bbox_to_anchor=(0.0, 0.10), frameon=True,
                 fontsize=7.5, handletextpad=0.4, borderpad=0.6)
leg1.get_title().set_fontsize(8); leg1.get_title().set_fontweight("bold")
leg1.get_frame().set_edgecolor("#bdbdbd"); leg1.get_frame().set_linewidth(0.6)
ax.add_artist(leg1)
leg2 = ax.legend(handles=cty_handles, title="Map features",
                 loc="upper left", bbox_to_anchor=(0.0, 0.99), frameon=True,
                 fontsize=7.5, handletextpad=0.4, borderpad=0.6)
leg2.get_title().set_fontsize(8); leg2.get_title().set_fontweight("bold")
leg2.get_frame().set_edgecolor("#bdbdbd"); leg2.get_frame().set_linewidth(0.6)

fig.tight_layout(pad=0.3)
base = os.path.join(OUT, "Figure_2_Map")
if "--final" in sys.argv:
    fig.savefig(base + ".tif", dpi=600, format="tiff",
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(base + ".png", dpi=600)
    fig.savefig(base + ".pdf")
    fig.savefig(base + ".eps")
    from PIL import Image
    im = Image.open(base + ".tif")
    print("Final saved to", OUT)
    print("TIFF px:", im.size, "-> inches:",
          round(im.size[0] / 600, 2), "x", round(im.size[1] / 600, 2), "@600dpi")
else:
    fig.savefig(base + "_preview.png", dpi=DPI_PREVIEW)
    print("preview saved", base + "_preview.png")
