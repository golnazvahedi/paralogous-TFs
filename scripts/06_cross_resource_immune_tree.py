#!/usr/bin/env python3
"""
Section 2 (broader view): vertebrate species with single-cell RNA-seq of immune
cells across ALL public resources (not only SCEA), on one time-calibrated tree.

This extends scripts/05 (SCEA-only) with species whose immune-cell scRNA-seq
lives outside SCEA: rhesus macaque (RIRA), mouse lemur (Tabula Microcebus), pig
(porcine lymphoid-organ immune atlas), cattle (Cattle Cell Atlas), horse (equine
PBMC atlas), and the cat/dog/goat/deer/hamster/tiger/pigeon/chicken/sheep
profiled in cross-species PBMC / livestock studies. The jawless-vertebrate
outgroup is lamprey: a 14-tissue cell atlas of Lethenteron reissneri (604k
cells/nuclei, CNGB CNP0005120; Nat Commun 2025, s41467-025-56153-w) plus an
immune-focused gill/blood/intestine lymphocyte study in Lampetra morii
(65k cells; Nat Commun 2024, s41467-024-51763-2). Both profile VLR lymphocytes,
granulocytes, monocytes/macrophages and DCs in lamprey hematopoietic tissues
(blood, supraneural body, gill, kidney, intestine), rooting the tree ~500 Mya.

Each tip is annotated with its resource and immune tissues, coloured by whether
the data are in SCEA or an external resource, and sized by tissue breadth
(atlas-grade = blood + secondary lymphoid + marrow; PBMC = blood only).

Topology and divergence times (Mya) are the accepted vertebrate phylogeny,
approximate values from TimeTree (Kumar et al. 2017, www.timetree.org).

Outputs
-------
results/figures/cross_resource_vertebrate_immune_tree.pdf / .png
results/intermediate/cross_resource_immune_species.tsv
log/06_cross_resource_immune_tree.log
"""
from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
FIGDIR = ROOT / "results" / "figures"; FIGDIR.mkdir(parents=True, exist_ok=True)
INT = ROOT / "results" / "intermediate"
LOG = ROOT / "log" / "06_cross_resource_immune_tree.log"
LOG.parent.mkdir(parents=True, exist_ok=True)
_logf = open(LOG, "w")
def log(*a):
    m = " ".join(str(x) for x in a); print(m); _logf.write(m + "\n")

# species metadata: scientific -> (common, in_scea, breadth, resource, tissues)
#   breadth: 'atlas' (blood+spleen/LN/marrow), 'limited' (1-2 immune studies),
#            'pbmc' (blood/PBMC only)
SP = {
 "Homo sapiens":          ("Human",        True,  "atlas",   "SCEA / Human Cell Atlas",        "blood, spleen, lymph node, marrow, thymus"),
 "Macaca mulatta":        ("Rhesus macaque", False,"atlas",   "RIRA atlas",                      "blood, lymph node, bone marrow"),
 "Microcebus murinus":    ("Mouse lemur",  False, "atlas",   "Tabula Microcebus",               "blood, spleen, bone marrow"),
 "Mus musculus":          ("Mouse",        True,  "atlas",   "SCEA / ImmGen, Tabula Muris",     "spleen, blood, marrow, lymph node, thymus"),
 "Rattus norvegicus":     ("Rat",          True,  "limited", "SCEA",                            "bone-marrow phagocytes, microglia"),
 "Mesocricetus auratus":  ("Hamster",      False, "pbmc",    "cross-species PBMC atlas",        "blood (PBMC)"),
 "Oryctolagus cuniculus": ("Rabbit",       True,  "limited", "SCEA / cross-species PBMC",       "bone-marrow phagocytes, PBMC"),
 "Sus scrofa":            ("Pig",          False, "atlas",   "porcine lymphoid-organ atlas",    "marrow, thymus, spleen, lymph node, PBMC"),
 "Bos taurus":            ("Cattle",       False, "atlas",   "Cattle Cell Atlas",               "PBMC, spleen"),
 "Ovis aries":            ("Sheep",        False, "limited", "ruminant scRNA-seq",              "PBMC / immune"),
 "Capra hircus":          ("Goat",         False, "pbmc",    "cross-species / ruminant",        "blood (PBMC)"),
 "Cervus":                ("Deer",         False, "pbmc",    "cross-species PBMC atlas",        "blood (PBMC)"),
 "Equus caballus":        ("Horse",        False, "pbmc",    "equine PBMC atlas",               "blood (PBMC, deep)"),
 "Canis lupus familiaris":("Dog",          False, "pbmc",    "cross-species PBMC atlas",        "blood (PBMC)"),
 "Felis catus":           ("Cat",          False, "pbmc",    "cross-species PBMC atlas",        "blood (PBMC)"),
 "Panthera tigris":       ("Tiger",        False, "pbmc",    "cross-species PBMC atlas",        "blood (PBMC)"),
 "Gallus gallus":         ("Chicken",      False, "limited", "avian immune scRNA-seq",          "spleen / immune"),
 "Columba livia":         ("Pigeon",       False, "pbmc",    "cross-species PBMC atlas",        "blood (PBMC)"),
 "Danio rerio":           ("Zebrafish",    True,  "atlas",   "SCEA / zebrafish atlases",        "kidney marrow (blood), spleen, blood"),
 "Lethenteron reissneri": ("Lamprey",      False, "atlas",   "lamprey cell atlas (CNP0005120) + L. morii immune scRNA", "blood, supraneural body, gill, kidney, intestine"),
}

# --- time-calibrated tree (ages Mya, tips at 0) -------------------------------
class N:
    def __init__(self, age, name=None, children=None, clade=None):
        self.age, self.name, self.children, self.clade = age, name, children or [], clade
def L(sp): return N(0, name=sp)
def I(age, ch, clade=None): return N(age, children=ch, clade=clade)

# Primates
human, macaque, lemur = L("Homo sapiens"), L("Macaca mulatta"), L("Microcebus murinus")
catarrhine = I(29, [human, macaque])
primates   = I(74, [catarrhine, lemur], "Primates")
# Glires
mouse, rat, hamster = L("Mus musculus"), L("Rattus norvegicus"), L("Mesocricetus auratus")
murinae   = I(21, [mouse, rat])
muroidea  = I(25, [murinae, hamster])
rabbit    = L("Oryctolagus cuniculus")
glires    = I(82, [muroidea, rabbit], "Glires")
euarch    = I(90, [primates, glires], "Euarchontoglires")
# Laurasiatheria: Carnivora + (Perissodactyla + Artiodactyla)
dog, cat, tiger = L("Canis lupus familiaris"), L("Felis catus"), L("Panthera tigris")
felidae   = I(11, [cat, tiger])
carnivora = I(54, [dog, felidae], "Carnivora")
horse     = L("Equus caballus")
pig       = L("Sus scrofa")
cattle, sheep, goat, deer = L("Bos taurus"), L("Ovis aries"), L("Capra hircus"), L("Cervus")
caprini   = I(7, [sheep, goat])
bovidae   = I(24, [cattle, caprini])
pecora    = I(28, [bovidae, deer])
artiodactyla = I(60, [pig, pecora], "Artiodactyla")
euungulata = I(76, [horse, artiodactyla])
laurasia  = I(78, [carnivora, euungulata], "Laurasiatheria")
boreo     = I(96, [euarch, laurasia], "Boreoeutheria")
# Birds
chicken, pigeon = L("Gallus gallus"), L("Columba livia")
birds     = I(93, [chicken, pigeon], "Aves")
amniota   = I(319, [boreo, birds], "Amniota")
zebrafish = L("Danio rerio")
euteleostomi = I(429, [amniota, zebrafish], "Euteleostomi")
# Cyclostome outgroup: lamprey (jawless vertebrate) roots the whole tree
lamprey   = L("Lethenteron reissneri")
root      = I(500, [euteleostomi, lamprey], "Vertebrata")

# layout
order = []
def collect(n):
    if not n.children: order.append(n)
    else: [collect(c) for c in n.children]
collect(root)
yof = {}
for i, lf in enumerate(order):
    yof[lf] = len(order) - 1 - i
def yset(n):
    if not n.children: return yof[n]
    ys = [yset(c) for c in n.children]; yof[n] = sum(ys) / len(ys); return yof[n]
yset(root)

# write table
rows = [{"species": s, "common": m[0], "in_SCEA": m[1], "coverage": m[2],
         "resource": m[3], "immune_tissues": m[4]} for s, m in SP.items()]
with open(INT / "cross_resource_immune_species.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)
log("species total:", len(SP), "| in SCEA:", sum(m[1] for m in SP.values()),
    "| external:", sum(not m[1] for m in SP.values()))

# --- draw ---------------------------------------------------------------------
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
fig, ax = plt.subplots(figsize=(12.5, 9.2))
SCEA_C, EXT_C = "#1f7a8c", "#e07b39"
SIZE = {"atlas": 150, "limited": 80, "pbmc": 48}

def draw(n):
    for c in n.children:
        ax.plot([n.age, c.age], [yof[c], yof[c]], color="#666666", lw=1.4, zorder=1)
        draw(c)
    if n.children:
        ys = [yof[c] for c in n.children]
        ax.plot([n.age, n.age], [min(ys), max(ys)], color="#666666", lw=1.4, zorder=1)
draw(root)

LABEL_NODES = [root, euteleostomi, amniota, boreo, laurasia, euarch, glires,
               primates, carnivora, artiodactyla, birds]
for nd in LABEL_NODES:
    ax.plot(nd.age, yof[nd], "o", ms=3.5, color="#333333", zorder=3)
    ax.annotate(f"{int(nd.age)}", (nd.age, yof[nd]), textcoords="offset points",
                xytext=(2, 4), fontsize=7, color="#333333")
    if nd.clade:
        ax.annotate(nd.clade, (nd.age, yof[nd]), textcoords="offset points",
                    xytext=(3, -9), fontsize=6.8, style="italic", color="#999999")

for lf in order:
    common, in_scea, breadth, resource, tissues = SP[lf.name]
    col = SCEA_C if in_scea else EXT_C
    y = yof[lf]
    ax.scatter([0], [y], s=SIZE[breadth], color=col, zorder=4,
               edgecolor="white", linewidth=0.8)
    sci = lf.name.replace(" ", r"\ ")
    ax.annotate(r"%s  ($\it{%s}$)" % (common, sci), (0, y),
                textcoords="offset points", xytext=(13, 3.5), fontsize=9,
                fontweight="bold")
    ax.annotate(f"{resource}  -  {tissues}", (0, y), textcoords="offset points",
                xytext=(13, -8), fontsize=6.8, color=col, style="italic")

ax.set_xlim(525, -185)
ax.set_ylim(-0.8, len(order) - 0.2)
ax.set_yticks([])
for s in ("top", "left", "right"): ax.spines[s].set_visible(False)
ax.set_xlabel("Million years ago (approximate divergence time, TimeTree)", fontsize=9)
ax.set_xticks([500, 429, 319, 96, 90, 82, 78, 60, 0])
ax.tick_params(axis="x", labelsize=7.5)

fig.text(0.5, 0.965, "Vertebrates with single-cell RNA-seq of immune cells "
         "across public resources", ha="center", fontsize=13, fontweight="bold")
fig.text(0.5, 0.93, f"{len(SP)} species spanning {int(root.age)} Myr of evolution"
         "   |   teal = data in EBI SCEA, orange = external resource   |   "
         "bubble size = tissue breadth (atlas / limited / PBMC-only)",
         ha="center", fontsize=8.4, color="#555555")

legend = [
    Line2D([0],[0], marker="o", color="w", markerfacecolor=SCEA_C, markersize=11, label="in EBI SCEA"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor=EXT_C, markersize=11, label="external resource"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor="#888888", markersize=12, label="atlas: blood+spleen/LN/marrow"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor="#888888", markersize=8.5, label="limited (1-2 studies)"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor="#888888", markersize=6, label="PBMC / blood only"),
]
ax.legend(handles=legend, loc="lower left", frameon=False, fontsize=7.8,
          ncol=1, handletextpad=0.4, labelspacing=0.7)

fig.subplots_adjust(top=0.90, bottom=0.07, left=0.03, right=0.99)
for ext in ("pdf", "png"):
    fig.savefig(FIGDIR / f"cross_resource_vertebrate_immune_tree.{ext}", dpi=200,
                bbox_inches="tight")
log("wrote", FIGDIR / "cross_resource_vertebrate_immune_tree.pdf")
log("wrote", INT / "cross_resource_immune_species.tsv")
_logf.close()
