#!/usr/bin/env python3
"""
Section 2: Vertebrates in the EBI Single Cell Expression Atlas (SCEA) that have
single-cell RNA-seq of immune cells, drawn on a time-calibrated species tree.

Inputs
------
inputs/scea_experiments.json   (full SCEA experiment list, fetched from
                                https://www.ebi.ac.uk/gxa/sc/json/experiments)

What it does
------------
1. Reads every SCEA experiment, counts per-species experiments and the subset
   whose description matches an immune / haematopoietic keyword.
2. Restricts to vertebrate species present in SCEA and flags which have >=1
   immune-cell dataset (Human, Mouse, Rat, Rabbit, Zebrafish). The other SCEA
   vertebrates (Marmoset, Sheep, Chicken, Xenopus) have only non-immune datasets
   and are shown greyed for phylogenetic context.
3. Draws an ultrametric, time-calibrated species tree. Topology and divergence
   times (Mya) are the accepted vertebrate phylogeny (TimeTree / Kumar et al.
   2017, www.timetree.org): Murinae 21, Glires 82, Euarchontoglires 90,
   Boreoeutheria 96, Amniota 319, Tetrapoda 352, Euteleostomi (bony-vertebrate
   crown) 429 Mya.

Outputs
-------
results/figures/scea_vertebrate_immune_tree.pdf / .png
results/intermediate/scea_vertebrate_immune_experiments.tsv
log/05_scea_vertebrate_immune_tree.log
"""
import json, re
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
JSON = ROOT / "inputs" / "scea_experiments.json"
FIGDIR = ROOT / "results" / "figures"; FIGDIR.mkdir(parents=True, exist_ok=True)
INT = ROOT / "results" / "intermediate"
LOG = ROOT / "log" / "05_scea_vertebrate_immune_tree.log"
LOG.parent.mkdir(parents=True, exist_ok=True)
_logf = open(LOG, "w")
def log(*a):
    m = " ".join(str(x) for x in a); print(m); _logf.write(m+"\n")

# --- 1. count immune-related experiments per species --------------------------
IMM = re.compile(r"\b("
    r"immune|immuno|blood|pbmc|spleen|splenic|lymph node|lymphoid|lymphocyte|"
    r"leukocyte|leucocyte|haematopoi|hematopoi|bone marrow|kidney marrow|marrow|"
    r"thymus|thymic|T cell|B cell|T-cell|B-cell|NK cell|natural killer|"
    r"innate lymphoid|ILC|myeloid|monocyte|mononuclear phagocyte|macrophage|"
    r"microglia|neutrophil|granulocyte|erythroid|erythrocyte|thrombocyte|"
    r"megakaryocyte|dendritic|CD4|CD8|CD41|LCK)\b", re.I)

exps = json.load(open(JSON))["experiments"]
from collections import defaultdict
n_all = defaultdict(int); imm = defaultdict(list)
for e in exps:
    sp = e["species"].strip()
    sp = {"Mus Musculus": "Mus musculus"}.get(sp, sp)   # fix one casing dup
    n_all[sp] += 1
    if IMM.search(e["experimentDescription"] or ""):
        imm[sp].append((e["experimentAccession"], e["experimentDescription"]))

# --- 2. vertebrate species in SCEA + metadata ---------------------------------
# common name, example immune tissues (from dataset descriptions)
VERT = {
    "Homo sapiens":          ("Human",     "blood/PBMC, spleen, lymph node, bone marrow, thymus"),
    "Callithrix jacchus":    ("Marmoset",  "embryo only - no immune dataset"),
    "Mus musculus":          ("Mouse",     "spleen, blood, bone marrow, lymph node, thymus"),
    "Rattus norvegicus":     ("Rat",       "bone-marrow mononuclear phagocytes, microglia"),
    "Oryctolagus cuniculus": ("Rabbit",    "bone-marrow mononuclear phagocytes"),
    "Ovis aries":            ("Sheep",     "gut epithelium only - no immune dataset"),
    "Gallus gallus":         ("Chicken",   "limb / neural crest / embryo - no immune dataset"),
    "Xenopus tropicalis":    ("Frog (Xenopus)", "spinal cord only - no immune dataset"),
    "Danio rerio":           ("Zebrafish", "kidney marrow (blood), spleen, blood, ILC-like"),
}

rows = []
for sp, (common, tissue) in VERT.items():
    n_imm = len(imm.get(sp, []))
    rows.append({"species": sp, "common": common, "n_experiments": n_all.get(sp, 0),
                 "n_immune_experiments": n_imm, "has_immune": n_imm > 0,
                 "example_immune_tissues": tissue,
                 "immune_accessions": ";".join(a for a, _ in imm.get(sp, []))})
import csv
with open(INT / "scea_vertebrate_immune_experiments.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)
meta = {r["species"]: r for r in rows}
log("SCEA vertebrate species:", len(VERT),
    "| with immune-cell scRNA-seq:", sum(r["has_immune"] for r in rows))
for r in sorted(rows, key=lambda x: -x["n_immune_experiments"]):
    log(f"  {r['common']:14s} immune_exps={r['n_immune_experiments']:<3d} "
        f"total={r['n_experiments']:<4d} {'<-- immune' if r['has_immune'] else ''}")

# --- 3. time-calibrated species tree (ages in Mya; tips at 0) -----------------
class N:
    def __init__(self, age, name=None, children=None, clade=None):
        self.age, self.name, self.children, self.clade = age, name, children or [], clade

def leaf(sp): return N(0, name=sp)
def inode(age, ch, clade): return N(age, children=ch, clade=clade)

human, marmoset = leaf("Homo sapiens"), leaf("Callithrix jacchus")
primates = inode(43, [human, marmoset], "Primates")
mouse, rat = leaf("Mus musculus"), leaf("Rattus norvegicus")
murinae = inode(21, [mouse, rat], None)
rabbit = leaf("Oryctolagus cuniculus")
glires = inode(82, [murinae, rabbit], "Glires")
euarch = inode(90, [primates, glires], "Euarchontoglires")
sheep = leaf("Ovis aries")
boreo = inode(96, [euarch, sheep], "Boreoeutheria")
chicken = leaf("Gallus gallus")
amniota = inode(319, [boreo, chicken], "Amniota")
xenopus = leaf("Xenopus tropicalis")
tetrapoda = inode(352, [amniota, xenopus], "Tetrapoda")
zebrafish = leaf("Danio rerio")
root = inode(429, [tetrapoda, zebrafish], "Euteleostomi")

# assign y by leaf order (top->bottom), x = age
order = []
def collect(n):
    if not n.children: order.append(n)
    else: [collect(c) for c in n.children]
collect(root)
yof = {}
for i, lf in enumerate(order):
    yof[lf] = len(order) - 1 - i      # top species highest y
def yset(n):
    if not n.children: return yof[n]
    ys = [yset(c) for c in n.children]
    yof[n] = sum(ys) / len(ys); return yof[n]
yset(root)

# --- draw ---------------------------------------------------------------------
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
fig, ax = plt.subplots(figsize=(11.5, 6.2))
IMMUNE_C, GREY_C = "#1f7a8c", "#b0b0b0"

def draw(n):
    for c in n.children:
        # horizontal branch from parent age -> child age at child y
        col = "#555555"
        ax.plot([n.age, c.age], [yof[c], yof[c]], color=col, lw=1.6, zorder=1)
        draw(c)
    if n.children:
        ys = [yof[c] for c in n.children]
        ax.plot([n.age, n.age], [min(ys), max(ys)], color="#555555", lw=1.6, zorder=1)
draw(root)

# node age labels for key splits
NODE_LABELS = [euarch, glires, murinae, boreo, amniota, tetrapoda, root, primates]
for nd in NODE_LABELS:
    lbl = f"{int(nd.age)}"
    ax.plot(nd.age, yof[nd], "o", ms=4, color="#333333", zorder=3)
    ax.annotate(lbl, (nd.age, yof[nd]), textcoords="offset points",
                xytext=(2, 4), fontsize=7.5, color="#333333")
    if nd.clade:
        ax.annotate(nd.clade, (nd.age, yof[nd]), textcoords="offset points",
                    xytext=(4, -9), fontsize=7, style="italic", color="#777777")

# tips: marker sized by immune-exp count + label
for lf in order:
    r = meta[lf.name]
    is_imm = r["has_immune"]
    col = IMMUNE_C if is_imm else GREY_C
    y = yof[lf]
    if is_imm:
        # bubble scaled by sqrt of count
        size = 40 + 38 * (r["n_immune_experiments"] ** 0.5)
        ax.scatter([0], [y], s=size, color=col, zorder=4,
                   edgecolor="white", linewidth=0.8)
        ax.annotate(str(r["n_immune_experiments"]), (0, y), ha="center", va="center",
                    fontsize=7, color="white", fontweight="bold", zorder=5)
    else:
        ax.scatter([0], [y], s=30, facecolor="white", edgecolor=col,
                   linewidth=1.4, zorder=4)
    sci = lf.name.replace(" ", r"\ ")
    name = r"%s  ($\it{%s}$)" % (r["common"], sci)
    ax.annotate(name, (0, y), textcoords="offset points", xytext=(12, 3.5),
                fontsize=9.5, fontweight="bold" if is_imm else "normal",
                color="black" if is_imm else "#888888")
    tcol = "#1f7a8c" if is_imm else "#a8a8a8"
    ax.annotate(r["example_immune_tissues"], (0, y), textcoords="offset points",
                xytext=(12, -8.5), fontsize=7.2, color=tcol, style="italic")

ax.set_xlim(445, -120)             # past(left) -> present(right); room for labels
ax.set_ylim(-0.8, len(order) - 0.2)
ax.set_yticks([])
for s in ("top", "left", "right"): ax.spines[s].set_visible(False)
ax.set_xlabel("Million years ago (divergence time, TimeTree)", fontsize=9)
ax.set_xticks([429, 352, 319, 96, 90, 82, 43, 21, 0])
ax.tick_params(axis="x", labelsize=7.5)

n_imm_sp = sum(r["has_immune"] for r in rows)
fig.text(0.5, 0.965, "Vertebrates with single-cell RNA-seq of immune cells in "
         "the EBI Single Cell Expression Atlas",
         ha="center", fontsize=12.5, fontweight="bold")
fig.text(0.5, 0.925,
         f"{n_imm_sp} of {len(VERT)} SCEA vertebrate species carry an immune-cell "
         f"dataset   |   filled bubble = # immune scRNA-seq experiments;   "
         f"open grey = in SCEA but no immune dataset",
         ha="center", fontsize=8.3, color="#555555")

legend = [
    Line2D([0],[0], marker="o", color="w", markerfacecolor=IMMUNE_C, markersize=10,
           label="immune-cell scRNA-seq present"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor="white",
           markeredgecolor=GREY_C, markersize=9, label="in SCEA, no immune dataset"),
]
ax.legend(handles=legend, loc="lower left", frameon=False, fontsize=8)

fig.subplots_adjust(top=0.88, bottom=0.10, left=0.04, right=0.98)
for ext in ("pdf", "png"):
    fig.savefig(FIGDIR / f"scea_vertebrate_immune_tree.{ext}", dpi=200,
                bbox_inches="tight")
log("wrote", FIGDIR / "scea_vertebrate_immune_tree.pdf")
log("wrote", INT / "scea_vertebrate_immune_experiments.tsv")
_logf.close()
