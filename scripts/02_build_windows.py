#!/usr/bin/env python3
"""
02_build_windows.py -- Extend each clustered TF locus by +/-100 kb.

Reuses (do NOT recompute) the iteration16 four-group catalog:
  inputs/clustered_paralogs.noKRABHOX.tsv   (cluster coordinates)
  inputs/TF_dup_2x2_classification.tsv      (per-gene dup_class4 -> cluster_id)

Only the two CLUSTERED classes are in scope:
  clustered_cis_ohnolog  vs  clustered_SSD_tandem

Output: results/clustered_loci_windows.pm100kb.tsv
"""
from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]
PAD = 100_000

CLUSTERED_CLASSES = {"clustered_cis_ohnolog", "clustered_SSD_tandem"}


def load_chrom_sizes(path):
    sizes = {}
    with open(path) as f:
        for line in f:
            c, s = line.rstrip("\n").split("\t")[:2]
            sizes[c] = int(s)
    return sizes


def load_cluster_class(path):
    """Map cluster_id -> dup_class4 (validated single-valued per cluster upstream)."""
    m = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            cid = row["cluster_id"].strip()
            if cid:
                m.setdefault(cid, row["dup_class4"])
    return m


def main():
    sizes = load_chrom_sizes(ROOT / "inputs" / "hg38.chrom.sizes")
    cls = load_cluster_class(ROOT / "inputs" / "TF_dup_2x2_classification.tsv")

    out_rows = []
    n_clipped = 0
    n_dropped_chrY = 0
    with open(ROOT / "inputs" / "clustered_paralogs.noKRABHOX.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            cid = row["cluster_id"]
            dup = cls.get(cid)
            if dup not in CLUSTERED_CLASSES:
                continue
            chrom = row["chrom"]
            # GM12878 is female -> no chrY ChIP signal; a chrY locus would read a
            # technical zero, not biology. Exclude chrY (default CTCF cell line).
            if chrom == "chrY":
                n_dropped_chrY += 1
                continue
            cstart = int(row["cluster_start"])
            cend = int(row["cluster_end"])
            win_start = cstart - PAD
            win_end = cend + PAD
            clipped = []
            if win_start < 0:
                win_start = 0
                clipped.append("left")
            csize = sizes.get(chrom)
            if csize is not None and win_end > csize:
                win_end = csize
                clipped.append("right")
            if clipped:
                n_clipped += 1
            out_rows.append({
                "cluster_id": cid,
                "dup_class4": dup,
                "family": row["family"],
                "parent_DBD_family": row["parent_DBD_family"],
                "chrom": chrom,
                "n_members": row["n_members"],
                "members": row["members"],
                "member_strands": row["member_strands"],
                "member_tss": row["member_tss"],
                "cluster_start": cstart,
                "cluster_end": cend,
                "span_bp": row["span_bp"],
                "win_start": win_start,
                "win_end": win_end,
                "win_len_bp": win_end - win_start,
                "clipped": ",".join(clipped) if clipped else "none",
            })

    out_rows.sort(key=lambda r: (r["dup_class4"], r["chrom"], r["win_start"]))
    cols = ["cluster_id", "dup_class4", "family", "parent_DBD_family", "chrom",
            "n_members", "members", "member_strands", "member_tss",
            "cluster_start", "cluster_end", "span_bp",
            "win_start", "win_end", "win_len_bp", "clipped"]
    outp = ROOT / "results" / "clustered_loci_windows.pm100kb.tsv"
    with open(outp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    from collections import Counter
    cc = Counter(r["dup_class4"] for r in out_rows)
    print(f"wrote {outp}")
    print(f"total clustered loci: {len(out_rows)}")
    for k, v in sorted(cc.items()):
        print(f"  {k}: {v}")
    print(f"windows clipped to chrom bounds: {n_clipped}")
    print(f"chrY loci dropped (GM12878 female, no chrY signal): {n_dropped_chrY}")


if __name__ == "__main__":
    main()
