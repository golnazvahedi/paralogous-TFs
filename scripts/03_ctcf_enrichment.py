#!/usr/bin/env python3
"""
03_ctcf_enrichment.py -- Bin the extended clustered-locus windows and quantify
CTCF enrichment from the ENCODE GM12878 narrowPeak file.

Two binning representations (per iteration19 CLAUDE.md, caveat iii):
  (A) FIXED-WIDTH bins (default 5 kb) across each raw window  -> per-locus x per-bin matrix
      (loci have different native spans, so fixed-bp bins do NOT align across loci;
       used only for the per-locus heatmap / raw description).
  (B) LENGTH-NORMALIZED metagene: LEFT 100 kb flank (F fixed-bp bins) + cluster BODY
      (B scaled bins) + RIGHT 100 kb flank (F fixed-bp bins). Relative-position axis
      -1..0 (left flank), 0..1 (body), 1..2 (right flank). This is what is averaged
      and compared across loci of different sizes.

CTCF metrics per bin:
  - n_summit    : number of CTCF peak summits falling in the bin (primary density)
  - cov_frac    : fraction of bin bp covered by CTCF peak intervals

Per-locus metrics (summit-based; peaks/kb):
  - peaks_per_kb_window, peaks_per_kb_body, peaks_per_kb_flank
  - obs_exp_window : (peaks/bp in window) / (genome-wide peaks/bp)
  - nearest_peak_to_tss_min / _mean (bp, over member TSSs)
  - n_peaks_between_members : summits strictly between the outermost member TSSs

Outputs:
  results/ctcf/ctcf_bins.<binsize>.matrix.tsv     (fixed-width, long format)
  results/ctcf/ctcf_metagene.matrix.tsv           (length-normalized, long format)
  results/ctcf/ctcf_per_locus_summary.tsv         (per-locus metrics + dup_class4)
  log/03_ctcf_enrichment.log
"""
from pathlib import Path
import csv
import gzip
import argparse
from collections import defaultdict
from bisect import bisect_left, bisect_right

ROOT = Path(__file__).resolve().parents[1]
PAD = 100_000
PRIMARY = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def load_peaks(path):
    """Return dict chrom -> (starts[], ends[], summits[]) sorted by start."""
    raw = defaultdict(list)
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            chrom = p[0]
            if chrom not in PRIMARY:
                continue
            start, end = int(p[1]), int(p[2])
            # narrowPeak col10 (index 9) = summit offset from start; -1 if absent
            off = int(p[9]) if len(p) > 9 and p[9] not in (".", "-1") else (end - start) // 2
            if off < 0:
                off = (end - start) // 2
            summit = start + off
            raw[chrom].append((start, end, summit))
    peaks = {}
    for c, lst in raw.items():
        lst.sort()
        starts = [x[0] for x in lst]
        ends = [x[1] for x in lst]
        summits = sorted(x[2] for x in lst)
        peaks[c] = (starts, ends, summits, lst)
    return peaks


def count_summits(peaks, chrom, s, e):
    """Count peak summits in [s, e)."""
    if chrom not in peaks:
        return 0
    summits = peaks[chrom][2]
    return bisect_left(summits, e) - bisect_left(summits, s)


def coverage_bp(peaks, chrom, s, e):
    """Total bp of [s,e) covered by peak intervals."""
    if chrom not in peaks:
        return 0
    lst = peaks[chrom][3]
    starts = peaks[chrom][0]
    # peaks with start < e; iterate a bounded slice
    hi = bisect_right(starts, e)
    cov = 0
    for i in range(hi):
        ps, pe = lst[i][0], lst[i][1]
        if pe <= s:
            continue
        cov += min(pe, e) - max(ps, s)
    return max(cov, 0)


def nearest_peak_dist(peaks, chrom, pos):
    """Distance from pos to nearest peak summit (0 if inside a peak interval)."""
    if chrom not in peaks:
        return None
    starts, ends, summits, lst = peaks[chrom]
    # inside any interval?
    hi = bisect_right(starts, pos)
    for i in range(max(0, hi - 1), min(len(lst), hi + 1)):
        if lst[i][0] <= pos < lst[i][1]:
            return 0
    # nearest summit
    j = bisect_left(summits, pos)
    best = None
    for k in (j - 1, j):
        if 0 <= k < len(summits):
            d = abs(summits[k] - pos)
            best = d if best is None else min(best, d)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peaks", default=str(ROOT / "inputs" / "encode" / "ENCFF796WRU.bed.gz"))
    ap.add_argument("--binsize", type=int, default=5000)
    ap.add_argument("--flank-bins", type=int, default=20, help="metagene bins per 100kb flank")
    ap.add_argument("--body-bins", type=int, default=20, help="metagene scaled bins across body")
    ap.add_argument("--subdir", default="", help="write CTCF outputs to results/ctcf/<subdir>")
    args = ap.parse_args()

    peaks = load_peaks(args.peaks)
    sizes = {}
    with open(ROOT / "inputs" / "hg38.chrom.sizes") as f:
        for line in f:
            c, s = line.rstrip("\n").split("\t")[:2]
            sizes[c] = int(s)

    # genome-wide expected rate = total summits / total primary bp (chroms present in peaks)
    total_peaks = sum(len(peaks[c][2]) for c in peaks)
    total_bp = sum(sizes[c] for c in peaks if c in sizes)
    exp_rate = total_peaks / total_bp  # peaks per bp

    loci = []
    with open(ROOT / "results" / "clustered_loci_windows.pm100kb.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            loci.append(row)

    out_dir = ROOT / "results" / "ctcf" / args.subdir if args.subdir else ROOT / "results" / "ctcf"
    out_dir.mkdir(parents=True, exist_ok=True)

    fixed_rows = []
    meta_rows = []
    summ_rows = []

    F = args.flank_bins
    B = args.body_bins

    for row in loci:
        cid = row["cluster_id"]
        dup = row["dup_class4"]
        chrom = row["chrom"]
        ws, we = int(row["win_start"]), int(row["win_end"])
        cs, ce = int(row["cluster_start"]), int(row["cluster_end"])
        tss = [int(x) for x in row["member_tss"].split(",") if x != ""]

        # ---- (A) fixed-width bins ----
        b = ws
        idx = 0
        while b < we:
            be = min(b + args.binsize, we)
            n = count_summits(peaks, chrom, b, be)
            cov = coverage_bp(peaks, chrom, b, be)
            fixed_rows.append({
                "cluster_id": cid, "dup_class4": dup, "chrom": chrom,
                "bin_index": idx, "bin_start": b, "bin_end": be,
                "n_summit": n, "cov_frac": round(cov / (be - b), 6) if be > b else 0,
            })
            b = be
            idx += 1

        # ---- (B) length-normalized metagene ----
        # left flank: [cs-PAD, cs) split into F bins; but respect clipped win_start
        left0 = cs - PAD
        flank_bp = PAD / F
        # left flank bins
        for i in range(F):
            bs = left0 + i * flank_bp
            bend = left0 + (i + 1) * flank_bp
            bs_c, be_c = max(int(bs), ws), min(int(bend), we)
            if be_c <= bs_c:
                n, cov, width = 0, 0, 0
            else:
                n = count_summits(peaks, chrom, bs_c, be_c)
                cov = coverage_bp(peaks, chrom, bs_c, be_c)
                width = be_c - bs_c
            meta_rows.append({
                "cluster_id": cid, "dup_class4": dup,
                "region": "left_flank", "meta_index": i,
                "rel_pos": round(-1 + (i + 0.5) / F, 4),
                "n_summit": n, "width_bp": width,
                "cov_frac": round(cov / width, 6) if width else 0,
                "density_per_kb": round(n / (width / 1000), 6) if width else 0,
            })
        # body bins (scaled)
        body_len = max(ce - cs, 1)
        bstep = body_len / B
        for i in range(B):
            bs = cs + i * bstep
            bend = cs + (i + 1) * bstep
            bs_c, be_c = int(bs), int(bend)
            if be_c <= bs_c:
                be_c = bs_c + 1
            n = count_summits(peaks, chrom, bs_c, be_c)
            cov = coverage_bp(peaks, chrom, bs_c, be_c)
            width = be_c - bs_c
            meta_rows.append({
                "cluster_id": cid, "dup_class4": dup,
                "region": "body", "meta_index": F + i,
                "rel_pos": round((i + 0.5) / B, 4),
                "n_summit": n, "width_bp": width,
                "cov_frac": round(cov / width, 6) if width else 0,
                "density_per_kb": round(n / (width / 1000), 6) if width else 0,
            })
        # right flank
        right0 = ce
        for i in range(F):
            bs = right0 + i * flank_bp
            bend = right0 + (i + 1) * flank_bp
            bs_c, be_c = max(int(bs), ws), min(int(bend), we)
            if be_c <= bs_c:
                n, cov, width = 0, 0, 0
            else:
                n = count_summits(peaks, chrom, bs_c, be_c)
                cov = coverage_bp(peaks, chrom, bs_c, be_c)
                width = be_c - bs_c
            meta_rows.append({
                "cluster_id": cid, "dup_class4": dup,
                "region": "right_flank", "meta_index": F + B + i,
                "rel_pos": round(1 + (i + 0.5) / F, 4),
                "n_summit": n, "width_bp": width,
                "cov_frac": round(cov / width, 6) if width else 0,
                "density_per_kb": round(n / (width / 1000), 6) if width else 0,
            })

        # ---- per-locus summary metrics ----
        win_len = we - ws
        body_len_real = ce - cs
        flank_len = win_len - body_len_real  # both flanks combined
        n_win = count_summits(peaks, chrom, ws, we)
        n_body = count_summits(peaks, chrom, cs, ce)
        n_flank = n_win - n_body
        peaks_per_kb_window = n_win / (win_len / 1000)
        peaks_per_kb_body = n_body / (body_len_real / 1000) if body_len_real > 0 else 0
        peaks_per_kb_flank = n_flank / (flank_len / 1000) if flank_len > 0 else 0
        obs_rate = n_win / win_len
        obs_exp = obs_rate / exp_rate if exp_rate > 0 else 0

        dists = [nearest_peak_dist(peaks, chrom, t) for t in tss]
        dists = [d for d in dists if d is not None]
        nearest_min = min(dists) if dists else ""
        nearest_mean = round(sum(dists) / len(dists), 1) if dists else ""

        if len(tss) >= 2:
            lo, hi = min(tss), max(tss)
            n_between = count_summits(peaks, chrom, lo, hi)
            inter_member_dist = hi - lo
            between_per_kb = round(n_between / (inter_member_dist / 1000), 5) if inter_member_dist > 0 else ""
        else:
            n_between = ""
            inter_member_dist = ""
            between_per_kb = ""

        summ_rows.append({
            "cluster_id": cid, "dup_class4": dup, "family": row["family"],
            "chrom": chrom, "n_members": row["n_members"], "members": row["members"],
            "span_bp": row["span_bp"], "win_len_bp": win_len, "clipped": row["clipped"],
            "n_peaks_window": n_win, "n_peaks_body": n_body, "n_peaks_flank": n_flank,
            "peaks_per_kb_window": round(peaks_per_kb_window, 5),
            "peaks_per_kb_body": round(peaks_per_kb_body, 5),
            "peaks_per_kb_flank": round(peaks_per_kb_flank, 5),
            "obs_exp_window": round(obs_exp, 4),
            "nearest_peak_to_tss_min": nearest_min,
            "nearest_peak_to_tss_mean": nearest_mean,
            "n_peaks_between_members": n_between,
            "inter_member_tss_dist_bp": inter_member_dist,
            "peaks_between_per_kb": between_per_kb,
        })

    # write outputs
    fp = out_dir / f"ctcf_bins.{args.binsize}.matrix.tsv"
    with open(fp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(fixed_rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(fixed_rows)

    mp = out_dir / "ctcf_metagene.matrix.tsv"
    with open(mp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(meta_rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(meta_rows)

    sp = out_dir / "ctcf_per_locus_summary.tsv"
    with open(sp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summ_rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(summ_rows)

    logname = f"03_ctcf_enrichment.{args.subdir}.log" if args.subdir else "03_ctcf_enrichment.log"
    logp = ROOT / "log" / logname
    with open(logp, "w") as f:
        f.write(f"peaks file: {args.peaks}\n")
        f.write(f"total CTCF summits (primary chroms): {total_peaks}\n")
        f.write(f"total primary bp: {total_bp}\n")
        f.write(f"genome-wide expected rate: {exp_rate:.6e} peaks/bp "
                f"({exp_rate*1000:.5f} peaks/kb)\n")
        f.write(f"n loci: {len(loci)}\n")
        f.write(f"fixed binsize: {args.binsize}; flank_bins={F}; body_bins={B}\n")
        f.write(f"outputs: {fp.name}, {mp.name}, {sp.name}\n")

    print(f"genome-wide expected CTCF rate: {exp_rate*1000:.5f} peaks/kb")
    print(f"wrote {fp}\n      {mp}\n      {sp}\n      {logp}")


if __name__ == "__main__":
    main()
