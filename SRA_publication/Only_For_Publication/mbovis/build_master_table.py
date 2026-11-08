#!/usr/bin/env python3
"""
Build a single master sample table for the M. bovis / M. caprae study by
merging every analysis run so far.

Sources (all optional -- missing ones are skipped with a warning):

  1. Bactopia assembly-scan .............. assembly size, contigs, N50, GC
  2. TB-Profiler ......................... lineage, sublineage, depth,
                                           % mapped, resistance
  3. vSNP3 step 1 stats .................. spoligotype octal / binary / SB,
                                           read counts, zero-coverage,
                                           quality SNPs
  4. vSNP3 step 2 groups ................. clonal complex membership
                                           (read from the per-group FASTAs)
  5. RD typing ........................... RD9 / RD4 / La2-deletion / RD1
  6. Snippy SNP distance matrix .......... single-linkage cluster membership
                                           at a stated SNP threshold

Writes master_sample_table.csv and .xlsx.

Usage:
    python3 build_master_table.py [SNP_THRESHOLD]      # default 12
"""

import glob
import os
import sys

import pandas as pd

# --------------------------------------------------------------------------
# paths
# --------------------------------------------------------------------------
BACT = "/scratch/ma95362/Bio_project_publication/Mbovis"
VSNP = "/scratch/ma95362/vsnp3"
RDIR = "/scratch/ma95362/rd_typing"
OUTDIR = "/scratch/ma95362/Bio_project_publication/Mbovis/all_merged_results"

SNP_THRESHOLD = int(sys.argv[1]) if len(sys.argv) > 1 else 12


def first(pattern):
    """Most recent file matching a glob, or None."""
    hits = sorted(glob.glob(pattern))
    return hits[-1] if hits else None


def warn(msg):
    print(f"  [skip] {msg}")


frames = {}

# --------------------------------------------------------------------------
# 1. assembly-scan
# --------------------------------------------------------------------------
p = first(f"{BACT}/bactopia-runs/bactopia-*/merged-results/assembly-scan.tsv")
if p:
    d = pd.read_csv(p, sep="\t")
    keep = {
        "sample": "sample",
        "total_contig": "asm_contigs",
        "total_contig_length": "asm_length_bp",
        "n50_contig_length": "asm_n50",
        "contig_percent_a": "asm_pct_A",
        "contig_percent_t": "asm_pct_T",
    }
    d = d[[c for c in keep if c in d.columns]].rename(columns=keep)
    # GC from A/T, which is the contamination signal
    if "asm_pct_A" in d and "asm_pct_T" in d:
        d["asm_pct_GC"] = (100 - d.asm_pct_A - d.asm_pct_T).round(2)
    frames["assembly"] = d
    print(f"  assembly-scan: {len(d)} samples")
else:
    warn("assembly-scan.tsv not found")

# --------------------------------------------------------------------------
# 2. TB-Profiler
# --------------------------------------------------------------------------
p = first(f"{BACT}/bactopia-runs/tbprofiler-*/merged-results/tbprofiler.csv")
if p:
    d = pd.read_csv(p)
    keep = {
        "sample": "sample",
        "main_lineage": "tbp_lineage",
        "sub_lineage": "tbp_sublineage",
        "drtype": "tbp_drtype",
        "target_median_depth": "tbp_median_depth",
        "pct_reads_mapped": "tbp_pct_mapped",
        "num_dr_variants": "tbp_n_dr_variants",
        "num_other_variants": "tbp_n_other_variants",
        "pyrazinamide": "tbp_pyrazinamide",
    }
    d = d[[c for c in keep if c in d.columns]].rename(columns=keep)
    frames["tbprofiler"] = d
    print(f"  tbprofiler: {len(d)} samples")
else:
    warn("tbprofiler.csv not found")

# --------------------------------------------------------------------------
# 3. vSNP3 step 1 stats (spoligotypes)
# --------------------------------------------------------------------------
p = f"{VSNP}/step1_stats_all.csv"
if not os.path.exists(p):
    # rebuild from the per-sample spreadsheets
    sheets = sorted(glob.glob(f"{VSNP}/step1/*/*_stats.xlsx"))
    if sheets:
        pd.concat([pd.read_excel(f) for f in sheets], ignore_index=True).to_csv(p, index=False)
        print(f"  rebuilt {p} from {len(sheets)} spreadsheets")

if os.path.exists(p):
    d = pd.read_csv(p)
    keep = {
        "sample": "sample",
        "Spoligotype Octal Code": "spoligo_octal",
        "Spoligotype Binary Code": "spoligo_binary",
        "Spoligotype SB Number": "spoligo_SB",
        "Average Depth": "vsnp_mean_depth",
        "Unmapped Percent": "vsnp_pct_unmapped",
        "Percent Ref with Zero Coverage": "vsnp_pct_ref_zero_cov",
        "Quality SNPs": "vsnp_quality_snps",
        "Ambiguous SNPs": "vsnp_ambiguous_snps",
        "R1 Read Count": "reads_R1",
        "R1 Ave Length": "read_length",
    }
    d = d[[c for c in keep if c in d.columns]].rename(columns=keep)
    d["sample"] = d["sample"].astype(str)
    frames["spoligo"] = d
    print(f"  vsnp3 stats: {len(d)} samples")
else:
    warn("vSNP3 step1 stats not found")

# --------------------------------------------------------------------------
# 4. vSNP3 step 2 group membership
# --------------------------------------------------------------------------
rows = {}
for fa in glob.glob(f"{VSNP}/step2/*/*.fasta"):
    grp = os.path.basename(os.path.dirname(fa))
    for line in open(fa):
        if line.startswith(">"):
            s = line[1:].strip().replace(".vcf", "")
            if s.lower() == "root":
                continue
            rows.setdefault(s, []).append(grp)

if rows:
    def pick(groups):
        """Most specific group: drop the all-samples container."""
        g = [x for x in groups if x not in ("Mbovis-All",)]
        if any("Afr2CC" in x or "A2CC" in x for x in g):
            return "Af2"
        if any("caprae" in x.lower() for x in g):
            return "M. caprae (La2)"
        return "non-Af2" if g else "unassigned"

    d = pd.DataFrame(
        [{"sample": s,
          "vsnp_clonal_complex": pick(g),
          "vsnp_all_groups": ";".join(sorted(g))}
         for s, g in rows.items()]
    )
    frames["groups"] = d
    print(f"  vsnp3 groups: {len(d)} samples")
else:
    warn("vSNP3 step2 group FASTAs not found")

# --------------------------------------------------------------------------
# 5. RD typing
# --------------------------------------------------------------------------
p = f"{RDIR}/rd_typing_h37rv.tsv"
if os.path.exists(p):
    d = pd.read_csv(p, sep="\t")
    keep = {
        "sample": "sample",
        "RD9_pctdel": "RD9_pct_deleted",
        "BOVonly_pctdel": "RD4bov_pct_deleted",
        "CAPonly_pctdel": "La2del_pct_deleted",
        "RD1_pctdel": "RD1_pct_deleted",
        "call": "rd_species_call",
    }
    d = d[[c for c in keep if c in d.columns]].rename(columns=keep)
    d["sample"] = d["sample"].astype(str)
    frames["rd"] = d
    print(f"  RD typing: {len(d)} samples")
else:
    warn("rd_typing_h37rv.tsv not found")

# --------------------------------------------------------------------------
# 6. SNP distance -> single-linkage clusters
# --------------------------------------------------------------------------
p = first(f"{BACT}/bactopia-runs/snippy-*/core-snp.masked.distance.tsv")
if p:
    m = pd.read_csv(p, sep="\t", index_col=0)
    m.index = [str(i) for i in m.index]
    m.columns = [str(c) for c in m.columns]
    samples = [s for s in m.index if s.lower() not in ("reference", "root")]

    parent = {s: s for s in samples}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i, a in enumerate(samples):
        for b in samples[i + 1:]:
            try:
                if float(m.loc[a, b]) <= SNP_THRESHOLD:
                    union(a, b)
            except (KeyError, ValueError):
                pass

    comp = {}
    for s in samples:
        comp.setdefault(find(s), []).append(s)

    # name clusters by size, singletons left blank
    label, n = {}, 0
    for members in sorted(comp.values(), key=len, reverse=True):
        if len(members) > 1:
            n += 1
            for s in members:
                label[s] = f"C{n}"

    d = pd.DataFrame(
        [{"sample": s,
          f"snp_cluster_le{SNP_THRESHOLD}": label.get(s, "singleton"),
          "cluster_size": len(comp[find(s)]) if len(comp[find(s)]) > 1 else 1}
         for s in samples]
    )
    frames["clusters"] = d
    print(f"  SNP clusters (<= {SNP_THRESHOLD} SNPs): "
          f"{n} clusters among {len(samples)} samples")
else:
    warn("core-snp.masked.distance.tsv not found")

# --------------------------------------------------------------------------
# merge
# --------------------------------------------------------------------------
if not frames:
    sys.exit("No input files found -- check the paths at the top of this script.")

order = ["spoligo", "tbprofiler", "rd", "groups", "clusters", "assembly"]
master = None
for k in order:
    if k not in frames:
        continue
    d = frames[k]
    d["sample"] = d["sample"].astype(str)
    master = d if master is None else master.merge(d, on="sample", how="outer")

master = master.sort_values("sample").reset_index(drop=True)

# put the interpretive columns first
front = [c for c in [
    "sample", "rd_species_call", "tbp_lineage", "tbp_sublineage",
    "vsnp_clonal_complex", "spoligo_SB", "spoligo_octal",
    f"snp_cluster_le{SNP_THRESHOLD}", "cluster_size",
] if c in master.columns]
master = master[front + [c for c in master.columns if c not in front]]

csv_out = f"{OUTDIR}/master_sample_table.csv"
master.to_csv(csv_out, index=False)
print(f"\nWrote {csv_out}  ({len(master)} rows x {len(master.columns)} columns)")

try:
    xlsx_out = f"{OUTDIR}/master_sample_table.xlsx"
    with pd.ExcelWriter(xlsx_out, engine="openpyxl") as xl:
        master.to_excel(xl, sheet_name="samples", index=False)
    print(f"Wrote {xlsx_out}")
except Exception as e:  # openpyxl may be absent
    print(f"(xlsx not written: {e})")

print("\n--- preview ---")
with pd.option_context("display.max_columns", None, "display.width", 250):
    print(master[front].to_string(index=False))