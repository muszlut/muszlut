#!/usr/bin/env python3
"""
Confirm that the deletion detected in the Af2 isolates is RDAf2 as defined by
Berg et al. 2011 (J Bacteriol 193:670).

Berg et al. define RDAf2 as a 14,094 bp deletion that removes the complete
coding sequences of Mb0600c-Mb0609 and parts of Mb0599 and Mb0610.

This script reads the M. bovis AF2122/97 GenBank annotation (NC_002945.4),
pulls out the coordinates of BQ2027_MB0599 ... BQ2027_MB0610, and tests the
observed deletion interval against that definition.

No dependencies beyond the standard library.

Usage:
    python3 check_rdaf2.py [GBK] [START] [END]

Defaults to the vSNP3 AF2122/97 dependency file and the interval found in
this dataset.
"""

import re
import sys

GBK = sys.argv[1] if len(sys.argv) > 1 else (
    "/home/ma95362/vsnp3_reference_options/Mycobacterium_AF2122/NC_002945v4.gbk")
START = int(sys.argv[2]) if len(sys.argv) > 2 else 681580
END = int(sys.argv[3]) if len(sys.argv) > 3 else 695673

# MB0599 .. MB0610. Complement-strand genes carry a 'C' suffix in the 2017
# re-annotation (BQ2027_MB0600C = Berg's Mb0600c), so both forms are tried.
def resolve(n, genes):
    for tag in (f"BQ2027_MB{n:04d}", f"BQ2027_MB{n:04d}C"):
        if tag in genes:
            return tag
    return f"BQ2027_MB{n:04d}"

# ---------------------------------------------------------------- parse GBK
# Only CDS features are needed; a GenBank CDS block looks like
#      CDS             complement(681172..681975)
#                      /locus_tag="BQ2027_MB0599"
genes = {}
loc_re = re.compile(r"^\s{5}(CDS|gene)\s+(\S.*)$")
tag_re = re.compile(r'/locus_tag="([^"]+)"')
old_re = re.compile(r'/old_locus_tag="([^"]+)"')
gene_re = re.compile(r'/gene="([^"]+)"')

cur = None
with open(GBK) as fh:
    for line in fh:
        m = loc_re.match(line)
        if m:
            raw = m.group(2).strip()
            strand = "-" if raw.startswith("complement") else "+"
            nums = [int(x) for x in re.findall(r"\d+", raw)]
            cur = {"start": min(nums), "end": max(nums), "strand": strand,
                   "feature": m.group(1), "name": ""} if nums else None
            continue
        if cur is None:
            continue
        t = tag_re.search(line)
        if t:
            tag = t.group(1)
            # keep the CDS record in preference to the gene record
            if tag not in genes or cur["feature"] == "CDS":
                genes[tag] = dict(cur, tag=tag)
            cur = dict(cur, tag=tag)
            continue
        g = gene_re.search(line) or old_re.search(line)
        if g and cur.get("tag") in genes:
            genes[cur["tag"]]["name"] = g.group(1)

if not genes:
    sys.exit(f"No features parsed from {GBK} -- check the path.")

WANT = [resolve(n, genes) for n in range(599, 611)]

# ---------------------------------------------------------------- report
size = END - START + 1
print(f"annotation : {GBK}")
print(f"deletion   : NC_002945.4:{START:,}-{END:,}  ({size:,} bp)")
print(f"Berg et al.: 14,094 bp   -> {'MATCH' if size == 14094 else 'MISMATCH'}\n")

print(f"{'locus_tag':<16} {'gene':<10} {'start':>9} {'end':>9} {'str':>4} "
      f"{'bp in deletion':>15} {'status':>10}")
print("-" * 80)

full, partial, outside, missing = [], [], [], []
for tag in WANT:
    g = genes.get(tag)
    if not g:
        missing.append(tag)
        print(f"{tag:<16} {'?':<10} {'-':>9} {'-':>9} {'-':>4} "
              f"{'-':>15} {'NOT FOUND':>10}")
        continue
    ov = max(0, min(g["end"], END) - max(g["start"], START) + 1)
    glen = g["end"] - g["start"] + 1
    if ov == 0:
        status, bucket = "outside", outside
    elif ov == glen:
        status, bucket = "COMPLETE", full
    else:
        status, bucket = "partial", partial
    bucket.append(tag)
    print(f"{tag:<16} {g['name'][:10]:<10} {g['start']:>9,} {g['end']:>9,} "
          f"{g['strand']:>4} {ov:>10,}/{glen:<4,} {status:>10}")

# ---------------------------------------------------------------- verdict
print("\n" + "=" * 80)
exp_full = [resolve(n, genes) for n in range(600, 610)]   # MB0600c..MB0609
exp_part = [resolve(599, genes), resolve(610, genes)]

checks = [
    ("deletion size is 14,094 bp", size == 14094),
    ("10 genes MB0600-MB0609 completely deleted",
     sorted(full) == sorted(exp_full)),
    ("MB0599 and MB0610 partially deleted",
     sorted(partial) == sorted(exp_part)),
]
for label, ok in checks:
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}")

if all(ok for _, ok in checks):
    print("\nThe deletion matches RDAf2 as defined by Berg et al. 2011.")
    print("The name RDAf2 can be used in the manuscript.")
else:
    print("\nThe deletion does NOT fully match the published RDAf2 definition.")
    print("Report the coordinates and gene span, and do not use the name RDAf2")
    print("without explaining the difference.")
    if missing:
        print(f"\n(locus tags not found: {', '.join(missing)} -- this "
              f"annotation may use different tags; grep the GBK for 'MB059'.)")