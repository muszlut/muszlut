#!/usr/bin/env python3
"""List every annotated feature overlapping a coordinate window in a GenBank
file, with its locus_tag, old_locus_tag, gene name and product.

Usage: python3 show_gap.py <GBK> <START> <END>
"""
import re
import sys

GBK = sys.argv[1]
START = int(sys.argv[2]) if len(sys.argv) > 2 else 683106
END = int(sys.argv[3]) if len(sys.argv) > 3 else 685652

feat = None
out = []
with open(GBK) as fh:
    for line in fh:
        m = re.match(r"^\s{5}(\S+)\s+(\S.*)$", line)
        if m:
            if feat:
                out.append(feat)
            raw = m.group(2).strip()
            nums = [int(x) for x in re.findall(r"\d+", raw)]
            feat = {
                "type": m.group(1),
                "start": min(nums) if nums else 0,
                "end": max(nums) if nums else 0,
                "strand": "-" if raw.startswith("complement") else "+",
                "q": {},
            }
            continue
        if feat is None:
            continue
        q = re.search(r'/(\w+)="?([^"\n]+)"?', line)
        if q and q.group(1) in ("locus_tag", "old_locus_tag", "gene", "product"):
            feat["q"].setdefault(q.group(1), q.group(2).strip())
if feat:
    out.append(feat)

print(f"features overlapping {START:,}-{END:,}\n")
hits = [f for f in out
        if f["end"] >= START and f["start"] <= END and f["type"] != "source"]
if not hits:
    print("  none -- this interval is unannotated in this file.")
for f in hits:
    q = f["q"]
    print(f"{f['type']:<6} {f['start']:>9,}-{f['end']:<9,} {f['strand']}  "
          f"{q.get('locus_tag', '-'):<18} old={q.get('old_locus_tag', '-'):<12} "
          f"gene={q.get('gene', '-')}")
    if q.get("product"):
        print(f"       product: {q['product']}")