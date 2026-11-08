#!/usr/bin/env python3
"""
Generate a full set of iTOL annotation files for the eastern Ethiopia
M. bovis / M. caprae tree.

No dependencies beyond the standard library -- all assignments are recorded
here, taken from the analyses already run:

  species / RD typing ....... RD9, RD4 and La2-deletion coverage vs H37Rv
  clonal complex ............ vSNP3 3.31 defining-SNP groups
  sublineage ................ TB-Profiler 6.7.0 (Zwyer La barcodes)
  spoligotype ............... vSNP3 3.31 in silico, SB numbers
  SNP clusters .............. snp-dists, single linkage at <=12 SNPs

Produces, in OUTDIR:
  itol_species.txt            colour strip   M. bovis vs M. caprae
  itol_clonal_complex.txt     colour strip   Af2 / non-Af2 / M. caprae
  itol_sublineage.txt         colour strip   La1.2 / La1.3 / mixed / La2
  itol_spoligotype.txt        colour strip   SB number
  itol_spoligotype_text.txt   text labels    SB number beside each tip
  itol_labels.txt             tip relabels   "sample  SB number"
  itol_multibar_depth.txt     simple bar     mean depth per isolate

Drag any combination onto the tree in iTOL; they stack as separate rings
or strips.

Usage:  python3 make_itol_annotations.py [OUTDIR]
"""

import os
import sys

OUTDIR = sys.argv[1] if len(sys.argv) > 1 else \
    "/scratch/ma95362/Bio_project_publication/Mbovis/R_output/figures"
os.makedirs(OUTDIR, exist_ok=True)

# ---------------------------------------------------------------- assignments

AF2 = """17 151P 213P 262G 264G 264P 271G 281G 281P 283G 292G 300P 301G
         304G 304P 305G 308P 310P 342P 343P""".split()
NONAF2 = """120PG 231G 258G 259G 262P 288P 307G 318G 322G 323G 32_HK
            341P""".split()
CAPRAE = "097G_400 100GP 106GP 97G".split()
UNASSIGNED = "147 234GP".split()

SUBLINEAGE = {}
for s in AF2:
    SUBLINEAGE[s] = "La1.3"
for s in "120PG 231G 258G 262P 288P 307G 318G 322G 323G 32_HK 341P 234GP".split():
    SUBLINEAGE[s] = "La1.2"
for s in "147 259G".split():
    SUBLINEAGE[s] = "La1.2;La1.3"
for s in CAPRAE:
    SUBLINEAGE[s] = "La2"

SPOLIGO = {}
for grp, members in {
    "SB0134": "120PG 231G 258G 259G 262P 288P 307G 318G 322G 323G 32_HK 341P",
    "SB1176": "151P 17 262G 264G 264P 271G 292G 300P 301G 304G 304P 308P 310P 342P 343P",
    "SB0133": "281G 281P 283G",
    "SB0418": "097G_400 100GP 106GP",
    "SB1690": "97G",
    "unassigned": "147 213P 234GP 305G",
}.items():
    for s in members.split():
        SPOLIGO[s] = grp

# mean depth vs AF2122/97, from the vSNP3 step 1 stats
DEPTH = {
    "097G_400": 115.8, "100GP": 222.7, "106GP": 181.2, "120PG": 212.6,
    "147": 182.7, "151P": 197.8, "17": 217.7, "213P": 61.2, "231G": 83.8,
    "234GP": 34.6, "258G": 200.5, "259G": 99.3, "262G": 159.8, "262P": 186.7,
    "264G": 186.4, "264P": 148.6, "271G": 247.6, "281G": 158.6, "281P": 220.5,
    "283G": 92.8, "288P": 157.0, "292G": 251.3, "300P": 155.1, "301G": 169.7,
    "304G": 142.8, "304P": 230.7, "305G": 79.6, "307G": 166.8, "308P": 74.8,
    "310P": 157.7, "318G": 109.2, "322G": 132.0, "323G": 138.6,
    "32_HK": 148.7, "341P": 99.4, "342P": 91.9, "343P": 107.1, "97G": 65.4,
}

ALL = AF2 + NONAF2 + CAPRAE + UNASSIGNED
assert len(ALL) == 38, len(ALL)
assert len(set(ALL)) == 38

SPECIES = {s: ("M. caprae" if s in CAPRAE else "M. bovis") for s in ALL}
COMPLEX = {}
for s in ALL:
    COMPLEX[s] = ("Af2" if s in AF2 else
                  "non-Af2" if s in NONAF2 else
                  "M. caprae" if s in CAPRAE else "unassigned")

# ---------------------------------------------------------------- palettes

PAL_SPECIES = {"M. bovis": "#2a6f52", "M. caprae": "#6b4c9a"}
PAL_COMPLEX = {"Af2": "#1f6f8b", "non-Af2": "#c26a3d",
               "M. caprae": "#6b4c9a", "unassigned": "#999999"}
PAL_SUBLIN  = {"La1.2": "#c26a3d", "La1.3": "#1f6f8b",
               "La1.2;La1.3": "#d9b44a", "La2": "#6b4c9a"}
PAL_SPOLIGO = {"SB0134": "#c26a3d", "SB1176": "#1f6f8b", "SB0133": "#7fb3c8",
               "SB0418": "#6b4c9a", "SB1690": "#a88cc4", "unassigned": "#999999"}


def colorstrip(path, label, mapping, palette, width=25):
    keys = list(palette)
    lines = [
        "DATASET_COLORSTRIP", "SEPARATOR TAB",
        f"DATASET_LABEL\t{label}", f"COLOR\t{palette[keys[0]]}",
        f"STRIP_WIDTH\t{width}", "MARGIN\t3", "BORDER_WIDTH\t0.5",
        "SHOW_STRIP_LABELS\t0",
        f"LEGEND_TITLE\t{label}",
        "LEGEND_SHAPES\t" + "\t".join("1" for _ in keys),
        "LEGEND_COLORS\t" + "\t".join(palette[k] for k in keys),
        "LEGEND_LABELS\t" + "\t".join(keys),
        "DATA",
    ]
    for s in ALL:
        v = mapping.get(s)
        if v in palette:
            lines.append(f"{s}\t{palette[v]}\t{v}")
    open(path, "w").write("\n".join(lines) + "\n")
    print("wrote", path)


def textlabels(path, label, mapping, colour="#333333"):
    lines = [
        "DATASET_TEXT", "SEPARATOR TAB", f"DATASET_LABEL\t{label}",
        f"COLOR\t{colour}", "MARGIN\t8", "SHOW_INTERNAL\t0",
        "LABEL_ROTATION\t0", "ALIGN_TO_TREE\t1", "SIZE_FACTOR\t0.9",
        "DATA",
    ]
    for s in ALL:
        v = mapping.get(s)
        if v:
            #  node  label  position  colour  style  size_factor  rotation
            lines.append(f"{s}\t{v}\t-1\t{colour}\tnormal\t1\t0")
    open(path, "w").write("\n".join(lines) + "\n")
    print("wrote", path)


def relabels(path):
    lines = ["LABELS", "SEPARATOR TAB", "DATA"]
    for s in ALL:
        sb = SPOLIGO.get(s, "")
        sb = "" if sb == "unassigned" else sb
        lines.append(f"{s}\t{s}  {sb}".rstrip())
    open(path, "w").write("\n".join(lines) + "\n")
    print("wrote", path)


def simplebar(path, label, mapping, colour="#4b7f9b"):
    lines = [
        "DATASET_SIMPLEBAR", "SEPARATOR TAB", f"DATASET_LABEL\t{label}",
        f"COLOR\t{colour}", "WIDTH\t120", "MARGIN\t10",
        "DATASET_SCALE\t50\t100\t150\t200\t250", "DATA",
    ]
    for s in ALL:
        if s in mapping:
            lines.append(f"{s}\t{mapping[s]}")
    open(path, "w").write("\n".join(lines) + "\n")
    print("wrote", path)


# ---------------------------------------------------------------- write

colorstrip(os.path.join(OUTDIR, "itol_species.txt"),
           "Species", SPECIES, PAL_SPECIES)
colorstrip(os.path.join(OUTDIR, "itol_clonal_complex.txt"),
           "Clonal complex", COMPLEX, PAL_COMPLEX)
colorstrip(os.path.join(OUTDIR, "itol_sublineage.txt"),
           "Sublineage", SUBLINEAGE, PAL_SUBLIN)
colorstrip(os.path.join(OUTDIR, "itol_spoligotype.txt"),
           "Spoligotype", SPOLIGO, PAL_SPOLIGO)
textlabels(os.path.join(OUTDIR, "itol_spoligotype_text.txt"),
           "SB number", {k: ("" if v == "unassigned" else v)
                         for k, v in SPOLIGO.items()})
relabels(os.path.join(OUTDIR, "itol_labels.txt"))
simplebar(os.path.join(OUTDIR, "itol_depth.txt"),
          "Mean depth (AF2122/97)", DEPTH)

print(f"\n{len(ALL)} isolates annotated.")
print("Upload core-snp.treefile to itol.embl.de, then drag these files onto it.")
print("Note the Reference tip is intentionally unannotated.")