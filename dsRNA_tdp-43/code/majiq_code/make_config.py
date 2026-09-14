#!/usr/bin/env python3
"""
make_config.py -- write a MAJIQ builder config (.ini) from liu_sample_groups.csv.

The [experiments] section keys are build groups (TDPneg, TDPpos); the values are
the BAM basenames (Run IDs, no .bam) that belong to each group, comma-separated.
MAJIQ finds <basename>.bam and <basename>.bam.bai in one of the bamdirs.

Ref: MAJIQ builder config = [info] (bamdirs, genome, strandness) + [experiments].
"""
import argparse, csv, sys
from pathlib import Path

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-map", required=True, type=Path,
                    help="CSV with columns Run,patient,condition (TP/TN)")
    ap.add_argument("--bamdir", required=True,
                    help="absolute path to the dir with <Run>.bam + <Run>.bam.bai")
    ap.add_argument("--genome", default="hg38")
    ap.add_argument("--strandness", default="None",
                    choices=["None", "forward", "reverse"])
    ap.add_argument("--grp1-name", default="TDPneg")
    ap.add_argument("--grp2-name", default="TDPpos")
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    # condition code (TP/TN) -> build-group name
    cond2grp = {"TN": a.grp1_name, "TP": a.grp2_name}
    groups = {a.grp1_name: [], a.grp2_name: []}
    with open(a.sample_map) as fh:
        for row in csv.DictReader(fh):
            cond = row["condition"].strip()
            if cond not in cond2grp:            # skip Unsorted / anything else
                continue
            groups[cond2grp[cond]].append(row["Run"].strip())

    for g, members in groups.items():
        if not members:
            sys.exit(f"ERROR: build group {g!r} has no samples")

    lines = ["[info]",
             f"bamdirs={a.bamdir}",
             f"genome={a.genome}",
             f"strandness={a.strandness}",
             "",
             "[experiments]"]
    for g, members in groups.items():
        lines.append(f"{g}=" + ",".join(members))   # no spaces, comma-separated

    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text("\n".join(lines) + "\n")
    print(f"wrote {a.out}")
    for g, members in groups.items():
        print(f"  {g}: n={len(members)}  ({', '.join(members)})")
    return 0

if __name__ == "__main__":
    sys.exit(main())
