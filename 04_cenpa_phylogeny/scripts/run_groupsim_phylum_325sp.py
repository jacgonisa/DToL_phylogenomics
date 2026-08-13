#!/usr/bin/env python3
"""
Phylum-based GroupSim on the bnni alignment.

Groups all sequences (CENPA + H3, no archaea) by their broad taxonomic phylum:
  Viridiplantae, Fungi, Insects, Vertebrates, Invertebrates, Protists

CENPA sequences → phylum from iTOL annotation (already correct).
H3 sequences    → phylum inferred from DToL lowercase prefix.

Run from:  PhylogeneticProfiling/
  python3 19_curated_tree/run_groupsim_phylum.py
"""

import re
import subprocess
import sys
from pathlib import Path
from collections import Counter

sys.path.insert(0, str(Path("/home/jg2070/Desktop/dtol_review_August/PhylogeneticProfiling")))
from dtol_taxonomy import DToLTaxonomy, PHYLUM_MAP

# ── paths ──────────────────────────────────────────────────────────────────────
BASE      = Path(__file__).parent
GROUPSIM  = Path("/home/jg2070/Desktop/groupsim-py3/src/groupsim.py")

ALN          = BASE / "cenpa430_H3_archaea10.aligned.clipkit.325sp.fasta"
H3_ALL_FASTA = BASE / "H3_all.aligned.fasta"
ITOL         = BASE / "iTOL_bnni" / "iTOL_taxonomic_clade.txt"

OUT_DIR = BASE / "split_entropy" / "groupsim_phylum"
OUT_DIR.mkdir(parents=True, exist_ok=True)

ARCHAEA_IDS = {
    "OLS22332.1","OLS24873.1","OLS21974.1","KKK41979.1","KXH71038.1",
    "OLS18261.1","OLS16336.1","BAD86478.1","OIO61677.1","OIO41945.1",
}

GAP = {"-"}
# PHYLUM_MAP and DToLTaxonomy imported from dtol_taxonomy.py


# ── FASTA helpers ──────────────────────────────────────────────────────────────
def read_fasta(path: Path):
    seqs, cur, buf = {}, None, []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur is not None:
                seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


def trim_columns(seqs, gap_threshold=0.85):
    ids  = list(seqs.keys())
    arr  = [seqs[i] for i in ids]
    n, L = len(arr), len(arr[0])
    keep = [i for i in range(L)
            if sum(1 for s in arr if s[i] in GAP) / n <= gap_threshold]
    trimmed = {sid: "".join(arr[k][c] for c in keep)
               for k, sid in enumerate(ids)}
    return trimmed, keep


def write_fasta(seqs, path: Path):
    with open(path, "w") as fh:
        for sid, seq in seqs.items():
            fh.write(f">{sid}\n{seq}\n")


# ── assign phylum to every sequence ──────────────────────────────────────────
def assign_phyla(seqs, h3_pool, itol, tax: DToLTaxonomy):
    """
    Returns dict: sid → broad phylum label (or None if unresolvable).
    CENPA sequences: use iTOL label → PHYLUM_MAP (already phylum-annotated).
    H3 sequences:   use DToLTaxonomy (prefix → fine → broad).
    """
    assignments, unresolved = {}, []

    for sid in seqs:
        if sid in h3_pool:
            broad = tax.broad_group(sid)
        else:
            itol_label = itol.get(sid, "Unknown")
            broad = PHYLUM_MAP.get(itol_label) or tax.broad_group(sid)

        if broad:
            assignments[sid] = broad
        else:
            unresolved.append(sid)

    return assignments, unresolved


# ── main ───────────────────────────────────────────────────────────────────────
def run_for_gap(gap, seqs, assignments, groups_by_phylum):
    tag = str(gap).replace(".", "")
    print(f"\n── Gap threshold {gap} ──")

    # trim using ALL assigned sequences
    assigned_seqs = {k: seqs[k] for k in assignments}
    trimmed, keep = trim_columns(assigned_seqs, gap_threshold=gap)
    print(f"  columns kept: {len(keep)}")

    # order: groups alphabetically, sequences within each group
    group_order = sorted(groups_by_phylum.keys())
    ordered_ids = [sid for g in group_order
                   for sid in groups_by_phylum[g] if sid in trimmed]
    ordered_seqs = {sid: trimmed[sid] for sid in ordered_ids}

    # write trimmed FASTA
    aln_out = OUT_DIR / f"trimmed_phylum_gap{tag}.fasta"
    write_fasta(ordered_seqs, aln_out)

    # write group file
    grp_file = OUT_DIR / f"groups_phylum_gap{tag}.txt"
    with open(grp_file, "w") as fh:
        for g in group_order:
            ids_in_trim = [s for s in groups_by_phylum[g] if s in trimmed]
            if ids_in_trim:
                fh.write(g + ":" + ",".join(ids_in_trim) + "\n")
    print(f"  groups: " + ", ".join(
        f"{g}({len([s for s in groups_by_phylum[g] if s in trimmed])})"
        for g in group_order))

    # run GroupSim
    out_prefix = str(OUT_DIR / f"groupsim_phylum_gap{tag}")
    cmd = [
        sys.executable, str(GROUPSIM),
        "-k", str(grp_file),
        "-o", out_prefix,
        "-c", str(gap),
        "-g", "0.5",
        "-w", "3",
        "-l", "0.7",
        "-n",
        str(aln_out),
    ]
    print("  running GroupSim …")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("STDERR:", result.stderr[-2000:])
        raise RuntimeError(f"GroupSim failed (gap={gap})")
    print(result.stderr.strip())
    print(f"  done → {out_prefix}.*")


def main():
    print("Reading alignment …")
    seqs_all = read_fasta(ALN)
    seqs     = {k: v for k, v in seqs_all.items() if k not in ARCHAEA_IDS}
    print(f"  seqs (no archaea): {len(seqs)}")

    h3_pool = set(read_fasta(H3_ALL_FASTA).keys())
    tax     = DToLTaxonomy()   # loads 327_DTOL_clean.csv via dtol_taxonomy.py

    # parse iTOL (for CENPA phylum labels)
    itol = {}
    in_data = False
    for line in ITOL.read_text().splitlines():
        if line.strip() == "DATA":
            in_data = True
            continue
        if not in_data or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 3:
            itol[parts[0]] = parts[2]

    print("Assigning phyla …")
    assignments, unresolved = assign_phyla(seqs, h3_pool, itol, tax)
    print(f"  assigned: {len(assignments)}  unresolved: {len(unresolved)}")
    if unresolved:
        print(f"  unresolved sample: {unresolved[:5]}")

    # count CENPA vs H3 per phylum
    groups_by_phylum = {}
    for sid, phylum in assignments.items():
        groups_by_phylum.setdefault(phylum, []).append(sid)

    print("\nPhylum composition (CENPA / H3):")
    for g in sorted(groups_by_phylum):
        n_cenpa = sum(1 for s in groups_by_phylum[g] if s not in h3_pool)
        n_h3    = sum(1 for s in groups_by_phylum[g] if s in h3_pool)
        print(f"  {g}: {len(groups_by_phylum[g])} total  (CENPA={n_cenpa}, H3={n_h3})")

    for gap in [0.80, 0.85, 0.90]:
        run_for_gap(gap, seqs, assignments, groups_by_phylum)

    print("\nAll done. Outputs in:", OUT_DIR)


if __name__ == "__main__":
    main()
