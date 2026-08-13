#!/usr/bin/env python3
"""
Build per-class MAFFT alignments containing both vertebrate CENH3 and H3
sequences, ungapped from the full (untrimmed) 1681-col joint MSA.

Outputs: vertebrate_antibody/ntail/cenh3_h3_{cls}_unaligned.fasta
         vertebrate_antibody/ntail/cenh3_h3_{cls}_aligned.fasta
"""
import re, subprocess, sys
from pathlib import Path
import pandas as pd

BASE  = Path(__file__).parent
ALN   = BASE / "cenpa430_H3_archaea10.aligned.clipkit.325sp.fasta"
GRP_F = BASE / "split_entropy" / "groupsim" / "groups_gap085.txt"
TAX_F = Path("/home/jg2070/Desktop/dtol_review_August/2026_trees/"
             "annotation_centromeres/centromere_code_to_species.tsv")
OUT   = BASE / "vertebrate_antibody" / "ntail"
OUT.mkdir(parents=True, exist_ok=True)

VERT = {"Aves", "Mammalia", "Actinopterygii", "Reptilia"}

def read_fasta(p):
    seqs, cur, buf = {}, None, []
    for line in Path(p).read_text().splitlines():
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            if cur is not None: seqs[cur] = "".join(buf)
            cur = line[1:]; buf = []
        else: buf.append(line.upper())
    if cur is not None: seqs[cur] = "".join(buf)
    return seqs

def seq_to_code(sid):
    return re.sub(r"[0-9.].*$", "", sid.split()[0].split("_")[0].lower())

# ── load ──────────────────────────────────────────────────────────────────────
tax_df   = pd.read_csv(TAX_F, sep="\t")
code_map = {re.sub(r"[0-9.].*$", "", str(r.fasta_base).lower()): (r.taxa1, r.genus, r.species)
            for _, r in tax_df.iterrows()}

grps = {}
for line in GRP_F.read_text().splitlines():
    if ":" not in line: continue
    g, ids = line.strip().split(":", 1)
    grps[g] = ids.split(",")

seqs_all = read_fasta(ALN)   # full 1681-col MSA

# ── collect per-class ─────────────────────────────────────────────────────────
by_class = {"cenh3": {}, "h3": {}}

for sid in grps["CENPA"]:
    if sid not in seqs_all: continue
    info = code_map.get(seq_to_code(sid))
    if not info or info[0] not in VERT: continue
    cls = info[0]
    label = f"CENH3|{info[1]}_{info[2]}|{sid.split()[0]}"
    by_class["cenh3"].setdefault(cls, {})[label] = seqs_all[sid].replace("-", "")

for sid in grps["H3"]:
    if sid not in seqs_all: continue
    info = code_map.get(seq_to_code(sid))
    if not info or info[0] not in VERT: continue
    cls = info[0]
    label = f"H3|{info[1]}_{info[2]}|{sid.split()[0]}"
    by_class["h3"].setdefault(cls, {})[label] = seqs_all[sid].replace("-", "")

# ── write combined unaligned + run MAFFT ─────────────────────────────────────
for cls in sorted(set(by_class["cenh3"]) | set(by_class["h3"])):
    combined = {**by_class["cenh3"].get(cls, {}), **by_class["h3"].get(cls, {})}
    n_c = len(by_class["cenh3"].get(cls, {}))
    n_h = len(by_class["h3"].get(cls, {}))
    print(f"\n{cls}: {n_c} CENH3 + {n_h} H3 sequences")

    unaligned = OUT / f"cenh3_h3_{cls}_unaligned.fasta"
    with open(unaligned, "w") as fh:
        for label, seq in combined.items():
            fh.write(f">{label}\n{seq}\n")

    aligned = OUT / f"cenh3_h3_{cls}_aligned.fasta"
    print(f"  Running MAFFT...")
    result = subprocess.run(
        ["mafft", "--auto", "--quiet", str(unaligned)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"  MAFFT failed: {result.stderr[:200]}")
        continue
    aligned.write_text(result.stdout)
    n_cols = len(result.stdout.splitlines()[1]) if result.stdout else 0
    print(f"  → {aligned.name}  ({n_cols} cols)")

print("\nDone.")
