#!/usr/bin/env python3
"""
Extract vertebrate CENP-A sequences from the 325sp MSA.
Outputs per-class aligned FASTAs + metadata TSV.
"""
import re
from pathlib import Path
import pandas as pd

BASE    = Path(__file__).parent
ALN     = BASE / "cenpa430_H3_archaea10.aligned.clipkit.325sp.fasta"
GRP_F   = BASE / "split_entropy" / "groupsim" / "groups_gap085.txt"
TAX_F   = Path("/home/jg2070/Desktop/dtol_review_August/2026_trees/annotation_centromeres/centromere_code_to_species.tsv")
OUT_DIR = BASE / "vertebrate_antibody"
OUT_DIR.mkdir(exist_ok=True)

VERT_CLASSES = {"Aves", "Mammalia", "Actinopterygii", "Reptilia"}

def read_fasta(path):
    seqs, cur, buf = {}, None, []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            if cur is not None: seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]; buf = []
        else: buf.append(line.upper())
    if cur is not None: seqs[cur] = "".join(buf)
    return seqs

def seq_to_code(sid):
    return re.sub(r"[0-9.].*$", "", sid.split("_")[0].lower())

tax_df   = pd.read_csv(TAX_F, sep="\t")
code_map = {re.sub(r"[0-9.].*$", "", r.fasta_base.lower()): (r.taxa1, r.genus, r.species)
            for _, r in tax_df.iterrows()}

cenpa_ids = set()
for line in GRP_F.read_text().splitlines():
    if line.startswith("CENPA:"):
        cenpa_ids = set(line[6:].split(","))

seqs_all = read_fasta(ALN)

meta, vert_by_class = [], {}
for sid, seq in seqs_all.items():
    if sid not in cenpa_ids: continue
    code = seq_to_code(sid)
    if code not in code_map: continue
    taxa1, genus, species = code_map[code]
    if taxa1 not in VERT_CLASSES: continue
    meta.append({"seq_id": sid, "class": taxa1, "genus": genus, "species": species})
    vert_by_class.setdefault(taxa1, {})[sid] = seq

pd.DataFrame(meta).sort_values(["class", "genus", "species"]).to_csv(
    OUT_DIR / "vertebrate_cenpa_metadata.tsv", sep="\t", index=False)
print(f"Metadata: {len(meta)} sequences, {len({m['genus']+m['species'] for m in meta})} species")

NT_DIR = OUT_DIR / "ntail"
NT_DIR.mkdir(exist_ok=True)

all_verts = {}
for cls in sorted(vert_by_class):
    sdict = vert_by_class[cls]
    with open(OUT_DIR / f"cenpa_{cls}_aligned.fasta", "w") as fh:
        for sid, seq in sdict.items():
            fh.write(f">{sid}\n{seq}\n")
    # unaligned — ungap from the FULL MSA (not the trimmed one)
    with open(NT_DIR / f"cenpa_{cls}_unaligned.fasta", "w") as fh:
        for sid, seq in sdict.items():
            info = code_map.get(seq_to_code(sid))
            sp_label = f"{info[1]} {info[2]}" if info else sid
            fh.write(f">{sid} {sp_label}\n{seq.replace('-','')}\n")
    all_verts.update(sdict)
    n_sp = len({seq_to_code(s) for s in sdict})
    print(f"  {cls}: {len(sdict)} seqs, {n_sp} species → aligned + unaligned FASTAs")

with open(OUT_DIR / "cenpa_vertebrates_aligned.fasta", "w") as fh:
    for sid, seq in all_verts.items():
        fh.write(f">{sid}\n{seq}\n")
print(f"\nTotal vertebrate CENPA: {len(all_verts)} seqs → vertebrate_antibody/")
