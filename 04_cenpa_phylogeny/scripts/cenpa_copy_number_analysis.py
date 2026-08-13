#!/usr/bin/env python3
"""
CENPA copy-number analysis for the 325-sp DToL dataset.

For every species in the curated CENPA set:
  1. Parse curated sequences from cenpa_leaflabels_curated.fasta
  2. Resolve locus coordinates from sequence IDs
     → confirm each sequence comes from a DISTINCT genomic locus (no isoforms)
  3. BLAST all curated CENPA sequences against a reference panel
     (Human CENPA, Arabidopsis CENH3, Yeast CSE4, Drosophila CID, C. elegans HCP-3)
     using blastp from the of3_env conda environment
  4. Cross-reference with the 325sp count table (raw / curated / source)
  5. Output a comprehensive TSV table

Output: 04_cenpa_phylogeny/cenpa_copy_number_table.tsv
"""

import re, subprocess, os, sys, tempfile
from pathlib import Path
from collections import defaultdict

BASE   = Path("/home/jg2070/Desktop/dtol_review_August")
PP     = BASE / "PhylogeneticProfiling"
PUB    = BASE / "DToL_phylogenomics_publication_325genomes"
OUT    = PUB / "04_cenpa_phylogeny"

CURATED_FASTA  = PP / "19_curated_tree" / "cenpa_leaflabels_curated.fasta"
CENPA_TIPS_TXT = OUT / "cenpa_tips_from_contree.txt"   # authoritative: tips in the 325sp tree
REF_FASTA      = PP / "17_cenpa_phylogeny_with_rescue" / "reference_cenpa_sequences.fasta"
COUNT_TABLE    = PUB / "325sp_cenpa_count_table.tsv"
BLASTP         = Path("/home/jg2070/miniforge3/envs/of3_env/bin/blastp")
MAKEBLASTDB    = Path("/home/jg2070/miniforge3/envs/of3_env/bin/makeblastdb")

# Reference IDs we care about for the BLAST summary column
KEY_REFS = {
    "Homo_sapiens_HSAP011513_CENPA_CenpA":            "Human_CENPA",
    "Danio_rerio_DRER013080_CENPA_CenpA":             "Zebrafish_CENPA",
    "Arabidopsis_thaliana_ATHA004583_HTR12_CenpA":    "Arabidopsis_CENH3",
    "Oryza_sativa_OSAT034475_CenpA":                  "Rice_CENH3",
    "Saccharomyces_cerevisiae_SCER003239_CSE4_CenpA": "Yeast_CSE4",
    "Schizosaccharomyces_pombe_SPOM002782_cnp1_CenpA":"Spombe_CNP1",
    "Drosophila_melanogaster_DMEL028507_cid-PA_CenpA":"Drosophila_CID",
    "Caenorhabditis_elegans_CELE018214_hcp-3_CenpA":  "CElegans_HCP3",
}

# Species codes that are external references (skip from copy-number table)
REF_PREFIXES = ("Homo_", "Mus_", "Danio_", "Drosophila_", "Caenorhabditis_",
                "Arabidopsis_", "Saccharomyces_", "Schizosaccharomyces_",
                "Oryza_", "Aquilegia_", "Amborella_", "Selaginella_",
                "Physcomitrella_", "Chlamydomonas_", "Volvox_", "K11")


# ── 1. Parse curated CENPA fasta ─────────────────────────────────────────────
def read_fasta(path):
    seqs, cur, buf = {}, None, []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        if line.startswith(">"):
            if cur: seqs[cur] = "".join(buf)
            cur, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if cur: seqs[cur] = "".join(buf)
    return seqs


def is_ref(seq_id):
    return any(seq_id.startswith(p) for p in REF_PREFIXES)


def parse_locus(seq_id):
    """Return (species, scaffold, gene_idx, isoform_idx) from a sequence ID."""
    if "_rescue" in seq_id:
        sp = seq_id.split("_rescue")[0]
        n  = seq_id.split("_rescue")[1]
        return sp, "rescue", f"rescue_{n}", "1"
    for sep in ("_GRCg7b_chr_", "_chr_", "_SUPER_"):
        if sep in seq_id:
            sp   = seq_id.split(sep)[0]
            rest = seq_id.split(sep, 1)[1]
            parts = rest.rsplit(".", 1)
            gene_str = parts[0]
            iso      = parts[1] if len(parts) > 1 else "1"
            tokens   = gene_str.split("_")
            scaffold = "_".join(tokens[:-1])
            gene_idx = tokens[-1]
            return sp, scaffold, gene_idx, iso
    if "_ENA|" in seq_id:
        sp   = seq_id.split("_ENA|")[0]
        rest = seq_id.split("_ENA|", 1)[1]
        accession = rest.split("|")[1].split("_")[0]
        after_acc  = rest.split("_", 1)[1] if "_" in rest else rest
        parts = after_acc.rsplit(".", 1)
        gene_idx = parts[0]
        iso      = parts[1] if len(parts) > 1 else "1"
        return sp, accession, gene_idx, iso
    # fallback
    return seq_id.split("_")[0], "unknown", seq_id, "1"


print("Parsing CENPA sequences from contree tip list …")
cenpa_tips = set(CENPA_TIPS_TXT.read_text().split())
all_seqs   = read_fasta(CURATED_FASTA)
# sequence source: contree tips first, then alignment fasta for any missing
ALN_FASTA  = OUT / "cenpa430_H3_archaea10.aligned.clipkit.325sp.fasta"
if ALN_FASTA.exists():
    all_seqs.update({k: v for k, v in read_fasta(ALN_FASTA).items() if k not in all_seqs})

dtol_seqs = {k: v for k, v in all_seqs.items() if k in cenpa_tips}
print(f"  CENPA tips in tree: {len(cenpa_tips)}")
print(f"  Sequences resolved: {len(dtol_seqs)}")

# Map species → list of (seq_id, scaffold, gene_idx, isoform)
sp_loci = defaultdict(list)
for sid in dtol_seqs:
    sp, scaf, gene, iso = parse_locus(sid)
    sp_loci[sp].append({"seq_id": sid, "scaffold": scaf,
                         "gene_idx": gene, "isoform": iso})

# Check for isoforms (same sp + scaffold + gene_idx)
locus_keys = defaultdict(list)
for sp, entries in sp_loci.items():
    for e in entries:
        locus_keys[(sp, e["scaffold"], e["gene_idx"])].append(e["seq_id"])

isoform_pairs = {k: v for k, v in locus_keys.items() if len(v) > 1}
print(f"  Loci with multiple isoforms: {len(isoform_pairs)}")
if isoform_pairs:
    for k, v in list(isoform_pairs.items())[:5]:
        print(f"    {k}: {v}")


# ── 2. BLAST against references ───────────────────────────────────────────────
print("\nRunning BLAST against reference CENPA panel …")

BLAST_CACHE = OUT / "cenpa_blast_vs_refs_contree.tsv"  # rebuilt from tree tips
PIDENT_THRESHOLD = 40.0   # minimum % identity to count as confirmed CenpA

blast_hits = {}   # seq_id → best hit row (all results, no threshold yet)

if BLAST_CACHE.exists():
    print(f"  Loading cached BLAST results from {BLAST_CACHE.name} …")
    for line in BLAST_CACHE.read_text().splitlines()[1:]:   # skip header
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        qid = parts[0]
        if qid not in blast_hits:
            blast_hits[qid] = {
                "best_ref":  parts[1],
                "pct_id":    float(parts[2]),
                "aln_len":   int(parts[3]),
                "query_len": int(parts[4]),
                "evalue":    parts[6],
                "bitscore":  float(parts[7]),
            }
else:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        query_fa  = tmpdir / "query.fasta"
        db_path   = tmpdir / "refdb"
        blast_out = tmpdir / "blast.tsv"

        with open(query_fa, "w") as fh:
            for sid, seq in dtol_seqs.items():
                fh.write(f">{sid}\n{seq}\n")

        subprocess.run(
            [str(MAKEBLASTDB), "-in", str(REF_FASTA),
             "-dbtype", "prot", "-out", str(db_path), "-title", "cenpa_refs"],
            check=True, capture_output=True
        )
        subprocess.run(
            [str(BLASTP),
             "-query",  str(query_fa),
             "-db",     str(db_path),
             "-out",    str(blast_out),
             "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore",
             "-evalue", "1e-5",
             "-max_target_seqs", "3",
             "-num_threads", "4"],
            check=True, capture_output=True
        )

        raw_lines = []
        if blast_out.exists():
            for line in blast_out.read_text().splitlines():
                parts = line.split("\t")
                if len(parts) < 8:
                    continue
                raw_lines.append(line)
                qid = parts[0]
                if qid not in blast_hits:
                    blast_hits[qid] = {
                        "best_ref":  parts[1],
                        "pct_id":    float(parts[2]),
                        "aln_len":   int(parts[3]),
                        "query_len": int(parts[4]),
                        "evalue":    parts[6],
                        "bitscore":  float(parts[7]),
                    }

        # Save cache
        with open(BLAST_CACHE, "w") as fh:
            fh.write("qseqid\tsseqid\tpident\tlength\tqlen\tslen\tevalue\tbitscore\n")
            fh.write("\n".join(raw_lines) + "\n")
        print(f"  Saved BLAST cache: {BLAST_CACHE.name}")

print(f"  BLAST hits: {len(blast_hits)} / {len(dtol_seqs)} queries")

# Summarise best BLAST per key reference
def best_hit_for_ref(hits_for_sp, ref_id):
    """Best pct_id among all hits to a given reference for one species."""
    best = None
    for h in hits_for_sp:
        if h.get("best_ref") == ref_id:
            if best is None or h["pct_id"] > best:
                best = h["pct_id"]
    return best


# ── 3. Load 325sp count table ─────────────────────────────────────────────────
print("\nLoading 325sp count table …")
count_table = {}
if COUNT_TABLE.exists():
    import csv
    with open(COUNT_TABLE) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            count_table[row["species_code_full"]] = row
print(f"  Loaded {len(count_table)} species")


# ── 4. Build comprehensive table ──────────────────────────────────────────────
# Normalise species code for lookup (strip trailing digit)
def norm_sp(sp):
    m = re.match(r'^(.+?)(\d+)$', sp)
    return m.group(1) if m else sp

# Find species in count table matching a parsed sp code
def find_in_count_table(sp_parsed):
    """Try various normalised forms to match species in count_table."""
    for key in count_table:
        base = key.split(".")[0]
        if base == sp_parsed:
            return count_table[key]
        norm = norm_sp(base)
        if norm == sp_parsed or norm == norm_sp(sp_parsed):
            return count_table[key]
    return None

print("\nBuilding comprehensive copy-number table …")
rows = []

# Get all species from curated set (DToL only)
all_sp_parsed = sorted(sp_loci.keys())

for sp in all_sp_parsed:
    entries = sp_loci[sp]
    n_seqs  = len(entries)

    # Locus uniqueness (all sequences)
    unique_loci = set((e["scaffold"], e["gene_idx"]) for e in entries)
    n_loci      = len(unique_loci)
    n_isoforms  = n_seqs - n_loci   # should be 0

    # Filtered: only sequences with BLAST pct_id >= threshold
    entries_pass = [e for e in entries
                    if blast_hits.get(e["seq_id"], {}).get("pct_id", 0) >= PIDENT_THRESHOLD]
    n_seqs_pass   = len(entries_pass)
    unique_loci_f = set((e["scaffold"], e["gene_idx"]) for e in entries_pass)
    n_loci_pass   = len(unique_loci_f)

    # Scaffold summary
    scaffolds     = sorted(set(e["scaffold"] for e in entries))
    loci_per_scaf = defaultdict(set)
    for e in entries:
        loci_per_scaf[e["scaffold"]].add(e["gene_idx"])
    scaf_summary = "; ".join(
        f"{s}({len(g)} gene{'s' if len(g)>1 else ''})"
        for s, g in sorted(loci_per_scaf.items())
    )

    # BLAST summary for this species
    sp_hits = [blast_hits[e["seq_id"]]
               for e in entries if e["seq_id"] in blast_hits]
    n_with_hit   = len(sp_hits)
    best_overall_pct = max((h["pct_id"] for h in sp_hits), default=None)
    best_overall_ref = None
    if sp_hits:
        best_h = max(sp_hits, key=lambda h: h["pct_id"])
        best_overall_ref = KEY_REFS.get(best_h["best_ref"], best_h["best_ref"])

    # Per-key-ref best hits
    ref_hits = {}
    for ref_id, ref_label in KEY_REFS.items():
        best = None
        for h in sp_hits:
            if h["best_ref"] == ref_id:
                if best is None or h["pct_id"] > best:
                    best = h["pct_id"]
        ref_hits[ref_label] = f"{best:.1f}" if best is not None else "—"

    # Match count table — also try stripping haplotype/chunk suffixes
    ct = find_in_count_table(sp)
    if ct is None:
        sp_clean = re.sub(r'\.(hap\d+(\.\d+)?|chunk\d+)$', '', sp)
        ct = find_in_count_table(sp_clean)
    sp_code_full   = ct["species_code_full"]   if ct else sp
    cenpa_presence = ct["cenpa_presence"]       if ct else "?"
    cenpa_src      = ct["cenpa_source"]         if ct else "?"
    cenpa_raw      = ct["cenpa_count_raw"]      if ct else "?"
    cenpa_curated  = ct["cenpa_count_curated"]  if ct else "?"
    ctype          = ct["centromere_type"]       if ct else "?"

    # Duplication call — raw (all sequences)
    if n_loci == 1:
        dup_call = "single_copy"
    elif n_loci <= 3:
        dup_call = "low_copy"
    else:
        dup_call = "expanded"

    # Duplication call — filtered (>= 40% identity only)
    if n_loci_pass == 0:
        dup_call_filt = "no_confirmed_copy"
    elif n_loci_pass == 1:
        dup_call_filt = "single_copy"
    elif n_loci_pass <= 3:
        dup_call_filt = "low_copy"
    else:
        dup_call_filt = "expanded"

    # Confidence flag
    if n_with_hit == 0:
        confidence = "NO_BLAST_HIT"
    elif best_overall_pct is not None and best_overall_pct < PIDENT_THRESHOLD:
        confidence = f"BELOW_{int(PIDENT_THRESHOLD)}PCT"
    else:
        confidence = "PASS"

    # Known anomaly notes
    notes = ""
    if sp == "gyLeoLubr":
        notes = "CAUTION: 88 copies in a fungus — likely expanded H3-variant family, not true CenpA paralogs; cenpa_presence=Absent in 275-leaf tree"
    elif n_loci >= 10:
        notes = f"Highly expanded: {n_loci} distinct loci"
    elif cenpa_presence == "Absent" and n_loci > 0:
        notes = "Sequences present in curated set but species called Absent in phylogenetic profiling tree"

    rows.append({
        "species_code_full":    sp_code_full,
        "species_parsed":       sp,
        "centromere_type":      ctype,
        "cenpa_presence":       cenpa_presence,
        "cenpa_source":         cenpa_src,
        "cenpa_count_raw":      cenpa_raw,
        "cenpa_count_curated":  cenpa_curated,
        "n_curated_seqs":           n_seqs,
        "n_unique_loci":            n_loci,
        "n_isoforms_collapsed":     n_isoforms,
        "copy_number_call":         dup_call,
        "n_seqs_pass40pct":         n_seqs_pass,
        "n_loci_pass40pct":         n_loci_pass,
        "copy_number_call_40pct":   dup_call_filt,
        "n_scaffolds":          len(scaffolds),
        "loci_per_scaffold":    scaf_summary,
        "n_seqs_with_blast_hit":n_with_hit,
        "best_blast_pct_id":    f"{best_overall_pct:.1f}" if best_overall_pct is not None else "—",
        "best_blast_ref":       best_overall_ref or "—",
        "pct_id_Human_CENPA":   ref_hits["Human_CENPA"],
        "pct_id_Zebrafish_CENPA": ref_hits["Zebrafish_CENPA"],
        "pct_id_Arabidopsis_CENH3": ref_hits["Arabidopsis_CENH3"],
        "pct_id_Rice_CENH3":    ref_hits["Rice_CENH3"],
        "pct_id_Yeast_CSE4":    ref_hits["Yeast_CSE4"],
        "pct_id_Spombe_CNP1":   ref_hits["Spombe_CNP1"],
        "pct_id_Drosophila_CID":ref_hits["Drosophila_CID"],
        "pct_id_CElegans_HCP3": ref_hits["CElegans_HCP3"],
        "blast_confidence":     confidence,
        "notes":                notes,
        "seq_ids":              "|".join(e["seq_id"] for e in entries),
    })

# ── Add 325sp species with 0 curated CENPA ───────────────────────────────────
species_in_table = set()
for r in rows:
    species_in_table.add(r["species_code_full"])

for code, ct in count_table.items():
    if code not in species_in_table:
        rows.append({
            "species_code_full":    code,
            "species_parsed":       "",
            "centromere_type":      ct["centromere_type"],
            "cenpa_presence":       ct["cenpa_presence"],
            "cenpa_source":         ct["cenpa_source"],
            "cenpa_count_raw":      ct["cenpa_count_raw"],
            "cenpa_count_curated":  ct["cenpa_count_curated"],
            "n_curated_seqs":         0,
            "n_unique_loci":          0,
            "n_isoforms_collapsed":   0,
            "copy_number_call":       "absent_from_curated_tree",
            "n_seqs_pass40pct":       0,
            "n_loci_pass40pct":       0,
            "copy_number_call_40pct": "absent_from_curated_tree",
            "n_scaffolds":          0,
            "loci_per_scaffold":    "",
            "n_seqs_with_blast_hit":0,
            "best_blast_pct_id":    "—",
            "best_blast_ref":       "—",
            "pct_id_Human_CENPA":   "—",
            "pct_id_Zebrafish_CENPA":"—",
            "pct_id_Arabidopsis_CENH3":"—",
            "pct_id_Rice_CENH3":    "—",
            "pct_id_Yeast_CSE4":    "—",
            "pct_id_Spombe_CNP1":   "—",
            "pct_id_Drosophila_CID":"—",
            "pct_id_CElegans_HCP3": "—",
            "blast_confidence":     "NA",
            "notes":                "",
            "seq_ids":              "",
        })

# Sort: expanded first, then by n_loci desc, then alphabetically
rows.sort(key=lambda r: (-r["n_unique_loci"], r["species_code_full"]))

# ── 5. Write output ───────────────────────────────────────────────────────────
out_tsv = OUT / "cenpa_copy_number_table.tsv"
cols = [
    "species_code_full", "species_parsed", "centromere_type", "cenpa_presence",
    "cenpa_source", "cenpa_count_raw", "cenpa_count_curated",
    # raw counts (all sequences in curated tree)
    "n_curated_seqs", "n_unique_loci", "n_isoforms_collapsed", "copy_number_call",
    # filtered counts (sequences with >= 40% identity to any reference CenpA)
    "n_seqs_pass40pct", "n_loci_pass40pct", "copy_number_call_40pct",
    "n_scaffolds", "loci_per_scaffold",
    "n_seqs_with_blast_hit", "best_blast_pct_id", "best_blast_ref",
    "pct_id_Human_CENPA", "pct_id_Zebrafish_CENPA",
    "pct_id_Arabidopsis_CENH3", "pct_id_Rice_CENH3",
    "pct_id_Yeast_CSE4", "pct_id_Spombe_CNP1",
    "pct_id_Drosophila_CID", "pct_id_CElegans_HCP3",
    "blast_confidence", "notes", "seq_ids",
]
with open(out_tsv, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for row in rows:
        fh.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")

print(f"\nSaved: {out_tsv}")
print(f"  Total species in table: {len(rows)}")
print(f"  Single-copy:  {sum(1 for r in rows if r['copy_number_call']=='single_copy')}")
print(f"  Low-copy:     {sum(1 for r in rows if r['copy_number_call']=='low_copy')}")
print(f"  Expanded:     {sum(1 for r in rows if r['copy_number_call']=='expanded')}")
print(f"  NO_BLAST_HIT: {sum(1 for r in rows if r['blast_confidence']=='NO_BLAST_HIT')}")
print(f"  LOW_SIM:      {sum(1 for r in rows if r['blast_confidence']=='LOW_SIMILARITY')}")

print("\nTop expanded species:")
for r in rows[:15]:
    if r["copy_number_call"] == "expanded":
        print(f"  {r['species_code_full']:30s}  loci={r['n_unique_loci']:3d}  "
              f"type={r['centromere_type']:12s}  "
              f"best_blast={r['best_blast_pct_id']}% ({r['best_blast_ref']})")
