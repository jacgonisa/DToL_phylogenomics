We assume we already have all the BUSCO output, which was obtained running

```
busco -i genomes/ -l eukaryota_odb10 -m geno -o BUSCO_genomes -c 15 
```

We then parse the shared BUSCOs

```
python scripts/parse_shared_buscos.py
```


We can also detect absent BUSCOs. Here, absent BUSCOs would be those that are not labeled as "Complete". So "Duplicated" BUSCOs, which is the mejority of BUSCOs in polyploid species...

```
##Considering duplicated buscos as absent
python scripts/find_absent_busco.py


##Not considerating duplicated buscos as absent
python scripts/find_absent_busco_duplicatednotincluded.py

```



Retrieving BUSCOs fasta files
```
bash scripts/get_fasta_busco_bioython.sh
```
Trim alignments

```
bash scripts/trim_alignments.sh
```

And rename to change the name of the protein with only the species name, so we can concatenate

```
scripts/rename_alignments.sh

mv grouped_busco_fastas/review_allbuscos_442species/trimmed_msas/*.renamed.msa grouped_busco_fastas/review_allbuscos_442species/renamed_trimmed_msas/

```

And finally, concatenate alignment with iqtree built-in function


```
iqtree2 -p grouped_busco_fastas/review_allbuscos_442species/renamed_trimmed_msas/ --out-aln concatenated_alignment/concatenated_alignment_255busco_432species
```

Build the tree, with one model per partition

```
iqtree -s concatenated_alignment/concatenated_alignment_255busco_432species -p concatenated_alignment/concatenated_alignment_255busco_432species.nex -m MFP+MERGE -rcluster 10 -nt AUTO -B 1000
```

I am also building a faster tree, with one partition for the full supermatrix

```
iqtree \
  -s concatenated_alignment/concatenated_alignment_255busco_432species \
  -m MFP \
  -nt AUTO \
  -B 1000 \
  -pre trees/run_notpartitioned

```

