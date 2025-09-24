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



