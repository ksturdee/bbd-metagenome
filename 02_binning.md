BBD Metagenome: Binning
================
Kendra
2026-11-04

Binning of the assembled contigs into genomes with MetaBAT 2. Commands
for this step were run in `metabat/`.

To assess contig depth, contig names and their sequencing depths were
extracted from the `Flye` assembly information file to create the depth
file required for binning. The assembled contigs were also reorder to
match the order of contigs in the depth file.

``` bash
cat ../flye/bbd/assembly_info.txt  | grep -v '#' | 

awk '{OFS="\t";print $1,$3}' > depths_bbd.txt

cat depths_bbd.txt | awk '{print $1}' | xargs samtools faidx ../flye/bbd/assembly.fasta > assembly_bbd_ordered.fasta
```

Finally, `metaBAT 2` was used to group the reordered contigs into bins
based on sequence composition, and contig depth.
