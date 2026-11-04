BBD Metagenome: Bin Quality
================
Kendra
2026-11-04

<!--
Explain what this step does and why, and copy each command you ran into a
bash chunk. R chunks run from the top of the repository.
-->

Completeness and contamination of each bin were estimated with `CheckM`,
and the most complete bins were chosen for further analysis. Commands
for this step were run in `checkm/`.

Assembly contigs were sorted by sequence length to find the largest
contig:

``` bash
cat ../flye/bbd/assembly.fasta | bioawk -c fastx '{print length($seq),$name}' | sort -n
```

The largest contig in the bbd assembly was `contig_120`. This contig
must be manually saved as an independent bin, in the same folder as the
metabat bins, so that all the later steps treat it just like the other
bins. Then, must copy all of the bins, including `contig_120` into the
`bins_bbd` folder inside `CheckM`.

``` bash
samtools faidx ../flye/bbd/assembly.fasta contig_120 > ../metabat/bins_bbd/contig_120.fa
cp ../metabat/bins_bbd/*.fa bins_bbd/
```

To assess quality with `CheckM`, the slurm script ‘run_checkm_bbd.sh’
was created and run:

``` bash
sbatch run_checkm_bbd.sh
```

## Functional analysis of the BBD Metagenome

The purpose of binning is to sort the metagenome into species and strain
within the black band disease community. They can also be leveraged to
make inferences about metabolic processes present within a BBD lesion.
Overall, the purpose of this functional analysis is to investigate (1)
what taxonomic groups of microbes are present, and (2) what can they do
in terms of genes and biochemical pathways present.

The results were checked within `checkm_results_bbd`

``` bash
head checkm_results_bbd.txt 
```

Marker lineage for top few bins: o\_\_Cytophagales, p\_\_Cyanobacteria,
k\_\_Bacteria, and o\_\_Clostridiales.

Look further down the list:

``` bash
head -n 20 checkm_results_bbd.txt 
```

The analysis will focus on the most complete bins, i.e bins with a
completeness greater than 15% (more likely to contain the majority of
useful information on function within the bbd community). These bins
will likely represent the most abundant taxa (and best assembled).

Extract top 15% bins for completeness:

``` bash
cat checkm_results_bbd.txt | awk '$13>15'

cat checkm_results_bbd.txt | awk '$13>15{print $1}'

cat checkm_results_bbd.txt | awk '$13>15{print $1}' > bbd_topbins.txt
```

Saved output as `bbd_topbins.txt.`
