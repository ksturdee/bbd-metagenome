BBD Metagenome: Gene Annotation
================
Kendra
2026-11-19

<!--
Explain what this step does and why, and copy each command you ran into a
bash chunk. R chunks run from the top of the repository.
-->

Gene prediction and annotation of the most complete bins with Prokka.
Commands for this step were run in `prokka/`.

The goal here is to assess the genes present within the genomes of bbd
taxa within the mat to further our understanding of metabolic processes
present. Prokka will identify the gene sequences and identify homologous
genes with known functions.

To run `Prokka`, the bin `metabat.30` bin was chosen as a test:

``` bash
cat ../checkm/bbd_topbins.txt

bin=metabat.30
prokka --cpus 2 --outdir $bin --prefix $bin  ../metabat/bins_bbd/${bin}.fa

ls metabat.30
```

Next, `Prokka` was run separately on each genome listed in
`bbd_topbins.txt`. Each bin was given a seperate FASTA file, with its
bin identifier used to name any outputs. This was done using the slurm
script `run_prokka.sh` within `/prokka`. This produced predicted genes
and associated functional information for each genome bin for further
pathway analysis.
