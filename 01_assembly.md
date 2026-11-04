BBD Metagenome: Assembly
================
Kendra
2026-11-04

<!--
Explain what this step does and why, and copy each command you ran into a
bash chunk. R chunks run from the top of the repository.
-->

Nanopore reads from the black band disease (BBD) metagenome were
assembled using `Flye` in metagenome mode. Commands for this stage were
run from the `flye/` directory.

A symbolic link to the BBD nanopore read file was added to ensure that
the sequencing is accessed without duplicating the large input file:

``` bash
ln -s /pvol/data/metagenome_assembly/bbd_0.1.fastq.gz .
```

Then, the `flye` assembly was run using the slurm shell script
`run_flye_bbd.sh`, and the slurm job output was assessed:

``` bash
sbatch run_flye_bbd.sh

ls -l bbd/assembly.fasta
```

The completed assembly (written by `Flye`) was saved to as
`assembly.fasta`.
