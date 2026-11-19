# The Metagenome of Black Band Disease

```{=html}
<!--
This repository documents the assembly, binning and functional analysis of a nanopore metagenome from black band disease (BBD) lesion the coral Acropora hyacinthus, done for Workshop 4 of the MB5370 Marine Genomics module. The samples were collected at Orpheus Island in 2023 by Julia Hung. The data subset was provided for the MB5370 course; it is not the full sequencing dataset.
-->
```

## Background and aim

Black band disease is a microbial disease of corals. In this project, I assembled the nanopore reads, grouped contigs into bins, assessed bin completeness and contamination, and used taxonomic and gene annotations to investigate which microbial groups and candidate metabolic pathways were present. Pathway predictions are based on gene annotations.

## Data and reference

The raw reads and large intermediate files are not included in this repository. As well, the workshop describes and references the associated study, [“Microbial Metabolism and Disease Virulence Changes Across Day and Night in Coral Black Band Disease Lesions”](https://doi.org/10.1111/1462-2920.70219).

## Main results

CheckM estimated completeness over 15% for 15 bins. The two most complete entries were `contig_120` (97.62% completeness, 1.64% contamination) and `metabat.42` (97.50% completeness, 0.67% contamination). For the photosynthesis light-reaction pathway (`PWY-101`), four bins or contigs had at least one of the expected gene families. `metabat.42` had all four families and MinPath retained the pathway; MinPath also retained it for `metabat.15`, where one of four families was detected. For the tetrathionate-reduction pathway (`PWY-5360`), `metabat.12` had one of three expected families, but MinPath did not retain the pathway.

## Analysis

```{=html}
<!--
Each link opens the knitted version of one RMarkdown file. Knit each .Rmd,
then add both the .md file and its _files folder (if there is one) to git.
-->
```

1.  [Assembly](01_assembly.md)
2.  [Binning](02_binning.md)
3.  [Bin quality](03_bin_quality.md)
4.  [Taxonomy](04_taxonomy.md)
5.  [Gene annotation](05_gene_annotation.md)
6.  [Metabolic pathways](06_metabolic_pathways.md)
7.  [Data visualization source](07_Data_Visualization.Rmd)

- [Combined report](index.html)

## Repository layout

Each step has a numbered RMarkdown file at the top of the repository and a folder where its commands are run.

| Step               | Notes                       | Folder     | Tool      |
|--------------------|-----------------------------|------------|-----------|
| Assembly           | `01_assembly.Rmd`           | `flye/`    | Flye      |
| Binning            | `02_binning.Rmd`            | `metabat/` | MetaBAT 2 |
| Bin quality        | `03_bin_quality.Rmd`        | `checkm/`  | CheckM    |
| Taxonomy           | `04_taxonomy.Rmd`           | `gtdbtk/`  | GTDB-Tk   |
| Gene annotation    | `05_gene_annotation.Rmd`    | `prokka/`  | Prokka    |
| Metabolic pathways | `06_metabolic_pathways.Rmd` | `minpath/` | MinPath   |

`index.Rmd` joins all six steps into one web page. Knit it to make `index.html`.
