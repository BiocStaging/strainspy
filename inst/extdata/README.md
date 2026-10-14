# Example data provided in strainspy

All files are gzip-compressed. This is a random 200-sample subset of the
Parkinson's disease gut metagenome study by Wallen et al. 2022
(https://doi.org/10.1038/s41467-022-34667-x). There are 100 PD cases
and 100 controls. 

## example_metadata.csv.gz

Sample metadata with one row per sample. The first column is the
`run_accession`, the sample identifier that matches the profile files.
`Case_status` (`PD` or `Control`) is the outcome used in examples.
The remaining columns are host, diet, gastrointestinal symptom and medication
variables from the original study.

## example_sylph_profile.tsv.gz

Output of `sylph profile` for the 200 samples against the GTDB r220
species-representative genome database: 3,010 genomes detected, 48,494
sample-by-genome rows.

## example_taxonomy.tsv.gz

Two-column, headerless taxonomy for the GTDB r220 representative genomes:
assembly accession and GTDB lineage string (`d__;p__;c__;o__;f__;g__;s__`).
Read with `read_taxonomy()`. It covers every genome in the GTDB r220 representative database used for the bundled profiles.
Additional taxonomy files are available at https://doi.org/10.5281/zenodo.21796878.

## example_sourmash.csv.gz

`sourmash gather` output for 20 of the 200 samples against the same GTDB r220
genomes: 1,508 genomes, 7,032 rows, merged across samples into one CSV in the
layout produced by `merge_sourmash_files()`. This file can be read 
with `read_sourmash()`.

## metaphlan_merged.tsv.gz and metaphlan_taxonomy.tsv.gz

MetaPhlAn 4 species-level genome bin (SGB) relative abundances for the 200
samples (1,831 SGBs, 40,026 rows) and the matching SGB-to-lineage table. Both
are the outputs of `merge_metaphlan_files()` applied to the per-sample
MetaPhlAn profiles. These files can respectively be read using `read_metaphlan()` 
and `read_taxonomy()`.
