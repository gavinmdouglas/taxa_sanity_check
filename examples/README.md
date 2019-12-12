This folder contains the files required for the sanity checks run for the custom classifiers that are available as part of [Microbiome Helper](https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2019.7)).

The subfolders are:

* `test_infiles` - Test ASVs in FASTA format (and as QIIME 2 artifacts) corresponding to each amplicon region used for the Microbiome Helper custom classifiers.
* `simplistic_taxa` - The taxa assignments for the test ASVs based on the simplistic approach described in this repo.
* `usearch_global_out` - Output table of test ASV similarity matches against reference sequences (used for simplistic taxonomic assignment). Kept here in case looking at example outputs is useful.
* `taxa_diff` - Output comparison tables from `compare_simplistic_and_qiime2_taxa.py` contrasting taxonomy assignments by QIIME 2 classifiers with simplistic approach.

**These subfolders contain gzipped files that would need to be decompressed before using them as inputs.**

The raw commands to run QIIME 2 taxa assignment and compare the simplistic taxa assignments to these QIIME 2 outputs are found in `classifier_compare_example_cmds.txt`, which is very rough code just kept here in case so that it can be easily re-run in the future. 

