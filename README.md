### Summary

This repository contains two scripts for running a quick sanity check that a QIIME 2 taxonomic classifier (based on Naive Bayes) is performing roughly as expected.

The motivation for writing these scripts is that it is sometimes possible for these classifiers to break depending on the input dataset. In particular, one classifier I created seemed to work fine on a small test set, but on more realistic datasets it called most ASVs as the genus Alteromonas. This issue seems to have occurred for at least a few other people as well (see [here](https://forum.qiime2.org/t/wrong-taxonomy-qzv-file/5287/3)). This appears to only be an issue when creating custom classifiers and I believe it is mainly related to which reference sequences are included for training (i.e. using more conservative options for ```qiime feature-classifier extract-reads``` is likely all that is needed to resolve these kinds of problems).

Although it's easy to identify red flags like all ASVs being called the same genus, it would be easy to miss if a classifier was failing to classify only particular lineages. I wrote the basic scripts in this repository to address this issue, which can be used to check that no major misclassifications are occurring when classifying rarer lineages with QIIME 2.

### Simplistic taxa classifications

This sanity check is performed by comparing QIIME 2 classifications to simplistic classifications based on percent identity scores against the same database of reference sequences (but full-length). These identities are used to determine classifications at each taxonomic level based on >=90% of matching weights supporting a particular classification. The weight of each matching reference sequence to the query ASV is calculated as 2^id (where id is the percent identity). Only matches above 80% are retained and no lower levels are classified if a higher taxonomic level cannot be classified. There are also heuristics regarding whether it is possible to classify a given taxonomic rank - a minimum match of 99%, 97%, 95%, 93%, 91%, 89%, and 80% is required to make a classification at the species, genus, family, order, phylum, and kingom levels, respectively. There are many more sophisticated approaches for running taxonomic classification, but the advantage of this approach is that the general lineages should be correct (especially at higher levels) and it is very simple to back-track and figure out why classifications with this approach disagree with QIIME 2 classifiers.

### Steps

The first step is to run an alignment of the query ASVs against the reference sequences, which in my case was against the SILVA v132 99% 16S rRNA gene sequences (`SILVA_132_QIIME_release/rep_set/rep_set_all/99/silva132_99.fna`). This can be done using a number of methods, but one way is to use `vsearch`:

```
vsearch --usearch_global representative_sequences.fna \
        --db SILVA_132_QIIME_release/rep_set/rep_set_all/99/silva132_99.fna \
        --id 0.8 \
        --blast6out rep_seqs_vs_SILVA_blast6out.txt
```

You can then run the simplistic taxa classifier with the below command, which requires the `ray`, `numpy`, and `pandas` python packages:
```
python simplistic_taxa_classifier.py -i rep_seqs_vs_SILVA_blast6out.txt \
                                     -f representative_sequences.fna \
                                     -t SILVA_132_QIIME_release/taxonomy/16S_only/99/majority_taxonomy_7_levels_header.txt \
                                     -p 1 > simplistic_taxa.tsv
```

Note that you can speed this command up by specifying more processors (```-p```). Also, the reference taxonomy file differs from the default SILVA file, because I added a headerline. The head of the taxonomy file looks like this:

```
Feature ID      Taxon
KF494428.1.1396 D_0__Bacteria;D_1__Epsilonbacteraeota;D_2__Campylobacteria;D_3__Campylobacterales;D_4__Thiovulaceae;D_5__Sulfuricurvum;D_6__Sulfuricurvum sp. EW1
AF506248.1.1375 D_0__Bacteria;D_1__Cyanobacteria;D_2__Oxyphotobacteria;D_3__Nostocales;D_4__Nostocaceae;D_5__Nostoc PCC-73102;D_6__Nostoc sp. 'Nephroma expallidum cyanobiont 23'
EF603722.1.1487 D_0__Bacteria;D_1__Bacteroidetes;D_2__Bacteroidia;D_3__Bacteroidales;D_4__Muribaculaceae;D_5__uncultured bacterium;D_6__uncultured bacterium
KX826903.1.1438 D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;D_3__Alteromonadales;D_4__Pseudoalteromonadaceae;D_5__Pseudoalteromonas;D_6__Pseudoalteromonas sp.
...
```

The output taxonomy classifications should look like this:

```
Sequence        Kingdom Phylum  Class   Order   Family  Genus   Species
bacteriaV6V8_seq1       D_0__Bacteria   D_1__Bacteroidetes      D_2__Bacteroidia        D_3__Chitinophagales    D_4__uncultured Unclassified    Unclassified
bacteriaV6V8_seq2       D_0__Bacteria   D_1__Verrucomicrobia    D_2__Verrucomicrobiae   D_3__Verrucomicrobiales D_4__Rubritaleaceae     D_5__Rubritalea Unclassified
bacteriaV6V8_seq3       D_0__Bacteria   D_1__Gemmatimonadetes   D_2__BD2-11 terrestrial group   Unclassified    Unclassified    Unclassified    Unclassified
...
```

To compare these classifications to QIIME 2 classifications (that have been exported to TSV) you can use this command:
```
python compare_simplistic_and_qiime2_taxa.py -q taxonomy.tsv \
                                             -a simplistic_taxa.tsv \
                                             -o taxa_diff_comparison.tsv
```

Where ```taxonomy.tsv``` contains the QIIME 2 classificaitons and ```taxa_diff_comparison.tsv``` is the output table.

The output file looks like this:
```
asv	label_diff	any_diff	alt_Kingdom	alt_Phylum	alt_Class	alt_Order	alt_Family	alt_Genus	alt_Species	q2_Kingdom	q2_Phylum	q2_Class	q2_Order	q2_Family	q2_Genus	q2_Species
bacteriaV6V8_seq12296	5.0	5.0	D_0__Bacteria	D_1__Marinimicrobia (SAR406 clade)	D_2__uncultured bacterium	D_3__uncultured bacterium	D_4__uncultured bacterium	D_5__uncultured bacterium	D_6__uncultured bacterium	D_0__Bacteria	D_1__Marinimicrobia (SAR406 clade)	Ambiguous_taxa	Ambiguous_taxa	Ambiguous_taxa	Ambiguous_taxa	Ambiguous_taxa
bacteriaV6V8_seq538	4.0	4.0	D_0__Bacteria	D_1__Planctomycetes	D_2__OM190	D_3__uncultured bacterium	D_4__uncultured bacterium	D_5__uncultured bacterium	D_6__uncultured bacterium	D_0__Bacteria	D_1__Planctomycetes	D_2__OM190	D_3__uncultured deep-sea bacterium	D_4__uncultured deep-sea bacterium	D_5__uncultured deep-sea bacterium	D_6__uncultured deep-sea bacterium
b
...
...
...
bacteriaV6V8_seq4627	0.0	1.0	D_0__Bacteria	D_1__Proteobacteria	D_2__Deltaproteobacteria	D_3__Myxococcales	D_4__Sandaracinaceae	Unclassified	Unclassified	D_0__Bacteria	D_1__Proteobacteria	D_2__Deltaproteobacteria	D_3__Myxococcales	D_4__Sandaracinaceae	D_5__uncultured	Unclassified
bacteriaV6V8_seq4628	0.0	0.0	D_0__Bacteria	D_1__Planctomycetes	Unclassified	Unclassified	Unclassified	Unclassified	Unclassified	D_0__Bacteria	D_1__Planctomycetes	Unclassified	Unclassified	Unclassified	Unclassified	Unclassified
bacteriaV6V8_seq13805	0.0	1.0	D_0__Bacteria	D_1__Proteobacteria	D_2__Deltaproteobacteria	Unclassified	Unclassified	Unclassified	Unclassified	D_0__Bacteria	D_1__Proteobacteria	D_2__Deltaproteobacteria	D_3__PB19	Unclassified	Unclassified	Unclassified
```

Most of the columns contain the taxonomic labels output by QIIME 2 (```q2```) or the simplistic approach (```alt``` for "alternative"). Differences under ```label_diff``` correspond to differences in taxonomic labels where that rank is something besides "Unclassified" by each classifier (missing ranks are set to "Unclassified" by default when parsing the QIIME 2 output). These are the differences we're most interesed in since if any drastically different classifications between the two pipelines would have high ```label_diff``` values. The table is sorted by this column in descending order to help spot potentially problematic ASVs / lineages. The other column is ```any_diff```, which is the same thing but also includes "Unclassified" labels when calculating the number of differences in the taxonomic labels. This metric isn't as useful because many lineages can be classified with different precision based on a Naive Bayes approach vs simplistic percent identity.

