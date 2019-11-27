#!/usr/bin/python3

import argparse
import pandas as pd
import numpy as np
import sys

def main():

    parser = argparse.ArgumentParser(
        description=
"Compare taxa based on alternative approach with qiime2 taxa on the same set of ASVs. This is intended to be a quick sanity check on QIIME2 taxa classifications. Note that the alternative approach output table is tab-delimited for each taxonomic level.",

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-q", "--qiime2", metavar="PATH",
                        type=str, help="QIIME2 taxa table.",
                        required=True)

    parser.add_argument("-a", "--alternative", metavar="PATH",
                        type=str, help="Alternative (tab-delimited) taxa table.",
                        required=True)

    parser.add_argument("-o", "--output", metavar="PATH",
                        type=str, help="Output table.",
                        required=True)

    args = parser.parse_args()

    q2_taxa = pd.read_csv(args.qiime2, sep="\t", index_col=0)

    q2_taxa[['q2_Kingdom', 'q2_Phylum', 'q2_Class', 'q2_Order', 'q2_Family', 'q2_Genus', 'q2_Species']] = q2_taxa['Taxon'].str.split(pat=';', expand=True)

    q2_taxa.drop(labels=['Taxon', 'Confidence'], axis=1, inplace=True)

    q2_taxa.replace(to_replace=[None], value="Unclassified", inplace=True)

    alt_taxa = pd.read_csv(args.alternative, sep="\t", index_col=0)

    alt_taxa.rename(columns={"Kingdom": "alt_Kingdom",
                             "Phylum": "alt_Phylum",
                             "Class": "alt_Class",
                             "Order": "alt_Order",
                             "Family": "alt_Family",
                             "Genus": "alt_Genus",
                             "Species": "alt_Species"
                             }, inplace=True)

    all_asvs = list(set().union(list(alt_taxa.index), list(q2_taxa.index))) 

    # First check that all ASVs are found in each table:
    for asv in all_asvs:
        if asv in alt_taxa.index and asv not in q2_taxa.index:
            print("asv " + asv + " not found in qiime2 taxa table", file=sys.stderr)
        elif asv not in alt_taxa.index and asv in q2_taxa.index:
            print("asv " + asv + " not found in alternative taxa table", file=sys.stderr)


    # Merge taxa dfs:
    combined = alt_taxa.join(q2_taxa, how='outer')

    orig_cols = combined.columns.tolist()

    combined["any_diff"] = np.nan

    combined["label_diff"] = np.nan

    combined = combined.reindex(columns=["label_diff", "any_diff"] + orig_cols)

    # Figure out the number of differences in taxa levels for each ASV and
    # then sort table to be in decreasing order before printing out.
    for asv in all_asvs:
        any_diff = 0
        label_diff = 0

        taxa_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

        for level in taxa_levels:
            q2_level = 'q2_' + level
            alt_level = 'alt_' + level

            if combined.loc[asv, q2_level] != combined.loc[asv, alt_level]:
                any_diff += 1

                if combined.loc[asv, q2_level] != "Unclassified" and combined.loc[asv, alt_level] != "Unclassified":
                    label_diff += 1

        combined.loc[asv, "any_diff"] = any_diff
        combined.loc[asv, "label_diff"] = label_diff

    combined.sort_values(by="label_diff", axis=0, inplace=True, ascending=False)

    combined.to_csv(path_or_buf=args.output, sep="\t", index_label="asv")


if __name__ == '__main__':
    main()
