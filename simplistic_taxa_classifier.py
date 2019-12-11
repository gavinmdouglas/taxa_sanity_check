#!/usr/bin/python3

import argparse
import sys
import pandas as pd
import numpy as np
import gzip
import ray

def main():

    parser = argparse.ArgumentParser(
        description=
"Simplistic taxonomic classifier based on BLAST 6 output table.",

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="PATH",
                        type=str, help="Input table of BLAST6OUT.",
                        required=True)

    parser.add_argument("-f", "--fasta", metavar="PATH",
                        type=str, help="FASTA file of study ASVs.",
                        required=True)

    parser.add_argument("-t", "--taxa", metavar="PATH",
                        type=str, help="Reference taxonomy database.",
                        required=True)

    parser.add_argument("--min_prob", metavar="FLOAT",
                        type=float, help="Max probability for taxon to be unclassified (i.e. a taxon needs probability greater than this value to be called).",
                        required=False, default=0.9)

    parser.add_argument("-p", "--proc", metavar="INT",
                        type=int, help="Number of processes to run in parallel.",
                        required=False, default=1)

    args = parser.parse_args()

    ray.init(num_cpus=args.proc)

    print("Reading in files.", file=sys.stderr)

    seqs = read_fasta(args.fasta)

    # The taxonomic and match tables will be made into global variables to make them easier to run.
    matches = pd.read_csv(args.input, sep="\t", header=None,
                             names=['query', 'target', 'identity',
                             'align_length', 'num_mismatch', 'num_gap_opens',
                             'query_start', 'query_end', 'target_start',
                             'target_end', 'evalue', 'bitscore'])

    taxa = pd.read_csv(args.taxa, sep="\t", index_col=0)

    print("Started parsing taxa table.", file=sys.stderr)

    taxa[['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']] = taxa['Taxon'].str.split(pat=';', expand=True)

    print("Done parsing taxa table.", file=sys.stderr)

    print("Started transforming identity.", file=sys.stderr)

    matches['exp_identity'] = 2 ** matches['identity'] / 2 ** 80

    print("Done transforming identity.", file=sys.stderr)

    print("\t".join(['Sequence', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']))

    matches_mem = ray.put(matches)
    taxa_mem = ray.put(taxa)

    print("Classifying ASVs.", file=sys.stderr)

    assignments_raw = [assign_asv.remote(asv, matches_mem, taxa_mem, args.min_prob) for asv in seqs.keys()]

    for assignment in ray.get(assignments_raw):
        print(assignment)


@ray.remote
def assign_asv(asv_id, in_matches, in_taxa, prob_cutoff):
    '''Assign taxonomy to each ASV based hits to target sequences in BLAST6OUT
    table.'''

    taxa_max_cutoff = {}
    taxa_max_cutoff["Species"] = 99
    taxa_max_cutoff["Genus"] = 97
    taxa_max_cutoff["Family"] = 95
    taxa_max_cutoff["Order"] = 93
    taxa_max_cutoff["Class"] = 91
    taxa_max_cutoff["Phylum"] = 89
    taxa_max_cutoff["Kingdom"] = 80

    asv_matches = in_matches.loc[in_matches['query'] == asv_id]

    # If no matches then output as Unassigned and move on.
    if asv_matches.shape[0] == 0:
        return(asv_id + "\tUnassigned" + "\tUnclassified" * 6)

    # Otherwise weight all ASV matches by the sequence similarity, where
    # the weight value is given by 2^(similarity).
    # For each taxonomic level only accept a classification if >=90% of
    # the weight supports a call, otherwise leave that level blank.

    matching_ids = list(asv_matches['target'])

    matching_taxa = in_taxa.loc[matching_ids, :]

    max_identity = asv_matches['identity'].max()

    asv_taxonomy = ""

    last_unclassified = False

    for level in ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']:

        if max_identity < taxa_max_cutoff[level] or last_unclassified:
            asv_taxonomy = asv_taxonomy + "\t" + "Unclassified"
            continue

        asv_level_matches = asv_matches.copy()

        asv_level_matches[level] = list(matching_taxa.loc[asv_matches['target'], level])

        asv_level_matches['prob'] = asv_matches['exp_identity'] / asv_matches['exp_identity'].sum()

        asv_level_matches = asv_level_matches.loc[:, [level, 'prob']]

        asv_level_pivot = pd.pivot_table(asv_level_matches, index=level, values='prob', aggfunc=np.sum)

        max_prob = asv_level_pivot['prob'].max()

        if max_prob > prob_cutoff:
            asv_taxonomy = asv_taxonomy + "\t" + list(asv_level_pivot.loc[asv_level_pivot['prob'] == max_prob].index)[0]
        else:
            asv_taxonomy = asv_taxonomy + "\tUnclassified"
            last_unclassified = True

    return(asv_id + asv_taxonomy)


def read_fasta(filename):
    '''Read FASTA file into a dictionary - ids are keys and sequences are
    values.'''

    seq = {}

    name = None

    fasta_in = open(filename, "r")

    for line in fasta_in:

        if line[0] == ">":

         
            name = line[1:]

            name = name.rstrip("\r\n")

            seq[name] = ""

        else:

            line = line.rstrip("\r\n")

            seq[name] += line

    fasta_in.close()

    return seq


if __name__ == '__main__':
    main()
