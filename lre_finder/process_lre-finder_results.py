#!/usr/bin/env python3

import argparse
import logging
import os
import sys


# ---------------------------------------------------------------------------------------------------------------------
# Development notes
# ---------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------
# Usage notes
# ---------------------------------------------------------------------------------------------------------------------
# This script is used to call linezolid resistance from the output files of LRE-Finder

# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "This script is used to call linezolid resistance from the output files of LRE-Finder.\n" \
                  "This script outputs a file with one line and three columns: " \
                  "sample_id\tlinezolid_genotype\tlinezolid_phenotype.\n"

    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('I/O arguments')
    group.add_argument(
        "-i", "--sample_id", action="store", dest="sample_id",
        help="sample ID used as prefix of LRE-Finder output files",
        required=True, metavar="SAMPLE_ID"
    )
    group.add_argument(
        "-r", "--res_file", action="store", dest="res_file",
        help="LRE-Finder .res output file",
        required=True, metavar="RES"
    )
    group.add_argument(
        "-p", "--pos_file", action="store", dest="pos_file",
        help="LRE-Finder .pos output file",
        required=True, metavar="POS"
    )
    group.add_argument(
        "-s", "--output_suffix", action="store", dest="output_suffix",
        help="Output suffix for the LRE-Finder output generated by this script",
        required=False, metavar="SUFFIX", default=".lre_finder.txt"
    )
    return parser.parse_args()


# ------------------------------------------------------------------------------------
# Main program
# ------------------------------------------------------------------------------------

def _main():
    # Configure logging
    logging.basicConfig(
        format='%(asctime)s %(levelname)s: %(message)s',
        level=logging.INFO
    )
    # Get arguments
    args = parse_arguments()

    # Making sure expected input files exist
    if not os.path.exists(args.res_file):
        logging.error(f'Input file in --res_file {args.res_file} not found!')
        sys.exit(-1)
    if not os.path.exists(args.pos_file):
        logging.error(f'Input file in --pos_file {args.pos_file} not found!')
        sys.exit(-1)

    # Loading gene information from .res file
    # Only acquired genes are saved, i.e. 23S gene is ignored
    genes_found = list()
    logging.info(f"Loading gene information from {args.res_file}")
    with open(args.res_file, 'r') as res_file:
        for res_line in res_file:
            if not res_line.startswith('#'):
                gene_id = res_line.strip().split('\t')[0]
                if not gene_id.startswith('23S_'):
                    print('Acquired gene ' + gene_id + ' found')
                    genes_found.append(gene_id)

    # Loading gene information from .pos file
    # LRE-Finder only detects 23S rRNA mutations G2576T and G2505A
    # A cut-off of 10% mutations in position 2576 and/or position 2505 was therefore set in LRE-Finder for predicting
    # a linezolid resistance phenotype. This cut-off allows for detection of a single mutated 23S allele in both
    # E. faecalis and E. faecium, while ignoring low-level sequencing noise.
    mutations_found = list()
    per_cutoff = 10
    logging.info(f"Loading mutations information from {args.pos_file}")
    with open(args.pos_file, 'r') as pos_file:
        for pos_line in pos_file:
            if not pos_line.startswith('#'):
                items = pos_line.strip().split('\t')
                pos = int(items[1])
                # #Template	pos	ref	A	C	G	T	N	-	A[%]	C[%]	G[%]	T[%]	N[%]	-[%]
                if pos == 2576:
                    per_reads = round(float(items[12]))
                    if per_reads > per_cutoff:
                        mutation = "G2576T("+str(per_reads)+"%)"
                        mutations_found.append(mutation)
                if pos == 2505:
                    per_reads = round(float(items[9]))
                    if per_reads > per_cutoff:
                        mutation = "G2505A("+str(per_reads)+"%)"
                        mutations_found.append(mutation)

    # Calling resistance
    lin_phen = "S"
    if len(genes_found) > 0 or len(mutations_found):
        lin_phen = "R"

    # Writing output file
    print('All acquired genes found: ' + ';'.join(genes_found))
    print('All acquired mutations found: ' + ';'.join(mutations_found))
    header = 'sample_id\tpredicted_phenotype\tgenotype\n'
    tmp = mutations_found + genes_found
    genotype = '-'
    if len(tmp) > 0:
        genotype = ';'.join(tmp)
    output_line = args.sample_id + '\t' + lin_phen + '\t' + genotype + '\n'
    output_file = args.sample_id + args.output_suffix
    output = open(output_file, 'w')
    output.write(header + output_line)
    output.close()


if __name__ == "__main__":
    _main()
