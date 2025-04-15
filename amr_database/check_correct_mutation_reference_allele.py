#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
# from Bio.Alphabet import ThreeLetterProtein


# ---------------------------------------------------------------------------------------------------------------------
# Notes
# ---------------------------------------------------------------------------------------------------------------------

# This python script is used to check the correct reference allele of mutations in enterococci_amr_genes.xlsx


# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "Script is used to check the correct reference allele of mutations in enterococci_amr_genes"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-m", "--input_mutations", action="store", dest="input_mutations",
        help="mutations in 'mutation' column in enterococci_amr_genes.xlsx",
        required=True, metavar="MUTATIONS"
    )
    group.add_argument(
        "-g", "--input_genes", action="store", dest="input_genes",
        help="tab-delimited file containing gene_name and full path to gene's fasta file",
        required=True, metavar="GENES"
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

    # Making sure input files exist
    if not os.path.isfile(args.input_mutations):
        logging.error(f'Input mutations file {args.input_mutations} not found!')
        sys.exit(-1)
    if not os.path.isfile(args.input_genes):
        logging.error(f'Input genes file {args.input_genes} not found!')
        sys.exit(-1)

    # Saving the full path of genes
    gene_fasta_files = dict()
    with open(args.input_genes, 'r') as input_genes:
        for line in input_genes:
            (gene_name, sequence_file) = line.strip().split('\t')
            gene_fasta_files[gene_name] = sequence_file
            print("gene_fasta_files[" + gene_name + "] --> " + gene_fasta_files[gene_name])

    # Reading mutations
    with open(args.input_mutations, 'r') as input_mutations:
        for line in input_mutations:
            mutations = line.strip().split('+')
            for mutation in mutations:
                # extract mutation's gene name
                (gene_name, change) = mutation.strip().split('@')
                # if gene_name == 'rplK':
                if gene_name != 'na':
                    print(mutation)
                    # extract mutation's reference allele and position
                    ref = list(change)[0]
                    # determine whether nucleotide or amino acid change
                    type = 'unknown'
                    if ref.islower():
                        type = 'nucleotide'
                    if ref.isupper():
                        type = 'amino_acid'
                    if type == 'unknown':
                        print('Unknown reference type: ' + ref + ' in ' + change)
                    # extract positition
                    digits = [int(s) for s in list(change) if s.isdigit()]
                    pos = ''
                    for digit in digits:
                        pos += str(digit)
                    # manually changing pos coordinate
                    # if gene_name == 'liaF':
                    #     pos = int(pos) - 8
                    # if gene_name == '23SrRNA':
                    #     pos = int(pos) + 14
                    # make sure gene's corresponding FASTA file exists
                    if gene_name not in gene_fasta_files:
                        logging.error(f"{gene_name} (mutation {change}) not found in {args.input_genes}")
                    # read gene's FASTA file and translate to amino acid sequence
                    first_record = next(SeqIO.parse(gene_fasta_files[gene_name], "fasta"))
                    gene_BioSeq_nt = first_record.seq
                    gene_BioSeq_aa = gene_BioSeq_nt.translate()
                    # compare mutation's reference allele to gene's FASTA
                    ref_fasta = ''
                    if type == 'amino_acid':
                        ref_fasta = gene_BioSeq_aa[int(pos)-1]
                        # print('change ' + change + ' ref ' + ref + ' ref_fasta ' + ref_fasta)
                    if type == 'unknown':
                        ref_fasta = gene_BioSeq_aa[int(pos) - 1]
                        # print('change ' + change + ' ref ' + ref + ' ref_fasta ' + ref_fasta)
                    if type == 'nucleotide':
                        ref_fasta = gene_BioSeq_nt[int(pos)-1]
                        ref = ref.upper()
                        # print('change ' + change + ' ref ' + ref + ' ref_fasta ' + ref_fasta)
                    if ref != ref_fasta:
                        print('mutation ' + mutation + ' pos ' + str(pos) + ' change ' + change + ' ref ' + ref + ' ref_fasta ' + ref_fasta)


if __name__ == "__main__":
    _main()
