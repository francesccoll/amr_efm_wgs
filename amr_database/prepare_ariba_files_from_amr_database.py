#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from Bio import SeqIO

# ---------------------------------------------------------------------------------------------------------------------
# Notes
# ---------------------------------------------------------------------------------------------------------------------

# This python script is used to create ARIBA-compliant files from the AMR look-up table enterococci_amr_genes.xlsx
# First, save the 'look-up table' and 'gene_sequences' sheets in the xlsx file into separate tab-delimited CSV files
# NOTES:
#   - mutations labelled as undetermined_not_annotated had to be removed for ARIVBA prepareref not to fail


# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "Script is used to check the correct reference allele of mutations in enterococci_amr_genes"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-t", "--look_up_table", action="store", dest="look_up_table",
        help="'look-up table' sheet in enterococci_amr_genes.xlsx as a tab-delimited CSV file",
        required=True, metavar="TABLE"
    )
    group.add_argument(
        "-g", "--gene_sequences", action="store", dest="gene_sequences",
        help="'gene_sequences' sheet in enterococci_amr_genes.xlsx as a tab-delimited CSV file",
        required=True, metavar="GENES"
    )
    group.add_argument(
        "-d", "--gene_sequences_dir", action="store", dest="gene_sequences_dir",
        help="directory where gene FASTA sequences are stored",
        required=True, metavar="DIR"
    )
    group.add_argument(
        "-m", "--ariba_metadata_file", action="store", dest="ariba_metadata_file",
        help="name of ariba_metadata_file to save",
        required=True, metavar="META"
    )
    group.add_argument(
        "-f", "--ariba_fasta_file", action="store", dest="ariba_fasta_file",
        help="name of ariba_fasta_file to save",
        required=True, metavar="FASTA"
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
    if not os.path.isfile(args.look_up_table):
        logging.error(f'Input look_up_table file {args.look_up_table} not found!')
        sys.exit(-1)
    if not os.path.isfile(args.gene_sequences):
        logging.error(f'Input gene_sequences file {args.gene_sequences} not found!')
        sys.exit(-1)
    if not os.path.exists(args.gene_sequences_dir):
        logging.error(f'Input gene_sequences_dir directory {args.gene_sequences_dir} not found!')
        sys.exit(-1)

    # Saving paths to gene sequences
    gene_fasta_files = dict()
    with open(args.gene_sequences, 'r') as input_genes:
        input_genes.readline()
        for line in input_genes:
            (gene_name, _, _, _, sequence_file, *_) = line.strip().split('\t')
            gene_fasta_file = args.gene_sequences_dir + sequence_file
            gene_fasta_files[gene_name] = gene_fasta_file
            print("gene_fasta_files[" + gene_name + "] --> " + gene_fasta_files[gene_name])
            if not os.path.isfile(gene_fasta_file):
                logging.error(f'gene_fasta_file {gene_fasta_file} not found!')
                sys.exit(-1)
            # Making sure FASTA header matches gene_name
            first_record = next(SeqIO.parse(gene_fasta_files[gene_name], "fasta"))
            first_record_id = first_record.id
            if first_record_id != gene_name:
                logging.error(f'first_record_id {first_record_id} different from gene_name {gene_name}')
                sys.exit(-1)

    # Extracting mutations and genes --> creating ARIBA metadata file
    metadata_file = open(args.ariba_metadata_file, 'w')
    saved_genes = dict() # dictionary to avoid saving the same gene more than once
    saved_mutations = dict() # dictionary to avoid saving the same mutation more than once
    with open(args.look_up_table, 'r') as input_lines:
        input_lines.readline()
        for line in input_lines:
            #print(line)
            (_,_,_,effect,_,type,mutation,gene_names,*_) = line.strip().split('\t')
            # print(mutation)
            # Make sure gene_names are stored in gene_fasta_files
            gene_names = gene_names.split(';')
            for gene_name in gene_names:
                if gene_name not in gene_fasta_files:
                    logging.warning(f'gene_name {gene_name} not in gene_fasta_file')
            # saving individual acquired genes
            if mutation == 'na':
                for gene_name in gene_names:
                    if gene_name in gene_fasta_files:
                        if gene_name not in saved_genes:
                            # for each gene, whether acquired or core, will be added a as presence/absence gene
                            metadata_line = gene_name + '\t1\t0\t.\t.\t.' + '\n'
                            print('metadata_line ' + metadata_line)
                            metadata_file.write(metadata_line)
                        saved_genes[gene_name] = 1

            # saving individual mutations
            if mutation != 'na':
                mutations = mutation.split('+')
                for mut in mutations:
                    if mut not in saved_mutations:
                        (gene_name, change) = mut.split('@')
                        if gene_name not in gene_fasta_files:
                            logging.error(f'gene_name {gene_name} in {mut} not in gene_fasta_files')
                            sys.exit(-1)
                        # extract mutation's reference allele and position
                        ref = list(change)[0]
                        # determine whether nucleotide or amino acid change
                        type = 'unknown'
                        if ref.islower():
                            type = 'nucleotide'
                        if ref.isupper():
                             type = 'amino_acid'
                        if type == 'unknown':
                            logging.warning(f'Unknown reference type: {ref} in {mut} --> mutation not parsed')
                        if type == 'amino_acid':
                            metadata_line = gene_name + '\t1\t1\t' + change.upper() + '\t.\t.' + '\n'
                            print('metadata_line ' + metadata_line)
                            metadata_file.write(metadata_line)
                        if type == 'nucleotide':
                            metadata_line = gene_name + '\t0\t0\t' + change.upper() + '\t.\t.' + '\n'
                            print('metadata_line ' + metadata_line)
                            metadata_file.write(metadata_line)
                    saved_mutations[mut] = 1
                    saved_genes[gene_name] = 1
    metadata_file.close()

    # Creating ARIBA fasta file
    gene_sequences = []
    for gene_name in saved_genes:
        first_record = next(SeqIO.parse(gene_fasta_files[gene_name], "fasta"))
        gene_sequences.append(first_record)
    print("%i gene sequences saved" % len(gene_sequences))
    SeqIO.write(gene_sequences, args.ariba_fasta_file, "fasta")


if __name__ == "__main__":
    _main()