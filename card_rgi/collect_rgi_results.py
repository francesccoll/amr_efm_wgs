#!/usr/bin/env python3

import argparse
import logging
import os
import sys

# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "This script is used to collect WGS-predicted AMR phenotypes and genotypes by rgi. " \
                  "Tested with rgi v4.2.2 and CARD database v3.1.0\n " \
                  "It outputs two tab-delimited tables: one with WGS-predicted AMR phenotypes and another one with AMR " \
                  "genotypes, with as many samples as rows and as many columns as antibiotics.\n"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('Required: Input files')
    group.add_argument(
        "-i", "--list_of_samples", action="store", dest="list_of_samples",
        help="List of sample Ids",
        required=True, metavar="LIST_SAMPLES"
    )
    group.add_argument(
        "-a", "--list_of_antibiotics", action="store", dest="list_of_antibiotics",
        help="List of antibiotics to look for in ResFinder output files.",
        required=False, metavar="LIST_ABS"
    )
    group.add_argument(
        "-d", "--rgi_out_dir", action="store", dest="rgi_out_dir",
        help="Directory where RGI output files are stored",
        required=True, metavar="RGI_DIR"
    )
    group.add_argument(
        "-s", "--sample_out_suffix", action="store", dest="sample_out_suffix",
        help="RGI output suffix used to name sample's RGI output files (default '_rgi.txt') ",
        required=False, metavar="OUT_SUFFIX", default='_rgi.txt'
    )
    group.add_argument(
        "-p", "--output_phenotypes_file", action="store", dest="output_phenotypes_file",
        help="Output file with RGI WGS-predicted AMR phenotypes",
        required=True, metavar="PHENO_FILE"
    )
    group.add_argument(
        "-g", "--output_genotypes_file", action="store", dest="output_genotypes_file",
        help="Output file with RGI WGS-detected AMR genotypes",
        required=True, metavar="GENO_FILE"
    )
    group = parser.add_argument_group('Optional:')
    group.add_argument(
        "-l", "--percentage_length", action="store", dest="percentage_length",
        help="Percentage Length of Reference Sequence (default: 95)",
        required=False, metavar="PER_LENGTH", default=95
    )

    return parser.parse_args()


def check_file_exists(my_file):
    if not os.path.isfile(my_file):
        logging.error(f'File {my_file} not found!')
        sys.exit(-1)


def check_path_exists(my_dir):
    if not os.path.exists(my_dir):
        logging.error(f'Directory {my_dir} not found!')
        sys.exit(-1)


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
    logging.info('Making sure input files exist...')
    check_file_exists(args.list_of_samples)
    check_file_exists(args.list_of_antibiotics)
    check_path_exists(args.rgi_out_dir)
    logging.info('All required input files found.')

    # Saving samples ids
    logging.info('Saving sample ids.')
    samples = list()
    with open(args.list_of_samples, 'r') as lines:
        for line in lines:
            sample = line.strip()
            samples.append(sample)
    logging.info(f'{str(len(samples))} samples saved.')

    # Saving samples ids
    logging.info('Saving antibiotics ids.')
    antibiotics = list()
    with open(args.list_of_antibiotics, 'r') as lines:
        for line in lines:
            antibiotic = line.strip()
            antibiotics.append(antibiotic)
    logging.info(f'{str(len(antibiotics))} antibiotics saved.')

    # Extracting RGI phenotype and genotype results
    # NOTE: in RGI output files, multiple AMR genes/mutations can be found in different lines for the same antibiotic
    # Thus, the whole file must be processed first
    logging.info(f'Extracting AMR genes and mutations from RGI output files...')
    rgi_phenotypes = dict()
    rgi_genotypes = dict()
    for sample in samples:
        file = args.rgi_out_dir + sample + args.sample_out_suffix
        print(file)
        if not os.path.isfile(file):
            logging.warning(f'File {file} not found!')
            for antibiotic in antibiotics:
                pred_phen = 'NA'
                genotypes = 'NA'
                rgi_phenotypes[sample] = {}
                rgi_phenotypes[sample][antibiotic] = pred_phen
                rgi_genotypes[sample] = {}
                rgi_genotypes[sample][antibiotic] = genotypes
        if os.path.isfile(file):
            antibiotic_genotypes = dict()  # to save all observed AMR genes/mutations for the each antibiotic class
            with open(file, 'r') as lines:
                for line in lines:
                    if not line.startswith('ORF_ID'):
                        items = line.strip().split('\t')
                        model_type = items[11]
                        snps = items[12]
                        gene = items[8]
                        drug_class = items[14].strip().split('; ')
                        per_length = float(items[20])
                        if per_length >= float(args.percentage_length):
                            genotype = gene
                            if model_type == 'protein variant model':
                                genotype = gene + ' (' + snps + ')'
                            for drug in drug_class:
                                if antibiotic_genotypes.get(drug) is None:
                                    antibiotic_genotypes[drug] = genotype
                                else:
                                    antibiotic_genotypes[drug] = antibiotic_genotypes[drug] + ';' + genotype
            # for each antibiotic in input antibiotic list
            for antibiotic in antibiotics:
                pred_phen = 'S'
                genotypes = 'NA'
                if antibiotic in antibiotic_genotypes:
                    pred_phen = 'R'
                    genotypes = antibiotic_genotypes[antibiotic]
                if rgi_phenotypes.get(sample) is None:
                    rgi_phenotypes[sample] = {}
                    rgi_phenotypes[sample][antibiotic] = pred_phen
                else:
                    rgi_phenotypes[sample][antibiotic] = pred_phen
                if rgi_genotypes.get(sample) is None:
                    rgi_genotypes[sample] = {}
                    rgi_genotypes[sample][antibiotic] = genotypes
                else:
                    rgi_genotypes[sample][antibiotic] = genotypes

    # # Writing output files
    logging.info(f"Writing output files...")
    output_table = open(args.output_phenotypes_file, 'w')
    header = 'sample'
    for antibiotic in antibiotics:
        header += '\t' + antibiotic
    output_table.write(header + '\n')
    for sample in samples:
        newline = sample
        for antibiotic in antibiotics:
            if rgi_phenotypes[sample].get(antibiotic) is None:
                newline += '\t' + 'NA'
            else:
                newline += '\t' + rgi_phenotypes[sample][antibiotic]
        output_table.write(newline + '\n')

    output_table = open(args.output_genotypes_file, 'w')
    header = 'sample'
    for antibiotic in antibiotics:
        header += '\t' + antibiotic
    output_table.write(header + '\n')
    for sample in samples:
        newline = sample
        for antibiotic in antibiotics:
            if rgi_genotypes[sample].get(antibiotic) is None:
                newline += '\t' + 'NA'
            else:
                # Keeping only unique list of genotypes
                unique_genotypes = list(set(rgi_genotypes[sample][antibiotic].split(';')))
                newline += '\t' + ';'.join(unique_genotypes)
        output_table.write(newline + '\n')
    logging.info(f"Writing output files... DONE!")


if __name__ == "__main__":
    _main()
