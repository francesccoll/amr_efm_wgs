#!/usr/bin/env python3

import argparse
import logging
import os
import sys

# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "This script is used to collect WGS-predicted AMR phenotypes and genotypes generated by ResFinder. " \
                  "Tested with ResFinder 4.0.\n " \
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
        required=True, metavar="LIST_ABS"
    )
    group.add_argument(
        "-d", "--resfinder_out_dir", action="store", dest="resfinder_out_dir",
        help="Directory where ResFinder results are stored",
        required=True, metavar="RESFINDER_DIR"
    )
    group.add_argument(
        "-s", "--sample_out_dir_suffix", action="store", dest="sample_out_dir_suffix",
        help="Directory suffix used to name sample's ResFinder output directory (default '_resfinder') ",
        required=False, metavar="OUT_PREFIX", default='_resfinder'
    )
    group.add_argument(
        "-p", "--output_phenotypes_file", action="store", dest="output_phenotypes_file",
        help="Output file with ResFinder WGS-predicted AMR phenotypes",
        required=True, metavar="PHENO_FILE"
    )
    group.add_argument(
        "-g", "--output_genotypes_file", action="store", dest="output_genotypes_file",
        help="Output file with ResFinder WGS-detected AMR genotypes",
        required=True, metavar="GENO_FILE"
    )
    group.add_argument(
        "-o", "--organism", action="store", dest="organism",
        help="Species selected in ResFinder runs (e.g. enterococcus_faecium)",
        required=False, metavar="ORGANISM"
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
    check_path_exists(args.resfinder_out_dir)
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

    # Extracting ResFinder phenotype and genotype results
    resfinder_phenotypes = dict()
    resfinder_genotypes = dict()
    for sample in samples:
        file = args.resfinder_out_dir + sample + args.sample_out_dir_suffix + '/pheno_table.txt'
        if args.organism is not None:
            file = args.resfinder_out_dir + sample + args.sample_out_dir_suffix + '/pheno_table_' + args.organism + '.txt'
        print(file)
        if not os.path.isfile(file):
            logging.warning(f'File {file} not found!')
        else:
            with open(file, 'r') as lines:
                for line in lines:
                    for antibiotic in antibiotics:
                        pattern = antibiotic + '\t'
                        if line.startswith(pattern):
                            items = line.strip().split('\t')
                            ab, ab_class, pred_phen, match, genotype = 'Not_found','Not_found','Not_found','Not_found','Not_found'
                            if len(items) == 5:
                                ab, ab_class, pred_phen, match, genotype = items
                            if len(items) == 4:
                                ab, ab_class, pred_phen, match = items
                            if resfinder_phenotypes.get(sample) is None:
                                resfinder_phenotypes[sample] = {}
                                resfinder_phenotypes[sample][antibiotic] = pred_phen
                            else:
                                resfinder_phenotypes[sample][antibiotic] = pred_phen
                            if resfinder_genotypes.get(sample) is None:
                                resfinder_genotypes[sample] = {}
                                resfinder_genotypes[sample][antibiotic] = genotype
                            else:
                                resfinder_genotypes[sample][antibiotic] = genotype

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
            prediction = 'NA'
            if sample in resfinder_phenotypes:
                if antibiotic in resfinder_phenotypes[sample]:
                    prediction = resfinder_phenotypes[sample][antibiotic]
            newline += '\t' + prediction
        output_table.write(newline + '\n')

    output_table = open(args.output_genotypes_file, 'w')
    header = 'sample'
    for antibiotic in antibiotics:
        header += '\t' + antibiotic
    output_table.write(header + '\n')
    for sample in samples:
        newline = sample
        for antibiotic in antibiotics:
            prediction = 'NA'
            if sample in resfinder_genotypes:
                if antibiotic in resfinder_genotypes[sample]:
                    prediction = resfinder_genotypes[sample][antibiotic]
            newline += '\t' + prediction
        output_table.write(newline + '\n')
    logging.info(f"Writing output files... DONE!")


if __name__ == "__main__":
    _main()
