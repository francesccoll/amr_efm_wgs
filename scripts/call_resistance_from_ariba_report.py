#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from Bio import SeqIO
from Bio.Data import IUPACData
import gzip

# ---------------------------------------------------------------------------------------------------------------------
# Notes
# ---------------------------------------------------------------------------------------------------------------------

# This python script is used to call antibiotic susceptibility from ARIBA reports
# This script is currently under development, the following features still need to be implemented:
#   - if nucleotide indels are added to look-up table, from_garc_to_ariba_hgvs_ins need to be modified, atm expects aa
#   - Calculate number of copies for genes and call resistance based on number of genes
#   - Wildcard annotated mutations (e.g. pbp5_*! meaning any stop codon) not currently interpreted
#   - Hierarchical rules (that is when the effect of one determinant should over-rule that of another determinant)
#     has not been not implemented yet (e.g. a loss-of-function mutation in pbp5 conferring susceptibility should
#     over-rule any other pbp5 mutation conferring resistance)
#   - The detection of heterozygous calls only work for heterozygous nucleotide SNPs ATM. Detection of HET aa or indels
#     needs to be implemented
#   - Simplify how variant annotation is dealt with. At the moment, look-up table, ARIBA and SnpEff Snippy VCF each use
#     a different variant annotation
#   - To extract non-protein coding mutations (i.e. rRNA) from SnpEff-annotated Snippy VCF files

# Changes made with respect to previous version v0.0
#   - Added detection of heterozygous nucleotide SNPs
#   - Nucleotide changes from ARIBA report were converted to lower case to allow comparison with look-up table
# Changes made with respect to previous version v0.2
#   - ARIBA flags (https://github.com/sanger-pathogens/ariba/wiki/Task%3A-flag) are used to call gene presence,
#   as in some cases genes may be assembled in multiple contigs. If a gene has an 'assembled' ARIBA flag, then the
#   pc_assembled = str(args.ariba_assembled_threshold). This is to recover genes assembled into multiple contigs
#   - added an option to ignore ARIBA mutations TRUNC (truncated) and FSHIFT (frameshift) when labelling gene integrity
#   This is to test the effect of gene variants which may still be functional despite truncations or frameshifts.
# Changes made with respect to previous version v0.3
#   - added the option to read SnpEff-annotated Snippy VCF files to rescue mutations not detected by ARIBA
#   - changed expected format of input_genes file to extract: gene_id   locus_tag   path_to_fasta_file
#       locus_tag is needed to identify to identify genes in SnpEff-annotated Snippy VCF files
# added the option to open zipped vcf files

# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------


def parse_arguments():
    description = "This python script is used to call antibiotic susceptibility from ARIBA reports"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-t", "--look_up_table", action="store", dest="look_up_table",
        help="AMR 'look-up table' sheet in enterococci_amr_genes.xlsx as a tab-delimited CSV file",
        required=True, metavar="TABLE"
    )
    group.add_argument(
        "-r", "--ariba_report", action="store", dest="ariba_report",
        help="ARIBA report",
        required=True, metavar="ARIBA"
    )
    group.add_argument(
        "-g", "--input_genes", action="store", dest="input_genes",
        help="tab-delimited file containing gene_name, locus_tag and full path to gene's fasta file",
        required=True, metavar="GENES"
    )
    group.add_argument(
        "-s", "--sample_id", action="store", dest="sample_id",
        help="sample id used in saved output files",
        required=True, metavar="SAMPLE"
    )
    group.add_argument(
        "-p", "--output_prefix", action="store", dest="output_prefix",
        help="output prefix used to name output files",
        required=True, metavar="PREFIX"
    )
    group = parser.add_argument_group('Optional arguments: gene/mutation match')
    group.add_argument(
        "-v", "--vcf_file", action="store", dest="vcf_file",
        help="SnpEff-annotated VCF file to be used to scan for mutations at candidate genes "
             "(i.e. those specified in --input_genes)."
             "Recommended option, as ARIBA may miss detection of mnp and complex SNPs.",
        required=False, metavar="VCF"
    )
    group.add_argument(
        "-i", "--gene_pc_ident", action="store", dest="gene_pc_ident",
        help="threshold for percentage nucleotide identity of genes to call resistance (in percentage)",
        required=False, metavar="GENE_ID", default=90
    )
    group.add_argument(
        "-l", "--gene_pc_assembled", action="store", dest="gene_pc_assembled",
        help="threshold for minimum length of gene assembled to call resistance (in percentage)",
        required=False, metavar="GENE_LENGTH", default=60
    )
    group.add_argument(
        "-f", "--use_ariba_flag", action="store_true", dest="use_ariba_flag",
        help="whether to use ARIBA flag to determine if gene was assembled. See option --ariba_assembled_threshold",
        required=False
    )
    group.add_argument(
        "-a", "--ariba_assembled_threshold", action="store", dest="ariba_assembled_threshold",
        help="The percentage the reference gene sequence used by ARIBA to call a gene assembled."
             "This value is specified with ARIBA's option --assembled_threshold, which is 95 per cent by default.",
        required=False, metavar="ARIBA_AT", default=95
    )
    group.add_argument(
        "-x", "--ignore_lof_mutations", action="store_false", dest="ignore_lof_mutations",
        help="Ignore loss-of-function mutations (i.e. those labelled as TRUNC and FSHIFT by ARIBA) when determining "
             "gene integrity. By default a gene containing LOF mutations is not considered functional and will not be "
             "used to call resistance. By selecting this option, LOF mutations are not considered.",
        required=False, default=True
    )
    group.add_argument(
        "-m", "--minimum_ratio_het", action="store", dest="minimum_ratio_het",
        help="Minimum ratio of reads supporting alternative allele. "
             "Used to call resistance from heterozygous mutations (default: 0.10)",
        required=False, metavar="MIN_RATIO", default=0.10
    )

    return parser.parse_args()


def to_ariba_ref_name(gene_name):
    ref_name = str(gene_name)
    ref_name = ref_name.replace('\'', '_')
    ref_name = ref_name.replace('(', '_')
    ref_name = ref_name.replace(')', '_')
    ref_name = ref_name.replace('-', '_')
    ref_name = ref_name.replace('/', '_')
    return ref_name


def from_garc_to_ariba_hgvs_del(deletion):
    # In the look-up table, indels annotation is adapted from GARC rules to include both reference and alternative allele.
    # E.g. rpsJ@I51_del_IRATH. At aa position 51 the reference aa is I and aa 52-56 IRATH are deleted.
    # However, this differs from the notation used by ARIBA (assumed to be HGVS notation)
    # rpsJ@I52_del_RA --> R53_A54del
    # HDprotein@H35_del_S --> S36del
    (gene_name, change) = deletion.split('@')
    # extract reference position
    digits = [int(s) for s in list(change) if s.isdigit()]
    ref_pos = ''
    for digit in digits:
        ref_pos += str(digit)
    from_pos = int(ref_pos) + 1
    # change format
    new_change = 'na'
    new_deletion = 'na'
    parts = change.split('_')
    if len(parts) == 3:
        deleted_seq = list(parts[2])
        if len(deleted_seq) > 1:
            to_pos = from_pos + len(deleted_seq) - 1
            new_change = str(deleted_seq[0]) + str(from_pos) + '_' + deleted_seq[-1] + str(to_pos) + 'del'
        if len(deleted_seq) == 1:
            new_change = str(deleted_seq[0]) + str(from_pos) + 'del'
        new_deletion = gene_name + '@' + new_change
    return new_deletion


def from_garc_to_ariba_hgvs_ins(insertion, fasta_file):
    # In the look-up table, indels annotation is adapted from GARC rules to include both reference and alternative allele.
    # E.g. cls@L110_ins_MPL. At aa position 110 the reference aa is L and aa MPL are inserted after that.
    # However, this differs from the notation used by ARIBA (assumed to be HGVS notation)
    # cls@L110_ins_MPL --> L110_?111insMPL, where ? is the amino acid at position 111
    # pbp5@S466_ins_S --> S466_?467insS, where ? is the amino acid at position 467
    (gene_name, change) = insertion.split('@')
    # getting gene amino acid sequence
    first_record = next(SeqIO.parse(fasta_file, "fasta"))
    gene_BioSeq_nt = first_record.seq
    gene_BioSeq_aa = gene_BioSeq_nt.translate()
    # extract reference position
    digits = [int(s) for s in list(change) if s.isdigit()]
    ref_pos = ''
    for digit in digits:
        ref_pos += str(digit)
    # change format
    new_change = 'na'
    new_insertion = 'na'
    parts = change.split('_')
    ref_aa = str(list(parts[0])[0])
    next_aa = gene_BioSeq_aa[int(ref_pos) - 1 + 1]
    if len(parts) == 3:
        inserted_seq = str(parts[2])
        new_change = ref_aa + str(ref_pos) + '_' + next_aa + str(int(ref_pos)+1) + 'ins' + inserted_seq
        new_insertion = gene_name + '@' + new_change
    return new_insertion


def is_snp(mutation, variant_type):
    # This function is used to determine whether a SnpEff annotated variant is a SNP
    # NOTE: variant_type is extracted from SnpEff ANN field TYPE
    # NOTE: the are five types of variants annotated by SnpEff: snp, mnp, ins, del and complex
    # NOTE: in some occasions, complex variants may contain indels
    is_snp = False
    snp_types = ['snp', 'mnp', 'complex']
    indel_types = ['ins', 'del']
    if variant_type in indel_types:
        is_snp = False
    if variant_type in snp_types:
        is_snp = True
    if variant_type == "complex":
        if 'del' in mutation:
            is_snp = False
        if 'ins' in mutation:
            is_snp = False
        if 'dup' in mutation:
            is_snp = False
    return is_snp


# Adding new codings absent in IUPACData
# NOTE: these codings are added to deal with odd symbols in SnpEff-annotated protein changes
IUPACData.protein_letters_3to1["ext"] = ""
IUPACData.protein_letters_3to1["*?"] = ""
IUPACData.protein_letters_3to1["?"] = "?"
#   to encode for stop codons
IUPACData.protein_letters_3to1["Ter"] = "*"
IUPACData.protein_letters_3to1["*"] = "*"
#   or other hgvs symbols
IUPACData.protein_letters_3to1["dup"] = "dup"
IUPACData.protein_letters_3to1["fs"] = "fs"
IUPACData.protein_letters_3to1["del"] = "del"
IUPACData.protein_letters_3to1["ins"] = "ins"


def amino_acids_3to1(letters):
    n = 3
    aa3 = [letters[i:i + n] for i in range(0, len(letters), n)]
    aa1 = [IUPACData.protein_letters_3to1[aa3[i]] for i in range(0, len(aa3), 1) ]
    return ''.join(aa1)


def from_snpeff_aa_to_ariba_hgvs_snp(mutation):
    # SnpEff-annotated VCFs encode amino acid changes differently from ARIBA
    # AA changes like p.Thr19Ile or p.ValLeuThr66AlaLeuIle, need to be replaced to one-letter AA
    # extract reference position
    print(mutation)
    mutation = mutation.replace("p.", "")
    new_mutation = mutation  # if mutation annotation cannot be converted, original mutation will be returned
    digits = [int(s) for s in list(mutation) if s.isdigit()]
    ref_pos = ''
    for digit in digits:
        ref_pos += str(digit)
    print('ref_pos ' + str(ref_pos))
    parts = mutation.split(str(ref_pos))
    ref_aa = str(parts[0])
    alt_aa = str(parts[1])
    print("ref_aa ref_pos alt_aa " + ref_aa + ' ' + ref_pos + ' ' + alt_aa)
    if len(list(ref_aa)) > 0:
        ref_aa = amino_acids_3to1(ref_aa)
        alt_aa = amino_acids_3to1(alt_aa)
        print("ref_aa ref_pos alt_aa " + ref_aa + ' ' + ref_pos + ' ' + alt_aa)
        new_mutation = ref_aa + ref_pos + alt_aa
    return new_mutation


def from_snpeff_aa_to_ariba_hgvs_ins(mutation):
    # SnpEff-annotated VCFs encode amino acid changes differently from ARIBA
    # AA insertions such as p.Met119_Tyr120insSer or p.Ser308_Thr309insAla, need to be changed to one-letter HGVS format
    # This function is specifically design to change the annotation of conservative_inframe_insertion
    # A major group of conservative_inframe_insertion are annotated as dup
    # Frameshift and dup insertions can be annotated with from_snpeff_aa_to_ariba_hgvs_snp
    # see https://varnomen.hgvs.org/recommendations/protein/variant/insertion/
    # Format: “prefix”“amino_acids+positions_flanking”“ins”“inserted_sequence”, e.g. p.(Lys23_Leu24insArgSerGln)
    print(mutation)
    mutation = mutation.replace("p.", "")
    parts = mutation.split("_")
    new_mutation = mutation
    if len(parts) == 1:
        # e.g. p.Thr669delinsIleAla
        digits = [int(s) for s in list(mutation) if s.isdigit()]
        ref_pos = ''
        for digit in digits:
            ref_pos += str(digit)
        ref_aa = mutation.split(str(ref_pos))[0]
        ref_aa = amino_acids_3to1(ref_aa)
        ins_seq = mutation.split(str(ref_pos))[1]
        ins_seq = amino_acids_3to1(ins_seq)
        new_mutation = ref_aa + str(ref_pos) + ins_seq
    if len(parts) == 2:
        (left, right) = mutation.split("_")
        # extracting flanking positions
        digits = [int(s) for s in list(left) if s.isdigit()]
        ref_pos1 = ''
        for digit in digits:
            ref_pos1 += str(digit)
        digits = [int(s) for s in list(right) if s.isdigit()]
        ref_pos2 = ''
        for digit in digits:
            ref_pos2 += str(digit)
        # extracting flanking amino acids
        ref_aa1 = left.split(str(ref_pos1))[0]
        ref_aa2 = right.split(str(ref_pos2))[0]
        # extracting inserted sequence
        ins_seq = right.split(str(ref_pos2))[1]
        ins_seq = ins_seq.replace("ins", "")
        print("ref_aa1 ref_aa2 ins_seq " + ref_aa1 + ' ' + ref_aa2 + ' ' + ins_seq)
        # concerting 3 to 1 amino acid letters
        ref_aa1 = amino_acids_3to1(ref_aa1)
        ref_aa2 = amino_acids_3to1(ref_aa2)
        ins_seq = amino_acids_3to1(ins_seq)
        # re-creating new mutation
        new_mutation = ref_aa1 + str(ref_pos1) + '_' + ref_aa2 + str(ref_pos2) + 'ins' + ins_seq
    print(new_mutation)
    return new_mutation


def from_snpeff_aa_to_ariba_hgvs_del(mutation):
    # SnpEff-annotated VCFs encode amino acid changes differently from ARIBA
    # AA deletions such as p.Tyr370_Phe371del, p.His144del or p.Val19_Ser23delinsGly, need to be changed to one-letter
    # HGVS format
    # This function is specifically design to change the annotation of conservative_inframe_deletion
    # Frameshift deletions can be annotated with from_snpeff_aa_to_ariba_hgvs_snp
    # see https://varnomen.hgvs.org/recommendations/protein/variant/deletion/
    # Format: “prefix”“amino_acid(s)+position(s)_deleted”“del”, e.g. p.(Cys76_Glu79del)
    mutation = mutation.replace("p.", "")
    parts = mutation.split("_")
    new_mutation = mutation
    if len(parts) == 1:
        # e.g. p.His144del
        digits = [int(s) for s in list(mutation) if s.isdigit()]
        ref_pos = ''
        for digit in digits:
            ref_pos += str(digit)
        ref_aa = mutation.split(str(ref_pos))[0]
        ref_aa = amino_acids_3to1(ref_aa)
        new_mutation = ref_aa + str(ref_pos) + 'del'
    if len(parts) == 2:
        # e.g. p.Tyr370_Phe371del
        (left, right) = mutation.split("_")
        # extracting flanking positions
        digits = [int(s) for s in list(left) if s.isdigit()]
        ref_pos1 = ''
        for digit in digits:
            ref_pos1 += str(digit)
        digits = [int(s) for s in list(right) if s.isdigit()]
        ref_pos2 = ''
        for digit in digits:
            ref_pos2 += str(digit)
        # extracting flanking amino acids
        ref_aa1 = left.split(str(ref_pos1))[0]
        ref_aa2 = right.split(str(ref_pos2))[0]
        # extracting inserted sequence
        # NOTE: in most cases the deleted sequence will be just 'del'
        # that's why IUPACData.protein_letters_3to1["del"] = "del" was added
        del_seq = right.split(str(ref_pos2))[1]
        # concerting 3 to 1 amino acid letters
        ref_aa1 = amino_acids_3to1(ref_aa1)
        ref_aa2 = amino_acids_3to1(ref_aa2)
        del_seq = amino_acids_3to1(del_seq)
        # re-creating new mutation
        new_mutation = ref_aa1 + str(ref_pos1) + '_' + ref_aa2 + str(ref_pos2) + del_seq
    return new_mutation


def from_snpeff_aa_to_ariba_hgvs_dup(mutation):
    # SnpEff-annotated VCFs encode amino acid changes differently from ARIBA
    # AA duplications can take two forms: p.Ala3dup (one amino acid) or p.Ala3_Ser5dup (several amino acids)
    # This function is specifically design to change the annotation of AA duplications
    # see https://varnomen.hgvs.org/recommendations/protein/variant/duplication/
    print(mutation)
    mutation = mutation.replace("p.", "")
    parts = mutation.split("_")
    new_mutation = mutation
    if len(parts) == 1:
        # e.g. p.Ala3dup
        new_mutation = from_snpeff_aa_to_ariba_hgvs_snp(mutation)
    if len(parts) == 2:
        # e.g. p.Ala3_Ser5dup
        new_mutation = from_snpeff_aa_to_ariba_hgvs_del(mutation)
    print(new_mutation)
    return new_mutation


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
    logging.info(f'Making sure input files exist...')
    if not os.path.isfile(args.look_up_table):
        logging.error(f'Input look_up_table file {args.look_up_table} not found!')
        sys.exit(-1)
    if not os.path.isfile(args.ariba_report):
        logging.error(f'Input ariba_report file {args.ariba_report} not found!')
        sys.exit(-1)
    if not os.path.isfile(args.input_genes):
        logging.error(f'Input genes file {args.input_genes} not found!')
        sys.exit(-1)
    logging.info(f'Making sure input files exist. DONE.')

    # Saving the full path of genes
    gene_fasta_files = dict()
    gene_locus_tags = dict()  # dict needed to detect genes in VCF file
    with open(args.input_genes, 'r') as input_genes:
        for line in input_genes:
            (gene_name, locus_tag, sequence_file) = line.strip().split('\t')
            gene_fasta_files[gene_name] = sequence_file
            gene_locus_tags[locus_tag] = gene_name
            print("gene_fasta_files[" + gene_name + "] --> " + gene_fasta_files[gene_name])

    # Variables to store ARIBA report gene and mutation information
    gene_ariba_ids = dict()  # dict to save all genes detected by ARIBA
    gene_ariba_pc_ident = dict()  # dict to save percentage of nucleotide identity of genes
    gene_ariba_pc_assembled = dict()  # dict to save percentage of gene reference sequence assembled by ARIBA
    gene_ariba_all_var = dict()  # dict to save ALL variants found by ARIBA on gene sequence
    gene_ariba_all_var_effect = dict()  # dict to save ALL variants found by ARIBA on gene sequence
    gene_ariba_integrity = dict()  # dict to save the integrity of genes
    gene_ariba_copies = dict()  # dict to save the number of gene copies
    # Integrity labels: fully_assembled, partially_assembled, truncated and frameshift (or combinations of these)
    mutations_known = dict()  # dict to save known ARIBA mutations
    mutations_ariba_all = dict()  # dict to save ALL ARIBA mutations
    mutations_vcf = dict()  # dict to save mutations from VCF at candidate genes

    # Saving mutations at candidate genes from SnpEff annotated VCF file if option --vcf_file chosen
    if args.vcf_file:
        vcf_file_file = ""
        if '.gz' in args.vcf_file:
            vcf_file_file = gzip.open(args.vcf_file, 'rt')
        else:
            vcf_file_file = open(args.vcf_file, 'r')
        # with gzip.open(args.vcf_file, 'rt') as input_vcf_lines:
        with vcf_file_file as input_vcf_lines:
            input_vcf_lines.readline()
            for vcf_line in input_vcf_lines:
                if not vcf_line.startswith("#"):
                    variant_type = ''
                    items = vcf_line.split('\t')[7].split(';')
                    for item in items:
                        if item.startswith("TYPE"):
                            variant_type = item.replace("TYPE=", "")
                        if item.startswith("ANN"):
                            ann_items = item.split('|')
                            coding = ann_items[5]
                            mutation = ann_items[10]
                            locus_tag = ann_items[3]
                            # coding categories considered: transcript, intergenic_region and gene_variant
                            if coding != 'transcript':
                                mutation = ann_items[9]
                            # NOTE: mutations annotated at the DNA level (n.) are not changed
                            # dealing with AA changes: changing annotation to match 3-to-1 AA ARIBA annotation
                            if mutation.startswith("p."):
                                # print(mutation)
                                if is_snp(mutation, variant_type):
                                    mutation = from_snpeff_aa_to_ariba_hgvs_snp(mutation)
                                if variant_type == "ins":
                                    if mutation.endswith('fs'):
                                        mutation = from_snpeff_aa_to_ariba_hgvs_snp(mutation)
                                    elif mutation.endswith('dup'):
                                        mutation = from_snpeff_aa_to_ariba_hgvs_dup(mutation)
                                    else:
                                        mutation = from_snpeff_aa_to_ariba_hgvs_ins(mutation)
                                if variant_type == "del":
                                    if mutation.endswith('fs'):
                                        mutation = from_snpeff_aa_to_ariba_hgvs_snp(mutation)
                                    else:
                                        mutation = from_snpeff_aa_to_ariba_hgvs_del(mutation)
                            # Saving mutation if in candidate gene
                            if locus_tag in gene_locus_tags:
                                gene_name = gene_locus_tags[locus_tag]
                                mutation_id = gene_name + '@' + mutation
                                mutations_vcf[mutation_id] = 1

    # Parsing ARIBA report
    logging.info(f'Parsing ARIBA report {args.ariba_report}...')
    with open(args.ariba_report, 'r') as input_lines:
        input_lines.readline()
        for line in input_lines:
            # Expected fields:
            # #ariba_ref_name	ref_name	gene	var_only	flag	reads	cluster	ref_len	ref_base_assembled	pc_ident
            # ctg	ctg_len	ctg_cov	known_var	var_type	var_seq_type	known_var_change	has_known_var	ref_ctg_change
            # ref_ctg_effect	ref_start	ref_end	ref_nt	ctg_start	ctg_end	ctg_nt	smtls_total_depth	smtls_nts
            # smtls_nts_depth	var_description	free_text
            (_, ref_name, gene, var_only, flag, _, _, ref_len, ref_base_asse, pc_ident, *rest) = line.strip().split('\t')
            (_, _, _, known_var, var_type, var_seq_type, known_var_cg, has_known_var, ref_change, ref_effect, *rest2) = rest
            # Saving gene information
            gene_name = ref_name
            pc_assembled = str(float(ref_base_asse)/float(ref_len) * 100)
            # If ARIBA flag is used, then pc_assembled is known to be at least ariba_assembled_threshold
            if args.use_ariba_flag:
                flags = flag.split(',')
                if 'assembled' in flags:
                    pc_assembled = str(args.ariba_assembled_threshold)
            if gene_name not in gene_ariba_ids:
                gene_ariba_ids[gene_name] = 1
                gene_ariba_pc_ident[gene_name] = pc_ident
                gene_ariba_pc_assembled[gene_name] = pc_assembled
                gene_ariba_all_var[gene_name] = ref_change
                gene_ariba_all_var_effect[gene_name] = ref_effect
            else:
                # if percentage of assembled gene is greater than that saved in a previous ARIBA report line,
                # then the highest value of pc_assembled is saved, along with its corresponding pc_ident
                if float(pc_assembled) > float(gene_ariba_pc_assembled[gene_name]):
                    gene_ariba_pc_assembled[gene_name] = pc_assembled
                    gene_ariba_pc_ident[gene_name] = pc_ident
                if ref_change != '.':
                    gene_ariba_all_var[gene_name] = gene_ariba_all_var[gene_name] + ';' + ref_change
                    gene_ariba_all_var_effect[gene_name] = gene_ariba_all_var_effect[gene_name] + ';' + ref_effect
            # Saving mutation information
            mutation_id = 'na'
            # NOTE: nucleotide changes need to be changed to lowercase to be matched with look-up table
            if var_seq_type == 'n':
                known_var_cg = known_var_cg.lower()
            # NOTE: ARIBA report field ref_change is used to check if line corresponds to a genetic change
            if ref_change != '.':
                mutation_id = gene_name + '@' + ref_change
            if known_var == '1':
                # NOTE: known_var_cg saved instead of ref_change to save known mutations where alternative allele
                # matches gene reference amino acid
                mutation_id = gene_name + '@' + known_var_cg
                if has_known_var == '1':
                    mutations_known[mutation_id] = 1
                else:
                    # if known variant not called, check for heterozygous call
                    # NOTE: the code below only works to detect heterozygous nucleotide SNPs ATM, not aa or indels
                    (_,	_, ref_nt, _, _, ctg_nt, smtls_total_depth, smtls_nts, smtls_nts_depth, *_) = rest2
                    known_var_cg_alt = list(known_var_cg)[-1]
                    known_var_cg_alt = known_var_cg_alt.upper()
                    if len(ctg_nt) == 1:
                        if ',' in smtls_nts:
                            print('Heterozygous call detected at ' + gene_name + ' ' + known_var_cg + ': ' + smtls_nts)
                            het_alleles = smtls_nts.split(',')
                            het_alleles_depths = smtls_nts_depth.split(',')
                            for het_allele, het_allele_depth in zip(het_alleles, het_alleles_depths):
                                if het_allele == known_var_cg_alt:
                                    ratio = float(int(het_allele_depth)/int(smtls_total_depth))
                                    print('\tAlternative allele ' + het_allele + ' found at ' + str(ratio) + ' ratio.')
                                    if ratio > float(args.minimum_ratio_het):
                                        mutations_known[mutation_id] = ratio
            if mutation_id != 'na':
                mutations_ariba_all[mutation_id] = 1
    logging.info(f'Parsing ARIBA report {args.ariba_report}. DONE.')

    # Determine ARIBA gene integrity
    for gene_name in gene_ariba_ids:
        integrity_label = ''
        pc_assembled = gene_ariba_pc_assembled[gene_name]
        if float(pc_assembled) > float(args.gene_pc_assembled):
            integrity_label += 'fully_assembled'
        else:
            integrity_label = 'partially_assembled'
        if args.ignore_lof_mutations:
            if 'TRUNC' in gene_ariba_all_var_effect[gene_name]:
                integrity_label += '+truncated'
            if 'FSHIFT' in gene_ariba_all_var_effect[gene_name]:
                integrity_label += '+frameshift'
        gene_ariba_integrity[gene_name] = integrity_label
        print('gene_ariba_integrity[' + gene_name + '] --> ' + gene_ariba_integrity[gene_name])

    # Determine ARIBA gene copies --> to do

    # Printing all genes and mutations saved
    logging.info(f'Printing all extracted genes from ARIBA report.')
    for gene_name in gene_ariba_ids:
        print('\tgene_name ' + gene_name + ' pc_ident ' + gene_ariba_pc_ident[gene_name] + ' pc_assembled ' +
              gene_ariba_pc_assembled[gene_name] + ' all gene variants ' + gene_ariba_all_var[gene_name])
    logging.info(f'Printing all extracted known mutations from ARIBA report.')
    for mutation_id in mutations_known:
        print('\tmutation_id ' + mutation_id)

    # Variables to store look-up table gene and mutation information
    lu_table_abs = dict()  # dict to save all unique antibiotics from look-up table
    abs_res_genotype = dict()  # dict to save genotype ids responsible for predicted phenotypes
    abs_res_genotype_all = dict()  # dict to save ALL genotypes, including those not meeting all inclusion criteria
    abs_res_phenotype = dict()  # dict to save predicted antibiotic susceptibility phenotypes
    # the list ignore_effects is used to know what genetic determinants to ignore, based on their effect
    ignore_effects = ['none', 'undetermined', 'undetermined_not_annotated']  # list of effects to ignore
    # the two lists below are used to indicate what broad types of genetic determinants exist
    gene_types = ['acquired gene', 'acquired multiple genes', 'chromosomal gene']
    mutation_types = ['single mutation', 'mutation combination']

    # Parsing look-up table
    logging.info(f'Parsing look-up table {args.look_up_table}')
    with open(args.look_up_table, 'r') as input_lines:
        input_lines.readline()
        for line in input_lines:
            # mechanism	class	antibiotic	effect	id	type	mutation	gene_name
            (_, _, antibiotic, effect, id, type, mut, gn, *_) = line.strip().split('\t')
            lu_table_abs[antibiotic.strip()] = 1
            if effect not in ignore_effects:
                # check presence and integrity of genes
                if type in gene_types:
                    gene_names = gn.split(';')
                    gene_names_found = list()  # to save ganes meeting all inclusion criteria
                    gene_names_not_found = list()  # to save genes not meeting all inclusion criteria
                    for gene_name in gene_names:
                        gene_name_not_found = gene_name + ':'  # label for genes not meeting all inclusion criteria
                        gene_name_ariba = to_ariba_ref_name(gene_name)
                        if gene_name_ariba in gene_ariba_ids:
                            pc_ident = gene_ariba_pc_ident[gene_name_ariba]
                            if float(pc_ident) > float(args.gene_pc_ident):
                                if gene_ariba_integrity[gene_name_ariba] == 'fully_assembled':
                                    gene_names_found.append(gene_name)
                                else:
                                    gene_name_not_found += gene_ariba_integrity[gene_name_ariba]
                            else:
                                gene_name_not_found += 'low_sequence_identity'
                                gene_name_not_found += '+' + gene_ariba_integrity[gene_name_ariba]
                        else:
                            gene_name_not_found += 'not_found'
                        print(gene_name_not_found)
                        gene_names_not_found.append(gene_name_not_found)
                    # if ALL genes in genetic determinant are found by ARIBA > save effect
                    if gn == ';'.join(gene_names_found):
                        print('\tgene(s) ' + gn + ' (all) found: ' + effect + ' to ' + antibiotic + ' called')
                        if antibiotic in abs_res_phenotype:
                            abs_res_phenotype[antibiotic] = abs_res_phenotype[antibiotic] + ';' + effect
                        else:
                            abs_res_phenotype[antibiotic] = effect
                        if antibiotic in abs_res_genotype:
                            abs_res_genotype[antibiotic] = abs_res_genotype[antibiotic] + ';' + id
                        else:
                            abs_res_genotype[antibiotic] = id
                        if antibiotic in abs_res_genotype_all:
                            abs_res_genotype_all[antibiotic] = abs_res_genotype_all[antibiotic] + ';' + gn
                        else:
                            abs_res_genotype_all[antibiotic] = gn
                    else:
                        # if not > save why
                        ann = ';'.join(gene_names_not_found)
                        print('genetic determinant NOT found by ARIBA: ' + ann)
                        if antibiotic in abs_res_genotype_all:
                            abs_res_genotype_all[antibiotic] = abs_res_genotype_all[antibiotic] + ';' + ann
                        else:
                            abs_res_genotype_all[antibiotic] = ann

                # check presence of mutations
                if type in mutation_types:
                    mutations = mut.split('+')
                    mutations_found = list()  # to save mutations
                    for mutation in mutations:
                        # change mutation format for indels
                        if 'del' in mutation:
                            mutation = from_garc_to_ariba_hgvs_del(mutation)
                        if 'ins' in mutation:
                            (gene_name, *_) = mutation.split('@')
                            if gene_name not in gene_fasta_files:
                                logging.error(f'{gene_name} from {mutation} not in {args.input_genes}!')
                                sys.exit(-1)
                            mutation = from_garc_to_ariba_hgvs_ins(mutation, gene_fasta_files[gene_name])
                        # if mutation found in ARIBA report
                        if mutation in mutations_known:
                            mutations_found.append(mutation)
                        # If VCF is selected, then check if mutation found in VCF
                        if args.vcf_file:
                            if mutation in mutations_vcf:
                                if mutation not in mutations_found:
                                    mutations_found.append(mutation)
                    if mut == '+'.join(mutations_found):
                        print('\tmutation(s) ' + mut + ' (all) found: ' + effect + ' to ' + antibiotic + ' called')
                        # NOTE: because multiple genetic determinants may be found, effects will need to be resolved
                        if antibiotic in abs_res_phenotype:
                            abs_res_phenotype[antibiotic] = abs_res_phenotype[antibiotic] + ';' + effect
                        else:
                            abs_res_phenotype[antibiotic] = effect
                        if antibiotic in abs_res_genotype:
                            abs_res_genotype[antibiotic] = abs_res_genotype[antibiotic] + ';' + id
                        else:
                            abs_res_genotype[antibiotic] = id
                        if antibiotic in abs_res_genotype_all:
                            abs_res_genotype_all[antibiotic] = abs_res_genotype_all[antibiotic] + ';' + mut
                        else:
                            abs_res_genotype_all[antibiotic] = mut

    # Resolving effects from multiple genetic determinants.
    # NOTE: this step is important when hierarchical rules are included.
    # NOTE: ATM only removes redundant calls, e.g. when multiple 'resistance' effects were called
    for antibiotic in abs_res_phenotype:
        # get unique list of effects
        abs_res_phenotype[antibiotic] = ';'.join(list(set(abs_res_phenotype[antibiotic].split(';'))))
    for antibiotic in abs_res_genotype:
        # get unique list of genetic determinants
        abs_res_genotype[antibiotic] = ';'.join(list(set(abs_res_genotype[antibiotic].split(';'))))

    # Printing look-up table
    logging.info(f'All antibiotics susceptibilities predicted from ARIBA report and look-up table')
    logging.info(f'{str(len(lu_table_abs))} antibiotics extracted from look-up table')
    header = 'sample'
    newlineG1 = args.sample_id
    newlineG2 = args.sample_id
    newlinePh = args.sample_id
    output_genotype_file1 = open(args.output_prefix + '.geno.txt', 'w')
    output_genotype_file2 = open(args.output_prefix + '.geno_all.txt', 'w')
    output_phenotype_file = open(args.output_prefix + '.pheno.txt', 'w')
    for antibiotic in lu_table_abs:
        print(antibiotic)
        if antibiotic != 'unknown':
            header += '\t' + antibiotic
            susceptibility = 'susceptible'
            gen_determinants1 = 'not_found'
            gen_determinants2 = 'not_found'
            if antibiotic in abs_res_phenotype:
                susceptibility = abs_res_phenotype[antibiotic]
            if antibiotic in abs_res_genotype:
                gen_determinants1 = abs_res_genotype[antibiotic]
            if antibiotic in abs_res_genotype_all:
                gen_determinants2 = abs_res_genotype_all[antibiotic]
            newlineG1 += '\t' + gen_determinants1
            newlineG2 += '\t' + gen_determinants2
            newlinePh += '\t' + susceptibility
            print('\tantibiotic ' + antibiotic + ' susceptibility ' + susceptibility +
                  ' gen_determinants ' + gen_determinants1)
    logging.info(f'Saving output files {output_phenotype_file}, {output_genotype_file1} and {output_genotype_file2}')
    output_genotype_file1.write(header + '\n' + newlineG1 + '\n')
    output_genotype_file2.write(header + '\n' + newlineG2 + '\n')
    output_phenotype_file.write(header + '\n' + newlinePh + '\n')
    output_genotype_file1.close()
    output_genotype_file2.close()
    output_phenotype_file.close()


if __name__ == "__main__":
    _main()
