#!/bin/bash

set -e

# Positional arguments expected

fasta_file=$1
sample_name=$2
fastq_dir=$3

fastq_file=$fastq_dir$sample_name".fastq";
fastq1_file=$fastq_dir$sample_name"_1.fastq";
fastq2_file=$fastq_dir$sample_name"_2.fastq";

echo "Running fastaq to_perfect_reads "$fasta_file" "$fastq_file" 500 50 100 150"

fastaq to_perfect_reads $fasta_file $fastq_file 500 50 100 150

echo "Running fastaq deinterleave "$fastq_file" "$fastq1_file" "$fastq2_file

fastaq deinterleave $fastq_file $fastq1_file $fastq2_file

echo "Deleting "$fastq_file

rm $fastq_file

echo "Gzipping "$fastq1_file" and "$fastq2_file

gzip $fastq1_file

gzip $fastq2_file

echo "Done"
