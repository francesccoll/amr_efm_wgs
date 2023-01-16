# Antimicrobial resistance (amr) detection from *Enterococcus faecium* (efm) whole-genome sequences (wgs)
ARIBA-based pipeline for genotypic prediction of AMR from whole-genome sequences of *Enterococcus faecium*

This GitHub project contains the data and code necessary to reproduce the findings of the study 'Improved accuracy of antibiotic resistance determination from *Enterococcus faecium* whole-genome sequences', and includes the following directories:
* amr_database: directory with database (catalogue) of AMR genetic determinants for *E. faecium* and scripts needed to transform it into an ARIBA-compliant format.
* ariba_amr_database: directory with input files and scripts to run ARIBA and call resistance from *E. faecium* genomes based on our curated catalogue of AMR genetic determinants for *E. faecium*.
* amrfinder: directory with scripts and files needed to run AMRFinderPlus from *E. faecium* assemblies, and parse AMRFinderPlus results.
* card_rgi: directory with input files and scripts to run CARD RGI from *E. faecium* genomes, and parse CARD RGI results.
* lre_finder: directory with files and scripts needed to run LRE-Finder on paired-end fastq files, and parse LRE-Finder reports.
* resfinder: directory with scripts and files needed to run ResFinder (and PointFinder) from *E. faecium* genomes (assemblies and paired-end fastq.gz files) and parse ResFinder reports.
* mapping: directory with the bash scripts used to run Snippy on paired-end fastq files and assemblies.
* phenotype_genotype: directory with the R scripts used to compare phenotypic AST results with WGS-derived genotypic predictions by various tool (described above), and our ARIBA-based predictions (based on our curated catalogue of AMR genetic determinants for *E. faecium*).
* docker: directory with the Dockerfile needed to create a Docker image to run our ARIBA-based pipeline (under construction).

# Required dependencies

## ARIBA-based pipeline

* [ARIBA](https://github.com/sanger-pathogens/ariba) version >= 2.14.6
* [ARIBA](https://github.com/sanger-pathogens/ariba#required-dependencies) dependencies
* [Snippy](https://github.com/tseemann/snippy) version >= 4.6.0
* [fastaq](https://github.com/sanger-pathogens/Fastaq) version >= 3.17.0

## Other tools used
* [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) version 3.10.18; amrfinder database version: 2021-12-21.1
* [CARD RGI](https://github.com/arpcard/rgi) rgi version v5.2.1; CARD database v3.1.4
* [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) version 4.1.10; database downloaded 03/03/2022
* [LRE-Finder](https://bitbucket.org/genomicepidemiology/lre-finder/src/master/) version 1.0.0

# Citation

Coll F, *et al.* Improved accuracy of antibiotic resistance determination from *Enterococcus faecium* whole-genome sequences. (submitted for publication)
