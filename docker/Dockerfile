# Use a LTS Ubuntu version as parent image
FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive

WORKDIR /opt

# metadata
LABEL base.image="ubuntu:22.04"
LABEL container.version="0.1"
LABEL software="amrwgsefm"
LABEL software.version="1.0.0"
LABEL description="Docker image of ARIBA-based pipeline for genotypic prediction of AMR from E. faecium genomes"
LABEL website="XXXX"
LABEL lisence="XXXX"
LABEL maintainer="Francesc Coll"
LABEL maintainer.email="francesc.coll@lshtm.ac.uk"

# NOTE: Docker commands to install ARIBA adapted from: https://github.com/sanger-pathogens/ariba/blob/master/Dockerfile
ARG BOWTIE2_VERSION=2.4.5
ARG SPADES_VERSION=3.13.1
ARG ARIBA_VERSION=2.14.6
ARG ARIBA_TAG=master
ARG ARIBA_BUILD_DIR=/ariba
ARG MUMMER_BUILD_DIR=/mummer
ARG MUMMER_VER="4.0.0"

ARG LOCALE_COUNTRY=en_GB

# NOTE: libsimde-dev needed to install bowtie2 from source. ubuntu:22.04 needed to install libsimde-dev, libsimde-dev not available on 20.04
# NOTE: automake added to avoid compilation error: "fatal error: emmintrin.h: No such file or directory" > to remove

RUN apt-get update && apt-get install --no-install-recommends -y \
    build-essential \
    autoconf \
    automake \
    cd-hit \
    curl \
    git \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libcurl4-gnutls-dev \
    libsimde-dev \
    libssl-dev \
    python3-dev \
    python3-setuptools \
    python3-pip \
    python3-tk \
    python3-matplotlib \
    unzip \
    wget \
    zlib1g-dev


# Install locales
RUN apt-get update && apt-get install -y locales-all && rm -rf /var/lib/apt/lists/* 
# Set a default locale.
ENV LANG=${LOCALE_COUNTRY}.UTF-8 \
    LANGUAGE=${LOCALE_COUNTRY}:en 

# Install bowtie2
# RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
#  unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
#  rm -r bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip

# bowtie2 installed from source to avoid error:
# /lib64/ld-linux-x86-64.so.2: No such file or directory
RUN wget -O bowtie2-v${BOWTIE2_VERSION}.tar.gz https://github.com/BenLangmead/bowtie2/archive/refs/tags/v${BOWTIE2_VERSION}.tar.gz && \
  tar -xvf bowtie2-v${BOWTIE2_VERSION}.tar.gz && \
  rm bowtie2-v${BOWTIE2_VERSION}.tar.gz && \
  cd bowtie2-${BOWTIE2_VERSION} && \
  make && \
  make install

# install SPAdes
RUN wget https://github.com/ablab/spades/releases/download/v${SPADES_VERSION}/SPAdes-${SPADES_VERSION}-Linux.tar.gz && \
  tar -xzf SPAdes-${SPADES_VERSION}-Linux.tar.gz && \
  rm -r SPAdes-${SPADES_VERSION}-Linux.tar.gz


# Need MPLBACKEND="agg" to make matplotlib work without X11, otherwise get the error
# _tkinter.TclError: no display name and no $DISPLAY environment variable
ENV ARIBA_BOWTIE2=/opt/bowtie2-${BOWTIE2_VERSION}-linux-x86_64/bowtie2 ARIBA_CDHIT=cdhit-est MPLBACKEND="agg"
ENV PATH=$PATH:/opt/SPAdes-${SPADES_VERSION}-Linux/bin
# ENV PATH=$PATH:/opt/bowtie2-${BOWTIE2_VERSION}-linux-x86_64

RUN ln -s -f /usr/bin/python3 /usr/local/bin/python

# installing mummer
# Docker commands to install mummer extracted from:
# https://github.com/StaPH-B/docker-builds/blob/master/mummer/4.0.0/Dockerfile
RUN wget https://github.com/mummer4/mummer/releases/download/v${MUMMER_VER}rc1/mummer-${MUMMER_VER}rc1.tar.gz && \
  tar -xvf mummer-${MUMMER_VER}rc1.tar.gz && \
  rm mummer-${MUMMER_VER}rc1.tar.gz && \
  cd mummer-${MUMMER_VER}rc1 && \
  ./configure && \
  make && \
  make install && \
  ldconfig && \
  cd .. && \
  mkdir /data


# Downloading Ariba files
RUN wget -O ariba-v${ARIBA_VERSION}.tar.gz -q https://github.com/sanger-pathogens/ariba/archive/refs/tags/v${ARIBA_VERSION}.tar.gz \
  && tar -zxf ariba-v${ARIBA_VERSION}.tar.gz


# Installing Ariba
RUN pip3 install ariba
# RUN mkdir -p $ARIBA_BUILD_DIR
# COPY /opt/ariba-${ARIBA_VERSION}/* $ARIBA_BUILD_DIR
# RUN cd $ARIBA_BUILD_DIR \
# RUN cd /opt/ariba-${ARIBA_VERSION}/ \
#  && python3 setup.py test \
#  && python3 setup.py install
# && rm -rf $ARIBA_BUILD_DIR


