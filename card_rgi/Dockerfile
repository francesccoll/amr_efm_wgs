# Use a LTS Ubuntu version as parent image
FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

WORKDIR /opt

# metadata
LABEL base.image="ubuntu:20.04"
LABEL container.version="0.1"
LABEL software="rgi_card"
LABEL software.version="5.2.1"
LABEL description="Docker image of RGI CARD not using conda"
LABEL website="https://github.com/arpcard/rgi"
LABEL lisence="https://github.com/arpcard/rgi/blob/master/LICENSE"
LABEL maintainer="Francesc Coll"
LABEL maintainer.email="francesc.coll@lshtm.ac.uk"

# get some system essentials
RUN apt-get update && apt-get install --no-install-recommends -y \
    curl \
    build-essential \
    automake \
    cmake \
    pkg-config \
    zlib1g-dev \
    unzip \
    autoconf \
    check \
    libtool \
    libsubunit-dev \
    git \
    wget \
    python3.8 \
    python3.8-dev \
    python3-setuptools \
    python3-pip

# NOTE: Building HTSlib requires libbzip2 development 
# NOTE: HTSlib requires liblzma development files
RUN apt-get install -y libbz2-dev liblzma-dev

# getting RGI python dependencies
RUN pip3 install six biopython filetype pytest mock pandas matplotlib seaborn pyfaidx pyahocorasick pysam beautifulsoup4 requests lxml

# Installing BLAST
# from source as indicated here: https://www.ncbi.nlm.nih.gov/books/NBK569861/
# NOTE: for Mac, BLAST is installed from source, see Dockerfile_Mac
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz
RUN tar -xvf ncbi-blast-2.12.0+-x64-linux.tar.gz
RUN cp /opt/ncbi-blast-2.12.0+/bin/* /usr/local/bin/
RUN rm -r /opt/ncbi-blast-2.12.0+/

# Building diamond from source, adapted from https://github.com/bbuchfink/diamond/blob/master/Dockerfile
# RUN pip3 install diamond
WORKDIR /opt
RUN wget https://github.com/bbuchfink/diamond/archive/refs/tags/v2.0.14.tar.gz
RUN tar -xvf v2.0.14.tar.gz 
RUN cd diamond-2.0.14/
RUN mkdir build 
RUN cd build/
RUN cmake /opt/diamond-2.0.14
RUN make -j4 && make install
RUN rm -r /opt/diamond-2.0.14

WORKDIR /opt

# Installing prodigal from source
RUN git clone https://github.com/hyattpd/prodigal/
RUN cd /opt/prodigal/ && \
   make install

# get latest version of the repo
RUN pip3 install git+https://github.com/arpcard/rgi.git

# get and load latest version of CARD
RUN rgi auto_load




 
