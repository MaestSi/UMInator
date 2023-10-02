FROM continuumio/miniconda3

MAINTAINER Simone Maestri <simone.maestri@unimi.it>

RUN apt-get update  && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install -y \
    build-essential \
    wget \
    unzip \
    bzip2 \
    git \
    libidn11* \
    nano \
    less \
    wget \
    r-base \
    bc \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    gfortran \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN conda config --add channels r && \
conda config --add channels anaconda && \
conda config --add channels conda-forge && \
conda config --add channels bioconda

RUN conda create -n UMInator_env pip bioconductor-biostrings r-kernsmooth r-factoextra python==3.10
RUN /opt/conda/envs/UMInator_env/bin/python -m pip install medaka
RUN conda install -n UMInator_env emboss
RUN conda install -n UMInator_env vsearch
RUN conda install -n UMInator_env seqtk
RUN conda install -n UMInator_env mafft
RUN conda install -n UMInator_env minimap2
RUN conda install -n UMInator_env bwa
RUN conda install -n UMInator_env samtools
RUN conda install -n UMInator_env bcftools

RUN conda create -n UMInator_supp_env cutadapt pip racon
RUN /opt/conda/envs/UMInator_supp_env/bin/python -m pip install NanoPlot
RUN /opt/conda/envs/UMInator_supp_env/bin/python -m pip install NanoFilt
RUN ln -s /opt/conda/envs/UMInator_supp_env/bin/NanoPlot /opt/conda/envs/UMInator_env/bin
RUN ln -s /opt/conda/envs/UMInator_supp_env/bin/NanoFilt /opt/conda/envs/UMInator_env/bin
RUN ln -s /opt/conda/envs/UMInator_supp_env/bin/racon /opt/conda/envs/UMInator_env/bin
RUN ln -s /opt/conda/envs/UMInator_supp_env/bin/cutadapt /opt/conda/envs/UMInator_env/bin
ENV PATH=$PATH:/opt/conda/envs/UMInator_env/bin
