#!/bin/bash

nextflow -c UMInator.conf run UMInator.nf \
--FW_adapter="CAAGCAGAAGACGGCATACGAGAT" \
--RV_adapter="AATGATACGGCGACCACCGAGATC" \
--FW_primer="AGRGTTYGATYMTGGCTCAG" \
--RV_primer="CGACATCGAGGTGCCAAAC" \
--fastq_files=/home/simone/pipelines/UMInator/test_reads.fastq \
--results_dir=/home/simone/pipelines/UMInator/test_output \
-profile docker
