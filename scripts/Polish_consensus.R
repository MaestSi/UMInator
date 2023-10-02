#
# Copyright 2023 Simone Maestri. All rights reserved.
# Simone Maestri <simone.maestri@unimi.it>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

args = commandArgs(trailingOnly=TRUE)

for(v in args)
{
  vTmp <- strsplit(v,"=")[[1]]
  assign(vTmp[[1]],vTmp[[2]])
}

#load BioStrings package
suppressMessages(library(Biostrings))

Polish_consensus <- function(draft_consensus, fastq_file, num_threads, TRP, medaka_model) {
  if (file.exists(draft_consensus)) {
    num_threads <- as.numeric(num_threads)
    TRP <- as.numeric(TRP)
    target_reads_polishing <- TRP
    fasta_file <- gsub(pattern = "\\.fastq", replacement = ".fasta", x = fastq_file)
    UMI_name <- gsub(pattern = "\\.fastq", replacement = "", x = basename(fastq_file))
    sample_name <- basename(dirname(fastq_file))
    sample_dir <- gsub(pattern = "readsUMIsAssignment", replacement = paste0("consensusPolishing/", sample_name), x = dirname(dirname(fastq_file)))
    logfile <- paste0(sample_dir, "/logfile.txt")
    UMI_dir <- paste0(sample_dir, "/", UMI_name)
    dir.create(UMI_dir)
    logfile <- paste0(dirname(UMI_dir), "/logfile.txt")
    num_reads_UMI <- as.double(system(paste0("cat ", fasta_file, " | grep \"^>\" | wc -l"), intern=TRUE))
    num_threads_medaka <- min(num_threads, 8)
    paf_file <- paste0(UMI_dir, "/", UMI_name, ".paf")
    racon_consensus <- paste0(UMI_dir, "/", UMI_name, "_polished_consensus_tmp.fasta")
    seed <- 2
    consensus_polished <- paste0(UMI_dir, "/", UMI_name, "_polished_consensus.fasta")
    if (num_reads_UMI < target_reads_polishing) {
      target_reads_polishing <- num_reads_UMI
    }
    polishing_reads_fq <- paste0(UMI_dir, "/", UMI_name, "_polishing_", target_reads_polishing, "_reads.fastq")
    system(paste0("seqtk sample -s ", seed , " ", fastq_file, " ",  target_reads_polishing, " > ", polishing_reads_fq))
    cat(text = paste0("Running Racon for consensus polishing of sample ", sample_name, " - UMI ", UMI_name), sep = "\n")
    system(paste0("minimap2 -x ava-ont ", draft_consensus, " ", polishing_reads_fq, " > ", paf_file))
    system(paste0("racon -t ", num_threads, " -m 8 -x -6 -g -8 -w 500 --no-trimming ", polishing_reads_fq, " ", paf_file, " ", draft_consensus, " > ", racon_consensus))
    if (length(which(readLines(racon_consensus) != "")) == 0) {
      cat(text = paste0("WARNING: Failed to run Racon for sample ", sample_name, " - UMI ", UMI_name), sep = "\n")
      cat(text = paste0("WARNING: Failed to run Racon for sample ", sample_name, " - UMI ", UMI_name),  file = logfile, sep = "\n", append = TRUE)
      system(paste0("cp ", draft_consensus, " ", racon_consensus))
    }
    cat(text = paste0("Running Medaka for consensus polishing of sample ", sample_name, " - UMI ", UMI_name), sep = "\n")
    system(paste0("medaka_consensus -i ", polishing_reads_fq, " -d ", racon_consensus, " -m ", medaka_model, " -t ", num_threads_medaka, " -o ", UMI_dir, "/medaka_consensus"))
    system(paste0("cp ", UMI_dir, "/medaka_consensus/consensus.fasta ", consensus_polished))  
  }
}

Polish_consensus(draft_consensus, fastq_file, num_threads, TRP, medaka_model)