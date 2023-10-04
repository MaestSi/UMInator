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

Obtain_draft_consensus <- function(fastq_file, TRC, PLUR, num_threads, fast_alignment_flag, fast_consensus_flag, min_UMI_freq) {
  TRC <- as.numeric(TRC)
  target_reads_consensus <- TRC
  PLUR <- as.numeric(PLUR)
  fast_consensus_flag <- as.numeric(fast_consensus_flag)
  min_UMI_freq <- as.numeric(min_UMI_freq)
  num_threads <- as.numeric(num_threads)
  fast_alignment_flag <- as.numeric(fast_alignment_flag)
  fasta_file <- gsub(pattern = "\\.fastq", replacement = ".fasta", x = fastq_file)
  sample_name <- basename(dirname(fastq_file))
  UMI_name <- gsub(pattern = "\\.fastq", replacement = "", x = basename(fastq_file))
  sample_dir <- gsub(pattern = "readsUMIsAssignment", replacement = paste0("draftConsensusCalling/", sample_name), x = dirname(dirname(fastq_file)))
  logfile <- paste0(sample_dir, "/logfile.txt")
  UMI_dir <- paste0(sample_dir, "/", UMI_name)
  draft_consensus_tmp1 <- paste0(sample_dir, "/", UMI_name, "/", UMI_name, "_draft_consensus_tmp1.fasta")
  draft_consensus_tmp2 <- paste0(sample_dir, "/", UMI_name, "/", UMI_name, "_draft_consensus_tmp2.fasta")
  draft_consensus <- paste0(sample_dir, "/", UMI_name, "/", UMI_name, "_draft_consensus.fasta")
  num_reads_UMI <- as.double(system(paste0("cat ", fasta_file, " | grep \"^>\" | wc -l"), intern=TRUE))
  #create consensus sequence
  #if not at least min_UMI_freq reads are assigned to the UMI, consensus calling is skipped
  if (num_reads_UMI < min_UMI_freq) {
    cat(text = paste0("WARNING: Only ", num_reads_UMI, " reads available for sample ", sample_name, " - ", UMI_name, "; skipping consensus sequence generation"), sep = "\n")
    cat(text = paste0("WARNING: Only ", num_reads_UMI, " reads available for sample ", sample_name, " - ", UMI_name, "; skipping consensus sequence generation"),  file = logfile, sep = "\n", append = TRUE)
  } else {
    dir.create(UMI_dir)
    if (num_reads_UMI >= min_UMI_freq && num_reads_UMI < target_reads_consensus) {
      target_reads_consensus <- num_reads_UMI
      cat(text = paste0("WARNING: Only ", num_reads_UMI, " reads available for sample ", sample_name, " - ", UMI_name), sep = "\n")
      cat(text = paste0("WARNING: Only ", num_reads_UMI, " reads available for sample ", sample_name, " - ", UMI_name),  file = logfile, sep = "\n", append = TRUE)
    } 
    plurality_value <- PLUR*target_reads_consensus
    sequences <- readDNAStringSet(fasta_file, "fasta")
    ws <- width(sequences)
    amplicon_length <- ceiling(mean(ws))
    draft_reads_fq <- paste0(UMI_dir, "/", UMI_name, "_draft_", target_reads_consensus, "_reads.fastq")
    draft_reads_fa <- paste0(UMI_dir, "/", UMI_name, "_draft_", target_reads_consensus, "_reads.fasta")
    seed <- 1
    system(paste0("seqtk sample -s ", seed , " ", fastq_file, " ",  target_reads_consensus, " > ", draft_reads_fq))
    system(paste0("seqtk seq -A ", draft_reads_fq, " > ", draft_reads_fa))
    
    #produce accurate draft consensus sequence with mafft + EMBOss cons
    if(fast_consensus_flag != 1) {
      mfa_file <- gsub(pattern = "\\.fasta$", replacement = ".mfa", x = draft_reads_fa)
      if (fast_alignment_flag == 1) {
        system(paste0("mafft --auto --thread ", num_threads, " --adjustdirectionaccurately ", draft_reads_fa, " > ", mfa_file))
      } else {
        system(paste0("mafft -linsi --thread ", num_threads, " --adjustdirectionaccurately ", draft_reads_fa, " > ", mfa_file))
      }
      system(paste0("cons -sequence ", mfa_file, " -plurality ", plurality_value, " -outseq ", draft_consensus_tmp1))
    #produce fast draft consensus sequence with VSEARCH
    } else {
       system(paste0("vsearch --threads ", num_threads, " --cluster_smallmem ", draft_reads_fa, " --usersort --id 0.75 --iddef 1 --strand both --clusterout_sort --consout ", draft_consensus_tmp1, " --minwordmatches 0"))
    }
    system(paste0("cat ", draft_consensus_tmp1, " | head -n2 | sed 's/[nN]//g' > ", draft_consensus_tmp2))
    DNAStringSet_obj <- readDNAStringSet(draft_consensus_tmp2, "fasta")
    DNAStringSet_obj_renamed <- DNAStringSet_obj
    original_headers <- names(DNAStringSet_obj)
    names(DNAStringSet_obj_renamed) <- paste0("consensus_", sample_name, "_", UMI_name, "_", num_reads_UMI, "_reads")
    sequences <- seq(DNAStringSet_obj)
    writeXStringSet(x = DNAStringSet_obj_renamed, filepath = draft_consensus, format = "fasta", width = 20000)
  }
}

Obtain_draft_consensus(fastq_file, TRC, PLUR, num_threads, fast_alignment_flag, fast_consensus_flag, min_UMI_freq)