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

Bin_reads <- function(fastq_file, map_file, outdir) {
  #read fastq file
  fastq_reads <- readDNAStringSet(filepath = fastq_file, format = "fastq", with.qualities = TRUE)
  names(fastq_reads) <- gsub(x = names(fastq_reads), pattern = " .*", replacement = "")
  #read alignments file
  num_cols <- max(count.fields(map_file, sep = "\t"))
  map <- read.table(map_file, header = FALSE, sep = "\t")
  reads <- map[, 1]
  UMIs <- gsub(x = gsub(x = map[, 2], pattern = "centroid=", replacement = ""), pattern = ";.*", replacement = "")
  #assign reads to UMIs
  data <- data.frame(reads, UMIs)
  data_split <- split(data$reads, data$UMIs)
  maxLen <- 20000
  #cycle across UMIs and split reads to files
  for (u in 1:length(data_split)) {
    fastq_reads_curr_UMI <- fastq_reads[which(names(fastq_reads) %in% data_split[[u]])]
    tbd <- names(fastq_reads_curr_UMI[which(width(fastq_reads_curr_UMI) >= maxLen)])
    if (length(tbd) > 0) {
      cat(sprintf("Discarding the following reads, because they are longer than %d bp: %s\n", maxLen, paste0(tbd, collapse = ", ")))
      fastq_reads_curr_UMI <- fastq_reads_curr_UMI[which(width(fastq_reads_curr_UMI) < maxLen)]
    }
    writeXStringSet(x = fastq_reads_curr_UMI, filepath = paste0(outdir, "/", names(data_split[u]), "_chunk_", basename(fastq_file)), append = TRUE, format = "fastq")
  }
  #write unbinned reads to file
  unbinned_reads_ids <- setdiff(names(fastq_reads), reads)
  if (length(unbinned_reads_ids) > 0) {
    unbinned_reads <- fastq_reads[which(names(fastq_reads) %in% unbinned_reads_ids)]
    tbd <- names(unbinned_reads[which(width(unbinned_reads) >= maxLen)])
    if (length(tbd) > 0) {
      cat(sprintf("Discarding the following reads, because they are longer than %d bp: %s\n", maxLen, paste0(tbd, collapse = ", ")))
      unbinned_reads <- unbinned_reads[which(width(unbinned_reads) < maxLen)]
    }
    writeXStringSet(x = unbinned_reads, filepath = paste0(outdir, "/unbinned_chunk_", basename(fastq_file)), append = TRUE, format = "fastq")
  }
}

Bin_reads(fastq_file, map_file, outdir)