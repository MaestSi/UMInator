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

#from UMI -> reads to read -> UMIs
Reformat_alignment <- function(alignment, max_NM_mean, max_NM_sd) {
  matches <- c() 
  #cycle across alignments
  for (i in 1:nrow(alignment)) {
    #extract current alignment
    UMI_curr_alignment <- alignment[i, ]
    #extract UMI name, size, reads alignments and NM
    UMI_name <- UMI_curr_alignment[[1]]
    UMI_size <- as.numeric(gsub(x = UMI_name, pattern = ".*;size=", replacement = ""))
    prim_al_read_name <- UMI_curr_alignment[[2]]
    prim_al_NM <- as.numeric(gsub(pattern = "NM:i:", replacement = "", x = UMI_curr_alignment[[3]]))
    sec_al <- strsplit(x = UMI_curr_alignment[[4]], split = ";")[[1]]
    sec_al_read_name <- gsub(x = gsub(x = sec_al, pattern = "XA:Z:", replacement = ""), pattern = ",.*", replacement = "")
    sec_al_NM <- as.numeric(gsub(x = sec_al, pattern = ".*,", replacement = ""))
    NM_mean <- mean(c(prim_al_NM, sec_al_NM))
    NM_sd <- sd(c(prim_al_NM, sec_al_NM))
    #if the average number of NM per UMI is < max_NM_mean and the sd of the number of NM per UMI is < max_NM_sd, keep the UMI
    if (NM_mean < max_NM_mean && NM_sd < max_NM_sd) {
      maxDiff_flag <- "keep"
    } else {
      maxDiff_flag <- "discard"
    }
    #create a dataframe with the current set of matches
    matches_curr <- data.frame(read_name = c(prim_al_read_name, sec_al_read_name), UMI_name, UMI_size, NM = c(prim_al_NM, sec_al_NM), maxDiff_flag)
    matches <- rbind(matches, matches_curr)
  }
  #sort by read name, NM and UMI_size
  matches <- matches[with(matches, order(read_name, NM, -UMI_size)), ]
  #for each read, keep only the UMI with the lowest NM or the higher UMI size in case of ties
  tmp <- lapply(split(matches, matches$read_name), function(x) data.frame(x[1, ]))
  matches_best_hit <- do.call("rbind", tmp)
  matches_best_hit <- matches_best_hit[matches_best_hit$maxDiff_flag == "keep", ]
  return(matches_best_hit)
}

Filter_UMIs <- function(alignment_file_1, alignment_file_2, map_file, max_NM_mean, max_NM_sd) {
  #read alignment UMI 1
  num_cols_1 <- max(count.fields(alignment_file_1, sep = "\t"))
  alignment_1 <- read.table(alignment_file_1, fill = TRUE, quote = "", comment.char = "", colClasses = c("character", "NULL", "character", rep("NULL", 9), "character", rep("NULL", num_cols_1 - 14), "character"))
  #from UMI -> reads to read -> UMI for alignment_1
  matches_1 <- Reformat_alignment(alignment_1, max_NM_mean, max_NM_sd)
  #read alignment UMI 2
  num_cols_2 <- max(count.fields(alignment_file_2, sep = "\t"))
  alignment_2 <- read.table(alignment_file_2, fill = TRUE, quote = "", comment.char = "", colClasses = c("character", "NULL", "character", rep("NULL", 9), "character", rep("NULL", num_cols_2 - 14), "character"))
  #from UMI -> reads to read -> UMI for alignment_2
  matches_2 <- Reformat_alignment(alignment_2, max_NM_mean, max_NM_sd)
  #merge the two dataframes
  reads_both_matches <- intersect(matches_1$read_name, matches_2$read_name)
  matches_1_2 <- data.frame(read_name = reads_both_matches, UMI_name_1 = matches_1$UMI_name, UMI_name_2 = matches_2$UMI_name)
  #retain the read only in case both start and end map to the same UMI
  matches_1_2 <- matches_1_2[which(matches_1_2$UMI_name_1 == matches_1_2$UMI_name_2), ]
  map <- matches_1_2[, c(1, 2)]
  #write to file
  write.table(x = map, file = map_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

Filter_UMIs(alignment_file_1, alignment_file_2, map_file, max_NM_mean, max_NM_sd)