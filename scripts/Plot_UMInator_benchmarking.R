library("ggplot2")
library("taxize")
library("ggpubr")

#MetaBlast files with single-read alignment results
longread_umi_file <- "/path/to/longread_umi_blast_hits_unique*.txt"
UMInator_file <- "/path/to/UMInator_blast_hits_unique*.txt"
raw_reads_file <- "/path/to/Raw_reads_blast_hits_unique*.txt"
#MetaBlast files with summary alignment results
longread_umi_summary_file <- "/path/to/longread_umi_summary_blast_hits_unique*.txt"
UMInator_summary_file <- "/path/to/UMInator_summary_blast_hits_unique*.txt"
raw_reads_summary_file <- "/path/to/Raw_reads_summary_blast_hits_unique*.txt"
#Set of target species expected to be in the mock community
target_species <- c("Bacillus subtilis", "Staphylococcus aureus", "Limosilactobacillus fermentum", "Listeria monocytogenes", "Escherichia coli", "Pseudomonas aeruginosa", "Enterococcus faecalis", "Salmonella enterica")

#define functions
Retrieve_taxa = function(blast_hits_file) {
  #read file
  blast_hits <- read.table(blast_hits_file,sep = "\t", comment.char = "", header = TRUE, quote = "")
  #retrieve alignment identity and taxid
  blast_hits_alid <- blast_hits[, "Alignment.identity.perc."]
  blast_hits_taxid <- as.character(blast_hits[, "Taxonomy.id"])
  #retrieve species name
  max_attempts <- 10
  num_uids_chunk <- 100
  #process NBCI taxonomy UID in chunks
  if (length(blast_hits_taxid) < num_uids_chunk) {
    chunks_list <- list(1:length(num_uids_chunk))
  } else {
    chunks_list <- split(1:length(blast_hits_taxid), ceiling(seq(from = 1, to = length(blast_hits_taxid))/num_uids_chunk))
  }
  raw_classification <- c()
  for (i in 1:length(chunks_list)) {
    cat(sprintf("Processing NCBI taxonomy UID chunk %d out of %d (%.2f%%)\n", i, length(chunks_list), i*100/length(chunks_list)))
    raw_classification_curr <- NULL
    attempt <- 1
    while(length(which(!is.na(raw_classification_curr))) != lapply(chunks_list, length)[[i]] && attempt <= max_attempts) {
      raw_classification_curr <- try(suppressMessages(classification(blast_hits_taxid[chunks_list[[i]]], db = 'ncbi')))
      if (attempt > 1)  cat(sprintf("Running attempt %d\n", attempt))
      attempt <- attempt + 1
    }
    raw_classification <- c(raw_classification, raw_classification_curr)
  }
  cat("Reformatting taxonomy\n")
  species_name <- c()
  full_taxonomy <- c()
  for (i in 1:length(raw_classification)) {
    if (i %% 100 == 0) {
      cat(sprintf("Reformatted taxonomy for %d entries out of %d (%.2f%%)\n", i, length(blast_hits_taxid), 100*i/length(blast_hits_taxid)))
    }
    if (!is.null(ncol(raw_classification[[i]]))) {
      if (nrow(raw_classification[[i]]) > 1) {
        subspecies_name_curr <- raw_classification[[i]][which(raw_classification[[i]][, 2] == "subspecies"), 1]
        species_name_curr <- raw_classification[[i]][which(raw_classification[[i]][, 2] == "species"), 1]
        genus_name_curr <- raw_classification[[i]][which(raw_classification[[i]][, 2] == "genus"), 1]
        family_name_curr <- raw_classification[[i]][which(raw_classification[[i]][, 2] == "family"), 1]
        order_name_curr <- raw_classification[[i]][which(raw_classification[[i]][, 2] == "order"), 1]
        class_name_curr <- raw_classification[[i]][which(raw_classification[[i]][, 2] == "class"), 1]
        phylum_name_curr <- raw_classification[[i]][which(raw_classification[[i]][, 2] == "phylum"), 1]
        kingdom_name_curr <- raw_classification[[i]][which(raw_classification[[i]][, 2] == "superkingdom"), 1]
        species_name[i] <- species_name_curr
        full_taxonomy[i] <- paste(kingdom_name_curr, phylum_name_curr, class_name_curr, order_name_curr, family_name_curr, genus_name_curr, species_name_curr, subspecies_name_curr, sep = ";")
      } else {
        species_name[i] <- NA
        full_taxonomy[i] <- "Unclassified;;;;;;;"
      }
    } else {
      species_name[i] <- NA
      full_taxonomy[i] <- "Unclassified;;;;;;;"
    }
  }
  return(list(species_name, blast_hits_alid))
}
Reformat_freq = function(blast_hits_summary_file, target_species) {
  #read file
  blast_hits_summary <- read.table(blast_hits_summary_file, sep = "\t", comment.char = "", header = TRUE, quote = "")
  #Assign non-target species to "Other"
  blast_hits_summary[!blast_hits_summary$Species %in% target_species, "Species"] <- "Other"
  #Merge counts with the same species
  blast_hits_species_counts <- unlist(lapply(split(blast_hits_summary[, "Read.Counts"], blast_hits_summary[, "Species"]), function(x) sum(x)))
  #Compute relative frequency
  blast_hits_species_counts <- blast_hits_species_counts/sum(blast_hits_species_counts)
  return(blast_hits_species_counts)
}
lm_eqn <- function(df){
  m <- lm(UMInator ~ longread_umi, df);
  eq <- substitute(~~italic(r)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

#Retrieve species name and alignment identity for each read and pipeline
raw_reads_taxid_species_alid <- Retrieve_taxa(raw_reads_file)
raw_reads_taxid_species <- raw_reads_taxid_species_alid[[1]]
raw_reads_alid <- raw_reads_taxid_species_alid[[2]]
longread_umi_taxid_species_alid <- Retrieve_taxa(longread_umi_file)
longread_umi_taxid_species <- longread_umi_taxid_species_alid[[1]]
longread_umi_alid <- longread_umi_taxid_species_alid[[2]]
UMInator_taxid_species_alid <- Retrieve_taxa(UMInator_file)
UMInator_taxid_species <- UMInator_taxid_species_alid[[1]]
UMInator_alid <- UMInator_taxid_species_alid[[2]]

#plot
df <- data.frame(alignment_id = c(raw_reads_alid, UMInator_alid, longread_umi_alid), Pipeline = factor(c(rep("Raw reads", length(raw_reads_alid)), rep("UMInator", length(UMInator_alid)), rep("longread_umi", length(longread_umi_alid))), levels = c("Raw reads", "longread_umi", "UMInator")))
df_raw <- df[df$Pipeline == "Raw reads", ]
df_consensus <- df[df$Pipeline != "Raw reads", ]

p1 <- ggplot(df_raw, aes(x = alignment_id, color = Pipeline, fill = Pipeline)) +
  geom_density(alpha = 0.3, adjust = 2) +
  ylab("Probability density") +
  xlab("Alignment identity%") +
  scale_color_manual(values=c("#F8766D")) +
  scale_fill_manual(values=c("#F8766D")) +
  xlim(c(80, 100)) +
  theme(axis.text=element_text(size=8)) + 
  theme(axis.title=element_text(size=8)) 

p2 <- ggplot(df_consensus, aes(x = alignment_id, color = Pipeline, fill = Pipeline)) +
  geom_density(alpha = 0.3, adjust = 2) +
  ylab("Probability density") +
  xlab("Alignment identity%") +
  scale_color_manual(values=c("#00BA38", "#619CFF")) +
  scale_fill_manual(values=c("#00BA38", "#619CFF")) +
  theme(axis.text=element_text(size=8)) + 
  theme(axis.title=element_text(size=8)) 

#Evaluate relative abundance for each target species and pipeline
raw_reads_species_counts <- Reformat_freq(raw_reads_summary_file, target_species)
longread_umi_species_counts <- Reformat_freq(longread_umi_summary_file, target_species)
UMInator_species_counts <- Reformat_freq(UMInator_summary_file, target_species)

#plot
df2 <- data.frame(Species = factor(c(names(raw_reads_species_counts), names(UMInator_species_counts), names(longread_umi_species_counts)), levels = c(target_species, "Other")), Counts = c(unname(raw_reads_species_counts), unname(UMInator_species_counts), unname(longread_umi_species_counts)), Pipeline = factor(c(rep("Raw reads", length(raw_reads_species_counts)), rep("UMInator", length(UMInator_species_counts)), rep("longread_umi", length(longread_umi_species_counts))), levels = c("Raw reads", "longread_umi", "UMInator")))

p3 <- ggplot(df2, aes(x = Species, y = Counts, fill = Pipeline)) +
  geom_bar(position="dodge", stat="identity", alpha = 0.8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Rel. frequency") +
  theme(axis.text=element_text(size=8)) + 
  theme(axis.title=element_text(size=8)) 

df2_consensus_tmp <- df2[df2$Pipeline != "Raw reads", ]
df2_consensus <- data.frame(Species = df2_consensus_tmp$Species[which(df2_consensus_tmp$Pipeline == "longread_umi")],
                            UMInator =  df2_consensus_tmp$Counts[which(df2_consensus_tmp$Pipeline == "UMInator")],
                            longread_umi = df2_consensus_tmp$Counts[which(df2_consensus_tmp$Pipeline == "longread_umi")])

p4 <- ggplot(df2_consensus, aes(x = UMInator, y = longread_umi, label = Species)) +
  geom_point(size = 5, aes(colour = Species)) +
  geom_text(hjust = 0.5, vjust = 1, angle = 0, size = 1.7) +
  geom_smooth(method = lm, color = "black") +
  geom_text(x = -0.1, y = 0.27, label = lm_eqn(df2_consensus[, c(2, 3)]), parse = TRUE) +
  theme(legend.position = "none") +
  xlim(c(-0.2, 0.5))


ggarrange(p1, p2, p3, p4, 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

ggsave("UMInator_benchmarking.png", width = 6, height = 6)  
