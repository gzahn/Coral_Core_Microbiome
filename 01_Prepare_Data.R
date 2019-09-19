# Combine P acuta and P lutea phyloseq objects and prepare for analyses ####

# Load Packages
library(phyloseq)
library(tidyverse)
library(vegan)

# Load data ####
pluteus <- readRDS("./data/P_lutea_clean_phyloseq.RDS")
pacuta <- readRDS("./data/P_acuta_clean_phyloseq.RDS")


# fix column names
pacuta@sam_data$Location <- as.character(pacuta@sam_data$Island.Collected.From)
pacuta@sam_data$Library.ID <-  as.character(pacuta@sam_data$SampleName)

# Parse components for merging
pacuta.otu <- otu_table(pacuta)
pacuta.tax <- tax_table(pacuta)
pacuta.sam <- as(sample_data(pacuta), "data.frame")
pluteus.otu <- otu_table(pluteus) 
pluteus.tax <- tax_table(pluteus)
pluteus.sam <- as(sample_data(pluteus),"data.frame")

# Merge phyloseq objects ####
ps <- merge_phyloseq(pluteus.otu,sample_data(pluteus.sam),pluteus.tax,pacuta.otu,sample_data(pacuta.sam),pacuta.tax)

# Fix metadata
ps@sam_data$Species[ps@sam_data$Species=="1"] <- "P acuta"

# look at read coverage for all merged samples
ggplot(mapping=aes(x=1:nrow(otu_table(ps)),y=rowSums(otu_table(ps)))) + geom_point(aes(color=ps@sam_data$Species)) + 
  labs(x="Sample number",y="Read count",color="Species") + theme_bw()

# Add read count to sampledata
sample_data(ps)$ReadDepth <- rowSums(otu_table(ps))

# Investigate taxonomic classification rates of reads and ESVs ####
# fraction of reads with tax assignment at each tax rank
Sys.time()
reads_classified <- psmelt(ps) %>%
  gather("Rank","Name",rank_names(ps)) %>%
  group_by(Rank) %>%
  summarize(Reads_classified = sum(Abundance * !is.na(Name))) %>%
  mutate(Frac_classified = Reads_classified / sum(sample_sums(ps))) %>%
  mutate(Rank = factor(Rank, rank_names(ps))) %>%
  arrange(Rank)
Sys.time()

# fraction of ASVs with tax assignment at each tax rank
Sys.time()
ESVs_classified <- tax_table(ps) %>%
  as("matrix") %>%
  as_tibble(rownames = "OTU") %>%
  gather("Rank","Name",rank_names(ps)) %>%
  group_by(Rank) %>%
  summarize(OTUs_classified = sum(!is.na(Name))) %>%
  mutate(Frac_classified = OTUs_classified / ntaxa(ps)) %>%
  mutate(Rank = factor(Rank, rank_names(ps))) %>%
  arrange(Rank)
Sys.time()

sink("./output/Taxonomic_Classification_Stats.txt")
reads_classified
ESVs_classified
sink(NULL)


# Collapse based on taxonomy at genus level ####
ps2 = tax_glom(ps,taxrank = rank_names(ps)[6])


# Remove low-abundance OTUs ####

# Look at taxa sums 
summary(taxa_sums(ps2))

# Keep OTUs with at least 10 occurrances
ps.min <- subset_taxa(ps2, taxa_sums(ps2) >= 10)

saveRDS(ps.min,"./data/combined_and_cleaned_ps_object.RDS")


