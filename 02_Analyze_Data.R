# Core microbiome analyses
# Depends on 01_Prepare_Data.R

# Load packages
library(phyloseq)
library(tidyverse)
library(vegan)

# Load cleaned data
ps.min <- readRDS("./data/combined_and_cleaned_ps_object.RDS")

# species name punctuation
ps.min@sam_data$Species[ps.min@sam_data$Species == "P lutea"] <- "P. lutea"
ps.min@sam_data$Species[ps.min@sam_data$Species == "P acuta"] <- "P. acuta"

# Ordination ####

ord1=ordinate(ps.min,method = "NMDS",color="Species")
plot_ordination(ps.min,ord1,color="Species") + theme_minimal() +
  theme(legend.text = element_text(face = "italic"))
ggsave("./output/figs/Pacuta_vs_Plutea_NMDSPlot.png",dpi = 300,height = 8,width = 10)


mds <- metaMDS(otu_table(ps.min))
png("./output/figs/NMDS_Stressplot.png")
stressplot(mds)
dev.off()


# Calculate distance matrix
comm.dist <- vegdist(otu_table(ps.min),method = "bray")

plot(comm.dist)

# Convert to presence-absence ####
ps_pa <- transform_sample_counts(ps.min, function(abund) 1*(abund>0))


# prepare data for heatmap and plot ####
pa = as(t(otu_table(ps_pa)),"matrix")

cols = plyr::mapvalues(ps_pa@sam_data$Species,from=unique(ps_pa@sam_data$Species),to=c("Blue","Red"))

heatmap(pa,Rowv = NA,ColSideColors = cols,Colv = NA,col=gray.colors(10))

# Blue = P lutea, Red = P acuta
heatmap(t(as.matrix(otu_table(ps_pa))),Rowv = NA,ColSideColors = cols,Colv = NA, labRow = NA, col = gray.colors(2))

# maybe invert colors for above heatmap???


# Find genera tha overlap between species of corals ####

Plutea.cols = grep("ABB",x = colnames(pa))
Pacuta.cols = grep("AOO",x = colnames(pa))
Plutea.matrix = pa[,Plutea.cols]
Pacuta.matrix = pa[,Pacuta.cols]

shared.genera = rowSums(Plutea.matrix) > 0 & rowSums(Pacuta.matrix) > 0
shared.genera = which(shared.genera == TRUE)
sg <- ps_pa@tax_table[shared.genera,"Genus"]
sg <- as.data.frame(sg)

sink("./output/shared_Genera_Pactua-Plutea.txt")
sg$Genus
sink(NULL)


# Prep shared taxa data frame
shared.genera.names = c(ps_pa@tax_table[shared.genera,"Genus"])
shared.family.names = c(ps_pa@tax_table[shared.genera,"Family"])
shared.order.names = c(ps_pa@tax_table[shared.genera,"Order"])
shared.phylum.names = c(ps_pa@tax_table[shared.genera,"Phylum"])

df.shared.taxa <- data.frame(Genus = shared.genera.names,Family = shared.family.names,Order=shared.order.names,Phylum=shared.phylum.names)

# Bar plot, colored by phylum
ggplot(df.shared.taxa) + geom_bar(aes(x=reorder(Order,Order,function(x)-length(x)),fill=Phylum)) +
  theme(axis.text.x = element_text(angle=90),
        axis.title = element_text(face="bold")) + 
  labs(x="Order",y="Shared Count")
ggsave("./output/figs/Shared_Genera_Pacuta-Plutea_within_each-order.png",dpi=300)


# PermANOVA ####
mod1 <- adonis(otu_table(ps) ~ sample_data(ps)$Species + sample_data(ps)$ReadDepth)
mod1
