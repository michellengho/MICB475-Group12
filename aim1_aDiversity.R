library("phyloseq")
library("tidyverse")
library("ape")
library("vegan")
library("picante")

# Load in the metadata, OTU table, taxonomy file, and phylogenetic tree. 
meta <- read_delim("project2_data/ryan_metadata_2.0.tsv", delim="\t")
otu <- read_delim(file = "project2_data/feature-table.txt", delim="\t", skip=1)
tax <- read_delim("project2_data/taxonomy.tsv", delim="\t")
phylotree <- read.tree("project2_data/exported-tree/tree.nwk")

# Adjust files to be read into a phyloseq object. Make the phyloseq object.
samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$'sample-id'
SAMP <- sample_data(samp_df)

otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

taxa_mat <- tax |> 
  select(-Confidence) |>
  separate(col = Taxon, sep = "; ",
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) |>
  as.matrix()
taxa_mat <- taxa_mat[,-1]
rownames(taxa_mat) <- tax$`Feature ID`
TAXA <- tax_table(taxa_mat)

TREE <- phy_tree(phylotree)

ryan <- phyloseq(OTU, SAMP, TAXA, TREE)

sample_data(ryan)
otu_table(ryan)
tax_table(ryan)
phy_tree(ryan)

# Remove non-bacterial sequences, if any
ryan_filt <- subset_taxa(ryan,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

# Rarefy samples
rarecurve(t(as.data.frame(otu_table(ryan_filt))), cex=0.1)
ryan_rare <- rarefy_even_depth(ryan_filt, rngseed = 1, sample.size = 13)
#sample size is decided based on table.qzv and a-rarefaction.qzv

#### Alpha diversity ######
plot_richness(ryan_rare) 

plot_richness(ryan_rare, measures = c("Shannon","Chao1")) 

gg_richness <- plot_richness(ryan_rare, x = "inflammation.and.smoking.combined", measures = c("Shannon","Chao1")) +
  xlab("inflammation and smoking combined") +
  geom_boxplot()
gg_richness

ggsave(filename = "plot_richness.png"
       , gg_richness
       , height=4, width=6)

richness_metrics <- estimate_richness(ryan_rare)

# Add observed features to sample data (the code didn't work specifically for this line)
sample_data(ryan_rare)$Observed <- richness_metrics$Observed

# Plot observed features
plot_observed <- ggplot(sample_data(ryan_rare), aes(inflammation.and.smoking.combined, Observed)) + 
  geom_boxplot() +
  xlab("Inflammation and Smoking Combined") +
  ylab("Number of Observed Features (OTUs/ASVs)")
ggsave(filename = "plot_observed_features_inflammation_smoking.png", plot_observed, height = 4, width = 6)


# phylogenetic diversity
# calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(ryan_rare)), phy_tree(ryan_rare),
                 include.root=F) 

# add PD to metadata table
sample_data(ryan_rare)$PD <- phylo_dist$PD

# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(ryan_rare), aes(inflammation.and.smoking.combined, PD)) + 
  geom_boxplot() +
  xlab("Inflammation and Smoking Status") +
  ylab("Phylogenetic Diversity")

# view plot
plot.pd


