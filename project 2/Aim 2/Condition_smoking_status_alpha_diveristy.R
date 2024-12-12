library("phyloseq")
library("tidyverse")
library("ape")
library("picante")
library("vegan")

# Load metadata
metafp <- "ryan_metadata_new.tsv"
meta <- read_delim(metafp, delim="\t")

# Filter metadata for IBD and smoker/non-smoker
meta <- meta %>%
  filter(Condition %in% c("Crohn's Disease", "Ulcerative Colitis") & 
           Smoking_status %in% c("Smoker", "Non-smoker"))

# Load OTU table
otufp <- "feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

# Load taxonomy table
taxfp <- "taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

# Load phylogenetic tree
phylotreefp <- "tree.nwk"
phylotree <- read.tree(phylotreefp)
class(phylotree) # Check tree class

# Create OTU matrix
# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

# Prepare sample data
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$`sample-id`
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

# Prepare taxonomy table
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX) # Verify taxonomy table class

# Merge all into a phyloseq object
ryan <- phyloseq(OTU, SAMP, TAX, phylotree)

# Remove non-bacterial sequences
ryan_filt <- subset_taxa(ryan,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

# Prune phyloseq object to retain only filtered metadata samples
ryan_filt <- prune_samples(sample_names(ryan_filt) %in% meta$`sample-id`, ryan_filt)

# Rarefy samples
ryan_rare_2 <- rarefy_even_depth(ryan_filt, rngseed = 1, sample.size = 2198)

# Retained and removed samples analysis
retained_samples <- sample_names(ryan_rare_2)
removed_samples <- setdiff(meta$`sample-id`, retained_samples)

# Count the number of retained samples per group (Condition + Smoking_status)
retained_counts <- meta %>%
  filter(`sample-id` %in% retained_samples) %>%
  group_by(Condition, Smoking_status) %>%
  summarise(count = n(), .groups = "drop")
print("Retained Samples Group Counts:")
print(retained_counts)

# Count the number of removed samples per group
removed_samples <- setdiff(meta$`sample-id`, retained_samples)
removed_counts <- meta %>%
  filter(`sample-id` %in% removed_samples) %>%
  group_by(Condition, Smoking_status) %>%
  summarise(count = n(), .groups = "drop")
print("Removed Samples Group Counts:")
print(removed_counts)

# Aim 2
# Calculate alpha diversity metrics for the rarefied data
alpha_div_2 <- estimate_richness(ryan_rare_2, measures = c("Shannon", "Observed"))

# Add sample metadata for plotting
alpha_div_2 <- cbind(alpha_div_2, sample_data(ryan_rare_2))

# Plot Shannon Diversity Index by Condition, colored by Smoking_status
shannon_plot_2 <- ggplot(alpha_div_2, aes(x = Condition, y = Shannon, color = Smoking_status)) +
  geom_boxplot() +
  labs(
    title = "Shannon Diversity Index by Condition and Smoking Status",
    x = "Condition",
    y = "Shannon Diversity Index",
    color = "Smoking Status"
  )

ggsave("shannon_diversity_2.png", plot = shannon_plot_2, height = 5, width = 6)

# Plot Observed Features by Condition, colored by Smoking_status
observed_plot_2 <- ggplot(alpha_div_2, aes(x = Condition, y = Observed, color = Smoking_status)) +
  geom_boxplot() +
  labs(
    title = "Observed Features by Condition and Smoking Status",
    x = "Condition",
    y = "Observed Features",
    color = "Smoking Status"
  )

ggsave("observed_features_2.png", plot = observed_plot_2, height = 5, width = 6)

# Calculate Chao1 richness for the rarefied data
chao1_div <- estimate_richness(ryan_rare_2, measures = c("Chao1"))

# Add Chao1 data to the alpha diversity table
alpha_div_2 <- cbind(alpha_div_2, Chao1 = chao1_div$Chao1)

# Visualizations for Chao1
chao1_plot <- ggplot(alpha_div_2, aes(x = Condition, y = Chao1, color = Smoking_status)) +
  geom_boxplot() +
  labs(
    title = "Chao1 Richness by Condition and Smoking Status",
    x = "Condition",
    y = "Chao1 Richness",
    color = "Smoking Status"
  )

ggsave("chao1_richness.png", plot = chao1_plot, height = 5, width = 6)

# Calculate Faith's Phylogenetic Diversity (PD)
phylo_dist <- pd(t(otu_table(ryan_rare_2)), phy_tree(ryan_rare_2), include.root = FALSE)

# Add Faith's PD to sample metadata
sample_data(ryan_rare_2)$Faith_PD <- phylo_dist$PD

# Confirm that Faith's PD has been added
head(sample_data(ryan_rare_2)$Faith_PD)

# Convert sample_data to a dataframe for ggplot
faith_pd_data <- as.data.frame(sample_data(ryan_rare_2))

# Boxplot for Faith's PD
faith_pd_plot <- ggplot(faith_pd_data, aes(x = Condition, y = Faith_PD, color = Smoking_status)) +
  geom_boxplot() +
  labs(
    title = "Faith's Phylogenetic Diversity by Condition and Smoking Status",
    x = "Condition",
    y = "Faith's PD",
    color = "Smoking Status"
  )

# Save the plot
ggsave("faith_pd_plot.png", plot = faith_pd_plot, height = 5, width = 6)

# Statistical Tests
# Combine sample data with alpha diversity estimates
samp_dat_wdiv <- data.frame(sample_data(ryan_rare_2), estimate_richness(ryan_rare_2))

# Kruskal-Wallis Tests
kruskal_shannon <- kruskal.test(Shannon ~ combined_condition, data = samp_dat_wdiv)
print(kruskal_shannon)

kruskal_chao1 <- kruskal.test(Chao1 ~ combined_condition, data = samp_dat_wdiv)
print(kruskal_chao1)

kruskal_faith_pd <- kruskal.test(Faith_PD ~ combined_condition, data = samp_dat_wdiv)
print(kruskal_faith_pd)

kruskal_observed <- kruskal.test(Observed ~ combined_condition, data = samp_dat_wdiv)
print(kruskal_observed)

# Perform pairwise Wilcoxon test for Chao1
if (kruskal_chao1$p.value < 0.05) {
  pairwise_chao1 <- pairwise.wilcox.test(samp_dat_wdiv$Chao1, samp_dat_wdiv$combined_condition, p.adjust.method = "BH")
  print("Pairwise comparisons for Chao1:")
  print(pairwise_chao1)
}