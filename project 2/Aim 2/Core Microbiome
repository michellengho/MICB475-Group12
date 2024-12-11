# Load necessary packages for analysis
library("tidyverse")
library("phyloseq")
library("microbiome")
library("ape")
library("vegan")

# Load necessary packages for plotting a venn diagram
library("ggVennDiagram")
library("ggplot2")
library("UpSetR")

# Organize the meta, otu, taxonomy, and tree files
metafp <- "ryan_metadata_new.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "tree.nwk"
phylotree <- read.tree(phylotreefp)
class(phylotree)

# Save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
class(OTU)

# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$`sample-id`
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into =
             c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

# Merge all into a phyloseq object
ryan <- phyloseq(OTU, SAMP, TAX, phylotree)

# Subset samples by combined condition of smoking status and disease condition
cd_smoker <- subset_samples(ryan, combined_condition == "Crohn's Disease/Smoker")
cd_non_smoker <- subset_samples(ryan, combined_condition == "Crohn's Disease/Non-smoker")
uc_smoker <- subset_samples(ryan, combined_condition == "Ulcerative Colitis/Smoker")
uc_non_smoker <- subset_samples(ryan, combined_condition == "Ulcerative Colitis/Non-smoker")

# Identify core members with a detection of 0.01 and 
cd_smoker_ASVs <- core_members(cd_smoker, detection = 0.01, prevalence = 0.65)
cd_non_smoker_ASVs <- core_members(cd_non_smoker, detection = 0.01, prevalence = 0.65)
uc_smoker_ASVs <- core_members(uc_smoker, detection = 0.01, prevalence = 0.65)
uc_non_smoker_ASVs <- core_members(uc_non_smoker, detection = 0.01, prevalence = 0.65)
# The detection is set to 0.01 (1% relative abundance) in order to filter out rare things
# and only consider those that are somewhat abundant 
# The prevalence threshold is set to 0.65 in order to see differences between groups
# which means it  has to be present in 65% of samples in that group to be considered

# View core microbes
view(cd_smoker_ASVs)
view(uc_smoker_ASVs)
view(cd_non_smoker_ASVs)
view(uc_non_smoker_ASVs)

# Plot venn diagram
venn_plot <- ggVennDiagram(x = list(CD_Smokers = cd_smoker_ASVs, CD_Non_Smokers = cd_non_smoker_ASVs, UC_Smokers = uc_smoker_ASVs, UC_Non_Smokers = uc_non_smoker_ASVs)) +
  scale_fill_gradient(low = "#B6D0E2", high = "#6082B6") +
  scale_x_continuous(expand = expansion(mult = .1))

# Save venn diagram
ggsave(filename = "core_microbiome.png", plot = venn_plot, width = 9, height = 7)
