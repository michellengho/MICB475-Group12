library("phyloseq")
library("tidyverse")
library("ape")
library("picante")
library("vegan")

metafp <- "ryan_metadata_new.tsv"
meta <- read_delim(metafp, delim="\t")

# Replace "T colon" with "Transverse" in Biopsy_location
meta <- meta %>%
  mutate(Biopsy_location = ifelse(Biopsy_location == "T colon", "Transverse", Biopsy_location))

# Filter metadata for IBD and smoker/non-smoker
meta <- meta %>%
  filter(Condition %in% c("Crohn's Disease", "Ulcerative Colitis") & 
           Smoking_status %in% c("Smoker", "Non-smoker"))

otufp <- "feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "tree.nwk"
phylotree <- read.tree(phylotreefp)
class(phylotree)

# save everything except first column (OTU ID) into a matrix
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
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
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

# Remove non-bacterial sequences, if any
ryan_filt <- subset_taxa(ryan,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

# Prune phyloseq object to retain only filtered metadata samples
ryan_filt <- prune_samples(sample_names(ryan_filt) %in% meta$`sample-id`, ryan_filt)

# Rarefy samples
ryan_filt <- rarefy_even_depth(ryan_filt, rngseed = 1, sample.size = 2198)

# Extract the sample names from the rarefied phyloseq object
retained_samples <- sample_names(ryan_filt)

# Calculate Shannon for the rarefied data
alpha <- estimate_richness(ryan_filt, measures = c("Shannon"))

# Add sample metadata for plotting
alpha <- cbind(alpha, as.data.frame(sample_data(ryan_filt)))

# Summarize Shannon Diversity by Group
alpha <- alpha %>%
  # Group data by Condition, Smoking_status, and Biopsy_location
  group_by(Condition, Smoking_status, Biopsy_location) %>%
  # Calculate mean and standard deviation for Shannon index within each group
  summarise(
    mean_shannon = mean(Shannon, na.rm = TRUE), # Mean Shannon diversity
    sd_shannon = sd(Shannon, na.rm = TRUE),     # Standard deviation of Shannon diversity
    .groups = "drop"                            # Drop grouping after summarisation
  )

# Filter out "Right colon undefined"
alpha <- alpha %>%
  filter(Biopsy_location != "Right colon undefined") %>%
  mutate(
    Biopsy_location = factor(
      Biopsy_location,
      levels = c(
        "Caecum", "Ascending", "Transverse", "Splenic flexture",
        "Descending", "Sigmoid", "Rectosigmoid", "Rectum", 
        "Right colon"
      )
    )
  )

# Heatmap with Mean Shannon Diversity
Heatmap <- ggplot(alpha, aes(x = Biopsy_location, y = interaction(Condition, Smoking_status, sep = "/"))) +
  geom_tile(aes(fill = mean_shannon), color = "white") +  # Add tiles with color based on mean_shannon
  geom_text(aes(label = round(mean_shannon, 2)), size = 3, color = "black") +  # Add text for mean values
  scale_fill_gradient(low = "#B6D0E2", high = "#6082B6", name = "Mean Shannon") +  # Color gradient for the heatmap
  labs(
    x = "Biopsy Location",
    y = "Condition / Smoking Status",
    title = "Heatmap of Shannon Diversity (Mean)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold and rotate x-axis text
    axis.text.y = element_text(size = 10, face = "bold"),  # Bold y-axis text
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),   # Center and style the title
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white", colour = NA),  # Set white background for panel
    plot.background = element_rect(fill = "white", colour = NA),   # Set white background for plot area
    axis.line = element_line(colour = "black"),  # Keep axis lines black
    legend.background = element_rect(fill = "white", colour = NA)  # Set white background for legend
  )



ggsave(
  filename = "shannon_heatmap.png",
  plot = Heatmap,
  width = 10,
  height = 7,
  dpi = 1200
)