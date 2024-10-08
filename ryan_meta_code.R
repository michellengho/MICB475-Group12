library(tidyverse)

# Load metadata
ryan_meta <- read_delim(file = "ryan_metadata.tsv", delim = "\t")
ryan_meta

#unite condition and smoking status column
ryanmeta_new <- unite(ryan_meta, col="combined condition", `Condition`, `Smoking.status`, sep = "/", remove = FALSE)
ryanmeta_new

write.table(ryanmeta_new, file = "new_ryan_metadata.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

