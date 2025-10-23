
library(GenomicRanges)
library(ggplot2)
library(dplyr)


de_exons <- read.delim("/home/roger/Downloads/compare/Sig_exons/LOW_vs_HIGH_c", header = TRUE, sep = "\t")
all_exons <-  read.delim("/home/roger/Downloads/INCLUSION_LEVELS_FULL-hg38-9.tab", header = TRUE, sep = "\t")

# Load your DE exon results and background
# Calculate basic statistics
de_stats <- data.frame(
  mean_LENGTH = mean(de_exons$LENGTH),
  median_LENGTH = median(de_exons$LENGTH),
  sd_LENGTH = sd(de_exons$LENGTH)
)


bg_stats <- data.frame(
  mean_LENGTH = mean(all_exons$LENGTH),
  median_LENGTH = median(all_exons$LENGTH),
  sd_LENGTH = sd(all_exons$LENGTH)
)

# Statistical test
LENGTH_test <- wilcox.test(de_exons$LENGTH, all_exons$LENGTH)

# Create comprehensive visualization
p1 <- ggplot() +
  geom_density(aes(x = log10(all_exons$LENGTH), color = "Background"), 
               alpha = 0.7) +
  geom_density(aes(x = log10(de_exons$LENGTH), color = "DE Exons"), 
               alpha = 0.7) +
  labs(title = "Exon LENGTH Distribution",
       x = "log10(Exon LENGTH)", y = "Density") +
  theme_minimal()

print(p1)
print(LENGTH_test)


######## exon coordinates ############
library(clusterProfiler)
library(org.Hs.eg.db)

coord <-  read.delim("/home/roger/Downloads/exons_coordinates.csv", header = FALSE, sep = ",")
alt_exp_exons <-  read.delim("/home/roger/Downloads/compare/Alt_Exons/AltEx-LOW_vs_HIGH_c.txt", header = FALSE, sep = "\t")
detected_exons <- alt_exp_exons[which(alt_exp_exons$V1 %in% coord$V1),]

