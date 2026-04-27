# alpha_diversity analysis -------------------------------------------------------

library(readr)
library(dplyr)
library(tibble)
library(vegan)
library(ggplot2)
library(tidyverse)
library(readxl)
library(ggpubr)

# Reading the abundance matrix
mat_df <- read_tsv("/Users/yqhhh/Desktop/species_adversity/bracken_species_matrix.tsv", show_col_types = FALSE)

# The first column contains species names, which are then converted into a matrix.
mat <- mat_df %>%
  column_to_rownames("name") %>%
  as.matrix()

# Transpose to Sample × Species
mat <- t(mat)

# read metadata
meta <- read_excel("/Users/yqhhh/Desktop/metadata.xlsx")


# Ensure consistent sample order
meta <- meta %>% filter(Sample %in% rownames(mat))
mat <- mat[meta$Sample, , drop = FALSE]

# alpha diversity
alpha_df <- data.frame(
  Sample = rownames(mat),
  Shannon = diversity(mat, index = "shannon"),
  Simpson = diversity(mat, index = "simpson"),
  Richness = specnumber(mat)
)

alpha_df <- left_join(alpha_df, meta, by = "Sample")

write_tsv(alpha_df, "/Users/yqhhh/Desktop/alpha_diversity_results.tsv")

# Shannon box plot (by continent)）
my_colors <- c("#ffa88d", "#FFD84F", "#36CFFD", "#5CBC71")

comparisons_list <- list(
  c("Asia", "Europe"),
  c("Asia", "North America"),
  c("Asia", "Africa"),
  c("Europe", "North America"),
  c("Europe", "Africa"),
  c("North America", "Africa")
)

alpha_df$Continent <- factor(
  alpha_df$Continent,
  levels = c("Asia", "Europe", "North America", "Africa")
)

p1 <- ggplot(alpha_df, aes(x = Continent, y = Shannon, color = Continent)) +
  geom_boxplot(outlier.shape = NA,           
               fill = "white",               
               linewidth = 0.8) +            
  geom_jitter(width = 0.2,                  
              size = 1.5,                  
              alpha = 0.8,                 
              aes(color = Continent)) +
  scale_color_manual(values = my_colors) +  
  theme_classic() +
  labs(title = "Shannon diversity", 
       x = "Continent", 
       y = "Shannon") +
  stat_compare_means(
    comparisons = comparisons_list,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = FALSE,
    step.increase = 0.08
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )



print(p1)

ggsave("/Users/yqhhh/Desktop/alpha_shannon_boxplot.pdf", p1, width = 6, height = 5)



# Simpson box plot
p2 <- ggplot(alpha_df, aes(x = Continent, y = Simpson, color = Continent)) +
  geom_boxplot(
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.8
  ) +
  geom_jitter(
    width = 0.2,
    size = 1.5,
    alpha = 0.8
  ) +
  scale_color_manual(values = my_colors) +
  theme_classic() +
  labs(
    title = "Simpson diversity",
    x = "Continent",
    y = "Simpson",
    color = "Continent"
  ) +
  stat_compare_means(
    comparisons = comparisons_list,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = FALSE,
    step.increase = 0.06
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  ) +
  coord_cartesian(
    ylim = c(
      min(alpha_df$Simpson, na.rm = TRUE),
      max(alpha_df$Simpson, na.rm = TRUE) * 1.35
    )
  )

p2


ggsave("/Users/yqhhh/Desktop/alpha_simpson_boxplot.pdf", p2, width = 6, height = 5)

# Richness box plot
p3 <- ggplot(alpha_df, aes(x = Continent, y = Richness, color = Continent)) +
  geom_boxplot(
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.8
  ) +
  geom_jitter(
    width = 0.2,
    size = 1.5,
    alpha = 0.8
  ) +
  scale_color_manual(values = my_colors) +
  theme_classic() +
  labs(
    title = "Richness",
    x = "Continent",
    y = "Richness",
    color = "Continent"
  ) +
  stat_compare_means(
    comparisons = comparisons_list,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = FALSE,
    step.increase = 0.06
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  ) +
  coord_cartesian(
    ylim = c(
      min(alpha_df$Richness, na.rm = TRUE),
      max(alpha_df$Richness, na.rm = TRUE) * 1.35
    )
  )

p3

ggsave("/Users/yqhhh/Desktop/alpha_richness_boxplot.pdf", p3, width = 6, height = 5)

# Between-group test
if ("Continent" %in% colnames(alpha_df)) {
  kw_shannon <- kruskal.test(Shannon ~ Continent, data = alpha_df)
  kw_simpson <- kruskal.test(Simpson ~ Continent, data = alpha_df)
  kw_richness <- kruskal.test(Richness ~ Continent, data = alpha_df)
  
  capture.output(kw_shannon, file = "3.30species_adversity/kruskal_shannon.txt")
  capture.output(kw_simpson, file = "3.30species_adversity/kruskal_simpson.txt")
  capture.output(kw_richness, file = "3.30species_adversity/kruskal_richness.txt")
}





# beta_diversity analysis -----------------------------------------------------------------

library(readr)
library(dplyr)
library(tibble)
library(vegan)
library(ggplot2)

# Reading the abundance matrix
mat_df <- read_tsv("/Users/yqhhh/Desktop/species_adversity/bracken_species_matrix.tsv", show_col_types = FALSE)

mat <- mat_df %>%
  column_to_rownames("name") %>%
  as.matrix()

# Convert to sample × species
mat <- t(mat)

# Relative abundance (recommended)
mat_rel <- decostand(mat, method = "total")

# Read metadata
meta <- read_excel("/Users/yqhhh/Desktop/metadata.xlsx")
meta <- meta %>% filter(Sample %in% rownames(mat_rel))
mat_rel <- mat_rel[meta$Sample, , drop = FALSE]

# Bray-Curtis distance
dist_mat <- vegdist(mat_rel, method = "bray")

# PCoA
pcoa <- cmdscale(dist_mat, eig = TRUE, k = 2)

pcoa_df <- data.frame(
  Sample = rownames(mat_rel),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2]
) %>%
  left_join(meta, by = "Sample")

write_tsv(pcoa_df, "/Users/yqhhh/Desktop/pcoa_results.tsv")

# Variance Explained Proportion
var_explained <- round(100 * pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]), 2)

ad <- adonis2(dist_mat ~ Continent, data = meta, permutations = 999)

# Variance Explained Proportion
r2 <- ad$R2[1]
pval <- ad$`Pr(>F)`[1]


label_text <- paste0(
  "PERMANOVA\n",
  "R² = ", round(r2, 3), "\n",
  "p = ", signif(pval, 3)
)



my_colors <- c(
  "Asia" = "#ffa88d",
  "Europe" = "#FFD84F",
  "North America" = "#36CFFD",
  "Africa" = "#5CBC71"
)





p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Continent)) +
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  stat_ellipse(
    aes(group = Continent, color = Continent), 
    linetype = 2,
    linewidth = 0.8
  ) +
  scale_color_manual(values = my_colors) +     
  theme_classic() +
  labs(
    title = "PCoA based on Bray-Curtis",
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)"),
    color = "Continent"
  ) +
  annotate(
    "text",
    x = max(pcoa_df$PC1),
    y = max(pcoa_df$PC2),
    label = label_text,
    hjust = 0.2,
    vjust = 1,
    size = 4
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )



p


ggsave("/Users/yqhhh/Desktop/beta_pcoa_bray.pdf", p, width = 6, height = 5)

# PERMANOVA
if ("Continent" %in% colnames(meta)) {
  ad <- adonis2(dist_mat ~ Continent, data = meta, permutations = 999)
  capture.output(ad, file = "3.30species_adversity/permanova_continent.txt")
}