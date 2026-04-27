# KEGG Data Processing and Visualization --------------------------------

library(tidyverse)
library(stringr)
library(RColorBrewer)

# Read the source file
df <- read.table(
  "/Users/yqhhh/Desktop/egg_vis/KO_relative_abundance_with_level_info_cleaned.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

df <- df %>%
  mutate(
    ko_id = str_trim(ko_id),
    level1 = str_trim(level1),
    level2 = str_trim(level2),
    pathway_name = str_trim(pathway_name)
  )

df_long <- df %>%
  pivot_longer(
    cols = 5:ncol(df),
    names_to = "sample",
    values_to = "abundance"
  )

# Filter Invalid Categories
plot_df <- df_long %>%
  filter(
    !is.na(level1), !is.na(level2),
    level1 != "", level2 != "",
    level1 != "Brite Hierarchies",
    level1 != "Not Included in Pathway or Brite"
  ) %>%
  group_by(level1, level2) %>%
  summarise(
    Abundance = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  )

class_order <- plot_df %>%
  group_by(level1) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(level1)

plot_df <- plot_df %>%
  mutate(level1 = factor(level1, levels = class_order)) %>%
  arrange(level1, desc(Abundance))

plot_df$level2 <- factor(plot_df$level2, levels = plot_df$level2)

n_class <- length(unique(plot_df$level1))
pal <- brewer.pal(max(3, min(8, n_class)), "Set2")

# Draw
p <- ggplot(plot_df, aes(x = level2, y = Abundance, fill = level1)) +
  geom_col(color = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = pal) +
  labs(
    x = "Subclass",
    y = "Total abundance",
    fill = "Class"
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(face = "bold", size = 15),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(t = 20, r = 20, b = 100, l = 40)
  )

print(p)

# Save
ggsave(
  "/Users/yqhhh/Desktop/egg_vis/KEGG_pathway_classification_final_sorted-1.pdf",
  p, width = 16, height = 10
)

ggsave(
  "/Users/yqhhh/Desktop/egg_vis/KEGG_pathway_classification_final_sorted.png",
  p, width = 22, height = 10, dpi = 300
)





# COG Data Processing and Visualization--------------------------------------

library(tidyverse)
library(RColorBrewer)

# Reading data
cog <- read.table(
  "/Users/yqhhh/Desktop/egg_vis/COG_relative_abundance.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

mapping <- read.csv(
  "/Users/yqhhh/Desktop/egg_vis/COG_function_mapping.csv",
  stringsAsFactors = FALSE
)

colnames(mapping) <- c("COG_code", "COG_function")


cog_long <- cog %>%
  pivot_longer(
    cols = -sample,
    names_to = "COG_code",
    values_to = "abundance"
  )

# Calculate the total relative abundance of each COG letter.
plot_df <- cog_long %>%
  group_by(COG_code) %>%
  summarise(
    Total_abundance = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(mapping, by = "COG_code") %>%
  arrange(desc(Total_abundance))


plot_df$COG_code <- factor(plot_df$COG_code, levels = plot_df$COG_code)

# Generate Colors
n_cog <- nrow(plot_df)


if (n_cog <= 12) {
  my_cols <- brewer.pal(n_cog, "Set3")
} else {
  my_cols <- colorRampPalette(brewer.pal(12, "Set3"))(n_cog)
}

names(my_cols) <- as.character(plot_df$COG_code)


plot_df$legend_label <- paste0(plot_df$COG_code, ": ", plot_df$COG_function)
legend_labels <- plot_df$legend_label
names(legend_labels) <- plot_df$COG_code

# Draw
p <- ggplot(plot_df, aes(x = COG_code, y = Total_abundance, fill = COG_code)) +
  geom_col(color = "grey30", linewidth = 0.4) +
  scale_fill_manual(
    values = my_cols,
    labels = legend_labels,
    name = "COG category"
  ) +
  labs(
    x = "COG category",
    y = "Total relative abundance"
  ) +
  theme_bw(base_size = 14) +
  scale_y_continuous(expand = c(0.01, 0))+
  theme(
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(t = 20, r = 30, b = 20, l = 20)
  )


print(p)

# Save
ggsave(
  "/Users/yqhhh/Desktop/4.18egg_vis/COG_barplot.pdf",
  p, width = 14, height = 8
)

ggsave(
  "/Users/yqhhh/Desktop/4.18egg_vis/COG_barplot.png",
  p, width = 14, height = 8, dpi = 300
)






# CAZy Data Processing and Visualization --------------------------------------------

library(tidyverse)
library(stringr)

# Read the CAZy relative abundance matrix.
cazy <- read.table(
  "/Users/yqhhh/Desktop/egg_vis/CAZy_relative_abundance.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

# Convert to Long Format
cazy_long <- cazy %>%
  pivot_longer(
    cols = -sample,
    names_to = "CAZy_family",
    values_to = "abundance"
  )

# Classified into six major categories
cazy_long <- cazy_long %>%
  mutate(
    CAZy_class = case_when(
      str_detect(CAZy_family, "^GH")  ~ "GHs",
      str_detect(CAZy_family, "^GT")  ~ "GTs",
      str_detect(CAZy_family, "^PL")  ~ "PLs",
      str_detect(CAZy_family, "^CE")  ~ "CEs",
      str_detect(CAZy_family, "^CBM") ~ "CBMs",
      str_detect(CAZy_family, "^AA")  ~ "AAs",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(CAZy_class))

# Aggregate Overall Relative Abundance
plot_df <- cazy_long %>%
  group_by(CAZy_class) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    percentage = total_abundance / sum(total_abundance) * 100
  ) %>%
  arrange(desc(percentage))

# Full-name labels.
plot_df <- plot_df %>%
  mutate(
    CAZy_label = case_when(
      CAZy_class == "GHs"  ~ "GHs: Glycoside Hydrolases",
      CAZy_class == "GTs"  ~ "GTs: GlycosylTransferases",
      CAZy_class == "PLs"  ~ "PLs: Polysaccharide Lyases",
      CAZy_class == "CEs"  ~ "CEs: Carbohydrate Esterases",
      CAZy_class == "CBMs" ~ "CBMs: Carbohydrate-Binding Modules",
      CAZy_class == "AAs"  ~ "AAs: Auxiliary Activities"
    )
  )

# Manually Calculate Label Position
plot_df <- plot_df %>%
  mutate(
    ymax = cumsum(percentage),
    ymin = lag(ymax, default = 0),
    label_pos = (ymax + ymin) / 2
  )

# Color
my_cols <- c(
  "GTs: GlycosylTransferases" = "#E6E49C",
  "GHs: Glycoside Hydrolases" = "#85C4BB",
  "CEs: Carbohydrate Esterases" = "#F07D6D",
  "CBMs: Carbohydrate-Binding Modules" = "#7AA7C7",
  "PLs: Polysaccharide Lyases" = "#B8B3D9",
  "AAs: Auxiliary Activities" = "#F0A857"
)

# Donut Chart
p <- ggplot(plot_df, aes(x = 2, y = percentage, fill = CAZy_label)) +
  geom_col(color = "white", linewidth = 0.2, width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.8) +   
  geom_text(
    aes(y = label_pos, label = paste0(round(percentage, 2), "%")),
    x = 2.55,
    size = 5
  ) +
  
  scale_fill_manual(values = my_cols) +
  theme_void(base_size = 14) +
  labs(fill = "CAZy class") +
  theme(
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 12)
  )

print(p)

# Save
ggsave(
  "/Users/yqhhh/Desktop/egg_vis/CAZy_class_donut_outside_labels.pdf",
  p, width = 10, height = 8
)

ggsave(
  "/Users/yqhhh/Desktop/egg_vis/CAZy_class_donut_outside_labels.png",
  p, width = 10, height = 8, dpi = 300
)





# EC Data Processing and Visualization -------------------------------------------------------------

library(tidyverse)
library(stringr)
library(RColorBrewer)

# Read in the EC relative abundance matrix.
ec <- read.table(
  "/Users/yqhhh/Desktop/egg_vis/EC_relative_abundance.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

# Convert to Long Format
ec_long <- ec %>%
  pivot_longer(
    cols = -sample,
    names_to = "EC",
    values_to = "abundance"
  )

# Extract EC Level 1 Category (First Digit)
ec_long <- ec_long %>%
  mutate(
    EC_class_num = str_extract(EC, "^[0-9]+"),
    EC_class = case_when(
      EC_class_num == "1" ~ "1. Oxidoreductases",
      EC_class_num == "2" ~ "2. Transferases",
      EC_class_num == "3" ~ "3. Hydrolases",
      EC_class_num == "4" ~ "4. Lyases",
      EC_class_num == "5" ~ "5. Isomerases",
      EC_class_num == "6" ~ "6. Ligases",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(EC_class))

# Summary of Total Relative Abundance for Six Major Categories
plot_df <- ec_long %>%
  group_by(EC_class) %>%
  summarise(
    Total_abundance = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_abundance))


plot_df$EC_class <- factor(
  plot_df$EC_class,
  levels = c(
    "1. Oxidoreductases",
    "2. Transferases",
    "3. Hydrolases",
    "4. Lyases",
    "5. Isomerases",
    "6. Ligases"
  )
)

my_cols <- c(
  "1. Oxidoreductases" = "#8DD3C7",
  "2. Transferases"    = "#FFFFB3",
  "3. Hydrolases"      = "#BEBADA",
  "4. Lyases"          = "#FB8072",
  "5. Isomerases"      = "#80B1D3",
  "6. Ligases"         = "#FDB462"
)

# Draw
p <- ggplot(plot_df, aes(x = EC_class, y = Total_abundance, fill = EC_class)) +
  geom_col(color = "grey30", linewidth = 0.3, width = 0.75) +
  scale_fill_manual(values = my_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(
    x = "EC class",
    y = "Total relative abundance",
    fill = "EC class"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 11),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

print(p)

# Save
ggsave(
  "/Users/yqhhh/Desktop/egg_vis/EC_6class_barplot.pdf",
  p, width = 10, height = 6
)

ggsave(
  "/Users/yqhhh/Desktop/egg_vis/EC_6class_barplot.png",
  p, width = 10, height = 6, dpi = 300
)
