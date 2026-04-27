# Species abundance results visualization-----------------------------------------------

# =========================
# 1.Load package
# =========================
library(readr)
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(RColorBrewer)


world_file <- "/Users/yqhhh/Desktop/bracken_species_World_processed.tsv"
asia_file <- "/Users/yqhhh/Desktop/bracken_species_Asia_processed.tsv"
europe_file <- "/Users/yqhhh/Desktop/bracken_species_Europe_processed.tsv"
north_america_file <- "/Users/yqhhh/Desktop/bracken_species_North_America_processed.tsv"
africa_file <- "/Users/yqhhh/Desktop/bracken_species_Africa_processed.tsv"

output_dir <- "/Users/yqhhh/Desktop"

read_species_data <- function(file_path) {
  df <- read_tsv(file_path, show_col_types = FALSE)
  
  cat("\n=============================\n")
  cat("Reading file:\n", file_path, "\n")
  cat("Column names:\n")
  print(names(df))
  
  if (!"Species" %in% names(df)) {
    stop(paste0("The Species column was not found in the file.", file_path))
  }
  
  if (!"Sum" %in% names(df)) {
    stop(paste0("The Sum column was not found in the file：", file_path))
  }
  
  df |>
    slice(-1) |>
    mutate(
      Species = trimws(as.character(Species)),
      Sum = suppressWarnings(as.numeric(Sum))
    ) |>
    filter(
      !is.na(Species),
      Species != "",
      !Species %in% c("Homo sapiens", "Homo_sapiens"),
      !is.na(Sum),
      Sum > 0
    )
}

world_df <- read_species_data(world_file)
asia_df <- read_species_data(asia_file)
europe_df <- read_species_data(europe_file)
north_america_df <- read_species_data(north_america_file)
africa_df <- read_species_data(africa_file)

top_n <- 20

world_top_species <- world_df |>
  group_by(Species) |>
  summarise(
    Sum = sum(Sum, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(Sum)) |>
  slice_head(n = top_n) |>
  pull(Species)

cat("\nTop species from World:\n")
print(world_top_species)

summarise_species <- function(data, dataset_name, top_species) {
  data |>
    mutate(
      category = if_else(Species %in% top_species, Species, "Others")
    ) |>
    group_by(category) |>
    summarise(
      Sum = sum(Sum, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(dataset = dataset_name)
}

plot_df <- bind_rows(
  summarise_species(world_df, "World", world_top_species),
  summarise_species(asia_df, "Asia", world_top_species),
  summarise_species(europe_df, "Europe", world_top_species),
  summarise_species(north_america_df, "North America", world_top_species),
  summarise_species(africa_df, "Africa", world_top_species)
)

category_levels <- c(world_top_species, "Others")
dataset_levels <- c("World", "Asia", "Europe", "North America", "Africa")

plot_df <- plot_df |>
  mutate(
    category = factor(category, levels = category_levels),
    dataset = factor(dataset, levels = dataset_levels)
  ) |>
  complete(dataset, category, fill = list(Sum = 0)) |>
  group_by(dataset) |>
  mutate(
    total_sum = sum(Sum, na.rm = TRUE),
    proportion = if_else(total_sum > 0, Sum / total_sum, 0),
    label = if_else(
      proportion >= 0.03,
      percent(proportion, accuracy = 0.1),
      ""
    )
  ) |>
  ungroup()

cat("\nCheck proportion sum for each dataset:\n")
plot_df |>
  group_by(dataset) |>
  summarise(prop_sum = sum(proportion)) |>
  print()

n_cat <- length(category_levels)
base_cols <- brewer.pal(min(12, max(3, n_cat - 1)), "Set3")

if (n_cat > length(base_cols)) {
  palette_values <- colorRampPalette(base_cols)(n_cat)
} else {
  palette_values <- base_cols[seq_len(n_cat)]
}

palette_values <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#A6CEE3", "#FBB4AE", "#B2DF8A",
  "#CAB2D6", "#FDBF6F", "#E0ECF4", "#F1EEF6", "#E5F5E0",
  "#F2E5D7" ) 

palette_values <- palette_values[seq_len(n_cat)]

p <- ggplot(plot_df, aes(x = dataset, y = proportion, fill = category)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  scale_fill_manual(values = palette_values, drop = FALSE) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Species composition across World and regions",
    x = NULL,
    y = "Relative abundance (%)",
    fill = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p)

ggsave(
  filename = file.path(output_dir, "species_world_region_stacked_barplot_fixed.png"),
  plot = p,
  width = 11,
  height = 6.5,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "species_world_region_stacked_barplot_fixed.pdf"),
  plot = p,
  width = 11,
  height = 6.5
)



# bracken Visualization-----Genus --------Other classification levels follow the same principle

library(readr)
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)

# =========================
# =========================
world_file <- "/Users/yqhhh/Desktop/bracken_species_World_processed.tsv"
asia_file <- "/Users/yqhhh/Desktop/bracken_species_Asia_processed.tsv"
europe_file <- "/Users/yqhhh/Desktop/bracken_species_Europe_processed.tsv"
north_america_file <- "/Users/yqhhh/Desktop/bracken_species_North_America_processed.tsv"
africa_file <- "/Users/yqhhh/Desktop/bracken_species_Africa_processed.tsv"

output_dir <- "/Users/yqhhh/Desktop/"

# =========================
# =========================
read_genus_data <- function(file_path) {
  df <- read_tsv(file_path, show_col_types = FALSE)
  
  required_cols <- c("Genus", "Species", "Sum")
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "The file is missing the following:",
        paste(missing_cols, collapse = ", "),
        "\n file：",
        file_path
      )
    )
  }
  
  df |>
    slice(-1) |>
    mutate(
      Genus = trimws(as.character(Genus)),
      Species = trimws(as.character(Species)),
      Sum = suppressWarnings(as.numeric(Sum))
    ) |>
    filter(
      !is.na(Genus),
      Genus != "",
      !Genus %in% c("Homo", "Homo sapiens", "Homo_sapiens"),
      !is.na(Species),
      Species != "",
      !Species %in% c("Homo sapiens", "Homo_sapiens"),
      !is.na(Sum),
      Sum > 0
    )
}

# =========================
# Import data from various regions
# =========================
world_df <- read_genus_data(world_file)
asia_df <- read_genus_data(asia_file)
europe_df <- read_genus_data(europe_file)
north_america_df <- read_genus_data(north_america_file)
africa_df <- read_genus_data(africa_file)

# =========================
# First summarize by Genus
# Sum of sums for all Species under a Genus
# =========================
aggregate_to_genus <- function(data) {
  data |>
    group_by(Genus) |>
    summarise(
      Sum = sum(Sum, na.rm = TRUE),
      .groups = "drop"
    )
}

world_genus <- aggregate_to_genus(world_df)
asia_genus <- aggregate_to_genus(asia_df)
europe_genus <- aggregate_to_genus(europe_df)
north_america_genus <- aggregate_to_genus(north_america_df)
africa_genus <- aggregate_to_genus(africa_df)

# =========================
# Calculate the top 20 Genus using the World
# =========================
top_n <- 20

world_top_genus <- world_genus |>
  arrange(desc(Sum)) |>
  slice_head(n = top_n) |>
  pull(Genus)

print("Top genus from World:")
print(world_top_genus)

# =========================
# Each region is categorized according to the World's top 20 Genus.
# =========================
summarise_genus <- function(genus_data, dataset_name, top_genus) {
  genus_data |>
    mutate(
      category = if_else(Genus %in% top_genus, Genus, "Others")
    ) |>
    group_by(category) |>
    summarise(
      Sum = sum(Sum, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(dataset = dataset_name)
}

plot_df <- bind_rows(
  summarise_genus(world_genus, "World", world_top_genus),
  summarise_genus(asia_genus, "Asia", world_top_genus),
  summarise_genus(europe_genus, "Europe", world_top_genus),
  summarise_genus(north_america_genus, "North America", world_top_genus),
  summarise_genus(africa_genus, "Africa", world_top_genus)
)

# =========================
# Complete the missing categories and convert them to proportions
# =========================
category_levels <- c(world_top_genus, "Others")
dataset_levels <- c("World", "Asia", "Europe", "North America", "Africa")

plot_df <- plot_df |>
  mutate(
    category = factor(category, levels = category_levels),
    dataset = factor(dataset, levels = dataset_levels)
  ) |>
  complete(dataset, category, fill = list(Sum = 0)) |>
  group_by(dataset) |>
  mutate(
    total_sum = sum(Sum, na.rm = TRUE),
    proportion = if_else(total_sum > 0, Sum / total_sum, 0)
  ) |>
  ungroup() |>
  mutate(
    label = if_else(
      proportion >= 0.03,
      percent(proportion, accuracy = 0.1),
      ""
    )
  )

# Check if the ratio for each region is 1
print(
  plot_df |>
    group_by(dataset) |>
    summarise(prop_sum = sum(proportion), .groups = "drop")
)

# =========================
# Color scheme
# =========================
n_cat <- length(category_levels)

palette_values <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#A6CEE3", "#FBB4AE", "#B2DF8A",
  "#CAB2D6", "#FDBF6F", "#E0ECF4", "#F1EEF6", "#E5F5E0",
  "#F2E5D7"
)

palette_values <- palette_values[seq_len(n_cat)]

# =========================
# Drawing
# =========================
p <- ggplot(plot_df, aes(x = dataset, y = proportion, fill = category)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  scale_fill_manual(
    values = setNames(palette_values, category_levels),
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Genus composition across World and regions",
    x = NULL,
    y = "Relative abundance (%)",
    fill = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p)

# =========================
# Save picture
# =========================
ggsave(
  filename = file.path(output_dir, "genus_world_region_stacked_barplot.png"),
  plot = p,
  width = 11,
  height = 6.5,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "genus_world_region_stacked_barplot.pdf"),
  plot = p,
  width = 11,
  height = 6.5
)
