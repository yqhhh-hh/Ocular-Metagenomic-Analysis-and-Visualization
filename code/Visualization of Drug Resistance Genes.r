# 耐药基因的相关可视化 --------------------------------------------------------------


library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(circlize)
library(ComplexHeatmap)
library(grid)

# 读入数据 ---------------------------

input_file <- "/Users/yqhhh/Desktop/all_merged_AMR_with_species_time_continent.tsv"

amr_df <- readr::read_tsv(input_file, show_col_types = FALSE)
names(amr_df) <- trimws(names(amr_df))

required_cols <- c("AMR class", "Time", "Continent", "Strain", "Pathogen")
missing_cols <- setdiff(required_cols, names(amr_df))

if (length(missing_cols) > 0) {
  stop(paste0("The following are missing：", paste(missing_cols, collapse = ", ")))
}

# Data Cleaning ---------------------------

plot_df <- amr_df |>
  transmute(
    amr_class = as.character(`AMR class`),
    time = as.character(Time),
    continent = as.character(Continent),
    strain = as.character(Strain),
    pathogen = tolower(str_squish(as.character(Pathogen)))
  ) |>
  mutate(
    amr_class = str_squish(amr_class),
    time = str_squish(time),
    continent = str_squish(continent),
    strain = str_squish(strain),
    pathogen = if_else(pathogen == "yes", "yes", "no")
  ) |>
  filter(
    !is.na(amr_class), amr_class != "",
    !is.na(time), time != "",
    !is.na(continent), continent != "",
    !is.na(strain), strain != ""
  ) |>
  mutate(
    time_num = suppressWarnings(as.numeric(time)),
    col_id = paste(time, continent, sep = "___")
  )

# The top 20 strains are retained, and the rest are Others. ---------------------------

top_n_strain <- 20

top_strains <- plot_df |>
  count(strain, sort = TRUE) |>
  slice_head(n = top_n_strain) |>
  pull(strain)

plot_df <- plot_df |>
  mutate(
    strain = if_else(strain %in% top_strains, strain, "Others")
  )

# The top 10 AMR classes are retained, and the rest are classified as Others. ---------------------------

top_n_class <- 10

class_count_raw <- plot_df |>
  count(amr_class, sort = TRUE)

top_classes <- class_count_raw |>
  filter(amr_class != "Others") |>
  slice_head(n = top_n_class) |>
  pull(amr_class)

plot_df <- plot_df |>
  mutate(
    amr_class = if_else(amr_class %in% top_classes, amr_class, "Others")
  )

# Strain order: top 20 sorted by abundance, others last.---------------------------

strain_order <- plot_df |>
  count(strain, sort = TRUE) |>
  mutate(
    sort_key = if_else(strain == "Others", Inf, row_number())
  ) |>
  arrange(sort_key) |>
  pull(strain)

strain_order <- c(setdiff(strain_order, "Others"), "Others")
strain_order <- unique(strain_order)

# Pathogenicity markers: Mark "yes" whenever "yes" appears; "others" is always marked "no". ---------------------------

strain_pathogen <- plot_df |>
  group_by(strain) |>
  summarise(
    pathogen_flag = if_else(any(pathogen == "yes") & strain[1] != "Others", "yes", "no"),
    .groups = "drop"
  ) |>
  mutate(
    strain = factor(strain, levels = strain_order)
  ) |>
  arrange(strain)

# Column order: Time × Continent---------------------------

col_info <- plot_df |>
  distinct(col_id, time, time_num, continent) |>
  arrange(time_num, time, continent)

col_order <- col_info$col_id

# Heatmap Matrix: Specific Quantity ---------------------------

heat_df <- plot_df |>
  count(strain, col_id, name = "n") |>
  complete(
    strain = strain_order,
    col_id = col_order,
    fill = list(n = 0)
  )

heat_mat <- heat_df |>
  pivot_wider(names_from = col_id, values_from = n) |>
  as.data.frame()

rownames(heat_mat) <- heat_mat$strain
heat_mat$strain <- NULL
heat_mat <- as.matrix(heat_mat)

# Total number on the left ---------------------------

row_count <- plot_df |>
  count(strain) |>
  mutate(strain = factor(strain, levels = strain_order)) |>
  arrange(strain) |>
  pull(n)

names(row_count) <- strain_order

# Right side AMR class percentage matrix ---------------------------

class_df <- plot_df |>
  count(strain, amr_class, name = "n") |>
  group_by(strain) |>
  mutate(
    pct = n / sum(n)
  ) |>
  ungroup()

class_mat <- class_df |>
  select(strain, amr_class, pct) |>
  pivot_wider(
    names_from = amr_class,
    values_from = pct,
    values_fill = 0
  ) |>
  as.data.frame()

rownames(class_mat) <- class_mat$strain
class_mat$strain <- NULL
class_mat <- as.matrix(class_mat)

# Adjust the order of the class column: after top 10, other columns. ---------------------------

class_order <- plot_df |>
  count(amr_class, sort = TRUE) |>
  pull(amr_class)

class_order <- c(setdiff(class_order, "Others"), "Others")
class_order <- unique(class_order)

class_mat <- class_mat[, class_order[class_order %in% colnames(class_mat)], drop = FALSE]

# Color scheme ---------------------------

max_count <- max(heat_mat, na.rm = TRUE)

col_fun <- circlize::colorRamp2(
  c(0, 1, max(2, round(max_count / 4)), max(2, round(max_count / 2)), max_count),
  c("#FFFFFF", "#FEE0D2", "#FC9272", "#FB6A4A", "#CB181D")
)

time_levels <- unique(as.character(col_info$time))
time_colors <- structure(
  grDevices::hcl.colors(length(time_levels), "YlOrRd"),
  names = time_levels
)

continent_levels <- unique(as.character(col_info$continent))
continent_colors <- c(
  "Asia" = "#ffa88d",
  "Europe" = "#FFD84F",
  "North America" = "#36CFFD",
  "Africa" = "#5CBC71")

missing_cont <- setdiff(continent_levels, names(continent_colors))
if (length(missing_cont) > 0) {
  extra_cols <- structure(
    grDevices::hcl.colors(length(missing_cont), "Set 2"),
    names = missing_cont
  )
  continent_colors <- c(continent_colors, extra_cols)
}
continent_colors <- continent_colors[continent_levels]

class_colors <- c(
  "Multidrug" = "#E64B35",           
  "Aminoglycoside" = "#4DBBD5",      
  "Penicillin_beta-Lactam" = "#3C5488", 
  "Diaminopyrimidine" = "#00A087",   
  "Tetracycline" = "#F39B7F",         
  "Phenicol" = "#8491B4",            
  "Fluoroquinolone" = "#91D1C2",     
  "Macrolide" = "#DC0000",            
  "Peptide" = "#7E6148",            
  "Lincosamide" = "#A1C181",
  "Others" = "#B09C85"            
)



pathogen_map <- strain_pathogen |>
  mutate(
    strain = as.character(strain)
  ) |>
  select(strain, pathogen_flag)

row_label_col <- pathogen_map$pathogen_flag[
  match(rownames(heat_mat), pathogen_map$strain)
]

row_label_col <- ifelse(row_label_col == "yes", "#d7301f", "#222222")

# Top two lines of comments ---------------------------

time_anno <- as.character(col_info$time)
continent_anno <- as.character(col_info$continent)

top_ha <- HeatmapAnnotation(
  Time = time_anno,
  Continent = continent_anno,
  col = list(
    Time = time_colors,
    Continent = continent_colors
  ),
  simple_anno_size = unit(4, "mm"),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)

# Left-side total bar chart ---------------------------

left_counts <- row_count[rownames(heat_mat)]

left_ha <- rowAnnotation(
  Counts = anno_barplot(
    left_counts,
    gp = gpar(fill = "#9ECAE1", col = NA),
    border = FALSE,
    width = unit(30, "mm"),
    axis_param = list(
      side = "bottom",
      gp = gpar(fontsize = 8)
    )
  ),
  N = anno_text(
    left_counts,
    just = "left",
    location = 0,
    gp = gpar(fontsize = 8, col = "#444444"),
    width = unit(12, "mm")
  ),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)

# Right side AMR class percentage stacked bar chart ---------------------------

right_ha <- rowAnnotation(
  `AMR class (%)` = anno_barplot(
    class_mat[rownames(heat_mat), , drop = FALSE],
    gp = gpar(fill = class_colors, col = "black", lwd = 0.4),
    bar_width = 0.8,
    border = FALSE,
    beside = FALSE,
    baseline = 0,
    axis_param = list(
      at = c(0, 0.5, 1),
      labels = c("0%", "50%", "100%")
    ),
    width = unit(42, "mm")
  ),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)

# Main heatmap ---------------------------

ht <- Heatmap(
  heat_mat,
  name = "Count",
  col = col_fun,
  top_annotation = top_ha,
  left_annotation = left_ha,
  right_annotation = right_ha,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  width = unit(90, "mm"), 
  row_names_side = "left",
  row_names_gp = gpar(
    fontsize = 10,
    fontface = "italic",
    col = row_label_col
  ),
  rect_gp = gpar(col = "#BDBDBD", lwd = 0.8),   # 灰色边框
  heatmap_legend_param = list(
    title = "AMR count"
  )
)

amr_class_lgd <- Legend(
  title = "AMR class",
  labels = names(class_colors),
  legend_gp = gpar(fill = class_colors, col = NA),
  ncol = 1
)

# preview ---------------------------

grid::grid.newpage()

ht_drawn <- ComplexHeatmap::draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "bottom",
  annotation_legend_list = list(amr_class_lgd),
  merge_legends = FALSE,
  padding = grid::unit(c(6, 6, 16, 6), "mm")
)


pdf("amr_heatmap_red_3.31.pdf", width = 14, height = 9, useDingbats = FALSE)

ht_drawn <- ComplexHeatmap::draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "bottom",
  annotation_legend_list = list(amr_class_lgd),
  merge_legends = FALSE,
  padding = grid::unit(c(6, 6, 16, 6), "mm"),
  newpage = TRUE
)

dev.off()