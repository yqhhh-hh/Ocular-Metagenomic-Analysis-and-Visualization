# MGE---Drug resistance genes combined with visualization --------------------------------

 # Load packages -----------------------------------------------------------
 
 library(readr)
 library(dplyr)
 library(tidyr)
 library(stringr)
 library(ggplot2)
 library(scales)
 library(ggtext)
 
 # Paths -------------------------------------------------------------------
 
 input_file <- "/Users/yqhhh/Desktop/AMR_with_MGE_plasmid_yesno.tsv"
 output_dir <- "/Users/yqhhh/Desktop"
 
 dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
 
 # Settings ----------------------------------------------------------------
 
 mge_cols <- c("Tns", "IS", "MITE", "ICE", "IME", "ComTn", "Plasmid")
 top_n <- 20
 
 # Helper functions --------------------------------------------------------
 
 clean_yes_no <- function(x) {
   x |>
     as.character() |>
     str_trim() |>
     str_to_title() |>
     recode(
       "Yes" = "Yes",
       "No" = "No",
       .default = "No"
     )
 }
 
 wrap_gene_label <- function(x, width = 24) {
   x |>
     stringr::str_replace_all("_", "_ ") |>
     stringr::str_wrap(width = width) |>
     stringr::str_replace_all("\n", "<br>")
 }
 
 # Read data ---------------------------------------------------------------
 
 dat <- read_tsv(input_file, show_col_types = FALSE)
 
 required_cols <- c("GENE", mge_cols)
 missing_cols <- setdiff(required_cols, names(dat))
 
 if (length(missing_cols) > 0) {
   stop(paste0("Missing Required List: ", paste(missing_cols, collapse = ", ")))
 }
 
 dat <- dat |>
   mutate(
     GENE = as.character(GENE),
     GENE = str_replace_all(GENE, "’", "'"),
     GENE = str_replace_all(GENE, "“|”", "\""),
     GENE = str_trim(GENE),
     GENE = na_if(GENE, ""),
     across(all_of(mge_cols), clean_yes_no)
   ) |>
   filter(!is.na(GENE))
 
 # Keep top 20 genes -------------------------------------------------------
 
 top_gene <- dat |>
   count(GENE, sort = TRUE, name = "gene_total") |>
   slice_head(n = top_n) |>
   pull(GENE)
 
 dat <- dat |>
   mutate(
     gene_group = if_else(GENE %in% top_gene, GENE, "Others")
   )
 
 # Classify MGE / No MGE ---------------------------------------------------
 
 dat <- dat |>
   mutate(
     mge_status = if_else(
       if_any(all_of(mge_cols), ~ .x == "Yes"),
       "MGE",
       "No MGE"
     )
   )
 
 # Summarise ---------------------------------------------------------------
 
 plot_dat <- dat |>
   count(gene_group, mge_status, name = "n") |>
   complete(
     gene_group,
     mge_status = c("MGE", "No MGE"),
     fill = list(n = 0)
   )
 
 gene_total <- dat |>
   count(gene_group, name = "gene_total")
 
 plot_dat <- plot_dat |>
   left_join(gene_total, by = "gene_group")
 
 # Order genes -------------------------------------------------------------
 
 gene_order_df <- gene_total |>
   mutate(is_others = gene_group == "Others") |>
   arrange(is_others, desc(gene_total), gene_group)
 
 gene_order <- gene_order_df$gene_group
 
 plot_dat <- plot_dat |>
   mutate(
     gene_group = factor(gene_group, levels = rev(gene_order)),
     mge_status = factor(mge_status, levels = c("No MGE", "MGE"))
   )
 
 # Labels ------------------------------------------------------------------
 
 gene_levels <- levels(plot_dat$gene_group)
 
 gene_palette <- setNames(
   c(
     "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628",
     "#F781BF", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
     "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77",
     "#D95F02", "#7570B3", "#66A61E", "#A6761D", "#E7298A",
     "#BDBDBD"
   ),
   gene_levels
 )
 
 gene_palette["Others"] <- "#BDBDBD"
 
 gene_label_df <- gene_total |>
   mutate(
     gene_name_wrapped = wrap_gene_label(gene_group, width = 24),
     label = paste0(
       "<span style='color:", gene_palette[gene_group], ";'>",
       gene_name_wrapped, " (", gene_total, ")",
       "</span>"
     )
   )
 
 gene_labels <- setNames(gene_label_df$label, gene_label_df$gene_group)
 
 # Build fill colors for dark/light version --------------------------------
 
 plot_dat <- plot_dat |>
   mutate(
     fill_group = paste(gene_group, mge_status, sep = "__")
   )
 
 fill_values <- c()
 
 for (g in gene_levels) {
   fill_values[paste0(g, "__MGE")] <- gene_palette[g]
   fill_values[paste0(g, "__No MGE")] <- alpha(gene_palette[g], 0.35)
 }
 
 # Percentage labels -------------------------------------------------------
 
 label_dat <- plot_dat |>
   filter(mge_status == "MGE") |>
   mutate(
     pct = if_else(gene_total > 0, n / gene_total, 0),
     pct_label = percent(pct, accuracy = 0.1),
     x_pos = gene_total + max(gene_total) * 0.02
   )
 
 # Plot --------------------------------------------------------------------
 
 p <- ggplot(
   plot_dat,
   aes(x = n, y = gene_group, fill = fill_group)
 ) +
   geom_col(width = 0.82, color = "white", linewidth = 0.25) +
   geom_text(
     data = label_dat,
     aes(x = x_pos, y = gene_group, label = pct_label),
     inherit.aes = FALSE,
     hjust = 0,
     size = 3.5
   ) +
   annotate(
     "text",
     x = max(plot_dat$n) * 0.85,
     y = length(levels(plot_dat$gene_group)) + 0.2,
     label = "MGE (dark)   No MGE (light)",
     hjust = 0,
     size = 5,
     fontface = "bold"
   ) +
   scale_fill_manual(values = fill_values, guide = "none") +
   scale_y_discrete(labels = gene_labels) +
   scale_x_continuous(
     expand = expansion(mult = c(0, 0.25))
   ) +
   coord_cartesian(clip = "off") +
   labs(
     title = "Top 20 genes with MGE status",
     x = "Count",
     y = NULL
   ) +
   theme_bw(base_size = 12) +
   theme(
     panel.grid.major.y = element_blank(),
     panel.grid.minor = element_blank(),
     legend.position = "none",
     plot.title = element_text(hjust = 0.5, face = "bold"),
     axis.title.y = element_blank(),
     axis.text.y = ggtext::element_markdown(size = 10),
     axis.text.x = element_text(size = 10)
   )
 
 # Preview -----------------------------------------------------------------
 
 print(p)
 
 # Save --------------------------------------------------------------------
 
 ggsave(
   filename = file.path(output_dir, "gene_top20_mge_binary_barplot_1.pdf"),
   plot = p,
   width = 7.5,
   height = 8,
   bg = "white"
 )
 
 

# Relationship between strain and MGE -----------------------------------------------------------

 # Load packages -----------------------------------------------------------
 
 library(readr)
 library(dplyr)
 library(tidyr)
 library(stringr)
 library(ggplot2)
 library(scales)
 library(ggtext)
 
 # Paths -------------------------------------------------------------------
 
 input_file <- "/Users/yqhhh/Desktop/AMR_with_MGE_plasmid_yesno.tsv"
 output_dir <- "/Users/yqhhh/Desktop"
 
 dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
 
 # Settings ----------------------------------------------------------------
 
 mge_cols <- c("Tns", "IS", "MITE", "ICE", "IME", "ComTn", "Plasmid")
 top_n <- 20
 
 # Helper functions --------------------------------------------------------
 
 clean_yes_no <- function(x) {
   x |>
     as.character() |>
     str_trim() |>
     str_to_title() |>
     recode(
       "Yes" = "Yes",
       "No" = "No",
       .default = "No"
     )
 }
 
 wrap_strain_label <- function(x, width = 64) {
   x |>
     stringr::str_replace_all("_", "_ ") |>
     stringr::str_wrap(width = width) |>
     stringr::str_replace_all("\n", "<br>")
 }
 
 # Read data ---------------------------------------------------------------
 
 dat <- read_tsv(input_file, show_col_types = FALSE)
 
 required_cols <- c("Strain", mge_cols)
 missing_cols <- setdiff(required_cols, names(dat))
 
 if (length(missing_cols) > 0) {
   stop(paste0("缺少必须列: ", paste(missing_cols, collapse = ", ")))
 }
 
 dat <- dat |>
   mutate(
     Strain = as.character(Strain),
     Strain = str_replace_all(Strain, "’", "'"),
     Strain = str_replace_all(Strain, "“|”", "\""),
     Strain = str_trim(Strain),
     Strain = na_if(Strain, ""),
     across(all_of(mge_cols), clean_yes_no)
   ) |>
   filter(!is.na(Strain))
 
 # Keep top 20 strains -----------------------------------------------------
 
 top_strain <- dat |>
   count(Strain, sort = TRUE, name = "strain_total") |>
   slice_head(n = top_n) |>
   pull(Strain)
 
 dat <- dat |>
   mutate(
     strain_group = if_else(Strain %in% top_strain, Strain, "Others")
   )
 
 # Classify MGE / No MGE ---------------------------------------------------
 
 dat <- dat |>
   mutate(
     mge_status = if_else(
       if_any(all_of(mge_cols), ~ .x == "Yes"),
       "MGE",
       "No MGE"
     )
   )
 
 # Summarise ---------------------------------------------------------------
 
 plot_dat <- dat |>
   count(strain_group, mge_status, name = "n") |>
   complete(
     strain_group,
     mge_status = c("MGE", "No MGE"),
     fill = list(n = 0)
   )
 
 strain_total <- dat |>
   count(strain_group, name = "strain_total")
 
 plot_dat <- plot_dat |>
   left_join(strain_total, by = "strain_group")
 
 # Order strains -----------------------------------------------------------
 
 strain_order_df <- strain_total |>
   mutate(is_others = strain_group == "Others") |>
   arrange(is_others, desc(strain_total), strain_group)
 
 strain_order <- strain_order_df$strain_group
 
 plot_dat <- plot_dat |>
   mutate(
     strain_group = factor(strain_group, levels = rev(strain_order)),
     mge_status = factor(mge_status, levels = c("No MGE", "MGE"))
   )
 
 # Labels ------------------------------------------------------------------
 
 strain_levels <- levels(plot_dat$strain_group)
 
 strain_palette <- setNames(
   c(
     "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628",
     "#F781BF", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
     "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77",
     "#D95F02", "#7570B3", "#66A61E", "#A6761D", "#E7298A",
     "#FF724D"
   ),
   strain_levels
 )
 
 strain_palette["Others"] <- "#BDBDBD"
 
 strain_label_df <- strain_total |>
   mutate(
     strain_name_wrapped = wrap_strain_label(strain_group, width = 24),
     label = paste0(
       "<span style='color:", strain_palette[strain_group], ";'>",
       strain_name_wrapped, " (", strain_total, ")",
       "</span>"
     )
   )
 
 strain_labels <- setNames(strain_label_df$label, strain_label_df$strain_group)
 
 # Build fill colors for dark/light version --------------------------------
 
 plot_dat <- plot_dat |>
   mutate(
     fill_group = paste(strain_group, mge_status, sep = "__")
   )
 
 fill_values <- c()
 
 for (s in strain_levels) {
   fill_values[paste0(s, "__MGE")] <- strain_palette[s]
   fill_values[paste0(s, "__No MGE")] <- alpha(strain_palette[s], 0.35)
 }
 
 # Percentage labels -------------------------------------------------------
 
 label_dat <- plot_dat |>
   filter(mge_status == "MGE") |>
   mutate(
     pct = if_else(strain_total > 0, n / strain_total, 0),
     pct_label = percent(pct, accuracy = 0.1),
     x_pos = strain_total + max(strain_total) * 0.02
   )
 
 # Plot --------------------------------------------------------------------
 
 p <- ggplot(
   plot_dat,
   aes(x = n, y = strain_group, fill = fill_group)
 ) +
   geom_col(width = 0.82, color = "white", linewidth = 0.25) +
   geom_text(
     data = label_dat,
     aes(x = x_pos, y = strain_group, label = pct_label),
     inherit.aes = FALSE,
     hjust = 0,
     size = 3.5
   ) +
   scale_fill_manual(values = fill_values, guide = "none") +
   scale_y_discrete(labels = strain_labels) +
   scale_x_continuous(
     expand = expansion(mult = c(0, 0.1))
   ) +
   coord_cartesian(clip = "off") +
   labs(
     title = "Top 20 strains with MGE status",
     x = "Count",
     y = NULL
   ) +
   theme_bw(base_size = 12) +
   theme(
     panel.grid.major.y = element_blank(),
     panel.grid.minor = element_blank(),
     legend.position = "none",
     plot.title = element_text(hjust = 0.5, face = "bold"),
     axis.title.y = element_blank(),
     axis.text.y = ggtext::element_markdown(size = 10),
     axis.text.x = element_text(size = 10)
   )
 
 # Preview -----------------------------------------------------------------
 
 print(p)
 
 # Save --------------------------------------------------------------------
 
 ggsave(
   filename = file.path(output_dir, "strain_top20_mge_binary_barplot_1.pdf"),
   plot = p,
   width = 7.5,
   height = 8,
   bg = "white"
 )
 
 
 

# The relationship between different continents and MGE --------------------------------------
 

 # Load packages -----------------------------------------------------------
 
 library(readr)
 library(dplyr)
 library(tidyr)
 library(stringr)
 library(ggplot2)
 library(scales)
 
 # Paths -------------------------------------------------------------------
 
 input_file <- "/Users/yqhhh/Desktop/AMR_with_MGE_3.31_plasmid_yesno.tsv"
 output_dir <- "/Users/yqhhh/Desktop"
 
 dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
 
 # Settings ----------------------------------------------------------------
 
 mge_cols <- c("Tns", "IS", "MITE", "ICE", "IME", "ComTn", "Plasmid")
 
 # Helper functions --------------------------------------------------------
 
 clean_yes_no <- function(x) {
   x |>
     as.character() |>
     str_trim() |>
     str_to_title() |>
     recode(
       "Yes" = "Yes",
       "No" = "No",
       .default = "No"
     )
 }
 
 clean_text <- function(x) {
   x |>
     as.character() |>
     str_replace_all("’", "'") |>
     str_replace_all("“|”", "\"") |>
     str_trim() |>
     na_if("")
 }
 
 # Read data ---------------------------------------------------------------
 
 dat <- read_tsv(input_file, show_col_types = FALSE)
 
 required_cols <- c("Continent", mge_cols)
 missing_cols <- setdiff(required_cols, names(dat))
 
 if (length(missing_cols) > 0) {
   stop(paste0("Missing Required List: ", paste(missing_cols, collapse = ", ")))
 }
 
 dat <- dat |>
   mutate(
     Continent = clean_text(Continent),
     across(all_of(mge_cols), clean_yes_no)
   ) |>
   filter(!is.na(Continent))
 
 # Classify MGE / No MGE ---------------------------------------------------
 
 dat <- dat |>
   mutate(
     mge_status = if_else(
       if_any(all_of(mge_cols), ~ .x == "Yes"),
       "MGE",
       "No MGE"
     )
   )
 
 # Summarise ---------------------------------------------------------------
 
 plot_dat <- dat |>
   count(Continent, mge_status, name = "n") |>
   complete(
     Continent,
     mge_status = c("MGE", "No MGE"),
     fill = list(n = 0)
   )
 
 continent_total <- dat |>
   count(Continent, name = "continent_total")
 
 plot_dat <- plot_dat |>
   left_join(continent_total, by = "Continent")
 
 # Order continents --------------------------------------------------------
 
 continent_order <- continent_total |>
   arrange(desc(continent_total), Continent) |>
   pull(Continent)
 
 plot_dat <- plot_dat |>
   mutate(
     Continent = factor(Continent, levels = continent_order),
     mge_status = factor(mge_status, levels = c("No MGE", "MGE"))
   )
 
 # Colors ------------------------------------------------------------------
 
 continent_levels <- levels(plot_dat$Continent)
 
 continent_palette <- c(
   "Asia" = "#ffa88d",
   "Europe" = "#FFD84F",
   "North America" = "#36CFFD",
   "Africa" = "#5CBC71"
 )
 
 missing_continent <- setdiff(continent_levels, names(continent_palette))
 
 if (length(missing_continent) > 0) {
   extra_cols <- c("#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#7E6148")
   continent_palette[missing_continent] <- extra_cols[seq_along(missing_continent)]
 }
 
 continent_label_df <- continent_total |>
   mutate(
     label = paste0(Continent, " (", continent_total, ")")
   )
 
 continent_labels <- setNames(
   continent_label_df$label,
   continent_label_df$Continent
 )
 
 # Build fill colors for dark/light version --------------------------------
 
 plot_dat <- plot_dat |>
   mutate(
     fill_group = paste(Continent, mge_status, sep = "__")
   )
 
 fill_values <- c()
 
 for (ct in continent_levels) {
   fill_values[paste0(ct, "__MGE")] <- continent_palette[ct]
   fill_values[paste0(ct, "__No MGE")] <- alpha(continent_palette[ct], 0.35)
 }
 
 # Percentage labels -------------------------------------------------------
 
 label_dat <- plot_dat |>
   filter(mge_status == "MGE") |>
   mutate(
     pct = if_else(continent_total > 0, n / continent_total, 0),
     pct_label = percent(pct, accuracy = 0.1),
     y_pos = continent_total + max(continent_total) * 0.03
   )
 
 # Plot --------------------------------------------------------------------
 
 p <- ggplot(
   plot_dat,
   aes(x = Continent, y = n, fill = fill_group)
 ) +
   geom_col(width = 0.72, color = "white", linewidth = 0.3) +
   geom_text(
     data = label_dat,
     aes(x = Continent, y = y_pos, label = pct_label),
     inherit.aes = FALSE,
     vjust = 0,
     size = 4.2
   ) +
   scale_fill_manual(values = fill_values, guide = "none") +
   scale_x_discrete(labels = continent_labels) +
   scale_y_continuous(
     expand = expansion(mult = c(0, 0.16))
   ) +
   coord_cartesian(clip = "off") +
   labs(
     title = "Continent distribution with MGE status",
     x = NULL,
     y = "Count"
   ) +
   theme_bw(base_size = 13) +
   theme(
     panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     legend.position = "none",
     plot.title = element_text(hjust = 0.5, face = "bold"),
     axis.text.x = element_text(
       size = 12,
       face = "bold",
       colour = continent_palette[continent_levels]
     ),
     axis.text.y = element_text(size = 11),
     plot.margin = margin(20, 20, 20, 20)
   )
 
 # Preview -----------------------------------------------------------------
 
 print(p)
 
 # Save --------------------------------------------------------------------
 
 ggsave(
   filename = file.path(output_dir, "continent_mge_vertical_barplot.pdf"),
   plot = p,
   width = 8,
   height = 7,
   bg = "white"
 )
 

 

 
 
# Add Sankey diagram relationships between genes, strains, and continents. ---------------------
 
 
 # Load packages -----------------------------------------------------------
 
 library(readr)
 library(dplyr)
 library(stringr)
 library(ggplot2)
 library(ggalluvial)
 
 # Paths -------------------------------------------------------------------
 
 input_file <- "/Users/yqhhh/Desktop/AMR_with_MGE_3.31_plasmid_yesno.tsv"
 output_dir <- "/Users/yqhhh/Desktop"
 
 dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
 
 # Settings ----------------------------------------------------------------
 
 top_n_strain <- 10
 top_n_gene <- 10
 
 # Helper ------------------------------------------------------------------
 
 clean_text <- function(x) {
   x |>
     as.character() |>
     str_replace_all("’", "'") |>
     str_replace_all("“|”", "\"") |>
     str_trim() |>
     na_if("")
 }
 
 # Read data ---------------------------------------------------------------
 
 dat <- read_tsv(input_file, show_col_types = FALSE)
 
 required_cols <- c("Continent", "Strain", "GENE")
 missing_cols <- setdiff(required_cols, names(dat))
 
 if (length(missing_cols) > 0) {
   stop(paste0("Missing Required List: ", paste(missing_cols, collapse = ", ")))
 }
 
 dat <- dat |>
   mutate(
     Continent = clean_text(Continent),
     Strain = clean_text(Strain),
     GENE = clean_text(GENE)
   ) |>
   filter(
     !is.na(Continent),
     !is.na(Strain),
     !is.na(GENE)
   )
 
 # Keep top groups ---------------------------------------------------------
 
 top_strain <- dat |>
   count(Strain, sort = TRUE, name = "n") |>
   slice_head(n = top_n_strain) |>
   pull(Strain)
 
 top_gene <- dat |>
   count(GENE, sort = TRUE, name = "n") |>
   slice_head(n = top_n_gene) |>
   pull(GENE)
 
 dat <- dat |>
   mutate(
     strain_group = if_else(Strain %in% top_strain, Strain, "Other strains"),
     gene_group = if_else(GENE %in% top_gene, GENE, "Other genes")
   )
 
 # Order levels ------------------------------------------------------------
 
 continent_levels <- dat |>
   count(Continent, sort = TRUE, name = "n") |>
   pull(Continent)
 
 strain_levels <- dat |>
   count(strain_group, sort = TRUE, name = "n") |>
   mutate(is_other = strain_group == "Other strains") |>
   arrange(is_other, desc(n), strain_group) |>
   pull(strain_group)
 
 gene_levels <- dat |>
   count(gene_group, sort = TRUE, name = "n") |>
   mutate(is_other = gene_group == "Other genes") |>
   arrange(is_other, desc(n), gene_group) |>
   pull(gene_group)
 
 dat <- dat |>
   mutate(
     Continent = factor(Continent, levels = continent_levels),
     strain_group = factor(strain_group, levels = strain_levels),
     gene_group = factor(gene_group, levels = gene_levels)
   )
 
 # Build alluvial data -----------------------------------------------------
 
 plot_dat <- dat |>
   count(Continent, strain_group, gene_group, name = "Freq")
 
 # Colors ------------------------------------------------------------------
 
 # 
 continent_palette <- c(
   "Asia" = "#ffa88d",
   "Europe" = "#FFD84F",
   "North America" = "#36CFFD",
   "Africa" = "#5CBC71"
 )
 
 missing_continent <- setdiff(levels(dat$Continent), names(continent_palette))
 
 if (length(missing_continent) > 0) {
   extra_cols <- c("#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#7E6148")
   continent_palette[missing_continent] <- extra_cols[seq_along(missing_continent)]
 }
 
 #
 all_nodes <- c(
   levels(dat$Continent),
   levels(dat$strain_group),
   levels(dat$gene_group)
 )
 
 node_palette <- setNames(
   c(
     "#4DBBD5", "#E64B35", "#00A087", "#3C5488", "#F39B7F", "#8491B4",
     "#91D1C2", "#DC0000", "#7E6148", "#B09C85", "#F39B7F", "#D55E00",
     "#CC79A7", "#0072B2", "#E69F00", "#009E73", "#56B4E9", "#999999",
     "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
     "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928",
     "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
     "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
     "#BDBDBD", "#969696"
   )[seq_along(all_nodes)],
   all_nodes
 )
 
 #
 node_palette[names(continent_palette)] <- continent_palette[names(continent_palette)]
 
 # Plot --------------------------------------------------------------------
 
 p <- ggplot(
   plot_dat,
   aes(
     axis1 = Continent,
     axis2 = strain_group,
     axis3 = gene_group,
     y = Freq
   )
 ) +
   geom_alluvium(
     aes(fill = Continent),
     width = 0.18,
     alpha = 0.35,
     knot.pos = 0.4,
     color = "white",
     linewidth = 0.15
   ) +
   
   geom_stratum(
     aes(fill = after_stat(stratum)),
     width = 0.18,
     color = "grey35",
     linewidth = 0.3
   ) +
   scale_y_continuous(
     expand = expansion(mult = c(0.02, 0.05))
   ) +
   geom_text(
     stat = "stratum",
     aes(label = after_stat(stratum)),
     size = 3.2,
     color = "black"
   ) +
   scale_x_discrete(
     limits = c("Continent", "Strain", "GENE"),
     expand = c(0.05, 0.05)
   ) +
   scale_fill_manual(values = node_palette, guide = "none") +
   labs(
     title = "Relationships Among Continent, Strain, and Gene",
     x = NULL,
     y = "Count"
   ) +
   theme_minimal(base_size = 13) +
   theme(
     panel.grid = element_blank(),
     axis.text.x = element_text(face = "bold", size = 13),
     axis.text.y = element_blank(),
     axis.title.y = element_text(face = "bold"),
     plot.title = element_text(hjust = 0.5, face = "bold"),
     panel.background = element_rect(fill = "white", color = NA),
     plot.background = element_rect(fill = "white", color = NA)
   )
 
 # Preview -----------------------------------------------------------------
 
 print(p)
 
 # Save --------------------------------------------------------------------
 
 
 ggsave(
   filename = file.path(output_dir, "continent_strain_gene_colored_nodes_sankey.pdf"),
   plot = p,
   width = 14,
   height = 9,
   bg = "white"
 )
 

 
 
 

# Representative coexisting drug resistance annotations ----- Visualization 
 
 
 
 library(readr)
 library(dplyr)
 library(stringr)
 library(ggplot2)
 

 #=========================
 pattern_file <- "/Users/yqhhh/Desktop/05_top_patterns_reverse_merged.tsv"
 amr_file <- "/Users/yqhhh/Desktop/AMR_with_MGE_3.31_plasmid_yesno.tsv"
 
 output_pdf <- "/Users/yqhhh/Desktop/AMR_class_legend_macaron.pdf"
 
 #=========================
 #=========================
 pattern_df <- read_tsv(pattern_file, show_col_types = FALSE)
 amr_df <- read_tsv(amr_file, show_col_types = FALSE)
 
 names(pattern_df) <- trimws(names(pattern_df))
 names(amr_df) <- trimws(names(amr_df))
 
 #=========================
 # Only take the top 10 patterns
 #=========================
 top_patterns <- pattern_df %>%
   slice_head(n = 10)
 
 # Extract the contigs corresponding to these patterns
 top_contigs <- top_patterns %>%
   pull(contigs) %>%
   str_split(";") %>%
   unlist()
 
 #=========================
 # 4. Extract the corresponding AMR class
 #=========================
 legend_classes <- amr_df %>%
   filter(Contigs %in% top_contigs) %>%
   filter(!is.na(`AMR class`)) %>%
   distinct(`AMR class`) %>%
   arrange(`AMR class`) %>%
   pull(`AMR class`)
 

 #=========================
 macaron_palette <- c(
   "#F4A7B9",
   "#F6C1A4",
   "#F9E0AE",
   "#CDE7BE",
   "#A8D8EA",
   "#B8B5FF",
   "#D5AAFF",
   "#FFC8DD",
   "#BDE0FE",
   "#E2F0CB"
 )
 

 n_cls <- length(legend_classes)
 
 if (n_cls > length(macaron_palette)) {
   stop("There are more than 10 AMR classes, so we need to expand the color scheme.")
 }
 
 amr_class_cols <- macaron_palette[1:n_cls]
 names(amr_class_cols) <- legend_classes
 
 #=========================
 # 6.Building legend data
 #=========================
 legend_df <- data.frame(
   AMR_class = legend_classes
 ) %>%
   mutate(
     AMR_class = factor(AMR_class, levels = rev(AMR_class))
   )
 
 #=========================
 # 7. Drawing
 #=========================
 p <- ggplot(legend_df, aes(y = AMR_class, fill = AMR_class)) +
   
   geom_point(
     aes(x = 1),
     shape = 22,
     size = 8
   ) +
   
   geom_text(
     aes(x = 1.8, label = AMR_class),
     hjust = 0,
     size = 5
   ) +
   
   scale_fill_manual(values = amr_class_cols) +
   
   xlim(1, 3) +
   
   labs(title = "AMR Class") +
   
   theme_minimal(base_size = 14) +
   theme(
     legend.position = "none",
     axis.title = element_blank(),
     axis.text = element_blank(),
     axis.ticks = element_blank(),
     panel.grid = element_blank(),
     plot.title = element_text(face = "bold", hjust = 0)
   )
 #=========================
 # 8. Save
 #=========================
 ggsave(
   output_pdf,
   p,
   width = 6,
   height = 0.5 * n_cls + 2
 )
 
 print(p)
 
 
 
 
 

# The painting represents the coexistence of drug resistance and its corresponding pie chart 

 library(readr)
 library(dplyr)
 library(stringr)
 library(ggplot2)
 library(ggforce)
 library(scales)

 
 # =========================
 top_file <- "/Users/yqhhh/Desktop/05_top_patterns_reverse_merged.tsv"
 all_pattern_file <- "/Users/yqhhh/Desktop/
 amr_file <- "/Users/yqhhh/Desktop/AMR_with_MGE_plasmid_yesno.tsv"
 
 # Read data
 # =========================
 top_df <- read_tsv(top_file, show_col_types = FALSE)
 all_pattern_df <- read_tsv(all_pattern_file, show_col_types = FALSE)
 amr_df <- read_tsv(amr_file, show_col_types = FALSE)
 
 names(all_pattern_df)[1] <- "Contigs"
 names(all_pattern_df)[ncol(all_pattern_df)] <- "pattern_id"
 
 # =========================
 # Top10 pattern
 # =========================
 top10_ids <- top_df %>%
   slice_head(n = 10) %>%
   pull(pattern_id)
 
 pattern_levels <- rev(top10_ids)
 
 # =========================
 # contig -> pattern 
 # =========================
 pattern_map <- all_pattern_df %>%
   transmute(
     Contigs = as.character(Contigs),
     pattern_id = as.character(pattern_id)
   ) %>%
   filter(pattern_id %in% top10_ids) %>%
   distinct()
 
 # =========================
 # Merge AMR
 # =========================
 df_plot <- amr_df %>%
   mutate(
     Contigs = as.character(Contigs),
     Plasmid = as.character(Plasmid)
   ) %>%
   left_join(pattern_map, by = "Contigs") %>%
   filter(!is.na(pattern_id)) %>%
   mutate(
     pattern_id = factor(pattern_id, levels = pattern_levels)
   )
 
 # =========================
 # 6. Plasmid Classification
 # =========================
 df_plot <- df_plot %>%
   mutate(
     Plasmid_group = case_when(
       Plasmid %in% c("Yes","YES","yes") ~ "Plasmid",
       Plasmid %in% c("No","NO","no") ~ "Not plasmid",
       TRUE ~ NA_character_
     )
   )
 
 # =========================
 # 7. Building pie data
 # =========================
 make_pie_data <- function(df, col) {
   df %>%
     transmute(
       pattern_id,
       group = .data[[col]]
     ) %>%
     filter(!is.na(group)) %>%
     group_by(pattern_id, group) %>%
     summarise(n = n(), .groups = "drop") %>%
     group_by(pattern_id) %>%
     mutate(
       prop = n / sum(n),
       ymax = cumsum(prop),
       ymin = lag(ymax, default = 0)
     ) %>%
     ungroup() %>%
     mutate(
       y0 = as.numeric(pattern_id)
     )
 }
 
 pie_plasmid <- make_pie_data(df_plot, "Plasmid_group")
 pie_continent <- make_pie_data(df_plot, "Continent")
 pie_strain <- make_pie_data(df_plot, "Strain")
 
 # =========================
 # 8. Color
 # =========================
 continent_colors <- c(
   "Asia" = "#ffa88d",
   "Europe" = "#FFD84F"
 )
 
 plasmid_colors <- c(
   "Plasmid" = "#8B0000",
   "Not plasmid" = "#FF8C69"
 )
 
 strain_levels <- sort(unique(pie_strain$group))
 strain_colors <- setNames(
   colorRampPalette(c("#FFB6C1","#B5EAD7","#C7CEEA","#FFDAC1","#E2F0CB"))(length(strain_levels)),
   strain_levels
 )
 
 # =========================
 # 9. Drawing
 # =========================
 plot_pie <- function(df, colors, title) {
   
   df <- df %>%
     filter(group %in% names(colors)) %>%
     mutate(group = factor(group, levels = names(colors)))
   
   ggplot(df) +
     geom_arc_bar(
       aes(
         x0 = 1,
         y0 = y0,
         r0 = 0,
         r = 0.33,
         start = 2*pi*ymin,
         end = 2*pi*ymax,
         fill = group
       ),
       color = NA,
       linewidth = 0
     ) +
     scale_fill_manual(values = colors) +
     scale_y_continuous(
       breaks = seq_along(pattern_levels),
       labels = NULL   
     ) +
     coord_fixed(xlim = c(0.5, 2.2)) +
     labs(title = title, fill = NULL) +
     theme_void(base_size = 13) +
     theme(
       plot.title = element_text(face = "bold", hjust = 0),
       axis.text.y = element_blank(),
       axis.ticks.y = element_blank(),
       legend.position = "right",
       legend.text = element_text(size = 11),
       legend.key.width = unit(0.45, "cm"),
       legend.key.height = unit(0.45, "cm"),
       legend.spacing.y = unit(0.25, "cm"),
       legend.key = element_rect(fill = "white", colour = NA),
       legend.background = element_blank()
     )
 }
 # =========================
 # 10. Generated graph
 # =========================
 p_plasmid <- plot_pie(pie_plasmid, plasmid_colors, "Plasmid")
 p_continent <- plot_pie(pie_continent, continent_colors, "Continent")
 p_strain <- plot_pie(pie_strain, strain_colors, "Strain")
 
 print(p_plasmid)
 print(p_continent)
 print(p_strain)
 
 
 
 out_dir <- "/Users/yqhhh/Desktop/pie_outputs_final"
 dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
 
 # ======================
 # PDF
 # ======================
 ggsave(
   file.path(out_dir, "pie_plasmid.pdf"),
   p_plasmid,
   width = 5,
   height = 10
 )
 
 ggsave(
   file.path(out_dir, "pie_continent.pdf"),
   p_continent,
   width = 5,
   height = 10
 )
 
 ggsave(
   file.path(out_dir, "pie_strain.pdf"),
   p_strain,
   width = 6,
   height = 10
 )