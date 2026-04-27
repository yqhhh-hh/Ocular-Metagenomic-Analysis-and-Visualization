# Visualization of annotated virulence genes -------------------

 
 
 library(readr)
 library(dplyr)
 library(stringr)
 

 input_file <- "/Users/yqhhh/Desktop/all_vfdb_merged_with_vf_category.tsv"
 
 
 df <- read_tsv(input_file, show_col_types = FALSE)
 names(df) <- trimws(names(df))
 
 df <- df %>%
   mutate(
     VF_category = str_squish(VF_category),
     GENE = str_squish(GENE)
   ) %>%
   filter(!is.na(VF_category), VF_category != "")
 
 vf_count <- df %>%
   count(VF_category, sort = TRUE)
 
 vf_type_n <- n_distinct(df$GENE)
 vf_gene_n <- nrow(df)
 
 # Fixed color mapping
 vf_colors <- c(
   "Adherence" = "#F9B5A8",           
   "Regulation" = "#FAD4C2",          
   "Motility" = "#A8D8EA",            
   "Immune modulation" = "#A8D9CF",   
   "Invasion" = "#A6BDD9",            
   "Exotoxin" = "#C2D0E0",            
   "Effector delivery system" = "#C9E6DD", 
   "Exoenzyme" = "#D4C2AE",           
   "Biofilm" = "#E6CCD9",            
   "Antimicrobial activity/Competitive advantage" = "#C9E6DD",
   "Antimicrobial activity/\nCompetitive advantage" = "#C9E6DD",
   "Nutritional/Metabolic factor" = "#B8E0D0",   
   "Stress survival" = "#C5D5E8",                
   "Others" = "#F0F0F0"                         
 )
 

 vf_order <- vf_count$VF_category
 vf_count$VF_category <- factor(vf_count$VF_category, levels = vf_order)
 
 p1 <- ggplot(vf_count, aes(x = 2, y = n, fill = VF_category)) +
   geom_col(width = 1, color = "white", linewidth = 0.2) +
   coord_polar(theta = "y") +
   xlim(0.5, 2.5) +
   scale_fill_manual(
     values = vf_colors,
     breaks = vf_order,
     labels = paste0(vf_count$VF_category, " (", vf_count$n, ")")
   ) +
   theme_void() +
   theme(
     legend.position = "right",
     legend.title = element_blank(),
     legend.text = element_text(size = 12),
     plot.margin = margin(20, 10, 10, 10),
     plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
   ) +
   annotate(
     "text",
     x = 0.6,
     y = 0,
     label = paste0("VF gene type = ", vf_type_n, "\nTotal VF gene = ", vf_gene_n),
     size = 4
   ) +
   ggtitle("Overall VF")
 
 
 
 input_file2 <- "/Users/yqhhh/Desktop/co-vf-arg-30_meta.tsv"
 
 df2 <- read_tsv(input_file2, show_col_types = FALSE)
 names(df2) <- trimws(names(df2))
 
 df2 <- df2 %>%
   mutate(
     VF_category = str_squish(VF_category),
     GENE_vfdb = str_squish(GENE_vfdb)
   ) %>%
   filter(!is.na(VF_category), VF_category != "")
 
 vf_count2 <- df2 %>%
   count(VF_category, sort = TRUE)
 
 vf_type_n2 <- n_distinct(df2$GENE_vfdb)
 vf_gene_n2 <- nrow(df2)
 
 # Key point: The order should be consistent with overall.
 vf_count2$VF_category <- factor(vf_count2$VF_category, levels = vf_order)
 
 p2 <- ggplot(vf_count2, aes(x = 2, y = n, fill = VF_category)) +
   geom_col(width = 1, color = "white", linewidth = 0.2) +
   coord_polar(theta = "y") +
   xlim(0.5, 2.5) +
   scale_fill_manual(
     values = vf_colors,
     breaks = vf_count2$VF_category,   
     labels = paste0(vf_count2$VF_category, " (", vf_count2$n, ")")
   ) +
   theme_void() +
   theme(
     legend.position = "right",
     legend.title = element_blank(),
     legend.text = element_text(size = 12),
     plot.margin = margin(20, 10, 10, 10),
     plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
   ) +
   annotate(
     "text",
     x = 0.6,
     y = 0,
     label = paste0("VF gene type = ", vf_type_n2, "\nTotal VF gene = ", vf_gene_n2),
     size = 4
   ) +
   ggtitle("AMR-associated VF")
 
 print(p2)
 
 print(p1)
 
 
 final_plot <- p1 + p2 +
   plot_layout(ncol = 2)  
 
 
 print(final_plot)
 
 ggsave("VF_Horizontal light color.pdf", final_plot,
        width = 16, height = 8, bg = "white")
 
 
 
 

# Sankey diagram representation of coexisting drug-resistant virulence genes -------------
 
 
 
 library(readr)
 library(dplyr)
 library(stringr)
 library(ggplot2)
 library(ggalluvial)
 
 #========================
 # 1. Read data
 #========================
 input_file <- "/Users/yqhhh/Desktop/co-vf-arg-30_meta(end).tsv"
 
 df <- read_tsv(input_file, show_col_types = FALSE)
 names(df) <- trimws(names(df))
 
 #========================
 # 2. Clean up fields
 #========================
 df <- df %>%
   mutate(
     VF_category = str_squish(VF_category),
     GENE_vfdb   = str_squish(GENE_vfdb),
     Strain      = str_squish(Strain),
     Continent   = str_squish(Continent),
     MGE         = str_squish(as.character(MGE))
   ) %>%
   filter(
     !is.na(VF_category), VF_category != "",
     !is.na(GENE_vfdb),   GENE_vfdb != "",
     !is.na(Strain),      Strain != "",
     !is.na(Continent),   Continent != "",
     !is.na(MGE),         MGE != ""
   )
 
 #========================
 # 3. Select Top N
 #========================
 top_type <- df %>%
   count(VF_category, sort = TRUE) %>%
   slice_head(n = 8) %>%
   pull(VF_category)
 
 top_gene <- df %>%
   count(GENE_vfdb, sort = TRUE) %>%
   slice_head(n = 20) %>%
   pull(GENE_vfdb)
 
 top_strain <- df %>%
   count(Strain, sort = TRUE) %>%
   slice_head(n = 3) %>%
   pull(Strain)
 
 #========================
 # 4. The rest are classified as Others
 #========================
 df2 <- df %>%
   mutate(
     VF_category = if_else(VF_category %in% top_type, VF_category, "Others"),
     GENE_vfdb   = if_else(GENE_vfdb %in% top_gene, GENE_vfdb, "Others"),
     Strain      = if_else(Strain %in% top_strain, Strain, "Others")
   )
 
 #========================
 # 5. Long tag line break
 #========================

 
 #========================
 # 6. Summary
 #========================
 plot_df <- df2 %>%
   count(VF_category, GENE_vfdb, Strain, Continent, MGE, name = "Freq") %>%
   ungroup()
 
 #========================
 # 7. Order
 #========================
 vf_order <- c(top_type, "Others")
 gene_order <- c(top_gene, "Others")
 strain_order <- c(top_strain, "Others")
 
 continent_order <- plot_df %>%
   count(Continent, wt = Freq, sort = TRUE) %>%
   pull(Continent) %>%
   as.character()
 
 mge_order <- plot_df %>%
   count(MGE, wt = Freq, sort = TRUE) %>%
   pull(MGE) %>%
   as.character()
 
 plot_df <- plot_df %>%
   mutate(
     VF_category = factor(VF_category, levels = vf_order),
     GENE_vfdb   = factor(GENE_vfdb, levels = gene_order),
     Strain      = factor(Strain, levels = strain_order),
     Continent   = factor(Continent, levels = continent_order),
     MGE         = factor(MGE, levels = mge_order)
   )
 
 #========================
 # 8. Color scheme
 #========================
 nice_cols <- c(
   "#F9B5A8", "#FAD4C2", "#A8D8EA", "#A8D9CF",
   "#A6BDD9", "#C2D0E0", "#C9E6DD", "#D4C2AE"
 )
 
 vf_colors <- setNames(nice_cols[seq_along(top_type)], top_type)
 vf_colors <- c(vf_colors, "Others" = "#BDBDBD")
 vf_colors <- vf_colors[levels(plot_df$VF_category)]
 
 #========================
 # 9. Drawing
 #========================
 p <- ggplot(
   plot_df,
   aes(
     axis1 = VF_category,
     axis2 = GENE_vfdb,
     axis3 = Strain,
     axis4 = Continent,
     axis5 = MGE,
     y = Freq
   )
 ) +
   geom_alluvium(
     aes(fill = VF_category),
     width = 0.14,
     alpha = 0.85,
     knot.pos = 0.35
   ) +
   geom_stratum(
     width = 0.14,
     fill = "grey96",
     color = "grey40",
     linewidth = 0.35
   ) +
   geom_text(
     stat = "stratum",
     aes(label = after_stat(stratum)),
     size = 3,
     fontface = "bold",
     lineheight = 0.9
   ) +
   scale_x_discrete(
     limits = c("VF_category", "Gene", "Strain", "Continent", "MGE"),
     expand = c(0.12, 0.03)
   ) +
   scale_fill_manual(values = vf_colors, drop = FALSE) +
   coord_cartesian(clip = "off") +
   labs(
     title = "AMR-associated virulence genes",
     x = NULL,
     y = NULL,
     fill = "VF_category"
   ) +
   theme_minimal(base_size = 12) +
   theme(
     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
     axis.text.y = element_blank(),
     axis.text.x = element_text(face = "bold", size = 12),
     axis.title.y = element_text(face = "bold", size = 12),
     panel.grid = element_blank(),
     legend.position = "right",
     legend.title = element_text(face = "bold", size = 12),
     legend.text = element_text(size = 11),
     plot.margin = margin(t = 15, r = 25, b = 15, l = 95)
   )
 
 print(p)
 
 #========================
 # 10. Export
 #========================
 ggsave(
   filename = "/Users/yqhhh/Desktop/Sankey_VFcategory_Gene_Strain_Continent_MGE.pdf",
   plot = p,
   width = 22,
   height = 11,
   bg = "white"
 )