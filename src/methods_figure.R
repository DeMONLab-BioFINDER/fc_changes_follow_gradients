library(proxy)
library(RColorBrewer)
library(kableExtra)
library(patchwork)

zero_out_mat <- function(mat, thresh) {
  library(matrixStats)
  thresholds <- colQuantiles(abs(mat), probs = thresh)
  mat[abs(mat) <= thresholds[col(mat)]] <- 0
  return(mat)
}


figure_path <- "paper/figures/conceptual_plot"
pastel_colors <- colorRampPalette(brewer.pal(9, "Pastel1"))(1000)

average_connectome <- readRDS("data/atlas_data/average_connectome_normalyoung.rds")
average_connectome[average_connectome<0] <- 0

subject_connectome <- readRDS(list.files("data/bf_src_data/connectomes", full.names = TRUE)[70])
subject_connectome[subject_connectome<0] <- 0


png(file.path(figure_path, "subject_connectome.png"), width = 800, height = 800)
par(bg=NA, mar = c(0, 0, 0, 0))
heatmap(subject_connectome, Rowv = NA, Colv = "Rowv",
        symm = TRUE, ColSideColors = pastel_colors,
        #hclustfun = hclust(, method = "complete"),
        RowSideColors = pastel_colors, col = hcl.colors(50),
        labRow = NA,
        labCol = NA,
        margins = c(0, 0)
)
dev.off()

zeroed_mat <- zero_out_mat(subject_connectome, 0.85)
affinity <- proxy::simil(zeroed_mat, method = "cosine", by_rows = FALSE) %>% proxy::as.matrix()
diag(affinity) <- 1

png(file.path(figure_path, "subject_affinity.png"), width = 800, height = 800)
par(bg=NA, mar = c(0, 0, 0, 0))
heatmap(affinity, Rowv = NA, Colv = "Rowv",
        symm = TRUE, ColSideColors = pastel_colors,
        #hclustfun = hclust(, method = "complete"),
        RowSideColors = pastel_colors, col = hcl.colors(100),
        labRow = NA,
        labCol = NA,
        margins = c(0, 0)
)
dev.off()

png(file.path(figure_path, "subject_col_ave.png"), width = 800, height = 800)
par(bg=NA, mar = c(0, 0, 0, 0))
heatmap(base::as.matrix(rbind(colMeans(affinity),colMeans(affinity))), Rowv = NA, Colv = NA,
        symm = FALSE, ColSideColors = pastel_colors,
        #hclustfun = hclust(, method = "complete"),
        #RowSideColors = pastel_colors,
        col = hcl.colors(100),
        labRow = NA,
        labCol = NA,
        margins = c(0, 0)
)
dev.off()
p_name <- "subject_col_ave.png"
img <- image_read(file.path(figure_path, p_name))
trimmed_img <- image_trim(img)
cropped_img <- image_crop(trimmed_img, geometry = "100%x10%")
image_write(cropped_img, path = file.path(figure_path, p_name))


png(file.path(figure_path, "average_connectome.png"), width = 800, height = 800)
par(bg=NA, mar = c(0, 0, 0, 0))
heatmap(average_connectome, Rowv = NA, Colv = "Rowv",
        symm = TRUE, ColSideColors = pastel_colors,
        #hclustfun = hclust(, method = "complete"),
        RowSideColors = pastel_colors, col = hcl.colors(100),
        labRow = NA,
        labCol = NA,
        margins = c(0, 0)
)
dev.off()

g1 <- grad_df %>% filter(study=="biofinder") %>% pull(gradient1)
g2 <- grad_df %>% filter(study=="biofinder") %>% pull(gradient2)
g3 <- grad_df %>% filter(study=="biofinder") %>% pull(gradient3)

g_list <- list(g1 = g1, g2 = g2, g3 = g3)

for (g in names(g_list)) {
  p_name <- paste0(g, ".png")
  png(file.path(figure_path, p_name), width = 800, height = 800)
  par(bg=NA, mar = c(0, 0, 0, 0))
  heatmap(base::as.matrix(rbind(g_list[[g]], g_list[[g]])), Rowv = NA, Colv = NA,
          symm = FALSE, ColSideColors = pastel_colors,
          #hclustfun = hclust(, method = "complete"),
          #RowSideColors = pastel_colors,
          col = hcl.colors(100),
          labRow = NA,
          labCol = NA,
          margins = c(0, 0)
  )
  dev.off()
  img <- image_read(file.path(figure_path, p_name))
  trimmed_img <- image_trim(img)
  cropped_img <- image_crop(trimmed_img, geometry = "100%x10%")
  image_write(cropped_img, path = file.path(figure_path, p_name))
}

atlas_geometry = readRDS("data/atlas_data/Schaefer2018_1000Parcels_geometry.rds")

net_names <- data.frame(name = c('Vis', 'SomMot', 'DorsAttn','SalVentAttn','Limbic', 'Cont', 'Default'),
                        col = c("#781286", "#4682B4", "#00760E", "#C43AFA", "#c7cc7a", "#E69422", "#CD3E4E"), #"#DCF8A4"
                        label = c(1:7))

pastel_brain <- atlas_geometry %>% 
  inner_join(data.frame(region = rois, p_cols = pastel_colors)) %>% 
  ggplot() +
  geom_sf(aes(fill = p_cols,
              geometry = geometry), linewidth= 0.2,
          show.legend = FALSE)+
  theme_void() +
  scale_fill_identity()

p_name <- "pastel_brain.png"
ggsave(file.path(figure_path, p_name), pastel_brain, width = 180, height = 180, units = "mm", dpi = 600, device = "png")


fits <- nodal_regression_fits(biofinder_df %>% 
                                filter(fmri_bl) %>%
                                mutate(motion = rsqa__MeanFD) %>% 
                                mutate(diagnosis = ifelse(diagnosis %in% c("Normal", "SCD"), "CN/SCD", diagnosis)) %>% 
                                mutate(diagnosis=factor(diagnosis, levels = c("CN/SCD", "MCI", "AD"))), inter_sub_sim_matrix = fc_measures_bf$affinity, 
                              roi_names = rois,
                              dep_var = "Nodal_Affinity",
                              model_formula = formula("Nodal_Affinity ~ age + diagnosis + sex + motion"))

# You have to printscreen this or use webshot
tidy(summary(fits[[1000]])) %>% 
  select(Term = term, Estimate = estimate, "t-value" = statistic) %>% 
  slice(2:4) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  mutate(Term = str_replace_all(Term, "diagnosis", "CN vs. ")) %>% 
  kbl(caption = "<span style='font-size: 21px;'>Regression results Parcel<sub>1000</sub></span>") %>%
  kable_classic(full_width = F, html_font = "Cambria", font_size = 21) -> tbl



bf_dx <-  plot_gradient_relationships(biofinder_df %>% 
                                        filter(fmri_bl) %>%
                                        mutate(motion = rsqa__MeanFD
                                        ) %>% 
                                        mutate(diagnosis = ifelse(diagnosis %in% c("Normal", "SCD"), "CN/SCD", diagnosis)) %>% 
                                        mutate(diagnosis=factor(diagnosis, levels = c("CN/SCD", "MCI", "AD"))), 
                                      gradient_data = grad_df %>% filter(study=="biofinder"), 
                                      gradients = c(1, 2, 3),
                                      vect = TRUE,
                                      r2_size = rel(4.5),
                                      gradient_colors = gradient_cols,
                                      list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                      mod_formula = formula(paste0(" ~ age + diagnosis + sex + rsqa__MeanFD")),
                                      covariates = c("sex", "rsqa__MeanFD"),
                                      filter_criteria = quo(),
                                      show_networks = FALSE,
                                      tag_prefix = "",
                                      tag_sep = "",
                                      layout_construction = "horizontal",
                                      include_gradient_plots = TRUE,
                                      right_term_side = FALSE,
                                      plt_title = "",
                                      cache_runs = FALSE)

p_dx <- bf_dx$plot &
  theme(text = element_text(size = 15),
        plot.tag = element_blank())

p_dx[[5]] <- p_dx[[5]] + ggtitle("MCI vs. CN")
p_dx[[6]] <- p_dx[[6]] + ggtitle("AD vs. CN")
#p_dx[[7]] <- p_dx[[7]] + ggtitle("AD")

img_width <- 180/25.4/2

p_name <- "figure_dx.png"
ggsave(file.path(figure_path, p_name), p_dx, width = img_width*3, height = img_width*0.9*3, units = "in", dpi = 600, device = "png")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, "33%x33%")
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)


#### legend stuff

get_net_legend <- function(){
  scale_factor <- 5
  x <- grad_df %>% filter(study == "biofinder") %>% 
    ggplot(aes(gradient1, gradient3, color = name)) +
    geom_point(alpha = 0.5) +
    labs(color = "Yeo Network") +
    guides(color = guide_legend(
      nrow = 1, 
      override.aes = list(size = 1.5*scale_factor)
    )) +
    scale_color_manual(values = net_names %>% select(name, col) %>% deframe) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(0.75, "cm"),
      legend.spacing.x = unit(1.5, "cm"),
      legend.direction = "horizontal",
      legend.title.position = "top",
      legend.text.position = "right",
      # legend.title = element_text(size = 6*scale_factor),
      # legend.text = element_text(size = 5*scale_factor, hjust = 0.75),
      legend.background = element_blank()
    )
  
  leg <- ggpubr::get_legend(x)
  leg <- ggpubr::as_ggplot(leg)
  leg
}
x <- get_net_legend()

p_name <- "net_legend.png"
img_width <- 45
ggsave(file.path(figure_path, p_name), x, width = img_width*scale_factor+75, height = img_width*0.1*scale_factor, units = "mm", dpi = 300, device = "png")
