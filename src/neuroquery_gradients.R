# NEUROQUERY STUFF

schaef_ns_ <- read_delim("stuff_for_revisions/schaefer1000_neurosynth_terms.csv")

roi_labs <- read_delim("stuff_for_revisions/Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm.Centroid_RAS.csv") |> 
  select(1, 2) |> 
  rename(region = `ROI Name`,
         parcel = `ROI Label`)

r_filt <- 0.05
gradient_study <- "biofinder"


schaef_ns_g1 <- schaef_ns_ |> 
  mutate(feature = str_remove_all(feature, "LDA50_abstract_weight__"),
         feature = str_remove(feature, "^\\d+_")) |> 
  arrange(parcel, desc(r)) |> 
  #group_by(parcel) |> 
  #slice_head(n=5) |> 
  filter(r > r_filt, .by = parcel) |> 
  left_join(roi_labs) |> 
  left_join(grad_df |> filter(study== gradient_study) |> select(starts_with("grad"), region)) |> 
  ungroup() |> 
  group_by(feature) |> 
  summarise(mentions = n(),
            mean_g1 = mean(gradient1)) 

schaef_ns_g2 <- schaef_ns_ |> 
  mutate(feature = str_remove_all(feature, "LDA50_abstract_weight__"),
         feature = str_remove(feature, "^\\d+_")) |> 
  arrange(parcel, desc(r)) |> 
  #group_by(parcel) |> 
  #slice_head(n=5) |> 
  filter(r > r_filt, .by = parcel) |> 
  left_join(roi_labs) |> 
  left_join(grad_df |> filter(study== gradient_study) |> select(starts_with("grad"), region)) |> 
  ungroup() |> 
  group_by(feature) |> 
  summarise(mentions = n(),
            mean_g2 = mean(gradient2)) 


schaef_ns_g3 <- schaef_ns_ |> 
  mutate(feature = str_remove_all(feature, "LDA50_abstract_weight__"),
         feature = str_remove(feature, "^\\d+_")) |> 
  arrange(parcel, desc(r)) |> 
  #group_by(parcel) |> 
  #slice_head(n=5) |> 
  filter(r > r_filt, .by = parcel) |> 
  left_join(roi_labs) |> 
  left_join(grad_df |> filter(study== gradient_study) |> select(starts_with("grad"), region)) |> 
  ungroup() |> 
  group_by(feature) |> 
  summarise(mentions = n(),
            mean_g3 = mean(gradient3)) 


g3p <- ggplot(schaef_ns_g3, aes(x = reorder(feature, mean_g3), y = mean_g3, fill = mean_g3)) +
  #geom_point() +
  geom_col() +
  geom_text(aes(label = mentions), 
            hjust = -0.3,    
            size = 3) +
  theme_minimal() +
  scale_fill_gradient2(
    low = gradient_cols[1, 3],
    mid = "white",
    high = gradient_cols[2, 3]
  ) +
  coord_flip()+
  theme(
    legend.position = "top"
  )


g2p <- ggplot(schaef_ns_g2, aes(x = reorder(feature, mean_g2), y = mean_g2, fill = mean_g2)) +
  #geom_point() +
  geom_col() +
  geom_text(aes(label = mentions), 
            hjust = -0.3,    # position above the bar
            size = 3) +
  theme_minimal() +
  scale_fill_gradient2(
    low = gradient_cols[1, 2],
    mid = "white",
    high = gradient_cols[2, 2]
  ) +
  coord_flip() +
  theme(
    legend.position = "top"
  )

g1p <- ggplot(schaef_ns_g1, aes(x = reorder(feature, mean_g1), y = mean_g1, fill = mean_g1)) +
  #geom_point() +
  geom_col() +
  geom_text(aes(label = mentions), 
            hjust = -0.3,    # position above the bar
            size = 3) +
  theme_minimal() +
  scale_fill_gradient2(
    low = gradient_cols[1, 1],
    mid = "white",
    high = gradient_cols[2, 1]
  ) +
  coord_flip() +
  theme(
    legend.position = "top"
  )

g1p + g2p + g3p


schaef_ns_g1 |> 
  ggplot() +
  geom_text_wordcloud(aes(label=feature)) +
  theme_minimal()


nq_results <- read_rds("stuff_for_revisions/schaefer1000_NQ_results.rds")

# nq_results <- data.frame()
# 
# for (i in 1:1000) {
#   nq_i <- read_delim(paste0("stuff_for_revisions/schaefer1000_NQ/feat_parc_",i,".csv"), show_col_types = FALSE)
#   nq_results <- rbind(nq_results, nq_i)
# }

library(httr2)
library(jsonlite)

url <- "http://www.cognitiveatlas.org/api/search"
resp <- request(url) |> req_perform()
ca <- resp_body_json(resp)
ca <- tolower(sapply(ca, function(x) x$name))

nq_terms <- nq_results$feature |> unique() |> str_remove_all("neuroquery6308_combined_tfidf__")
length(nq_terms[nq_terms %in% ca])


idx <- stringdist::amatch(nq_terms, ca, method = "jw", maxDist = 0.1)
nq_match <- data.frame(
  nq_term = nq_terms,
  ca_match = ifelse(is.na(idx), NA, ca[idx])
) |> filter(!is.na(ca_match)) |> 
  mutate(dist = stringdist::stringdist(tolower(nq_term), ca_match, method = "jw"))

nq_filt <- nq_match$nq_term
nq_filt |> length()

gradient_study ="biofinder"




yeo_nq <- nq_results |> 
  arrange(parcel, desc(r)) |> 
  mutate(feature = str_remove_all(feature, "neuroquery6308_combined_tfidf__")) |> 
  filter(r > 0.05) |> 
  filter(feature %in% nq_filt) |> 
  left_join(roi_labs) |> 
  left_join(roi_data) |> 
  ungroup() |> 
  group_by(feature, name, col) |> 
  summarise(mentions = n(),
            mean_r = mean(r)) |> 
  filter(mentions > 10, mean_r>0.06
  )


nq_g1 <- nq_results |> 
  arrange(parcel, desc(r)) |> 
  mutate(feature = str_remove_all(feature, "neuroquery6308_combined_tfidf__")) |> 
  filter(r > 0.05) |> 
  filter(feature %in% nq_filt) |> 
  left_join(roi_labs) |> 
  left_join(grad_df |> filter(study== gradient_study) |> select(starts_with("grad"), region)) |> 
  ungroup() |> 
  group_by(feature) |> 
  summarise(mentions = n(),
            mean_g = mean(gradient1), 
            mean_r = mean(r)) |> 
  rowwise() |> 
  # mutate(
  #   closest_match = cog_terms[which.min(stringdist::stringdist(tolower(feature), cog_terms, method = "jw"))],
  #   match_dist = min(stringdist::stringdist(tolower(feature), cog_terms, method = "jw"))
  # ) |>
  ungroup() |> 
  #filter(match_dist < 0.05) |> 
  filter(mentions > 15)


nq_g3 <- nq_results |> 
  arrange(parcel, desc(r)) |> 
  mutate(feature = str_remove_all(feature, "neuroquery6308_combined_tfidf__")) |> 
  filter(r > 0.05) |> 
  filter(feature %in% nq_filt) |> 
  left_join(roi_labs) |> 
  left_join(grad_df |> filter(study== gradient_study) |> select(starts_with("grad"), region)) |> 
  ungroup() |> 
  group_by(feature) |> 
  summarise(mentions = n(),
            mean_g = mean(gradient3), 
            mean_r = mean(r)) |> 
  # rowwise() |> 
  # mutate(
  #   closest_match = cog_terms[which.min(stringdist::stringdist(tolower(feature), cog_terms, method = "jw"))],
  #   match_dist = min(stringdist::stringdist(tolower(feature), cog_terms, method = "jw"))
  # ) |> 
  ungroup() |> 
  # filter(match_dist < 0.05) |> 
  filter(mentions > 15)

nq_g2 <- nq_results |> 
  arrange(parcel, desc(r)) |> 
  mutate(feature = str_remove_all(feature, "neuroquery6308_combined_tfidf__")) |> 
  filter(r > 0.05) |> 
  filter(feature %in% nq_filt) |> 
  left_join(roi_labs) |> 
  left_join(grad_df |> filter(study== gradient_study) |> select(starts_with("grad"), region)) |> 
  ungroup() |> 
  group_by(feature) |> 
  summarise(mentions = n(),
            mean_g = mean(gradient2), 
            mean_r = mean(r)) |> 
  # rowwise() |> 
  # mutate(
  #   closest_match = cog_terms[which.min(stringdist::stringdist(tolower(feature), cog_terms, method = "jw"))],
  #   match_dist = min(stringdist::stringdist(tolower(feature), cog_terms, method = "jw"))
  # ) |> 
  ungroup() |> 
  # filter(match_dist < 0.05) |> 
  filter(mentions > 10)


library(ggwordcloud)

groups <- yeo_nq$name |> unique()
n_groups <- 7
yeo_nq$x <- unlist(lapply(seq_len(n_groups), function(i) runif(sum(yeo_nq$name == groups[i]), (i-1)/n_groups, i/n_groups)))
yeo_nq$y <- runif(nrow(yeo_nq))


yeo_wc <- yeo_nq |> 
  group_by(name) |> 
  arrange(mentions, .by_group = TRUE ) |> 
  slice_head(n = 12) |> 
  ggplot()+
  geom_text_wordcloud(aes(label = feature, #x = x, y = y,#size = abs(mean_g), 
                          color = col, #x = mean_g
                          angle_group = name,
                          size = 12
  ),
  #mask = png::readPNG("stuff_for_revisions/wc_mask_brain_trans.png"), 
  #grid_margin = 0
  ) +
  #facet_wrap(~name) +
  scale_color_identity() +
  theme_void(base_size = 24) +
  theme(strip.text = element_blank(),
        plot.margin = unit( c(0,0,0,0) , units = "lines" ))

ggsave("/media/strg/repos/presentations/halftime_seminar/images/heatmaps/yeo_wordcloud.png", yeo_wc, width=6.4, height=6.4, bg="transparent")



ggsave("/media/strg/repos/presentations/halftime_seminar/images/heatmaps/g3_wordcloud.png", g3_wc, width=5, height=5, bg="transparent")


ggsave("wordcloud.png", g1_wc, width=6, height=6, bg="transparent")
ggsave("/media/strg/repos/presentations/halftime_seminar/images/heatmaps/g1_wordcloud.png", g1_wc, width=5, height=5, bg="transparent")

g2_wc <- nq_g2 |> 
  mutate(grad_side = ifelse(mean_g>0, "Visual", "Motor")) |> 
  # group_by(grad_side) |> 
  # arrange(mentions, .by_group = TRUE ) |> 
  # slice_head(n = 30) |> 
  # group_by(grad_side) |> 
  # slice_head(n=10) |> 
  ggplot()+
  geom_text_wordcloud(aes(label = feature, size = abs(mean_g), 
                          color = mean_g, x = mean_g),
                      #mask = png::readPNG("stuff_for_revisions/wc_mask_brain_trans.png"), 
                      #grid_margin = 0
  ) +
  #facet_wrap(~grad_side) +
  scale_color_gradient2(
    low = gradient_cols[1, 2],
    mid = "white",
    high = gradient_cols[2, 2]
  ) +
  theme_void(base_size = 24) 






g1_wc <- nq_g1 |> 
  mutate(grad_side = ifelse(mean_g>0, "Association", "Sensory"),
         abs_g  = abs(mean_g)) |> 
  mutate(
    mean_g = tanh(mean_g/20),
    #mean_g = sign(mean_g) * log1p(abs(mean_g))
  ) |> 
  group_by(grad_side) |> 
  arrange(mentions, .by_group = TRUE ) |>
  slice_head(n = 25) |>
  ggplot()+
  geom_text_wordcloud(aes(label = feature, size = abs(mean_g), 
                          color = mean_g, x = mean_g
                          ),# nudge_x = -1
                      #mask = png::readPNG("stuff_for_revisions/wc_mask_brain_trans.png"), 
                      , grid_margin = 15, seed = 6
  ) +
  #facet_wrap(~grad_side) +
  scale_color_gradient2(
    low = gradient_cols[1, 1],
    mid = "white",
    high = gradient_cols[2, 1]
  ) +
  xlim(-1, 1.8) +
  scale_size_area(max_size = 6) +
  theme_void(base_size = 24) 

#ggsave("/media/strg/repos/presentations/ddls2025/images/heatmaps/g1_wordcloud.png", g1_wc, width=8, height=4, bg="transparent")


g3_wc <- nq_g3 |> 
  mutate(grad_side = ifelse(mean_g>0, "Execution", "Representation")) |> 
  mutate(
    mean_g = tanh(mean_g/20),
    #mean_g = sign(mean_g) * log1p(abs(mean_g))
  ) |> 
  group_by(grad_side) |> 
  arrange(mentions, .by_group = TRUE ) |>
  slice_head(n = 25) |>
  ggplot()+
  geom_text_wordcloud(aes(label = feature, size = abs(mean_g), color = mean_g, x = mean_g),
                      , grid_margin = 15, seed = 6
                      ) +
  #facet_wrap(~grad_side) +
  scale_color_gradient2(
    low = gradient_cols[1, 3],
    mid = "white",
    high = gradient_cols[2, 3]
  ) +
  scale_size_area(max_size = 6)+
  xlim(-1, 1.8) +
  theme_void(base_size = 24) 

ggsave("/media/strg/repos/presentations/ddls2025/images/heatmaps/g3_wordcloud.png", g3_wc, width=8, height=4, bg="transparent")

library(magick)
img <- image_read("/media/strg/repos/presentations/ddls2025/images/heatmaps/g1_wordcloud.png")
trimmed <- image_trim(img)
image_write(trimmed, "/media/strg/repos/presentations/ddls2025/images/heatmaps/g1_wordcloud.png")

img <- image_read("/media/strg/repos/presentations/ddls2025/images/heatmaps/g3_wordcloud.png")
trimmed <- image_trim(img)
image_write(trimmed, "/media/strg/repos/presentations/ddls2025/images/heatmaps/g3_wordcloud.png")


grad_plots <- make_gradient_plots(grad_df |> filter(study=="biofinder"), add_shade = TRUE, shade_alpha = 0.01, shade_size = 0.05)


# x_range_df <- tibble(
#   gradient1 = seq(min(grad_df$gradient1, na.rm = TRUE), max(grad_df$gradient1, na.rm = TRUE), length.out = 100),
#   gradient2 = seq(min(grad_df$gradient2, na.rm = TRUE), max(grad_df$gradient2, na.rm = TRUE), length.out = 100),
#   gradient3 = seq(min(grad_df$gradient3, na.rm = TRUE), max(grad_df$gradient3, na.rm = TRUE), length.out = 100)
# ) |> pivot_longer(starts_with("grad"), names_to = "gradient", values_to = "value")

net_sep <- grad_df |> filter(study=="biofinder") |> 
  pivot_longer(starts_with("grad"), names_to = "gradient", values_to = "value") |>
  filter(gradient != "gradient2") |> 
  mutate(gradient = case_when(
    gradient == "gradient1" ~ "Sensory-Association",
    gradient == "gradient3" ~ "Representational-Executive"
  )) |> 
  ggplot(aes(x = value, y = name, fill = name)) +
  geom_boxplot() +
  geom_violin(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~gradient, scales = "free_x") +
  scale_fill_manual(values = net_names |> select(name, col) |> deframe()) +
  labs(y = "", x = "Placement of Yeo Networks on Gradients", fill = "") +
  guides(fill = guide_legend(nrow = 1)) +
  theme(legend.position = "top",
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        legend.text = element_text(hjust = 0.75), 
        legend.margin = margin(t = 20))

g1 <- grad_plots[[1]] + ggtitle("Sensory-Association Axis") + theme(plot.title = element_text(size = 22))
g3 <- grad_plots[[3]] + ggtitle("Representational-Executive Axis") + theme(plot.title = element_text(size = 22))
gplots <- (g1 + g3)
wcs <- (g1_wc +  g3_wc) 


test <- ggdraw() +
  draw_plot(g1, width = 0.5, height = 0.38, x = 0, y = 0.62) +
  draw_plot(g3, width = 0.5, height = 0.38, x = 0.5, y = 0.62) +
  draw_plot(g1_wc, width = 0.5, height = 0.6, x = 0.05, y = 0.17) +
  draw_plot(g3_wc, width = 0.5, height = 0.6, x = 0.55, y = 0.19) +
  draw_plot(net_sep, width = 1, height = 0.35, x = 0, y = 0)

test <- test + theme(plot.background = element_rect(colour = "black", linewidth = 1))

ggsave("/home/jorittmo/repositories/fc_changes_follow_gradients/paper/revision_figures/gradient_coverage.png", test, width=14, height=12, bg="white")


ggplot() +
  geom_sf(data = schaef1k$atlas, fill="grey80") +
  theme_void() +
  annotation_raster(png::readPNG("wordcloud.png"),
                    xmin = -0.5695375 , xmax = 2.0711176, ymin = -1.3635968, ymax = 0.3865248) 

grad1 <- schaef1k$atlas |> 
  ggplot() +
  geom_sf(aes(geometry = geometry), linewidth = 2, fill = "white", color = "white") + 
  theme_void() +
  theme(
    plot.margin = margin(0, 3, 0, 3, "cm")
  )

mask_file <- "stuff_for_revisions/wc_mask_brain.png"
ggsave(mask_file, grad1, bg = "black")

mask_img <- image_read(mask_file) |>
  image_transparent("white") |>
  image_convert(colorspace = "rgb", type = "truecoloralpha") |> 
  image_fx(expression = "r", channel = "red") |>    # force RGB
  image_fx(expression = "g", channel = "green") |> 
  image_fx(expression = "b", channel = "blue")

mask_file_trans <- "stuff_for_revisions/wc_mask_brain_trans.png"
image_write(mask_img, mask_file_trans)

x <- png::readPNG("stuff_for_revisions/wc_mask_brain_trans.png")

library(freesurferformats)
library(tidyverse)

dk_l <- read.fs.annot("/usr/local/freesurfer/subjects/fsaverage/label/lh.aparc.a2009s.annot")
sch_1000_l <- read.fs.annot("~/repositories/ggseg2/data/schaefer2018/fsaverage/label/lh.Schaefer2018_1000Parcels_7Networks_order.annot")


dk_r <- read.fs.annot("/usr/local/freesurfer/subjects/fsaverage/label/rh.aparc.a2009s.annot")
sch_1000_r <- read.fs.annot("~/repositories/ggseg2/data/schaefer2018/fsaverage/label/rh.Schaefer2018_1000Parcels_7Networks_order.annot")


dk_sch <- tibble(dk = c(dk_l$label_names, dk_r$label_names) , 
                 sch_1000 = c(sch_1000_l$label_names, sch_1000_r$label_names),
                 hemi = rep(c("l", "r"), each = length(sch_1000_l$label_names) ))

atlas_comp <- dk_sch |> group_by(dk, sch_1000, hemi) |> 
  summarise(n_vert = n()) |> 
  ungroup() |> 
  group_by(sch_1000, hemi) |> 
  slice_max(n_vert, n = 1, with_ties = FALSE) |> 
  filter(dk != "")

atlas_comp |> 
  right_join(schaef1k$atlas, join_by(sch_1000 == region)) |> 
  ggplot() +
  geom_sf(aes(geometry = geometry, fill = dk)) +
  theme_void()

atlas_comp <- atlas_comp |> 
  rename(region = sch_1000) |> 
  left_join(grad_df |> filter(study=="biofinder") |> select(region, starts_with("grad"))) 

atlas_comp_g1 <- atlas_comp |> 
  group_by(dk) |> 
  summarise(mean_g = mean(gradient1))

atlas_comp_g3 <- atlas_comp |> 
  group_by(dk) |> 
  summarise(mean_g = mean(gradient3))


comp_wc_g1 <- atlas_comp_g1 |> 
  mutate(grad_side = ifelse(mean_g>0, "Associative", "Sensory")) |> 
  ggplot()+
  geom_text_wordcloud(aes(label = dk, size = abs(mean_g), color = mean_g, x = mean_g)) +
  #facet_wrap(~grad_side) +
  scale_color_gradient2(
    low = gradient_cols[1, 1],
    mid = "white",
    high = gradient_cols[2, 1]
  ) +
  theme_void()


comp_wc_g3 <- atlas_comp_g3 |> 
  mutate(grad_side = ifelse(mean_g>0, "Exec", "Non-exec")) |> 
  ggplot()+
  geom_text_wordcloud(aes(label = dk, size = abs(mean_g), color = mean_g, x = mean_g))+
  #facet_wrap(~grad_side) +
  scale_color_gradient2(
    low = gradient_cols[1, 3],
    mid = "white",
    high = gradient_cols[2, 3]
  ) +
  theme_void()

comp_wc_g1 + comp_wc_g3

g1_p <- ggplot(nq_g1, aes(x = reorder(feature, mean_g), y = mean_g, fill = mean_g)) +
  #geom_point() +
  geom_col() +
  geom_text(aes(label = mentions), 
            hjust = -0.3,    # position above the bar
            size = 3) +
  theme_minimal() +
  scale_fill_gradient2(
    low = gradient_cols[1, 1],
    mid = "white",
    high = gradient_cols[2, 1]
  ) +
  coord_flip() +
  theme(
    legend.position = "top"
  )

g3_p <- ggplot(nq_g3, aes(x = reorder(feature, mean_g), y = mean_g, fill = mean_g)) +
  #geom_point() +
  geom_col() +
  geom_text(aes(label = mentions), 
            hjust = -0.3,    # position above the bar
            size = 3) +
  theme_minimal() +
  scale_fill_gradient2(
    low = gradient_cols[1, 3],
    mid = "white",
    high = gradient_cols[2, 3]
  ) +
  coord_flip() +
  theme(
    legend.position = "top"
  )

g1_p + g3_p
