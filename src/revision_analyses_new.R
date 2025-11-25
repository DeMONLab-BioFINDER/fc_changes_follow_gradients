###########################
# Broken down by pathology
###########################


figure_path <- "paper/revision_figures"
img_width <-  180 / 25.4
scaling_factor <-  2
magick_geom_scaling <- paste0(100/scaling_factor, "%x", 100/scaling_factor, "%")

# tau_ab_full <- plot_gradient_relationships(biofinder_df_all %>% filter(fmri_bl, !is.na(pathology_ad)),
#                                            gradient_data = grad_df |> filter(study== "biofinder"),
#                                            gradients = c(1, 3),
#                                            gradient_colors = gradient_cols,
#                                            list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
#                                            #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
#                                            add_shade = TRUE,
#                                            shade_alpha = 0.0075,
#                                            shade_size = 0.001,
#                                            empty_row_height = -0.15,
#                                            base_size_ = 16,
#                                            vect = TRUE,
#                                            mod_formula = formula(paste0("~ age + tau_pathology + abnorm_ab + sex + rsqa__MeanFD")),
#                                            covariates = c("sex", "rsqa__MeanFD"),
#                                            r2_size = rel(6),
#                                            filter_criteria = quo(),
#                                            show_networks = FALSE,
#                                            tag_prefix = "",
#                                            tag_sep = "",
#                                            layout_construction = "horizontal",
#                                            include_gradient_plots = TRUE,
#                                            right_term_side = FALSE,
#                                            plt_title = "Full sample",
#                                            plt_subtitle = TRUE,
#                                            cache_runs = FALSE)

#########################

tau_ab_cu <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, !is.na(pathology_ad), diagnosis %in% c("Normal", "SCD")
                                                                 ) |> 
                                           mutate(diag = diagnosis %in% c("Normal", "SCD")), 
                                           gradient_data = grad_df |> filter(study== "biofinder"), 
                                           gradients = c(1, 3),
                                           gradient_colors = gradient_cols,
                                           list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                           #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
                                           add_shade = TRUE, 
                                           shade_alpha = 0.0075,
                                           shade_size = 0.001,
                                           empty_row_height = -0.15,
                                           base_size_ = 16,
                                           vect = TRUE,
                                           mod_formula = formula(paste0("~ age + tau_pathology + abnorm_ab + sex + rsqa__MeanFD")),
                                           covariates = c("sex", "rsqa__MeanFD"),
                                           plt_title = "Cognitively Unimpaired",
                                           plt_subtitle = TRUE,
                                         brain_names = c("Age", "Tau Pathology", "Αβ+"),
                                         rectangle = TRUE)

tau_ab_cu_adj <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, !is.na(pathology_ad)) |> 
                                            mutate(has_diagnosis = diagnosis %in% c("MCI", "AD")), 
                                          gradient_data = grad_df |> filter(study== "biofinder"), 
                                          gradients = c(1, 3),
                                          gradient_colors = gradient_cols,
                                          list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                          #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
                                          add_shade = TRUE, 
                                          shade_alpha = 0.0075,
                                          shade_size = 0.001,
                                          empty_row_height = -0.15,
                                          base_size_ = 16,
                                          vect = TRUE,
                                          mod_formula = formula(paste0("~ age + tau_pathology + abnorm_ab + has_diagnosis + sex + rsqa__MeanFD")),
                                          covariates = c("sex", "rsqa__MeanFD"),
                                          plt_title = "Full sample",
                                          plt_subtitle = TRUE,
                                          brain_names = c("Age", "Tau Pathology", "Αβ+", "MCI/AD"),
                                          rectangle = TRUE)


# a <- tau_ab_full$plot + plot_annotation(theme = theme(plot.background = element_rect(color = "black", fill = NA, size = 1)))
# 
# b <- tau_ab_cu$plot + plot_annotation(theme = theme(plot.background = element_rect(color = "black", fill = NA, size = 1)))
# 
# ab_plot <- ggdraw() +
#   draw_plot(a, x = 0, y = 0, width = 0.495, height = 1) +
#   draw_plot(b, x = 0.505, y = 0, width = 0.495, height = 1)

p_name <- "tau_and_ab.png"
ggsave(file.path(figure_path, p_name), tau_ab_cu$plot,
       width = img_width*scaling_factor, height = img_width*0.7*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)

###########################
# Nodal Tau
###########################

tau <- readRDS("/media/strg/repos/fmri_project1/data/tau_pet/tau_pet_extracted/tau_pet_schaefer_1000_r.rds")

nodal_tau <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, !is.na(pathology_ad)#, diagnosis %in% c("Normal", "SCD")
                                                                 ), 
                                         t_mat = tau,
                                         gradient_data = grad_df |> filter(study== "biofinder"), 
                                         gradients = c(1, 3),
                                         gradient_colors = gradient_cols,
                                         list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                         #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
                                         add_shade = TRUE, 
                                         shade_alpha = 0.0075,
                                         shade_size = 0.001,
                                         empty_row_height = -0.15,
                                         base_size_ = 18,
                                         vect = TRUE,
                                         mod_formula = formula(paste0("~ age + sex + rsqa__MeanFD")),
                                         covariates = c("sex", "rsqa__MeanFD"),
                                         plt_title = "Parcelwise tau on parcelwise FCS regression",
                                         plt_subtitle = TRUE,
                                         rectangle = TRUE)

# nodal_tau_cont <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, !is.na(pathology_ad)) |> mutate(tau_glob = (cho12+cho34+cho56)/3),
#                                               t_mat = tau,
#                                               gradient_data = grad_df |> filter(study== "biofinder"),
#                                               gradients = c(1, 3),
#                                               gradient_colors = gradient_cols,
#                                               list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
#                                               #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
#                                               add_shade = TRUE,
#                                               shade_alpha = 0.0075,
#                                               shade_size = 0.001,
#                                               empty_row_height = -0.15,
#                                               base_size_ = 16,
#                                               vect = TRUE,
#                                               mod_formula = formula(paste0("~ age + tau_glob + sex + rsqa__MeanFD")),
#                                               covariates = c("sex", "rsqa__MeanFD"),
#                                               r2_size = rel(6),
#                                               filter_criteria = quo(),
#                                               show_networks = FALSE,
#                                               tag_prefix = "",
#                                               tag_sep = "",
#                                               layout_construction = "horizontal",
#                                               include_gradient_plots = TRUE,
#                                               right_term_side = FALSE,
#                                               plt_title = "Parcelwise tau controlled for total tau",
#                                               plt_subtitle = TRUE,
#                                               cache_runs = FALSE)



#a <-  + plot_annotation(theme = theme(plot.background = element_rect(color = "black", fill = NA, size = 1)))

#b <- nodal_tau_cont$plot + plot_annotation(theme = theme(plot.background = element_rect(color = "black", fill = NA, size = 1)))

# nodal_tau_plot <- ggdraw() +
#   draw_plot(a, x = 0, y = 0, width = 0.4, height = 1) +
#   draw_plot(b, x = 0.41, y = 0, width = 0.59, height = 1)
scaling_factor <- 2
p_name <- "nodal_tau.png"
ggsave(file.path(figure_path, p_name), nodal_tau$plot,
       width = img_width*scaling_factor, height = img_width*0.9*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)



##################################
# APOE
#################################


#img_width <-  9#180 / 25.4
scaling_factor <-  3
magick_geom_scaling <- paste0(100/scaling_factor, "%x", 100/scaling_factor, "%")

apoe <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl#, abnorm_ab == 0, pathology_ad<0.15
                                                            ) |> 
                                      mutate(has_diag = diagnosis %in% c("MCI", "AD")), 
                                    gradient_data = grad_df |> filter(study== "biofinder"), 
                                    gradients = c(1, 3),
                                    gradient_colors = gradient_cols,
                                    list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                    #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
                                    add_shade = TRUE, 
                                    shade_alpha = 0.0075,
                                    shade_size = 0.001,
                                    empty_row_height = -0.15,
                                    base_size_ = 16,
                                    vect = TRUE,
                                    mod_formula = formula(paste0("~ age + pathology_ad + apoe4 + has_diag + sex + rsqa__MeanFD")),
                                    covariates = c("sex", "rsqa__MeanFD"),
                                    plt_title = "Full sample",
                                    #brain_names = c("Age", "Pathology AD", "ε4 carrier"),
                                    plt_subtitle = TRUE,
                                    rectangle = TRUE)


apoe_abneg <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, 
                                                                  diagnosis %in% c("Normal", "SCD"), 
                                                                 # abnorm_ab==0
                                                                  ), 
                                    gradient_data = grad_df |> filter(study== "biofinder"), 
                                    gradients = c(1, 3),
                                    gradient_colors = gradient_cols,
                                    list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                    #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
                                    add_shade = TRUE, 
                                    shade_alpha = 0.0075,
                                    shade_size = 0.001,
                                    empty_row_height = -0.15,
                                    base_size_ = 16,
                                    vect = TRUE,
                                    mod_formula = formula(paste0("~ age + tau_pathology + abnorm_ab + apoe4 + sex + rsqa__MeanFD")),
                                    covariates = c("sex", "rsqa__MeanFD"),
                                    plt_title = "Cognitively Unimpaired (Aβ-)",
                                    brain_names = c("Age", "Tau Pathology", "ε4 carrier"),
                                    plt_subtitle = TRUE,
                                    rectangle = TRUE)


apoe_int <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, 
                                                                #diagnosis %in% c("Normal", "SCD"), 
                                                                #abnorm_ab==0
                                                                ) |> 
                                          mutate(has_diag = factor(diagnosis %in% c("MCI", "AD"))), 
                                    gradient_data = grad_df |> filter(study== "biofinder"), 
                                    gradients = c(1, 3),
                                    gradient_colors = gradient_cols,
                                    list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                    #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
                                    add_shade = TRUE, 
                                    shade_alpha = 0.0075,
                                    shade_size = 0.001,
                                    empty_row_height = -0.15,
                                    base_size_ = 16,
                                    vect = TRUE,
                                    mod_formula = formula(paste0("~ age + tau_pathology + apoe4*has_diag  + sex + rsqa__MeanFD")),
                                    covariates = c("sex", "rsqa__MeanFD"),
                                    right_term_side = FALSE,
                                    #plt_title = "Cognitively Unimpaired (Aβ-)",
                                    #brain_names = c("Age", "Tau Pathology", "ε4 carrier", "Tau×ε4"),
                                    plt_subtitle = TRUE,
                                    rectangle = TRUE)
apoe_int$plot
biofinder_df %>% filter(fmri_bl, abnorm_ab == 0, tau_pathology > 0.15) |> 
  select(apoe4) |> table()

text_label <- "APOE ε4\nFALSE  TRUE\nn=219  n=135"
apoe_dist <- biofinder_df %>% filter(fmri_bl, abnorm_ab == 0, tau_pathology > 0.2) |> 
  select(apoe4, age) |> 
  ggplot(aes(age, fill = apoe4)) +
  geom_histogram(color = "black", position = "dodge") +
  annotate(
    "text",
    x = 20, y = 20,                # adjust position as needed
    label = text_label,
    hjust = 0, vjust = 1,
    family = "mono",
    size = 7
  ) +
  ggtitle("Αβ-, Tau Pathology < 0.15") +
  labs(y = "Count", x = "Age") +
  theme(
    legend.position = "bottom",
    plot.background = element_rect(color = "black", fill = NA))


apoe_p <- ggdraw() +
  draw_plot(apoe$plot, x = 0, y = 0.5, width = 0.495, height = 0.5) +
  draw_plot(apoe_abneg$plot, x = 0.505, y = 0.5, width = 0.495, height = 0.5) +
  draw_plot(apoe_int$plot, x = 0, y = 0.0, width = 0.7, height = 0.49) +
  draw_plot(apoe_dist, x = 0.71, y = 0.0, width = 0.29, height = 0.49)

scaling_factor <- 3
p_name <- "apoe_plot.png"
ggsave(file.path(figure_path, p_name), apoe_p,
       width = img_width*scaling_factor, height = img_width*0.8*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)



apoe_p_and_tauab <- ggdraw() +
  draw_plot(tau_ab_cu_adj$plot + plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=20)),
                                                               plot.subtitle = element_text(margin = margin(l=20)))),
            x = 0, y = 0.505, width = 0.495, height = 0.495) +
  draw_plot(tau_ab_cu$plot + plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=20)),
                                                           plot.subtitle = element_text(margin = margin(l=20))))
            , x = 0.505, y = 0.505, width = 0.495, height = 0.495) +
  draw_plot(apoe$plot+ plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=20)),
                                                     plot.subtitle = element_text(margin = margin(l=20))))
            , x = 0, y = 0, width = 0.495, height = 0.495) +
  draw_plot(apoe_abneg$plot+ plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=20)),
                                                           plot.subtitle = element_text(margin = margin(l=20))))
            , x = 0.505, y = 0, width = 0.495, height = 0.495) +
  draw_plot_label("A") +
  draw_plot_label("B", x = 0.505) +
  draw_plot_label("C", x = 0.0, y = 0.495) +
  draw_plot_label("D", x = 0.505, y = 0.495) 

scaling_factor <- 3
p_name <- "apoe_plot_v1_w_tau_ab.png"
ggsave(file.path(figure_path, p_name), apoe_p_and_tauab,
       width = img_width*scaling_factor, height = img_width*0.8*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)


###########################
# Atrophy + vascular
###########################

new_data <- read_delim("stuff_for_revisions/jonathan__20250926_104238.csv")

x <- new_data |> 
  mutate(mri_date = as.Date(str_extract(csv_ct_composites__index, "(?<=__).*"), format = "%Y%m%d")) |> 
  filter(!is.na(mri_date)) |> 
  mutate(
    mri_bl = mri_date == min(mri_date),
    has_longitudinal = n()>1, .by = "sid"
  ) |> 
  select(sid, mri_bl, visit = Visit, visit_date, mri_date,
                        contains("microbleeds"),
                        contains("lacunes"),
                        contains("vol"),
                        samseg_wmhs_WMH_total_mm3,
                        (contains("ct") & contains("lr"))) |> 
  #select(!contains("adsign")) |> 
  mutate(ct_mean_lr = rowMeans(across(contains("ct") & contains("lr") & contains("surfwt")), na.rm = TRUE),
         microbleeds = factor(as.logical(Cerebral_microbleeds_dich)),
         hippo_vol = -rowMeans(across(contains("Hippocampus") & contains("vols"))), 
         lacunes = factor(as.logical(Lacunes_dich)),
         wmhs_tot = samseg_wmhs_WMH_total_mm3) |> 
  filter(mri_bl)


biofinder_new <- biofinder_df |> filter(fmri_bl) |> 
  inner_join(x, by = c("sid"))

ct_adsign <- plot_gradient_relationships(biofinder_new %>% filter(fmri_bl, !is.na(pathology_ad)) |>
                                           mutate(CT_adsig = -ct_adsign_lr),
                                         gradient_data = grad_df |> filter(study== "biofinder"),
                                         gradients = c(1, 3),
                                         gradient_colors = gradient_cols,
                                         list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                         #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
                                         add_shade = TRUE,
                                         shade_alpha = 0.0075,
                                         shade_size = 0.001,
                                         empty_row_height = -0.15,
                                         base_size_ = 16,
                                         vect = TRUE,
                                         mod_formula = formula(paste0(" ~ age + tau_pathology + CT_adsig + wmhs_tot + sex + rsqa__MeanFD")),
                                         covariates = c("sex", "rsqa__MeanFD"),
                                         plt_title = "Full Sample ",
                                         brain_names = c("Age", "Tau Pathology", "-CT AD signature"),
                                         plt_subtitle = TRUE,
                                         rectangle = TRUE)

ct_adsign$plot


cerebro_vascular <- plot_gradient_relationships(biofinder_new %>% filter(fmri_bl) |> 
                                                  mutate(
                                                    has_diag = diagnosis %in% c("MCI", "AD"),
                                                    CT_adsig = -ct_adsign_lr,
                                                         a_syn = factor(a_syn)), 
                                                gradient_data = grad_df |> filter(study== "biofinder"), 
                                                gradients = c(1, 3),
                                                gradient_colors = gradient_cols,
                                                list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                                #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
                                                add_shade = TRUE, 
                                                shade_alpha = 0.0075,
                                                shade_size = 0.001,
                                                empty_row_height = -0.15,
                                                base_size_ = 16,
                                                vect = TRUE,
                                                mod_formula = formula(paste0(" ~ age + tau_pathology + CT_adsig + wmhs_tot + lacunes + microbleeds + a_syn + has_diag + sex + rsqa__MeanFD")),
                                                covariates = c("sex", "rsqa__MeanFD"),
                                                plt_title = "Full sample with available variables",
                                                brain_names = c("Age", "Tau Pathology", "-CT AD signature", "WM hyperint", "Has Lacunes", "Has Microbleeds", "αSyn+"),
                                                brain_title_vjust = -5,
                                                brain_subtitle_vjust = -7.5,
                                                plt_subtitle = TRUE,
                                                group_n_title = TRUE,
                                                rectangle = TRUE)


biofinder_new |> 
  ggplot(aes(wmhs_tot, ct_adsign_lr)) +
  geom_point()

# ct_p <- ggdraw() +
#   draw_plot(ct_adsign$plot, x = 0, y = 0.5, width = 0.75, height = 0.49) +
#   draw_plot(ct_full$plot, x = 0, y = 0.0, width = 0.99, height = 0.49) 

img_width <-  180 / 25.4
scaling_factor <-  3
magick_geom_scaling <- paste0(100/scaling_factor, "%x", 100/scaling_factor, "%")

p_name <- "atrophy_and_cerebrovasc.png"
ggsave(file.path(figure_path, p_name), cerebro_vascular$plot,
       width = img_width*scaling_factor, height = img_width*0.465*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)


###########################
# Cognition
###########################

###############
# Linear stuff
###############
new_cog <- read_delim("stuff_for_revisions/jonathan__20250926_104238.csv")

cogs <- new_cog |> select(aqt_form:`letter_s_45-60`)

#cog <- read_delim("stuff_for_revisions/cognitive_tests-2025-08-29.csv") 
nic <- read_delim("stuff_for_revisions/nicola__20240326_143046.csv") 
cog <- nic |> select(1:5, mPACC_v1:fcsrt_immediate,
                     # adas_delayed_word_recall, adas_immediate_word_recall_average,
                     # letter_s, trailmaking_a, trailmaking_b, symbol_digit,
                     # animal_fluency, vosp_cube, bnt_15_2,
                     cognitive_test_date
                     )

cog_bl <- cog |> 
  filter(!is.na(cognitive_test_date)) |> 
  filter(cognitive_test_date == min(cognitive_test_date), .by = "sid") |> 
  # drop_na(adas_delayed_word_recall, 
  #         adas_immediate_word_recall_average, 
  #         letter_s, trailmaking_a, #bnt_15_2, #trailmaking_b,
  #         symbol_digit, animal_fluency, vosp_cube
  #         ) |>
  filter(Visit==0)

biofinder_cog <- cog_bl |> 
  inner_join(biofinder_df |> filter(fmri_bl), by = "sid") 
  

cog_mat <- biofinder_cog |> 
  select(
    sid,
    adas_delayed_word_recall,
    adas_immediate_word_recall_average,
    symbol_digit,
    trailmaking_a, 
    trailmaking_b,
    letter_s,
    animal_fluency,
    vosp_cube,
    bnt_15_2,
    aqt_color_form,
    
  ) |>
  mutate(
    trailmaking_a = -trailmaking_a,
    trailmaking_b = -trailmaking_b,
    adas_delayed_word_recall = -adas_delayed_word_recall,
    adas_immediate_word_recall_average = -adas_immediate_word_recall_average,
    aqt_color_form = -aqt_color_form
  ) |> drop_na() |> column_to_rownames("sid")
  



cog_fac <- psych::omega(cog_mat, nfactors = 4)

om_factors <- cog_fac[["scores"]]
refs <- cog_mat |> select(adas_delayed_word_recall, vosp_cube, animal_fluency, trailmaking_b)


rename_factors_by_corr <- function(omega_factors, ref_vars) {
  global_cog <- omega_factors[, 1]
  omega_factors <- omega_factors[, -1]
  label_map <- c(
    adas_delayed_word_recall = "memory",
    trailmaking_b            = "executive",
    animal_fluency           = "language",
    vosp_cube                = "visuospatial"
  )
  cor_mat <- cor(ref_vars, omega_factors, use = "pairwise.complete.obs")
  strongest_ref <- apply(abs(cor_mat), 2, function(x) names(x)[which.max(x)])
  new_names <- label_map[strongest_ref]
  colnames(omega_factors) <- new_names
  omega_factors <- cbind(omega_factors, global_cog)
  return(omega_factors)
}

renamed_omega <- rename_factors_by_corr(om_factors, refs)


biofinder_cog_comp <- biofinder_cog |> 
  inner_join(renamed_omega |> as_tibble(rownames = NA) |> rownames_to_column("sid")) 


cog_diff_plot <- function() {
  fit <- biofinder_cog_comp |> 
    lm(memory ~ scale(age) + scale(pathology_ad) + sex , data = _) 
  
  # Extract coefficients
  coef_df <- broom::tidy(fit) |> 
    filter(term %in% c("scale(age)", "scale(pathology_ad)")) |> 
    select(term, estimate, std.error) |> 
    mutate(ci_low = estimate - 1.96*std.error,
           ci_high = estimate + 1.96*std.error,
           term = recode(term,
                         "scale(age)" = "Age",
                         "scale(pathology_ad)" = "Pathology"))
  
  # Plot
  mem_plot <- ggplot(coef_df, aes(x = term, y = estimate)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(y = "Memory", x = "") +
    theme_bw(base_size = 14) 
  
  
  
  fit_exec <- biofinder_cog_comp |> 
    lm(executive ~ scale(age) + scale(pathology_ad) + sex, data = _) 
  
  # Extract coefficients
  coef_df <- broom::tidy(fit_exec) |> 
    filter(term %in% c("scale(age)", "scale(pathology_ad)")) |> 
    select(term, estimate, std.error) |> 
    mutate(ci_low = estimate - 1.96*std.error,
           ci_high = estimate + 1.96*std.error,
           term = recode(term,
                         "scale(age)" = "Age",
                         "scale(pathology_ad)" = "Pathology"))
  
  # Plot
  exec_plot <- ggplot(coef_df, aes(x = term, y = estimate)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(y = "Executive", x = "") +
    theme_bw(base_size = 14) 
  
  library(patchwork)
  cog_eff_plot <- wrap_plots(list(mem_plot, exec_plot), ncol = 1, guides = "collect", axes = "collect")  +
    plot_annotation(
      title = "Coefficients of cognitive effects",
      subtitle = "Memory ~ scale(age) + scale(AD pathology) + sex \nExecutive ~ scale(age) + scale(AD pathology) + sex",
      theme = theme(plot.background = element_rect(color = "black", fill = NA))) &
    theme(axis.text.x = element_text(size = rel(1.5)),
          axis.title.y = element_text(size = rel(1.5)),
          plot.subtitle = element_text(family = "mono", size = rel(0.75)))
  cog_eff_plot
}

cog_diff <- cog_diff_plot()


biofinder_cog_comp |> 
  ggplot(aes(memory, fill = diagnosis)) +
  geom_density(alpha = 0.3)

biofinder_cog_comp |> drop_na(age) |> 
  ggplot(aes(memory, fill = cut(age, 4))) +
  geom_density(alpha = 0.3)

biofinder_cog_comp |> 
  ggplot(aes(global_cog, memory, color = diagnosis)) +
  geom_point(alpha = 0.8)

biofinder_cog_comp |> 
  ggplot(aes(pathology_ad, memory)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm")

biofinder_cog_comp |> drop_na(age) |> 
  ggplot(aes(executive, fill = cut(age, 4))) +
  geom_density(alpha = 0.3)

biofinder_cog_comp |> 
  ggplot(aes(age, executive)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm")

biofinder_cog_comp |> 
  ggplot(aes(global_cog, executive, color = diagnosis)) +
  geom_point(alpha = 0.8) 


cognitive_composites <- plot_gradient_relationships(biofinder_cog_comp |> 
                                                      filter(fmri_bl), 
                                                gradient_data = grad_df |> filter(study== "biofinder"), 
                                                gradients = c(1, 3),
                                                gradient_colors = gradient_cols,
                                                list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                                #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
                                                add_shade = TRUE, 
                                                shade_alpha = 0.0075,
                                                shade_size = 0.001,
                                                empty_row_height = -0.15,
                                                base_size_ = 18,
                                                vect = TRUE,
                                                mod_formula = formula(paste0(" ~ age + pathology_ad + global_cog + executive + memory + sex + rsqa__MeanFD")),
                                                covariates = c("age", "pathology_ad", "sex", "rsqa__MeanFD"),
                                                plt_title = "Full Sample, with available data",
                                                #brain_names = c("Executive Cog", "Memory Cog"),
                                                plt_subtitle = TRUE,
                                                group_n_title = TRUE,
                                                rectangle = TRUE)

source("src/mass_mediation_src.R")

df <- biofinder_cog_comp |> filter(fmri_bl, 
                                   ) |> 
  drop_na(age, pathology_ad, sex) |> 
  mutate(has_diag = factor(diagnosis %in% c("MCI", "AD")),
         memory = resid(lm(memory ~ global_cog)),
         executive = resid(lm(executive ~ global_cog))
         ) |> 
  #left_join(cog_edit) |> 
  inner_join(fc_measures_bf$affinity |> as_tibble(rownames = NA) |> 
               rownames_to_column("image_file") |> 
               janitor::clean_names()) 

parcs <- df |> select(starts_with("x7")) |> colnames()

pw_mediation <- plot_mediation_gradients(subject_data = df, 
                                      treatments = c("age", "pathology_ad"), 
                                      outcomes = c("executive", "memory"),
                                      covariates = list(c("pathology_ad", "sex", "rsqa__MeanFD"), 
                                                        c("age", "sex", "rsqa__MeanFD")),
                                      gradient_data = grad_df |> filter(study=="biofinder"), 
                                      gradients = c(1, 3),
                                      parcels = parcs,
                                      gradient_colors = gradient_cols,
                                      plt_title = "Parcelwise Mediation Analysis",
                                      terms_nice = c("age", "AD path"),
                                      plot_spacing = 0.3,
                                      empty_row_height = -0.1,
                                      base_size_ = 16,
                                      add_shade = TRUE,
                                      shade_size = 0.01, 
                                      shade_alpha = 0.01,
                                      rectangle = TRUE)


omega_plot <- ggplotify::as.ggplot(~omega.diagram(cog_fac, size = 500, marg = c(0, 0, 0, 0), main = "", 
                                                  labels = c(adas_delayed_word_recall = "ADAS Delayed Recall",
                                                             adas_immediate_word_recall_average = "ADAS Immediate Recall",
                                                             bnt_15_2 = "Boston Naming",
                                                             vosp_cube = "VOSP Cube",
                                                             animal_fluency = "Animal Fluency",
                                                             letter_s = "Letter Fluency",
                                                             trailmaking_a = "Trail Making A",
                                                             trailmaking_b = "Trail Making B",
                                                             symbol_digit = "Symbol Digit",
                                                             aqt_color_form = "AQT Color–Form"),
                                                  flabels = c("pl", "Mem", "Visuo", "Lang", "Exec"))) + 
  ggtitle("Bifactor model of cognition (Schmid–Leiman Decomposition)") +
  theme(plot.background = element_rect(color = "black", linewidth = 1),
        plot.title = element_text(size = 22, hjust = 0.65, margin = margin(10, 0, -20, 0)),
        plot.margin = margin(0, 0, 0, -150))

ggdraw() +
  draw_plot(omega_plot, x = 0, y = 0.52, width = 0.5, height = 0.48) +
  draw_plot(cog_diff, x = 0.51, y = 0.52, width = 0.49, height = 0.48)



cognition_full <- 
  ggdraw() +
  draw_plot(omega_plot, x = 0, y = 0.52, width = 0.5, height = 0.48) +
  draw_plot(cog_diff, x = 0.51, y = 0.52, width = 0.49, height = 0.48) +
  draw_plot(cognitive_composites$plot, x = 0, y = 0, width = 0.545, height = 0.51) +
  draw_plot(pw_mediation, x = 0.555, y = 0, width = 0.445, height = 0.51) + 
  draw_plot_label("A", x = 0.48, y = 1, size = 18) +
  draw_plot_label("B", x = 0.98, y = 1, size = 18) +
  draw_plot_label("C", x = 0.525, y = 0.51, size = 18) +
  draw_plot_label("D", x = 0.98, y = 0.51, size = 18) 

scaling_factor <- 3
p_name <- "cognition_full.png"
ggsave(file.path(figure_path, p_name), cognition_full,
       width = img_width*scaling_factor, height = img_width*0.8*scaling_factor, 
       units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)








###########################
# Subtypes
###########################

# one_vs_all_contrasts <- function(f) {
#   levs <- levels(f)
#   k <- length(levs)
#   contr <- matrix(0, nrow = k, ncol = k)
#   for (i in seq_len(k)) {
#     contr[, i] <- ifelse(seq_len(k) == i, 1, -1/(k - 1))
#   }
#   colnames(contr) <- paste0(levs, "_vs_all")
#   rownames(contr) <- levs
#   contr
# }
# 
# contrasts(bf_subs$BF_subtype) <- one_vs_all_contrasts(bf_subs$BF_subtype)
# 
# contrasts(bf_subs$BF_subtype) <- contr.sum(4)

tau_subs <- read_delim("stuff_for_revisions/Tau_subtypes_20230925.csv")

sub_type_plot <- list()

for(st in 1:4) {
  
  bf_subs <- biofinder_df |> filter(fmri_bl) |> 
    left_join(tau_subs |> select(sid, BF_subtype, BF_stage, ml_subtype, MT_subtype)) |> 
    filter(BF_subtype != 0) |> 
    mutate(
      BF_subtype = ifelse(BF_subtype == st, paste0("_", st), "Rest"),
      BF_subtype = factor(BF_subtype, levels = c("Rest", paste0("_", st))),
      ml_subtype = factor(ml_subtype),
      MT_subtype = factor(MT_subtype)) 
  
  sub_plot <- plot_gradient_relationships(bf_subs, 
                                          gradient_data = grad_df |> filter(study== "biofinder"), 
                                          gradients = c(1, 3),
                                          gradient_colors = gradient_cols,
                                          list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                          add_shade = FALSE, 
                                          shade_alpha = 0.0075,
                                          shade_size = 0.001,
                                          empty_row_height = -0.15,
                                          base_size_ = 16,
                                          vect = TRUE,
                                          mod_formula = formula(" ~ age + BF_subtype + sex + rsqa__MeanFD"),
                                          covariates = c("sex", "rsqa__MeanFD"),
                                          plt_title = paste0("Subtype ", st, "vs Rest"),
                                          plt_subtitle = FALSE, rectangle = TRUE)
  
  
  sub_type_plot[[st]] <- sub_plot$plot
}


subtypes_vs_rest <- 
  ggdraw() +
  draw_plot(sub_type_plot[[1]], x = 0, y = 0.51, width = 0.495, height = 0.49) +
  draw_plot(sub_type_plot[[2]], x = 0.505, y = 0.51, width = 0.495, height = 0.49) +
  draw_plot(sub_type_plot[[3]], x = 0, y = 0, width = 0.495, height = 0.49) +
  draw_plot(sub_type_plot[[4]], x = 0.505, y = 0, width = 0.495, height = 0.49) 


scaling_factor <- 3
p_name <- "subtypes_vs_rest.png"
ggsave(file.path(figure_path, p_name), subtypes_vs_rest,
       width = img_width*scaling_factor, height = img_width*0.8*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)




###########################
# Interaction
###########################


int_grp1 <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl) |> 
                                          mutate(age_grp  = cut(age, quantile(age, c(0, 0.33, 0.66, 1), na.rm = TRUE))) |> 
                                          rename(path = pathology_ad), 
                                        gradient_data = grad_df |> filter(study== "biofinder"), 
                                        gradients = c(1, 3),
                                        gradient_colors = gradient_cols,
                                        list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                        #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
                                        add_shade = TRUE, 
                                        shade_alpha = 0.0075,
                                        shade_size = 0.001,
                                        empty_row_height = -0.15,
                                        base_size_ = 14,
                                        vect = TRUE,
                                        mod_formula = formula(" ~ path * age_grp + sex + rsqa__MeanFD"),
                                        covariates = c("sex", "rsqa__MeanFD"),
                                        r2_size = rel(6),
                                        filter_criteria = quo(),
                                        show_networks = FALSE,
                                        tag_prefix = "",
                                        tag_sep = "",
                                        layout_construction = "horizontal",
                                        include_gradient_plots = TRUE,
                                        right_term_side = FALSE,
                                        plt_title = "Full sample ",
                                        plt_subtitle = TRUE,
                                        cache_runs = FALSE)

x <- biofinder_df %>% filter(fmri_bl) |> 
  mutate(path_ntile  = ntile(pathology_ad, 3)) |> 
  group_by(path_ntile) |> 
  mutate(path_grp = paste0("(",paste0(round(range(pathology_ad), 2), collapse = "-"), ")"), 
         path_grp = ifelse(path_grp == "(NA-NA)", NA, path_grp))

int_grp2 <- plot_gradient_relationships(x, 
                                        gradient_data = grad_df |> filter(study== "biofinder"), 
                                        gradients = c(1, 3),
                                        gradient_colors = gradient_cols,
                                        list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                        #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
                                        add_shade = TRUE, 
                                        shade_alpha = 0.0075,
                                        shade_size = 0.001,
                                        empty_row_height = -0.15,
                                        base_size_ = 14,
                                        vect = TRUE,
                                        mod_formula = formula(" ~ age * path_grp + sex + rsqa__MeanFD"),
                                        covariates = c("sex", "rsqa__MeanFD"),
                                        r2_size = rel(6),
                                        filter_criteria = quo(),
                                        show_networks = FALSE,
                                        tag_prefix = "",
                                        tag_sep = "",
                                        layout_construction = "horizontal",
                                        include_gradient_plots = TRUE,
                                        right_term_side = FALSE,
                                        plt_title = "Full sample ",
                                        plt_subtitle = TRUE,
                                        cache_runs = FALSE)


int_plot <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl), 
                                           gradient_data = grad_df |> filter(study== "biofinder"), 
                                           gradients = c(1, 3),
                                           gradient_colors = gradient_cols,
                                           list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                           #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
                                           add_shade = TRUE, 
                                           shade_alpha = 0.0075,
                                           shade_size = 0.001,
                                           empty_row_height = -0.15,
                                           base_size_ = 21,
                                           vect = TRUE,
                                           mod_formula = formula(" ~ scale(age) * pathology_ad + sex + rsqa__MeanFD"),
                                           covariates = c("sex", "rsqa__MeanFD"),
                                           r2_size = rel(6),
                                           filter_criteria = quo(),
                                           show_networks = FALSE,
                                           tag_prefix = "",
                                           tag_sep = "",
                                           layout_construction = "horizontal",
                                           include_gradient_plots = TRUE,
                                           right_term_side = FALSE,
                                           plt_title = "Full sample ",
                                           plt_subtitle = TRUE,
                                           cache_runs = FALSE)


int_plot2 <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, diagnosis %in% c("Normal", "SCD"), abnorm_ab == 0), 
                                         gradient_data = grad_df |> filter(study== "biofinder"), 
                                         gradients = c(1, 3),
                                         gradient_colors = gradient_cols,
                                         list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                         #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
                                         add_shade = TRUE, 
                                         shade_alpha = 0.0075,
                                         shade_size = 0.001,
                                         empty_row_height = -0.15,
                                         base_size_ = 21,
                                         vect = TRUE,
                                         mod_formula = formula(" ~ scale(age) * pathology_ad + sex + rsqa__MeanFD"),
                                         covariates = c("sex", "rsqa__MeanFD"),
                                         r2_size = rel(6),
                                         filter_criteria = quo(),
                                         show_networks = FALSE,
                                         tag_prefix = "",
                                         tag_sep = "",
                                         layout_construction = "horizontal",
                                         include_gradient_plots = TRUE,
                                         right_term_side = FALSE,
                                         plt_title = "CU, Ab- ",
                                         plt_subtitle = TRUE,
                                         cache_runs = FALSE)

pw <- ggdraw() +
  draw_plot(int_plot$plot, x = 0, y = 0.5, width = 1, height = 0.49) +
  draw_plot(int_plot2$plot, x = 0, y = 0, width = 1, height = 0.49) 


img_width <-  180 / 25.4
scaling_factor <-  3
magick_geom_scaling <- paste0(100/scaling_factor, "%x", 100/scaling_factor, "%")

p_name <- "figure_agexpath.png"
ggsave(file.path(figure_path, p_name), pw,
       width = img_width*scaling_factor, height = img_width*1.2*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")


path_x <- fc_measures_bf$affinity |> as_tibble(rownames=NA) |> 
  rownames_to_column("image_file") |> 
  pivot_longer(-image_file, names_to = "region", values_to = "affinity") |> 
  left_join(grad_df |> filter(study == "biofinder") |> select(region, starts_with("grad"))) |> 
  group_by(image_file) |> 
  arrange(image_file, desc(gradient1)) |> 
  mutate(grad_side = case_when(
    row_number() <= 100 ~ "associative",
    row_number() > n() - 100 ~ "sensory",
    TRUE ~ NA_character_
  )) |> 
  filter(!is.na(grad_side)) |> 
  ungroup() |> 
  group_by(image_file, grad_side) |> 
  summarise(affinity = mean(affinity)) |> 
  ungroup() |> 
  right_join(biofinder_df |> filter(fmri_bl))

path_plot <- path_x |> 
  mutate(age_grp = cut(age, quantile(age, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)),
         above_median_age = age > median(age, na.rm = TRUE)) |> 
  drop_na(above_median_age) |> 
  ggplot(aes(pathology_ad, affinity, linetype = above_median_age, color = grad_side)) +
  geom_smooth(fill = "lightgray") +
  geom_point(alpha = 0.1) +
  scale_color_manual(values = c("associative" = gradient_cols[2, 1], "sensory" = gradient_cols[1, 1])) +
  labs(x = "AD Pathology", y = "Affinity", linetype = "Old age", color = "Domain") +
  theme(legend.position = "right")

age_x <- fc_measures_bf$affinity |> as_tibble(rownames=NA) |> 
  rownames_to_column("image_file") |> 
  pivot_longer(-image_file, names_to = "region", values_to = "affinity") |> 
  left_join(grad_df |> filter(study == "biofinder") |> select(region, starts_with("grad"))) |> 
  group_by(image_file) |> 
  arrange(image_file, desc(gradient3)) |> 
  mutate(grad_side = case_when(
    row_number() <= 100 ~ "communicative",
    row_number() > n() - 100 ~ "executive",
    TRUE ~ NA_character_
  )) |> 
  filter(!is.na(grad_side)) |> 
  ungroup() |> 
  group_by(image_file, grad_side) |> 
  summarise(affinity = mean(affinity)) |> 
  ungroup() |> 
  right_join(biofinder_df |> filter(fmri_bl))

age_plot <- age_x |> 
  #mutate(pathology_ad = pathology_ad + 0.0000001) |> 
  mutate(above_median_path = pathology_ad > median(pathology_ad, na.rm = TRUE)) |> 
  drop_na(pathology_ad) |> 
  drop_na(above_median_path) |> 
  ggplot(aes(age, affinity, linetype = above_median_path, color = grad_side)) +
  geom_smooth(fill = "lightgray") +
  geom_point(alpha = 0.1) +
  scale_color_manual(values = c("communicative" = gradient_cols[2, 3], "executive" = gradient_cols[1, 3])) +
  labs(x = "Age", y = "Affinity", linetype = "High path", color = "Domain") +
  theme(legend.position = "right")

  
path_plot + age_plot




###########################
# For aff vs conn
###########################

ave_con <- readRDS(file.path(atlas_dir, "average_connectome_normalyoung.rds"))


zero_out_mat <- function(mat, thresh) {
  library(matrixStats)
  thresholds <- colQuantiles(abs(mat), probs = thresh)
  mat[abs(mat) <= thresholds[col(mat)]] <- 0
  return(mat)
}



connectome <- ave_con
L <- matrix(connectome[!diag(nrow(connectome))], nrow = nrow(connectome)-1) # mask it
L <- zero_out_mat(L, 0.75)
L[L<0] <- 0
L <- proxyC::simil(L, margin = 2, method = "cosine")#, use_nan = TRUE
L <- as.matrix(L)


con <- ave_con
con[con<0] <- 0


aff_con <- tibble(affinity = colMeans(L), con = colMeans(con)) |> 
  mutate(index = row_number()) |> 
  rowwise() |> 
  mutate(diff_ = affinity - con) |> 
  ungroup() |> 
  inner_join(roi_data |> mutate(index=row_number())) 
# |> 
#   pivot_longer(-c(index, diff_), names_to = "type", values_to = "value") |> 
#   filter(value != 1)

p1 <- aff_con |> 
  ggplot(aes(x = diff_, y = region, col = col)) +
  geom_col(show.legend = TRUE) +
  scale_color_identity() +
  labs(x = "Affinity - Connectivity (diff)", y = "Parcel") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) 

p2 <- aff_con |> 
  ggplot(aes(con, affinity, colour = col)) +
  geom_point(alpha = 0.7) +
  scale_color_identity() +
  labs(y = "Affinity", x = "Connectivity") +
  theme_minimal(base_size = 14)

p1 + p2


##################################
# Causality comment
##################################

n <- 1000
mat <- matrix(0, n, n)

diag(mat) <- 1
cors <- runif(999, 0.7, 0.9)
mat[cbind(1:(n-1), 2:n)] <- cors
mat[cbind(2:n, 1:(n-1))] <- cors

test_str <- colMeans(mat)
comp1 <- prcomp(mat, scale. = TRUE)$x[1, ]

cor(test_str, comp1)

e <- eigen(mat)
pc1 <- e$vectors[,1]

plot(strength, type = "l", main = "PC1 of Local-Only Correlation Matrix")

ave_con <- readRDS(file.path(atlas_dir, "average_connectome_normalyoung.rds"))

zero_out_mat <- function(mat, thresh) {
  library(matrixStats)
  thresholds <- colQuantiles(abs(mat), probs = thresh)
  mat[abs(mat) <= thresholds[col(mat)]] <- 0
  return(mat)
}


connectome <- ave_con
L <- matrix(connectome[!diag(nrow(connectome))], nrow = nrow(connectome)-1) # mask it
L <- zero_out_mat(L, 0.75)
L[L<0] <- 0


  L <- Matrix::Matrix(L, sparse = TRUE)
  L <- proxyC::simil(L, margin = 2, method = "cosine"#, use_nan = TRUE
  )
  L <- as.matrix(L)
  diag(L) <- NA


aff <- colMeans(L, na.rm = TRUE)

ave_con[ave_con<0] <- 0

strength <- colMeans(ave_con)

pca_result <- prcomp(ave_con, scale. = TRUE)
g1 <- pca_result$x[, 1]
g2 <- pca_result$x[, 2]
g3 <- pca_result$x[, 3]

g1 <- grad_df |> filter(study == "biofinder") |> pull(gradient1)
g2 <- grad_df |> filter(study == "biofinder") |> pull(gradient2)
g3 <- grad_df |> filter(study == "biofinder") |> pull(gradient3)

cor(aff, g1)
cor(aff, g2)
cor(aff, g3)

cor(tibble(g1, g2, g3, affinity = aff, strength = strength)) |> round(5)






lm_t <- lm(t_path ~ aff_mean); t_res <- resid(lm_t)
lm_g <- lm(grad1 ~ aff_mean); g_res <- resid(lm_g)
r_partial <- cor(t_res, g_res)





subj_corr_g1 <- apply(fc_measures_bf$affinity[biofinder_df %>% filter(fmri_bl) |> pull(image_file), ], 1 , function(x) cor(x, g1))
subj_corr_g3 <- apply(fc_measures_bf$affinity[biofinder_df %>% filter(fmri_bl) |> pull(image_file), ], 1 , function(x) cor(x, g3))



orig_analysis <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl), 
                                             gradient_data = grad_df |> filter(study== "biofinder"), 
                                             gradients = c(1, 2, 3),
                                             gradient_colors = gradient_cols,
                                             list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                             empty_row_height = -0.1,
                                             base_size_ = 21,
                                             vect = TRUE,
                                             mod_formula = formula(paste0(" ~ age + pathology_ad + sex  + rsqa__MeanFD")),
                                             covariates = c("rsqa__MeanFD"),
                                             r2_size = rel(6),
                                             filter_criteria = quo(),
                                             show_networks = FALSE,
                                             tag_prefix = "",
                                             tag_sep = "",
                                             layout_construction = "horizontal",
                                             include_gradient_plots = TRUE,
                                             right_term_side = FALSE,
                                             plt_title = "",
                                             cache_runs = FALSE)



mean_aff <- aff
mean_aff <- fc_measures_bf$affinity[biofinder_df %>% filter(fmri_bl) |> pull(image_file), ] |> colMeans()

# t_path <- orig_analysis$tmaps |> filter(term == "pathology_ad") |> pull(nodal_affinity)
# lm_t <- lm(t_path ~ mean_aff); t_res <- resid(lm_t)
# lm_g <- lm(g1 ~ mean_aff); g_res <- resid(lm_g)
# cor(t_res, g_res)
# 
# 
# t_age <- orig_analysis$tmaps |> filter(term == "age") |> pull(nodal_affinity)
# lm_ta <- lm(t_age ~ mean_aff); ta_res <- resid(lm_ta)
# lm_g3 <- lm(g3 ~ mean_aff); g3_res <- resid(lm_g3)
# cor(ta_res, g3_res)

cor_path <- ppcor::pcor(data.frame(Gradient1 = g1, tmap_path = t_path, Affinity = mean_aff))$estimate |> 
  as.data.frame() |> 
  rownames_to_column("Var1") |> 
  pivot_longer(-Var1, names_to = "Var2", values_to = "correlation") |> 
  mutate(
    Var1 = factor(Var1, levels = unique(Var1)),
    Var2 = factor(Var2, levels = rev(unique(Var1)))  # reverse order
  ) |> 
  ggplot(aes(x = Var1, y = Var2, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = muted("blue"), high = muted("red"), midpoint = 0,   limits = c(-1, 1)) +
  scale_x_discrete(position = "top") +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 6, color = "darkgray") +
  coord_fixed() +
  labs(title = "Partial correlation pathology",
       x = NULL, y = NULL, fill = "Correlation") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "",
    #axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())


x <- orig_analysis$plot[[5]] + ggtitle("")
library(grid)
library(ggplotify)

grob1 <- as.grob(x)

diag_vars = c("Path_tval")

cor_path <- cor_path +
  annotation_custom(grob1,
                    xmin = which(levels(cor_path$data$Var1) == diag_vars[1]) - 0.5,
                    xmax = which(levels(cor_path$data$Var1) == diag_vars[1]) + 0.5,
                    ymin = which(levels(cor_path$data$Var2) == diag_vars[1]) - 0.5,
                    ymax = which(levels(cor_path$data$Var2) == diag_vars[1]) + 0.8
  ) 

cow <- ggdraw() +
  draw_plot(cor_path) +
  draw_plot(x, x = 0.1, y = 0.8, width = )


library(ggplotify)


cor_age <- ppcor::pcor(data.frame(Gradient3 = g3, tmap_age = t_age, Affinity = mean_aff))$estimate |> 
  as.data.frame() |> 
  rownames_to_column("Var1") |> 
  pivot_longer(-Var1, names_to = "Var2", values_to = "correlation") |> 
  mutate(
    Var1 = factor(Var1, levels = unique(Var1)),
    Var2 = factor(Var2, levels = rev(unique(Var1)))  # reverse order
  ) |> 
  ggplot(aes(x = Var1, y = Var2, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = muted("blue"), high = muted("red"), midpoint = 0,   limits = c(-1, 1)) +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 6, color = "darkgray") +
  scale_x_discrete(position = "top") +
  coord_fixed() +
  labs(title = "Partial correlation age",
       x = NULL, y = NULL, fill = "Correlation") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "",
    #axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())


data.frame(G1_corr = subj_corr_g1, G3_corr = subj_corr_g3) |> 
  pivot_longer(everything(), names_to = "Gradient", values_to = "Correlation") |> 
  ggplot(aes(x = Correlation)) +
  geom_histogram(color = "black") +
  labs(title = "Individual correlation with gradients" , y = "") +
  facet_wrap(~Gradient) -> subject_corr
  
sexed <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl), 
                                     gradient_data = grad_df |> filter(study== "biofinder"), 
                                     gradients = c(1, 3),
                                     gradient_colors = gradient_cols,
                                     list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                     add_shade = TRUE, 
                                     shade_alpha = 0.0075,
                                     shade_size = 0.001,
                                     empty_row_height = -0.15,
                                     base_size_ = 16,
                                     vect = TRUE,
                                     mod_formula = formula(" ~ age + pathology_ad + education + sex + rsqa__MeanFD"),
                                     covariates = c("rsqa__MeanFD"),
                                     plt_title = "Full Sample",
                                     brain_names = c("Age", "AD Pathology", "Education", "Sex"),
                                     plt_subtitle = TRUE, rectangle = TRUE)



sexd <- sexed$plot + plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=20)),
                                           plot.subtitle = element_text(margin = margin(l=20))))

part_cor <- (cor_path + cor_age) / subject_corr


caus_plot <-
  ggdraw() +
  draw_plot(sexd, x = 0, y = 0, width = 0.6, height = 1) +
  draw_plot_label("A") +
  draw_plot(part_cor, x = 0.6, y = 0, width = 0.4, height = 1) +
  draw_plot_label("B", x = 0.6) +
  draw_plot_label("C", x = 0.6, y = 0.55) 

scaling_factor <- 3
p_name <- "caus_plot.png"
ggsave(file.path(figure_path, p_name), caus_plot,
       width = img_width*scaling_factor, height = img_width*0.425*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)

###########################
# Gradients from AD
###########################


# ad_con <- con_cube_bf[, , biofinder_df |> filter(fmri_bl, diagnosis %in% c("MCI", "AD")) |> pull(image_file)]
# ad_con <- apply(ad_con, c(1, 2), mean)
# write_rds(ad_con, file.path(atlas_dir, "average_connectome_AD.rds"))
# 
# con_all <- con_cube_bf[, , biofinder_df |> filter(fmri_bl) |> pull(image_file)]
# con_all <- apply(con_all, c(1, 2), mean)
# #write_rds(con_all, file.path(atlas_dir, "average_connectome_ALL.rds"))

con_all <- readRDS(file.path(atlas_dir, "average_connectome_ALL.rds"))

marg_gradients <- read_csv("data/gradients/margulies_gradients.csv")

grad_list_all <- get_gradients(connectome_ests = list(everyone = con_all),
                              reference_gradients = marg_gradients,
                              n_gradients = c(1,2,3),
                              threshold = 0.0,
                              similarity_method = "cosine",
                              on_affinity = FALSE,
                              method = "pca",
                              visualize = FALSE)

grad_all <- grad_list_all$gradients |> filter(study == "everyone") 

all_grad <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, diagnosis %in% c("MCI", "AD")) |> 
                                          mutate(`-mPACC` = -mPACC_v1), 
                                        gradient_data = grad_all, 
                                        gradients = c(1, 3),
                                        gradient_colors = gradient_cols,
                                        list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                        add_shade = TRUE, 
                                        shade_alpha = 0.0075,
                                        shade_size = 0.001,
                                        empty_row_height = -0.15,
                                        base_size_ = 16,
                                        vect = TRUE,
                                        mod_formula = formula(" ~ age + pathology_ad + `-mPACC` + sex + rsqa__MeanFD"),
                                        covariates = c("sex", "rsqa__MeanFD"),
                                        plt_title = "Diagnosed MCI/AD (Aβ+), gradients derived from full sample ",
                                        brain_names = c("Age", "AD Pathology", "-mPACC"),
                                        plt_subtitle = TRUE,
                                        rectangle = TRUE)



scaling_factor <- 2
p_name <- "all_grad_plot.png"
ggsave(file.path(figure_path, p_name), all_grad$plot,
       width = img_width*scaling_factor, height = img_width*0.7*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)



###############################
# Across cohort pathology score and switched gradients
###############################

adni_df_patscore <- adni_df |> filter(fmri_bl) |> 
  drop_na(ab_ratio,
          braak1,
          braak34,
          braak56)

adni_mat <- adni_df_patscore |> filter(fmri_bl) |> 
  select(
    ab_ratio,
    braak1,
    braak34,
    braak56
  ) |> 
  rename(cho12 = braak1,
         cho34 = braak34,
         cho56 = braak56) |> as.matrix()

biofinder_df_patscore <- biofinder_df |>  filter(fmri_bl) |> 
  drop_na(
         ab_ratio,
         cho12,
         cho34,
         cho56
  )

bf_mat <- biofinder_df_patscore |> 
  select(ab_ratio,
         cho12,
         cho34,
         cho56) |> as.matrix()


withr::with_seed(12345, {
  space <- reduce_dimensionality(rbind(adni_mat, bf_mat), "euclidean")
  trajectory <- infer_trajectory(space)
})
trajectory <- reverse_trajectory(trajectory)

traject <- trajectory$time

adni_pat <- traject[1:nrow(adni_mat)]

bf_pat <- traject[(nrow(adni_mat)+1):(nrow(adni_mat) + nrow(bf_mat))]

biofinder_df_patscore$path_all <- bf_pat

adni_df_patscore$path_all <- adni_pat




bf_adni_grad <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl), 
                                        gradient_data = grad_df |> filter(study== "adni"), 
                                        gradients = c(1, 3),
                                        gradient_colors = gradient_cols,
                                        list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                        add_shade = TRUE, 
                                        shade_alpha = 0.0075,
                                        shade_size = 0.001,
                                        empty_row_height = -0.15,
                                        base_size_ = 16,
                                        vect = TRUE,
                                        mod_formula = formula(" ~ age + pathology_ad + sex + rsqa__MeanFD"),
                                        covariates = c("sex", "rsqa__MeanFD"),
                                        plt_title = "Gradients derived from ADNI in BioFINDER sample",
                                        brain_names = c("Age", "AD Pathology"),
                                        plt_subtitle = TRUE,
                                        rectangle = TRUE)

adni_bf_grad <- plot_gradient_relationships(adni_df %>% filter(fmri_bl), 
                                            gradient_data = grad_df |> filter(study== "biofinder"), 
                                            gradients = c(1, 3),
                                            gradient_colors = gradient_cols,
                                            list_of_parcel_data = list(nodal_affinity = fc_measures_adni$affinity),
                                            add_shade = TRUE, 
                                            shade_alpha = 0.0075,
                                            shade_size = 0.001,
                                            empty_row_height = -0.15,
                                            base_size_ = 16,
                                            vect = TRUE,
                                            id_var = "file_func",
                                            mod_formula = formula(" ~ age + pathology_ad + sex + rsqa__MeanFD"),
                                            covariates = c("sex", "rsqa__MeanFD"),
                                            plt_title = "Gradients derived from BioFINDER in ADNI sample",
                                            brain_names = c("Age", "AD Pathology"),
                                            plt_subtitle = TRUE,
                                            rectangle = TRUE)



bf_path_all <- plot_gradient_relationships(biofinder_df_patscore %>% filter(fmri_bl), 
                                            gradient_data = grad_df |> filter(study== "biofinder"), 
                                            gradients = c(1, 3),
                                            gradient_colors = gradient_cols,
                                            list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                            add_shade = TRUE, 
                                            shade_alpha = 0.0075,
                                            shade_size = 0.001,
                                            empty_row_height = -0.15,
                                            base_size_ = 16,
                                            vect = TRUE,
                                            mod_formula = formula(" ~ age + path_all + sex + rsqa__MeanFD"),
                                            covariates = c("sex", "rsqa__MeanFD"),
                                            plt_title = "BioFINDER, pathology score derived from both cohorts",
                                            brain_names = c("Age", "AD Pathology (all)"),
                                            plt_subtitle = TRUE,
                                            rectangle = TRUE)

adni_path_all <- plot_gradient_relationships(adni_df_patscore %>% filter(fmri_bl), 
                                            gradient_data = grad_df |> filter(study== "adni"), 
                                            gradients = c(1, 3),
                                            gradient_colors = gradient_cols,
                                            list_of_parcel_data = list(nodal_affinity = fc_measures_adni$affinity),
                                            add_shade = TRUE, 
                                            shade_alpha = 0.0075,
                                            shade_size = 0.001,
                                            empty_row_height = -0.15,
                                            base_size_ = 16,
                                            vect = TRUE,
                                            id_var = "file_func",
                                            mod_formula = formula(" ~ age + path_all + sex + rsqa__MeanFD"),
                                            covariates = c("sex", "rsqa__MeanFD"),
                                            plt_title = "ADNI, pathology score derived from both cohorts",
                                            brain_names = c("Age", "AD Pathology (all)"),
                                            plt_subtitle = TRUE,
                                            rectangle = TRUE)


switched_grads_all_path <- ggdraw() +
  draw_plot(bf_adni_grad$plot + plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=20)),
                                                               plot.subtitle = element_text(margin = margin(l=20)))),
            x = 0, y = 0.505, width = 0.495, height = 0.495) +
  draw_plot(adni_bf_grad$plot + plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=20)),
                                                           plot.subtitle = element_text(margin = margin(l=20))))
            , x = 0.505, y = 0.505, width = 0.495, height = 0.495) +
  draw_plot(bf_path_all$plot+ plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=20)),
                                                     plot.subtitle = element_text(margin = margin(l=20))))
            , x = 0, y = 0, width = 0.495, height = 0.495) +
  draw_plot(adni_path_all$plot+ plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=20)),
                                                           plot.subtitle = element_text(margin = margin(l=20))))
            , x = 0.505, y = 0, width = 0.495, height = 0.495) +
  draw_plot_label("A") +
  draw_plot_label("B", x = 0.505) +
  draw_plot_label("C", x = 0.0, y = 0.495) +
  draw_plot_label("D", x = 0.505, y = 0.495) 

scaling_factor <- 2.75
p_name <- "switched_grad_all_path.png"
ggsave(file.path(figure_path, p_name), switched_grads_all_path,
       width = img_width*scaling_factor, height = img_width*0.85*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)


###########################
# Age and pathology collinearity
###########################

pal <- ggsci::pal_futurama()(3)

biofinder_df |> 
  filter(fmri_bl, !is.na(abnorm_ab)) |> 
  mutate(`Aβ+` = as.logical(abnorm_ab),
         Diagnosis = case_when(
           diagnosis %in% c("Normal", "SCD") ~ "CU",
           TRUE ~ diagnosis
         ) |> factor(levels = c("CU", "MCI", "AD"))) |> 
  ggplot(aes(age, pathology_ad, color = Diagnosis#`Aβ+`
  )) +
  geom_point(alpha=0.2, size = 2) +
  labs(x = "Age", y = "Pathology Score") +
  scale_color_manual(
    values = c(pal[3], pal[1], pal[2]),
    labels = c("CU", "MCI (Aβ+)", "AD (Aβ+)")
  )+
  theme_bw(base_size = 12) +
  theme(legend.position = "top") -> age_path

p_name <- "age_path_dist.png"
ggsave(file.path(figure_path, p_name), age_path,
       width = 6, height = 4, units = "in", dpi = 300, device = "png", bg = "white")

###########################
# Minor comments
###########################

long_bf_ |> filter(!fmri_bl) |> 
  mutate(yearly_path = pathΔ/as.numeric(time)) |> 
  select(path_bl, yearly_path, pathΔ) |> 
  drop_na() |> 
  ggplot(aes(path_bl, pathΔ)) +
  geom_point() +
  stat_poly_eq(use_label("R2"), size = 12) +
  stat_poly_line()


###########################
# Stuff
###########################







source("src/util.R")
source("src/util_vis.R")

print("Running non-linear models")

gam_preds_affinity <- gam_pred_nodes(biofinder_df |> filter(fmri_bl),
                                     fc_matrix = fc_measures_bf$affinity,
                                     roi_names = rois, 
                                     id_var = "image_file",
                                     print_ticker = TRUE,
                                     model_formula = formula("FC ~ s(age) + s(tau_pathology) + sex + rsqa__MeanFD"))

nonlin_p <- plot_gams_v1(gam_predictions = gam_preds_affinity, 
                         grad_df = grad_df |> filter(study == "biofinder"),
                         biofinder_data = biofinder_df |> filter(fmri_bl), scale_fac = 3)

p_name <- "figure2.png"
img_width <- 90/25.4
scale_factor = 3
magick_geom_scaling <- paste0(100*0.66, "%x", 100*0.66, "%")
ggsave(file.path(figure_path, p_name), nonlin_p, 
       width = img_width*scale_factor, 
       height = img_width*scale_factor*1.3, 
       units = "in", dpi = 300,
       bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)





