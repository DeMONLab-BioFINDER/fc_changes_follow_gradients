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
                                         subtit_lookup = c(abnorm_ab = "Αβ_pos", rsqa__MeanFD = "motion"),
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
                                          subtit_lookup = c(abnorm_ab = "Αβ_pos", rsqa__MeanFD = "motion", has_diagnosis = "MCI_or_AD"),
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
                                    mod_formula = formula(paste0("~ age + pathology_ad + apoe4 + sex + rsqa__MeanFD")),
                                    covariates = c("sex", "rsqa__MeanFD"),
                                    plt_title = "Full sample",
                                    brain_names = c("Age", "Pathology AD", "ε4 carrier"),
                                    plt_subtitle = TRUE,
                                    subtit_lookup = c(apoe4 = "ε4 carrier", rsqa__MeanFD = "motion", pathology_ad = "AD pathology"),
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
                                    mod_formula = formula(paste0("~ age + tau_pathology + apoe4 + sex + rsqa__MeanFD")),
                                    covariates = c("sex", "rsqa__MeanFD"),
                                    plt_title = "Cognitively Unimpaired (Aβ-)",
                                    brain_names = c("Age", "Tau Pathology", "ε4 carrier"),
                                    plt_subtitle = TRUE,
                                    subtit_lookup = c(apoe4 = "ε4 carrier", rsqa__MeanFD = "motion"),
                                    rectangle = TRUE)


# apoe_int <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, 
#                                                                 #diagnosis %in% c("Normal", "SCD"), 
#                                                                 #abnorm_ab==0
#                                                                 ) |> 
#                                           mutate(has_diag = factor(diagnosis %in% c("MCI", "AD"))), 
#                                     gradient_data = grad_df |> filter(study== "biofinder"), 
#                                     gradients = c(1, 3),
#                                     gradient_colors = gradient_cols,
#                                     list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
#                                     #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
#                                     add_shade = TRUE, 
#                                     shade_alpha = 0.0075,
#                                     shade_size = 0.001,
#                                     empty_row_height = -0.15,
#                                     base_size_ = 16,
#                                     vect = TRUE,
#                                     mod_formula = formula(paste0("~ age + tau_pathology + apoe4*has_diag  + sex + rsqa__MeanFD")),
#                                     covariates = c("sex", "rsqa__MeanFD"),
#                                     right_term_side = FALSE,
#                                     #plt_title = "Cognitively Unimpaired (Aβ-)",
#                                     #brain_names = c("Age", "Tau Pathology", "ε4 carrier", "Tau×ε4"),
#                                     plt_subtitle = TRUE,
#                                     rectangle = TRUE)
# apoe_int$plot
# biofinder_df %>% filter(fmri_bl, abnorm_ab == 0, tau_pathology > 0.15) |> 
#   select(apoe4) |> table()

# text_label <- "APOE ε4\nFALSE  TRUE\nn=219  n=135"
# apoe_dist <- biofinder_df %>% filter(fmri_bl, abnorm_ab == 0, tau_pathology > 0.2) |> 
#   select(apoe4, age) |> 
#   ggplot(aes(age, fill = apoe4)) +
#   geom_histogram(color = "black", position = "dodge") +
#   annotate(
#     "text",
#     x = 20, y = 20,                # adjust position as needed
#     label = text_label,
#     hjust = 0, vjust = 1,
#     family = "mono",
#     size = 7
#   ) +
#   ggtitle("Αβ-, Tau Pathology < 0.15") +
#   labs(y = "Count", x = "Age") +
#   theme(
#     legend.position = "bottom",
#     plot.background = element_rect(color = "black", fill = NA))


# apoe_p <- ggdraw() +
#   draw_plot(apoe$plot, x = 0, y = 0.5, width = 0.495, height = 0.5) +
#   draw_plot(apoe_abneg$plot, x = 0.505, y = 0.5, width = 0.495, height = 0.5) +
#   draw_plot(apoe_int$plot, x = 0, y = 0.0, width = 0.7, height = 0.49) +
#   draw_plot(apoe_dist, x = 0.71, y = 0.0, width = 0.29, height = 0.49)
# 
# scaling_factor <- 3
# p_name <- "apoe_plot.png"
# ggsave(file.path(figure_path, p_name), apoe_p,
#        width = img_width*scaling_factor, height = img_width*0.8*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
# img <- magick::image_read(file.path(figure_path, p_name))
# img_resized <- magick::image_resize(img, magick_geom_scaling)
# magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)



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
         wmhs_tot = samseg_wmhs_WMH_total_mm3,
         ct_tot = (ct_lateraltemporal_surfwt_lr + ct_medialtemporal_surfwt_lr + 
                     ct_lateralparietal_surfwt_lr + ct_medialparietal_surfwt_lr 
                   + ct_frontal_surfwt_lr + ct_occipital_surfwt_lr)/6) |> 
  filter(mri_bl)


biofinder_new <- biofinder_df |> filter(fmri_bl) |> 
  inner_join(x, by = c("sid"))

cerebro_vascular <- plot_gradient_relationships(biofinder_new %>% filter(fmri_bl) |>
                                                  mutate(
                                                    apoe4 = factor(apoe4),
                                                    has_diag = diagnosis %in% c("MCI", "AD"),
                                                    CT_adsig = -ct_adsign_lr,
                                                    ct_tot = -ct_tot,
                                                    a_syn = factor(a_syn)),
                                                gradient_data = grad_df |> filter(study== "biofinder"),
                                                gradients = c(1, 3),
                                                gradient_colors = gradient_cols,
                                                list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                                #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
                                                add_shade = TRUE,
                                                shade_alpha = 0.0075,
                                                shade_size = 0.001,
                                                point_size = 1,
                                                empty_row_height = -0.15,
                                                brain_title_vjust = -3,
                                                brain_subtitle_vjust = -5,
                                                l_width = 0.025,
                                                base_size_ = 16,
                                                vect = TRUE,
                                                mod_formula = formula(paste0(" ~ age + pathology_ad + ct_tot + wmhs_tot + lacunes + microbleeds + a_syn + sex + rsqa__MeanFD")),
                                                covariates = c("sex", "rsqa__MeanFD"),
                                                plt_title = "Effects of age and AD pathology on FCS are not explained by other age-related factors",
                                                brain_names = c("Age", "AD Pathology", "-CT total", "WMHS tot", "Lacunes+", "Microbleeds+", "αSyn+"),
                                                #brain_title_vjust = -5,
                                                #brain_subtitle_vjust = -7.5,
                                                plt_subtitle = TRUE,
                                                subtit_lookup = c(pathology_ad = "AD pathology", ct_tot = "-ct_tot", a_syn = "a_syn_pos", rsqa__MeanFD = "motion"),
                                                group_n_title = TRUE,
                                                rectangle = TRUE)

# biofinder_new |> 
#   ggplot(aes(wmhs_tot, ct_adsign_lr)) +
#   geom_point()

# ct_p <- ggdraw() +
#   draw_plot(ct_adsign$plot, x = 0, y = 0.5, width = 0.75, height = 0.49) +
#   draw_plot(ct_full$plot, x = 0, y = 0.0, width = 0.99, height = 0.49) 

img_width <-  180 / 25.4
scaling_factor <-  3
magick_geom_scaling <- paste0(100/scaling_factor, "%x", 100/scaling_factor, "%")

p_name <- "atrophy_and_cerebrovasc.png"
ggsave(file.path(figure_path, p_name), cerebro_vascular$plot,
       width = img_width*scaling_factor, height = img_width*0.45*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)


###########################
# Cognition
###########################

###############
# Linear stuff
###############
# new_cog <- read_delim("stuff_for_revisions/jonathan__20250926_104238.csv")
# 
# cogs <- new_cog |> select(aqt_form:`letter_s_45-60`)

#cog <- read_delim("stuff_for_revisions/cognitive_tests-2025-08-29.csv") 
#cog <- read_delim("stuff_for_revisions/cognitive_tests-2025-08-29.csv") 
nic <- read_delim("stuff_for_revisions/nicola__20240326_143046.csv") 
# nic <- read_delim("stuff_for_revisions/jonathan__20250926_104238.csv")
# cog <- nic |> select(1:5, aqt_form:`letter_s_45-60`,
#                      # adas_delayed_word_recall, adas_immediate_word_recall_average,
#                      # letter_s, trailmaking_a, trailmaking_b, symbol_digit,
#                      # animal_fluency, vosp_cube, bnt_15_2,
#                      cognitive_test_date
# )

cog <- nic |> select(1:5, mPACC_v1:fcsrt_immediate,
                     # adas_delayed_word_recall, adas_immediate_word_recall_average,
                     # letter_s, trailmaking_a, trailmaking_b, symbol_digit,
                     # animal_fluency, vosp_cube, bnt_15_2,
                     cognitive_test_date
)

cog_bl <- cog |> 
  filter(!is.na(cognitive_test_date)) |> 
  filter(Visit==0) |> 
  select(
    sid,
    cognitive_test_date,
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
  filter(cognitive_test_date == min(cognitive_test_date), .by = "sid") 
  # drop_na(adas_delayed_word_recall,
  #         adas_immediate_word_recall_average,
  #         letter_s, trailmaking_a, bnt_15_2, trailmaking_b,
  #         symbol_digit, animal_fluency, vosp_cube
  #         ) 

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
  ) |> drop_na() |> 
  column_to_rownames("sid")


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


cog_diff_plot <- function(b_size = 16, subt_size = rel(0.75)) {
  fit <- biofinder_cog_comp |> 
    mutate(
           memory = resid(lm(memory ~ global_cog)),
           executive = resid(lm(executive ~ global_cog))
    ) |> 
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
    theme_bw(base_size = b_size) 
  
  
  
  fit_exec <- biofinder_cog_comp |> 
    mutate(
      memory = resid(lm(memory ~ global_cog)),
      executive = resid(lm(executive ~ global_cog))
    ) |> 
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
    theme_bw(base_size = b_size) 
  
  n = length(fit_exec$residuals)
  
  library(patchwork)
  cog_eff_plot <- wrap_plots(list(mem_plot, exec_plot), ncol = 1, guides = "collect", axes = "collect")  +
    plot_annotation(
      title = paste0("Effects on memory/executive function (n = ", n, ")"),
      subtitle = "Memory ~ scale(age) + scale(AD pathology) + sex \nExecutive ~ scale(age) + scale(AD pathology) + sex",
      theme = theme(plot.background = element_rect(color = "black", fill = NA))) &
    theme(axis.text.x = element_text(size = rel(1.5)),
          axis.title.y = element_text(size = rel(1.15)),
          plot.subtitle = element_text(family = "mono", size = subt_size))
  cog_eff_plot
}

cog_diff <- cog_diff_plot(subt_size = rel(0.9))


# biofinder_cog_comp |>
#   ggplot(aes(memory, fill = diagnosis)) +
#   geom_density(alpha = 0.3)
# 
# biofinder_cog_comp |> drop_na(age) |>
#   ggplot(aes(memory, fill = cut(age, 4))) +
#   geom_density(alpha = 0.3)
# 
# biofinder_cog_comp |>
#   ggplot(aes(global_cog, memory, color = diagnosis)) +
#   geom_point(alpha = 0.8)
# 
# biofinder_cog_comp |>
#   ggplot(aes(pathology_ad, memory)) +
#   geom_point(alpha = 0.3) +
#   geom_smooth(method = "lm")
# 
# biofinder_cog_comp |> drop_na(age) |>
#   ggplot(aes(executive, fill = cut(age, 4))) +
#   geom_density(alpha = 0.3)
# 
# biofinder_cog_comp |>
#   ggplot(aes(age, executive)) +
#   geom_point(alpha = 0.3) +
#   geom_smooth(method = "lm")
# 
# biofinder_cog_comp |>
#   ggplot(aes(global_cog, executive, color = diagnosis)) +
#   geom_point(alpha = 0.8)


cognitive_composites <- plot_gradient_relationships(biofinder_cog_comp |> 
                                                      filter(fmri_bl, #diagnosis %in% c("Normal", "SCD")
                                                      ), 
                                                    gradient_data = grad_df |> filter(study== "biofinder"), 
                                                    gradients = c(1, 3),
                                                    gradient_colors = gradient_cols,
                                                    list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                                    #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
                                                    add_shade = TRUE, 
                                                    shade_alpha = 0.0075,
                                                    shade_size = 0.001,
                                                    l_width = 0.05,
                                                    point_size = 0.75,
                                                    empty_row_height = -0.225,
                                                    base_size_ = 16,
                                                    brain_title_vjust = -0.5,
                                                    r_spin_size = 1,
                                                    vect = TRUE,
                                                    mod_formula = formula(paste0(" ~ age + pathology_ad + global_cog + executive + memory + sex + rsqa__MeanFD")),
                                                    covariates = c("age", "pathology_ad", "sex", "rsqa__MeanFD"),
                                                    plt_title = "Higher executive or memory function show positive gradient alignment",
                                                    brain_names = c("Global Cognition (g)", "Executive", "Memory"),
                                                    plt_subtitle = TRUE,
                                                    group_n_title = FALSE,
                                                    rectangle = TRUE)

cognitive_composites$plot <- cognitive_composites$plot + plot_annotation(theme = theme(plot.subtitle = element_text(vjust = -6)))
source("src/mass_mediation_src.R")

df <- biofinder_cog_comp |> filter(fmri_bl
) |> 
  drop_na(age, pathology_ad, sex) |> 
  mutate(has_diag = factor(diagnosis %in% c("MCI", "AD")),
         memory = resid(lm(memory ~ global_cog)),
         executive = resid(lm(executive ~ global_cog))
  ) |> 
  #left_join(cog_edit) |> 
  inner_join(fc_measures_bf$affinity|> as_tibble(rownames = NA) |> 
               rownames_to_column("image_file") |> 
               janitor::clean_names()) 

parcs <- df |> select(starts_with("x7")) |> colnames()

# pw_mediation <- plot_mediation_gradients(subject_data = df, 
#                                       treatments = c("age", "pathology_ad"), 
#                                       outcomes = c("global_cog"),
#                                       covariates = list(c("pathology_ad", "sex", "rsqa__MeanFD"), 
#                                                         c("age", "sex", "rsqa__MeanFD")),
#                                       gradient_data = grad_df |> filter(study=="biofinder"), 
#                                       gradients = c(1, 3),
#                                       parcels = parcs,
#                                       gradient_colors = gradient_cols,
#                                       plt_title = "Parcelwise Mediation Analysis",
#                                       terms_nice = c("age", "AD path"),
#                                       plot_spacing = 0.3,
#                                       empty_row_height = -0.1,
#                                       base_size_ = 16,
#                                       add_shade = TRUE,
#                                       shade_size = 0.01, 
#                                       shade_alpha = 0.01,
#                                       rectangle = TRUE)

pw_mediation <- plot_mediation_gradients(subject_data = df |> 
                                           mutate(age_ = age,
                                                  pathology_ad_ = pathology_ad), 
                                         treatments = c("age", "age_", "pathology_ad", "pathology_ad_"), 
                                         outcomes = c("executive", "memory", "executive", "memory"),
                                         covariates = list(c("pathology_ad", "sex", "rsqa__MeanFD"), 
                                                           c("pathology_ad", "sex", "rsqa__MeanFD"),
                                                           c("age", "sex", "rsqa__MeanFD"), 
                                                           c("age", "sex", "rsqa__MeanFD")),
                                         gradient_data = grad_df |> filter(study=="biofinder"), 
                                         gradients = c(1, 3),
                                         parcels = parcs,
                                         gradient_colors = gradient_cols,
                                         plt_title = waiver(),
                                         terms_nice = c("Age", "Age", "ADpath", "ADpath"),
                                         point_size = 0.75,
                                         plot_spacing = 0.3,
                                         empty_row_height = -0.15,
                                         padding_med_str = 0.5,
                                         l_width = 0.05,
                                         med_lab_size = 5,
                                         base_size_ = 16,
                                         n_lab_size = 11,
                                         add_shade = TRUE,
                                         shade_alpha = 0.0075,
                                         shade_size = 0.001,
                                         rectangle = TRUE)

library(psych)
library(ggraph)
library(tidygraph)

sl <- cog_fac$schmid$sl

sl <- sl[, 1:5]
#sl[abs(sl)<0.2] <- NA

item_names <- c(adas_delayed_word_recall = "ADAS Delayed Recall",
                adas_immediate_word_recall_average = "ADAS Immediate Recall",
                bnt_15_2 = "Boston Naming",
                vosp_cube = "VOSP Cube",
                animal_fluency = "Animal Fluency",
                letter_s = "Letter Fluency",
                trailmaking_a = "Trail Making A",
                trailmaking_b = "Trail Making B",
                symbol_digit = "Symbol Digit",
                aqt_color_form = "AQT Color–Form")

sl <- sl[enframe(item_names)$name, ]
rownames(sl) <- item_names

latent_factors <- colnames(renamed_omega) #c("g", "Memory", "Visuosp", "Language", "Executive")

colnames(sl) <- latent_factors

df <- as_tibble(sl, rownames = NA) |> rownames_to_column("item") 

edges <- df %>%
  pivot_longer(cols = -item, names_to = "factor", values_to = "loading") %>%
  filter(!is.na(loading), loading > 0.2) %>%
  mutate(loading = round(loading, 1)) |> 
  rename(from = factor, to = item)

nodes <- tibble(
  name = c(latent_factors, df$item),
  type = c(
    rep("latent", length(latent_factors)),
    rep("item", nrow(df))
  ),
  x = c(
    0, rep(1, 4),   # left side (latent variables)
    rep(0.5, nrow(df))                  # middle (items)
  ),
  y = c(0.5,
        seq(0.8, 0.2, length.out = 4),  # spread latent factors vertically
        seq(1, 0, length.out = nrow(df))                 # spread items vertically
  )
)

graph <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)


make_push_vector <- function(edges, push_df, default = 0, unit = "mm") {
  merged <- edges %>%
    left_join(push_df, by = c("to", "from")) %>%
    mutate(push = ifelse(is.na(push), default, push))
  unit_list <- lapply(merged$push, function(x) unit(x, unit))
  push_vec <- do.call(grid::unit.c, unit_list)
  return(push_vec)
}

push_vec <- make_push_vector(edges = edges, push_df = tibble(
  to = c("VOSP Cube", "ADAS Delayed Recall", "ADAS Immediate Recall"),
  from = c("Executive", "Memory", "Memory"),
  push = c(15, 5, 5)
))

omega_plot <- ggraph(graph, layout = "manual", x = x, y = y) +
  geom_edge_link(
    aes(label = loading, 
        #start_cap = circle(radius = 1.075),#label_rect(node1.name),
        end_cap = label_rect(node2.name, padding = margin(3, 3, 3, 3, "mm"))),
    angle_calc = 'along',
    arrow = arrow(angle = 15, length = unit(3, "mm"), type = "closed"),
    label_dodge = unit(2, 'mm'),
    #check_overlap = TRUE,
    label_push = push_vec,
    color = "grey30"
  ) +
  geom_node_label(aes(label = name, filter = type == "item"), 
                  label.padding = unit(2, "mm"),
                  label.r = unit(3, "mm"), size = 5) +
  geom_node_circle(aes(filter = type == "latent", r = 0.095), fill = "white") +
  geom_node_text(aes(label = name, filter = type == "latent"), size = 5) +
  theme_void(base_size = 16) +
  ggtitle(" Bifactor model of cognition (Schmid–Leiman)") +
  theme(plot.background = element_rect(color = "black", linewidth = 1))


cognition_full_new <- 
  ggdraw() +
  draw_plot(omega_plot, x = 0, y = 0.505, width = 0.40, height = 0.495) +
  draw_plot(cognitive_composites$plot + plot_annotation(theme = theme(plot.title = element_text(margin = margin(b = -8)))),
            x = 0.41, y = 0.505, width = 0.59, height = 0.495) +
  draw_plot(cog_diff #+ plot_annotation(theme = theme(plot.margin = margin(10, 5, 0, 5)))
            , x = 0.0, y = 0.0, width = 0.35, height = 0.495) +
  draw_plot(pw_mediation, x = 0.36, y = 0, width = 0.64, height = 0.495) + 
  draw_plot_label("A", x = 0.382, y = 1, size = 18, color = "#4e4e4e") +
  draw_plot_label("B", x = 0.98, y = 1, size = 18, color = "#4e4e4e") +
  draw_plot_label("C", x = 0.33, y = 0.495, size = 18, color = "#4e4e4e") +
  draw_plot_label("D", x = 0.98, y = 0.495, size = 18, color = "#4e4e4e") +
  draw_text("Parcelwise mediation", x = 0.365, y = 0.495, hjust = 0, vjust = 1.5, size = 20) +
  draw_text("Covariates: \nAge/AD pathology \nSex \nMotion", x = 0.3675, y = 0.485, hjust = 0, vjust = 1.5, size = 12)

scaling_factor <- 3
magick_geom_scaling <- paste0(100/scaling_factor, "%x", 100/scaling_factor, "%")
p_name <- "cognition_full.png"
ggsave(file.path(figure_path, p_name), cognition_full_new,
       width = img_width*scaling_factor, height = img_width*0.8*scaling_factor, 
       units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)


###########################
# New Subtypes analysis
###########################

# wo_res <- apply(df_cog, 1, function(x) which(x == min(x)))
# w_res <- apply(df_cog, 1, function(x) which(x == min(x)))
# data.frame(wo_res, w_res, log = wo_res == w_res) |> filter(!log)
# sum(wo_res == 1)
# sum(w_res == 1)

domains <- c("memory", "executive", "language", "visuospatial")

ctrl_stats <- biofinder_cog_comp |>
  filter(diagnosis %in% c("Normal", "SCD")) |>
  summarise(
    across(
      all_of(domains),
      list(
        mean = \(x) mean(x, na.rm = TRUE),
        sd   = \(x) sd(x, na.rm = TRUE)
      )
    )
  )

df_cog <- biofinder_cog_comp |>
  mutate(
    across(
      all_of(domains),
      ~ (.x - ctrl_stats[[paste0(cur_column(), "_mean")]]) /
        ctrl_stats[[paste0(cur_column(), "_sd")]]
    )
  )

impair_cutoff <- -1 

df_cog <- df_cog |>
  mutate(
    across(
      all_of(domains),
      ~ .x < impair_cutoff,
      .names = "imp_{.col}"
    )
  ) |>
  rowwise() |>
  mutate(
    n_impairments = sum(c_across(starts_with("imp_"))),
    cog_subtype = case_when(
      imp_memory ~ "amnestic",
      !imp_memory & n_impairments > 0 ~ "non-amnestic",
      TRUE ~ NA_character_
    )
  ) |>
  ungroup() |> 
  mutate(cog_subtype = ifelse(diagnosis %in% c("MCI", "AD"), cog_subtype, NA))

df_cog$cog_subtype |> table()

source("src/util.R")
source("src/util_vis.R")

gam_preds_wo_amnestic <- gam_pred_nodes(df_cog |> filter(is.na(cog_subtype) | cog_subtype == "non-amnestic" ) |> 
                                          mutate(pathology_ad == tau_pathology),
                                        fc_matrix = fc_measures_bf$affinity,
                                        roi_names = rois, 
                                        id_var = "image_file",
                                        print_ticker = TRUE,
                                        model_formula = formula("FC ~ s(age) + s(pathology_ad, k = 5) + sex + rsqa__MeanFD"))

point_effects <- do.call(rbind, gam_preds_wo_amnestic$pat_derivs)
point_effects_pat <- point_effects[nrow(point_effects), ]
point_effects <- point_effects[-nrow(point_effects), ]

g_cors <- apply(point_effects, 2, function(x) cor(x, grad_df |> filter(study == "biofinder") |> pull(gradient1)))

plot(point_effects_pat, g_cors, type = "l")


gam_preds_only_amnestic <- gam_pred_nodes(df_cog |>  filter(is.na(cog_subtype) | cog_subtype == "amnestic" )|> 
                                            mutate(pathology_ad == tau_pathology),
                                        fc_matrix = fc_measures_bf$affinity,
                                        roi_names = rois, 
                                        id_var = "image_file",
                                        print_ticker = TRUE,
                                        model_formula = formula("FC ~ s(age) + s(pathology_ad, k = 5) + sex + rsqa__MeanFD"))

point_effects_amn <- do.call(rbind, gam_preds_only_amnestic$pat_derivs)
point_effects_pat_amn <- point_effects_amn[nrow(point_effects_amn), ]
point_effects_amn <- point_effects_amn[-nrow(point_effects_amn), ]

g_cors_amn <- apply(point_effects_amn, 2, function(x) cor(x, grad_df |> filter(study == "biofinder") |> pull(gradient1)))

plot(point_effects_pat_amn, g_cors_amn, type = "l")


g3_cors <- c(apply(point_effects, 2, function(x) cor(x, grad_df |> filter(study == "biofinder") |> pull(gradient3))),
 apply(point_effects_amn, 2, function(x) cor(x, grad_df |> filter(study == "biofinder") |> pull(gradient3))))


n_nonamn <- df_cog |> filter(is.na(cog_subtype) | cog_subtype == "non-amnestic" ) |> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()
n_amn <- df_cog |>  filter(is.na(cog_subtype) | cog_subtype == "amnestic" )|> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()

cog_sub_gam <- tibble(g1 = c(g_cors_amn, g_cors), 
       #g3 = g3_cors,
       pathology = rep(point_effects_pat_amn, 2),
       sample = factor(rep(c("Amnestic" |> paste0(" (n = ", n_amn, ")"), "Non-amnestic" |> paste0(" (n = ", n_nonamn, ")")
       ), each = length(g_cors)),
       levels = c(c("Amnestic" |> paste0(" (n = ", n_amn, ")"), "Non-amnestic" |> paste0(" (n = ", n_nonamn, ")")))
       ))   |> 
  #pivot_longer(starts_with("g"), values_to = "g_corr", names_to = "gradient") |> 
  ggplot(aes(x = pathology, y = g1, color = sample)) +
  geom_line(show.legend = F, linewidth = 1) +
  ggsci::scale_color_nejm() +
  facet_wrap(~sample) +
  labs(y = "Gradient Correlation (r)", x = "AD Pathology") +
  ggtitle("Full control sample, subsets of patients")


###########################
# Subtypes
###########################

# bf_subs <- read_csv("stuff_for_revisions/Tau_subtypes_20230925.csv")
# bf_subs_ <- bf_subs |> relocate(sid) |> 
#   select(sid, MT_subtype, cognitive_status_baseline_variable)
# 
# bf_subs_df <- biofinder_df |> filter(fmri_bl) |> 
#   left_join(bf_subs_) |> 
#   mutate(tau_subtype = case_when(
#     MT_subtype == 0 ~ NA, 
#     MT_subtype == 1 & (diagnosis %in% c("MCI", "AD")) ~ "limbic", 
#     (MT_subtype %in% c(2, 3, 4)) & (diagnosis %in% c("MCI", "AD")) ~ "non-limbic",
#     TRUE ~ NA
#     ))



# gam_preds_non_limbic <- gam_pred_nodes(bf_subs_df |> filter(is.na(tau_subtype) | tau_subtype == "non-limbic" ),
#                                         fc_matrix = fc_measures_bf$affinity,
#                                         roi_names = rois, 
#                                         id_var = "image_file",
#                                         print_ticker = TRUE,
#                                         model_formula = formula("FC ~ s(age) + s(pathology_ad, k = 5) + sex + rsqa__MeanFD"))
# 
# point_effects_non_limbic <- do.call(rbind, gam_preds_non_limbic$pat_derivs)
# point_effects_pat <- point_effects_non_limbic[nrow(point_effects_non_limbic), ]
# point_effects_non_limbic <- point_effects_non_limbic[-nrow(point_effects_non_limbic), ]
# 
# g_cors_non_limbic <- apply(point_effects_non_limbic, 2, function(x) cor(x, grad_df |> filter(study == "biofinder") |> pull(gradient1)))
# 
# plot(point_effects_pat, g_cors_non_limbic, type = "l")


# gam_preds_limbic <- gam_pred_nodes(bf_subs_df |>  filter(is.na(tau_subtype) | tau_subtype == "limbic" ),
#                                    fc_matrix = fc_measures_bf$affinity,
#                                    roi_names = rois, 
#                                    id_var = "image_file",
#                                    print_ticker = TRUE,
#                                    model_formula = formula("FC ~ s(age) + s(pathology_ad, k = 5) + sex + rsqa__MeanFD"))
# 
# 
# point_effects_limbic <- do.call(rbind, gam_preds_limbic$pat_derivs)
# point_effects_pat_amn <- point_effects_limbic[nrow(point_effects_limbic), ]
# point_effects_limbic <- point_effects_limbic[-nrow(point_effects_limbic), ]
# 
# g_cors_limbic <- apply(point_effects_limbic, 2, function(x) cor(x, grad_df |> filter(study == "biofinder") |> pull(gradient1)))
# 
# plot(point_effects_pat_amn, g_cors_limbic, type = "l")
# 
# g3_cors <- c(apply(point_effects_non_limbic, 2, function(x) cor(x, grad_df |> filter(study == "biofinder") |> pull(gradient3))),
#              apply(point_effects_limbic , 2, function(x) cor(x, grad_df |> filter(study == "biofinder") |> pull(gradient3))))


df_cog |> filter(cog_subtype == "non-amnestic" ) |> 
  ggplot(aes(x = pathology_ad )) +
  geom_histogram(col = "black")


n_nonamn <- df_cog |> filter(is.na(cog_subtype) | cog_subtype == "non-amnestic" ) |> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()
n_amn <- df_cog |>  filter(is.na(cog_subtype) | cog_subtype == "amnestic" )|> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()
# n_nonlimb <- bf_subs_df |> filter(is.na(tau_subtype) | tau_subtype == "non-limbic" )|> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()
# n_limb <- bf_subs_df |>  filter(is.na(tau_subtype) | tau_subtype == "limbic" )|> drop_na(age, pathology_ad, sex, rsqa__MeanFD)|> count()

# 
# tibble(g1 = c(g_cors_amn, g_cors), 
#        #g3 = g3_cors,
#        pathology = rep(point_effects_pat_only_eoad, 2),
#        sample = factor(rep(c("Amnestic" |> paste0(" (n = ", n_amn, ")"), "Non-amnestic" |> paste0(" (n = ", n_nonamn, ")")
#                              ), each = length(g_cors)),
#                        levels = c(c("Amnestic" |> paste0(" (n = ", n_amn, ")"), "Non-amnestic" |> paste0(" (n = ", n_nonamn, ")")))
#        ))   
#   # bind_rows(
# limb_gams <-   tibble(g1 = c(g_cors_limbic, g_cors_non_limbic),
#            g3 = g3_cors,
#            pathology = rep(point_effects_pat_amn, 2),
#            sample = factor(rep(c("Limbic" |> paste0(" (n = ", n_limb, ")"), "Non-limbic"  |> paste0(" (n = ", n_nonlimb, ")")
#                                  ), each = length(g_cors)),
#                            levels = c(c("Limbic" |> paste0(" (n = ", n_limb, ")"), "Non-limbic"  |> paste0(" (n = ", n_nonlimb, ")")))
#            )) |> 
#   #   ) |> 
#   #filter(sample %in% c("MCI/AD wo 'amnestic'", "MCI/AD w only 'amnestic'")) |> 
#   #pivot_longer(starts_with("g"), values_to = "g_corr", names_to = "gradient") |> 
#   ggplot(aes(x = pathology, y = g1,# color = sample
#              )) +
#   geom_line(show.legend = F, linewidth = 1) +
#   theme_bw(base_size = 12) +
#   ggsci::scale_color_nejm() +
#   facet_wrap(~sample) +
#   #ylim(-0.1, 0.7) +
#   labs(y = "Gradient Correlation (r)", x = "AD Pathology") +
#   ggtitle("Limbic")





legend_labs <-  c("Ab42/40", "Braak12", "Braak34", "Braak56")
biof_path_plot <- biofinder_df %>% filter(fmri_bl, !is.na(age), !is.na(pathology_ad)) %>% 
  select(pathology_ad, ab_ratio, 
         starts_with("braak"), starts_with("cho")) %>% 
  mutate(across(where(is.numeric) & !pathology_ad, scale))
ab_mean <- attr(biof_path_plot$ab_ratio, c("scaled:center"))
ab_sd <- attr(biof_path_plot$ab_ratio, c("scaled:scale"))
ab_pos <- (0.08 - ab_mean)/ab_sd
biof_path_plot |> 
  pivot_longer(-pathology_ad, names_to = "scaled_pat_measures", values_to = "value") %>% 
  ggplot(aes(pathology_ad, value, color = scaled_pat_measures)) +
  geom_smooth() +
  guides(color = guide_legend(
    ncol=1, byrow=TRUE,
    title.position="top", 
    title.hjust = 0)) +
  geom_hline(yintercept = ab_pos, linetype = 3) +
  geom_text(inherit.aes = FALSE, 
            data = data.frame(x = 0.8, y = ab_pos),
            aes(x = x, y = y),
            label = "Aβ+", nudge_y = 0.15, size = 6) +
  labs(x = "Pathology score", y = "Scaled value")+
  ggsci::scale_color_nejm(name = "Pathology", labels = legend_labs) +
  theme(legend.position = "right",
        legend.title.position = "top",
        #ggside.panel.scale = 0.2,
        legend.text = element_text(size = rel(0.6)),
        legend.title = element_text(size = rel(0.8)))  

biofinder_df |> filter(fmri_bl) |> pull(ab_ratio) |> scale()
biofinder_df |> filter(fmri_bl, abnorm_ab == 0) |> pull(ab_ratio) |> hist()
  

# biofinder_cog_comp |> filter(diagnosis %in% c("MCI", "AD")) |> 
#   select(memory, executive, language, visuospatial) -> df_cog
# 
# set.seed(4)
# km <- kmeans(scale(df_cog), centers = 4, nstart = 1000)
# clust <- km$cluster
# cluster_labels <- clust
# 
# df_clustered <- biofinder_cog_comp |> 
#   filter(diagnosis %in% c("MCI", "AD")) |> 
#   mutate(cluster = factor(cluster_labels))

# domain_scores <- df_clustered |>
#   pivot_longer(c(memory, executive, language, visuospatial#, global_cog
#                  ),
#                names_to = "domain", values_to = "cog_score") |>
#   group_by(cluster, domain) |>
#   summarise(
#     mean_cog = mean(cog_score),
#     sd = sd(cog_score),
#     n = dplyr::n(),
#     se = sd / sqrt(n),
#     ci95 = 1.96 * se,
#     .groups = "drop"
#   ) |>
#   ggplot(aes(x = cluster, y = mean_cog, fill = domain)) +
#   geom_col(position = position_dodge(width = 0.9)) +
#   geom_errorbar(
#     aes(ymin = mean_cog - ci95, ymax = mean_cog + ci95),
#     position = position_dodge(width = 0.9),
#     width = 0.2, color = "black"
#   ) +
#   ggsci::scale_fill_nejm() +
#   labs(x = "Cluster", y = "Mean Domain Score", fill = "") +
#   theme_bw(base_size = 22) +
#   theme(legend.position = "top",
#         legend.text = element_text(hjust = 0.75),
#         legend.box.margin = margin(20, 20, 0, 0),
#         plot.background = element_rect(colour = "black", linewidth = 1))
# 
# df_clustered_cu <- biofinder_cog_comp |> 
#   left_join(df_clustered |> select(sid, cluster)) |> 
#   mutate(cluster = factor(ifelse(is.na(cluster), "CU", cluster), levels = c("CU", "1", "2","3", "4")))
  

# contrasts(df_clustered$cluster) <- contr.sum
# 
# cluster_cog <- plot_gradient_relationships(df_clustered, 
#                                            gradient_data = grad_df |> filter(study== "biofinder"), 
#                                            gradients = c(1, 3),
#                                            gradient_colors = gradient_cols,
#                                            list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
#                                            #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
#                                            add_shade = TRUE, 
#                                            shade_alpha = 0.0075,
#                                            shade_size = 0.001,
#                                            empty_row_height = -0.2,
#                                            base_size_ = 18,
#                                            vect = TRUE,
#                                            rec_miss_contr = TRUE,
#                                            mod_formula = formula(paste0(" ~ age + cluster + sex + rsqa__MeanFD")),
#                                            #brain_names = c("Age", "path", "Clust1 vs All", "Clust2 vs All", "Clust3 vs All", "Clust4 vs All"),
#                                            covariates = c("sex", "rsqa__MeanFD"),
#                                            plt_title = "Cognitively Impaired (MCI/AD), sum contrasts",
#                                            #brain_names = c("Executive Cog", "Memory Cog"),
#                                            plt_subtitle = TRUE,
#                                            group_n_title = TRUE,
#                                            rectangle = TRUE)
# cluster_cog$plot
# variant_plot <- 
#   ggdraw() +
#   draw_plot(domain_scores, x = 0, y = 0.5, width = 0.35, height = 0.5) +
#   draw_plot(cluster_cog$plot, x = 0.36, y = 0, width = 0.64, height = 1) +
#   draw_plot_label("A", x = 0.33, size = 22) +
#   draw_plot_label("B", x = 0.98, size = 22) 
# 
# scaling_factor <- 3
# p_name <- "clinical_variants.png"
# ggsave(file.path(figure_path, p_name), variant_plot,
#        width = img_width*scaling_factor, height = img_width*0.5*scaling_factor, 
#        units = "in", dpi = 300, device = "png", bg = "white")
# img <- magick::image_read(file.path(figure_path, p_name))
# img_resized <- magick::image_resize(img, magick_geom_scaling)
# magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)



###########################
# Interaction
###########################


int_plot <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl), 
                                        gradient_data = grad_df |> filter(study== "biofinder"), 
                                        gradients = c(1, 3),
                                        gradient_colors = gradient_cols,
                                        list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                        add_shade = TRUE, 
                                        shade_alpha = 0.0075,
                                        shade_size = 0.001,
                                        empty_row_height = -0.15,
                                        base_size_ = 12,
                                        vect = TRUE,
                                        mod_formula = formula(" ~ scale(age) * pathology_ad + sex + rsqa__MeanFD"),
                                        covariates = c("sex", "rsqa__MeanFD"),
                                        brain_names = c("Age (scaled)", "AD Pathology", "Age × ADpath"),
                                        rectangle = TRUE,
                                        plt_title = "Full sample",
                                        plt_subtitle = TRUE)


# variant_plot_agexpath <- 
#   ggdraw() +
#   draw_plot(domain_scores, x = 0, y = 0.505, width = 0.35, height = 0.495) +
#   draw_plot(int_plot$plot, x = 0.36, y = 0.505, width = 0.64, height = 0.495) +
#   draw_plot(cluster_cog$plot, x = 0.0, y = 0, width = 1.0 , height = 0.495) +
#   draw_plot_label("A", x = 0.33, size = 22) +
#   draw_plot_label("B", x = 0.98, size = 22) +
#   draw_plot_label("C", x = 0.98, y = 0.49, size = 22) 


scaling_factor <- 3
p_name <- "interaction_plot.png"
ggsave(file.path(figure_path, p_name), int_plot$plot,
       width = img_width*scaling_factor, height = img_width*0.9*scaling_factor, 
       units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)





eoad_gam <- gam_pred_nodes(biofinder_df %>% filter(fmri_bl) |> filter(!(diagnosis %in% c("MCI", "AD") & age>=65)),
                           fc_matrix = fc_measures_bf$affinity,
                           roi_names = rois, 
                           id_var = "image_file",
                           print_ticker = TRUE,
                           model_formula = formula("FC ~ s(age) + s(pathology_ad, k = 5) + sex + rsqa__MeanFD"))


eff_eoad_gam <- apply(do.call(rbind, eoad_gam$pat_derivs)[-1001, ], 2, function(x) cor(x, grad_df |> filter(study == "biofinder") |> pull(gradient1)))
plot(point_effects_pat, eff_eoad_gam, type = "l")


load_gam  <- gam_pred_nodes(biofinder_df %>% filter(fmri_bl) |> filter(!(diagnosis %in% c("MCI", "AD") & age<=65)),
                            fc_matrix = fc_measures_bf$affinity,
                            roi_names = rois, 
                            id_var = "image_file",
                            print_ticker = TRUE,
                            model_formula = formula("FC ~ s(age) + s(pathology_ad, k = 5) + sex + rsqa__MeanFD"))


eff_load_gam <- apply(do.call(rbind, load_gam$pat_derivs)[-1001, ], 2, function(x) cor(x, grad_df |> filter(study == "biofinder") |> pull(gradient1)))
plot(point_effects_pat, eff_load_gam, type = "l")

n_eoad <- biofinder_df %>% filter(fmri_bl) |> filter(!(diagnosis %in% c("MCI", "AD") & age>=65)) |> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()
n_load <- biofinder_df %>% filter(fmri_bl) |> filter(!(diagnosis %in% c("MCI", "AD") & age<=65)) |> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()

age_gams <- tibble(g1 = c(eff_load_gam , eff_eoad_gam),
       #g3 = g3_cors,
       pathology = rep(point_effects_pat_amn, 2),
       sample = factor(rep(c("MCI/AD >65" |> paste0(" (n = ", n_load, ")"), "MCI/AD <65"  |> paste0(" (n = ", n_eoad, ")")
       ), each = length(g_cors)),
       levels = c(c("MCI/AD >65" |> paste0(" (n = ", n_load, ")"), "MCI/AD <65"  |> paste0(" (n = ", n_eoad, ")")
       ))
       )) |> 
  ggplot(aes(x = pathology, y = g1, #color = sample
             )) +
  geom_line(show.legend = F, linewidth = 1) +
  ggsci::scale_color_nejm() +
  #ylim(-0.1, 0.75) +
  theme_bw(base_size = 12) +
  facet_wrap(~sample) +
  #ylim(-0.1, 0.7) +
  labs(y = "Gradient Correlation (r)", x = "AD Pathology") +
  ggtitle("Full control sample, subsets of patients") +
  theme(#plot.background = element_rect(color = "black"),
        plot.title = element_text(size = rel(1)
                                  ))

cog_sub_gam <- tibble(g1 = c(g_cors_amn, g_cors), 
                      #g3 = g3_cors,
                      pathology = rep(point_effects_pat, 2),
                      sample = factor(rep(c("Amnestic" |> paste0(" (n = ", n_amn, ")"), "Non-amnestic" |> paste0(" (n = ", n_nonamn, ")")
                      ), each = length(g_cors)),
                      levels = c(c("Amnestic" |> paste0(" (n = ", n_amn, ")"), "Non-amnestic" |> paste0(" (n = ", n_nonamn, ")")))
                      ))   |> 
  #pivot_longer(starts_with("g"), values_to = "g_corr", names_to = "gradient") |> 
  ggplot(aes(x = pathology, y = g1, #color = sample
             )) +
  geom_line(show.legend = F, linewidth = 1) +
  ggsci::scale_color_nejm() +
  #ylim(-0.1, 0.75) +
  facet_wrap(~sample) +
  theme_bw(base_size = 12) +
  labs(y = "Gradient Correlation (r)", x = "AD Pathology") +
  #ggtitle("Full control sample, subsets of patients") +
  theme(#plot.background = element_rect(color = "black"),
        plot.title = element_text(size = rel(1)
        ))


pal <- ggsci::pal_futurama()(3)

biofinder_df |> 
  filter(fmri_bl, !is.na(abnorm_ab)) |> 
  mutate(`Aβ+` = as.logical(abnorm_ab),
         Diagnosis = case_when(
           diagnosis %in% c("Normal", "SCD") ~ "CU",
           TRUE ~ diagnosis
         ) |> factor(levels = c("CU", "MCI", "AD"))) |> 
  ggplot(aes(pathology_ad, age, color = Diagnosis#`Aβ+`
  )) +
  geom_point(alpha=0.2, size = 0.75) +
  geom_smooth(data = biofinder_df |> 
                filter(fmri_bl, !is.na(abnorm_ab), pathology_ad>0.5), 
              color = "black",
              method = "lm") +
  geom_smooth(data = biofinder_df |> 
                filter(fmri_bl, !is.na(abnorm_ab), pathology_ad<0.5), 
              color = "black",
              method = "lm") +
  labs(x = "Pathology Score", y = "Age") +
  scale_color_manual(
    values = c(pal[3], pal[1], pal[2]),
    labels = c("CU", "MCI (Aβ+)", "AD (Aβ+)")
  )+
  theme_bw(base_size = 12) +
  theme(legend.position = "top") +
  theme(plot.background = element_rect(color = "black")) -> age_path

gam_plots <- 
  wrap_plots(list(age_gams, cog_sub_gam), ncol = 1, guides = "collect", axis_titles = "collect", axes = "collect") +
  plot_annotation(theme = theme(plot.background = element_rect(color = "black")))
  
# (gam_plots[[2]] + ggtitle("Without imputation") + ylim(-0.1, 0.8)) / (cog_sub_gam+ ggtitle("With imputation") + ylim(-0.1, 0.8)) / 
#   (limb_gams+ ggtitle("Tau subtypes") + ylim(-0.1, 0.8))

age_plots <- ggdraw()+ 
  draw_plot(int_plot$plot, x = 0, y = 0, width = 0.60,  height = 1) + 
  draw_plot(gam_plots, x = 0.61, y = 0.0, width = 0.39, height = 1) 

  #draw_plot(cog_sub_gam, x = 0.61, y = 0, width = 0.39, height = 0.495) 
  

scaling_factor <- 1
img_width <- 12
p_name <- "age_plot.png"
ggsave(file.path(figure_path, p_name), age_plots,
       width = img_width*scaling_factor, height = img_width*0.5*scaling_factor, 
       units = "in", dpi = 300, device = "png", bg = "white")
# img <- magick::image_read(file.path(figure_path, p_name))
# img_resized <- magick::image_resize(img, magick_geom_scaling)
# magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)





# pw <- ggdraw() +
#   draw_plot(int_plot$plot, x = 0, y = 0.5, width = 1, height = 0.49) +
#   draw_plot(int_plot2$plot, x = 0, y = 0, width = 1, height = 0.49) 
# 
# 
# img_width <-  180 / 25.4
# scaling_factor <-  3
# magick_geom_scaling <- paste0(100/scaling_factor, "%x", 100/scaling_factor, "%")
# 
# p_name <- "figure_agexpath.png"
# ggsave(file.path(figure_path, p_name), pw,
#        width = img_width*scaling_factor, height = img_width*1.2*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
# 




# path_x <- fc_measures_bf$affinity |> as_tibble(rownames=NA) |> 
#   rownames_to_column("image_file") |> 
#   pivot_longer(-image_file, names_to = "region", values_to = "affinity") |> 
#   left_join(grad_df |> filter(study == "biofinder") |> select(region, starts_with("grad"))) |> 
#   group_by(image_file) |> 
#   arrange(image_file, desc(gradient1)) |> 
#   mutate(grad_side = case_when(
#     row_number() <= 100 ~ "associative",
#     row_number() > n() - 100 ~ "sensory",
#     TRUE ~ NA_character_
#   )) |> 
#   filter(!is.na(grad_side)) |> 
#   ungroup() |> 
#   group_by(image_file, grad_side) |> 
#   summarise(affinity = mean(affinity)) |> 
#   ungroup() |> 
#   right_join(biofinder_df |> filter(fmri_bl))
# 
# path_plot <- path_x |> 
#   mutate(age_grp = cut(age, quantile(age, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)),
#          above_median_age = age > median(age, na.rm = TRUE)) |> 
#   drop_na(above_median_age) |> 
#   ggplot(aes(pathology_ad, affinity, linetype = above_median_age, color = grad_side)) +
#   geom_smooth(fill = "lightgray") +
#   geom_point(alpha = 0.1) +
#   scale_color_manual(values = c("associative" = gradient_cols[2, 1], "sensory" = gradient_cols[1, 1])) +
#   labs(x = "AD Pathology", y = "Affinity", linetype = "Old age", color = "Domain") +
#   theme(legend.position = "right")
# 
# age_x <- fc_measures_bf$affinity |> as_tibble(rownames=NA) |> 
#   rownames_to_column("image_file") |> 
#   pivot_longer(-image_file, names_to = "region", values_to = "affinity") |> 
#   left_join(grad_df |> filter(study == "biofinder") |> select(region, starts_with("grad"))) |> 
#   group_by(image_file) |> 
#   arrange(image_file, desc(gradient3)) |> 
#   mutate(grad_side = case_when(
#     row_number() <= 100 ~ "communicative",
#     row_number() > n() - 100 ~ "executive",
#     TRUE ~ NA_character_
#   )) |> 
#   filter(!is.na(grad_side)) |> 
#   ungroup() |> 
#   group_by(image_file, grad_side) |> 
#   summarise(affinity = mean(affinity)) |> 
#   ungroup() |> 
#   right_join(biofinder_df |> filter(fmri_bl))

# age_plot <- age_x |> 
#   #mutate(pathology_ad = pathology_ad + 0.0000001) |> 
#   mutate(above_median_path = pathology_ad > median(pathology_ad, na.rm = TRUE)) |> 
#   drop_na(pathology_ad) |> 
#   drop_na(above_median_path) |> 
#   ggplot(aes(age, affinity, linetype = above_median_path, color = grad_side)) +
#   geom_smooth(fill = "lightgray") +
#   geom_point(alpha = 0.1) +
#   scale_color_manual(values = c("communicative" = gradient_cols[2, 3], "executive" = gradient_cols[1, 3])) +
#   labs(x = "Age", y = "Affinity", linetype = "High path", color = "Domain") +
#   theme(legend.position = "right")

#   
# path_plot + age_plot




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



#mean_aff <- aff
mean_aff <- fc_measures_bf$affinity[biofinder_df %>% filter(fmri_bl) |> pull(image_file), ] |> colMeans()

t_path <- orig_analysis$tmaps |> filter(term == "pathology_ad") |> pull(nodal_affinity)
lm_t <- lm(t_path ~ mean_aff); t_res <- resid(lm_t)
lm_g <- lm(g1 ~ mean_aff); g_res <- resid(lm_g)
cor(t_res, g_res)


t_age <- orig_analysis$tmaps |> filter(term == "age") |> pull(nodal_affinity)
lm_ta <- lm(t_age ~ mean_aff); ta_res <- resid(lm_ta)
lm_g3 <- lm(g3 ~ mean_aff); g3_res <- resid(lm_g3)
cor(ta_res, g3_res)

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

# grob1 <- as.grob(x)
# 
# diag_vars = c("Path_tval")
# 
# cor_path <- cor_path +
#   annotation_custom(grob1,
#                     xmin = which(levels(cor_path$data$Var1) == diag_vars[1]) - 0.5,
#                     xmax = which(levels(cor_path$data$Var1) == diag_vars[1]) + 0.5,
#                     ymin = which(levels(cor_path$data$Var2) == diag_vars[1]) - 0.5,
#                     ymax = which(levels(cor_path$data$Var2) == diag_vars[1]) + 0.8
#   ) 
# 
# cow <- ggdraw() +
#   draw_plot(cor_path) +
#   draw_plot(x, x = 0.1, y = 0.8, width = )


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

original_grads <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, diagnosis %in% c("MCI", "AD")) |> 
                                          mutate(`-mPACC` = -mPACC_v1), 
                                        gradient_data = grad_df |> filter(study=="biofinder"), 
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
                                        plt_title = "Diagnosed MCI/AD (Aβ+), original gradients",
                                        brain_names = c("Age", "AD Pathology", "-mPACC"),
                                        plt_subtitle = TRUE,
                                        rectangle = TRUE)

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


##############
# Eigen values
############## 

g1 <- grad_df |> filter(study == "biofinder") |> pull(gradient1)
g2 <- grad_df |> filter(study == "biofinder") |> pull(gradient2)
g3 <- grad_df |> filter(study == "biofinder") |> pull(gradient3)

df_clustered_cu <- biofinder_cog_comp |> 
    left_join(df_clustered |> select(sid, cluster)) |>
    mutate(cluster = factor(ifelse(is.na(cluster), "CU", cluster), levels = c("CU", "1", "2","3", "4")))

eigen_g1 <- apply(fc_measures_bf$strength[biofinder_cog_comp %>% filter(fmri_bl) |> pull(image_file), ], 1 , function(x) x %*% g1)
eigen_g3 <- apply(fc_measures_bf$strength[biofinder_cog_comp %>% filter(fmri_bl) |> pull(image_file), ], 1 , function(x) x %*% g3)

df_clustered_cu  |> filter(fmri_bl) |> mutate(
  eigen_g1 = eigen_g1,
  eigen_g3 = eigen_g3
) |> 
  ggplot(aes(pathology_ad, eigen_g1)) +
  geom_point(alpha = 0.7) +
  stat_poly_line() 
  
df_clustered_cu  |> filter(fmri_bl) |> mutate(
  eigen_g1 = eigen_g1,
  eigen_g3 = eigen_g3
) |> 
  ggplot(aes(eigen_g1)) +
  geom_density()

df_clustered_cu |> filter(fmri_bl) |> mutate(
  eigen_g1 = scale(eigen_g1),
  eigen_g3 = scale(eigen_g3)
) |> 
  filter(diagnosis %in% c("Normal", "SCD")) |> 
  pivot_longer(starts_with("eigen"), names_to = "loading", values_to = "values") |> 
  mutate(abnorm_ab = tau_pathology> 0.2) |> 
  group_by(abnorm_ab, loading) |> 
  summarise(
    mean_load = mean(values),
    sd = sd(values),
    n = dplyr::n(),
    se = sd / sqrt(n),
    ci95 = 1.96 * se,
    .groups = "drop"
  ) |>
  ggplot(aes(abnorm_ab, y = mean_load, fill = loading)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(
    aes(ymin = mean_load - ci95, ymax = mean_load + ci95),
    position = position_dodge(width = 0.9),
    width = 0.2, color = "black"
  ) 

df_clustered_cu |> filter(fmri_bl) |> mutate(
  eigen_g1 = scale(eigen_g1),
  eigen_g3 = scale(eigen_g3)
) |> 
  filter(diagnosis %in% c("Normal", "SCD")) |> 
  pivot_longer(starts_with("eigen"), names_to = "loading", values_to = "values") |> 
  mutate(abnorm_ab = tau_pathology> 0.2) |> 
  ggplot(aes(abnorm_ab, y = values, fill = loading)) +
  geom_jitter(alpha = 0.3) +
  geom_violin()



biofinder_cog_comp |> filter(fmri_bl) |> 
  mutate(
  eigen_g1 = eigen_g1,
  eigen_g3 = eigen_g3
) |> 
  lm(age ~  eigen_g1 + eigen_g3 + rsqa__MeanFD, data = _) |> 
  summary()



biofinder_cog_comp |> filter(fmri_bl) |> 
  mutate(
    eigen_g1 = eigen_g1,
    eigen_g3 = eigen_g3
  ) |> 
  filter(diagnosis %in% c("Normal", "SCD")) |> 
  ggplot(aes(age, eigen_g3)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm")


biofinder_cog_comp |> filter(fmri_bl) |> 
  mutate(
    eigen_g1 = eigen_g1,
    eigen_g3 = eigen_g3
  ) |> 
  ggplot(aes(cho12, eigen_g1)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm")


# library(mediation)
# m1 <- biofinder_cog_comp |> filter(fmri_bl) |> 
#   mutate(
#     eigen_g1 = eigen_g1,
#     eigen_g3 = eigen_g3
#   ) |> 
#   lm(eigen_g3 ~ age + pathology_ad + sex + rsqa__MeanFD, data = _) 
# 
# m2 <- biofinder_cog_comp |> filter(fmri_bl) |> 
#   mutate(
#     executive = resid(lm(executive~global_cog)),
#     memory = resid(lm(memory~global_cog)),
#     eigen_g1 = eigen_g1,
#     eigen_g3 = eigen_g3
#   ) |> 
#   lm(executive ~ age + eigen_g3 + pathology_ad + sex + rsqa__MeanFD, data = _) 
# 
# med <- mediation::mediate(m1, m2, treat = "age", mediator = "eigen_g3")
# summary(med)
