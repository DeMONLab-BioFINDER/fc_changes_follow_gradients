
figure_path <- "paper/revision_figures"
img_width <-  180 
scaling_factor <-  1

###########################
# Gradient coverage
###########################

source("src/plot_gradient_coverage.R")

test <- plot_gradient_coverage(grad_df = grad_df |> filter(study == "biofinder"),
                               base_size = 7)

gradient_cov <- plot_gradient_coverage(grad_df = grad_df |> filter(study == "biofinder"),
                                       nq_filter_terms = test$nq_filter_terms,
                                       rasterize = TRUE,
                                       wordcloud_max_size = 4.9,
                                       wordcloud_grid_margin = 1,
                                       network_legend_text_scale = 1,
                                       network_legend_key_lines = 1,
                                       base_size = 7)

p_name <- "gradient_coverage.pdf"
ggsave(file.path(figure_path, p_name), 
       gradient_cov$plot, 
       width=180, 
       height=180, 
       units = "mm",
       bg="white")

gradient_cov$data$g1_wordcloud |> mutate(panel = "SA wordcloud") |> 
  bind_rows(gradient_cov$data$g3_wordcloud |> mutate(panel = "RE wordcloud")) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_wordclouds.csv")))

gradient_cov$data$network_separation |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_boxplots.csv")))

###########################
# Broken down by pathology
###########################


tau_ab_cu <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, !is.na(pathology_ad), diagnosis %in% c("Normal", "SCD")
                                                                 ) |> 
                                           mutate(diag = diagnosis %in% c("Normal", "SCD")), 
                                           gradient_data = grad_df |> filter(study== "biofinder"), 
                                           gradients = c(1, 3),
                                           gradient_colors = gradient_cols,
                                           list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                         vect = TRUE,
                                         empty_row_height = -0.15,
                                         brain_title_size = 6,
                                         axis_text_size = 5,
                                         axis_title_size = 6,
                                         scatter_title_vjust = 0,
                                         base_size_ = 7,
                                         rasterize = TRUE,
                                         ggrastr_dpi = 300,
                                         add_shade = TRUE, 
                                         shade_alpha = 0.1,
                                         shade_size = 0.1,
                                         r_spin_size = 0.60,
                                         point_size = 0.05,
                                         point_alpha = 0.2,
                                         mod_formula = formula(paste0("~ age + tau_pathology + abnorm_ab + sex + rsqa__MeanFD")),
                                         covariates = c("sex", "rsqa__MeanFD"),
                                         plt_title = "Cognitively Unimpaired",
                                         plt_title_size = rel(0.9),
                                         plt_subtitle_size = rel(0.8),
                                         plt_subtitle = TRUE,
                                         subtit_lookup = c(abnorm_ab = "Ab_pos", rsqa__MeanFD = "motion"),
                                         brain_names = c("Age", "Tau Pathology", "Ab+"),
                                         rectangle = TRUE,
                                         plot_net_legend = FALSE,
                                         net_legend_x = 0.2,
                                         net_legend_y = 0.01)





tau_ab_cu_adj <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, !is.na(pathology_ad)) |> 
                                            mutate(has_diagnosis = diagnosis %in% c("MCI", "AD")), 
                                          gradient_data = grad_df |> filter(study== "biofinder"), 
                                          gradients = c(1, 3),
                                          gradient_colors = gradient_cols,
                                          list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                          vect = TRUE,
                                          empty_row_height = -0.15,
                                          brain_title_size = 6,
                                          axis_text_size = 5,
                                          axis_title_size = 6,
                                          scatter_title_vjust = 0,
                                          base_size_ = 7,
                                          rasterize = TRUE,
                                          ggrastr_dpi = 300,
                                          add_shade = TRUE, 
                                          shade_alpha = 0.1,
                                          shade_size = 0.1,
                                          r_spin_size = 0.65,
                                          point_size = 0.05,
                                          point_alpha = 0.3,
                                          mod_formula = formula(paste0("~ age + tau_pathology + abnorm_ab + has_diagnosis + sex + rsqa__MeanFD")),
                                          covariates = c("sex", "rsqa__MeanFD"),
                                          plt_title = "Full sample",
                                          plt_subtitle = TRUE,
                                          subtit_lookup = c(abnorm_ab = "Ab_pos", rsqa__MeanFD = "motion", has_diagnosis = "MCI_or_AD"),
                                          brain_names = c("Age", "Tau Pathology", "Ab+", "MCI/AD"),
                                          rectangle = TRUE,
                                          plot_net_legend = FALSE,
                                          net_legend_x = 0.1,
                                          net_legend_y = 0.01)



##################################
# APOE
#################################


apoe <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl#, abnorm_ab == 0, pathology_ad<0.15
                                                            ) |> 
                                      mutate(has_diag = diagnosis %in% c("MCI", "AD")), 
                                    gradient_data = grad_df |> filter(study== "biofinder"), 
                                    gradients = c(1, 3),
                                    gradient_colors = gradient_cols,
                                    list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                    empty_row_height = -0.15,
                                    brain_title_size = 6,
                                    axis_text_size = 5,
                                    axis_title_size = 6,
                                    scatter_title_vjust = 0,
                                    base_size_ = 7,
                                    vect = TRUE,
                                    rasterize = TRUE,
                                    ggrastr_dpi = 300,
                                    add_shade = TRUE, 
                                    shade_alpha = 0.1,
                                    shade_size = 0.1,
                                    r_spin_size = 0.75,
                                    point_size = 0.05,
                                    point_alpha = 0.3,
                                    plot_net_legend = FALSE,
                                    net_legend_x = 0.2,
                                    net_legend_y = 0.01,
                                    mod_formula = formula(paste0("~ age + pathology_ad + apoe4 + sex + rsqa__MeanFD")),
                                    covariates = c("sex", "rsqa__MeanFD"),
                                    plt_title = "Full sample",
                                    brain_names = c("Age", "AD Pathology", "e4 carrier"),
                                    plt_subtitle = TRUE,
                                    subtit_lookup = c(apoe4 = "e4 carrier", rsqa__MeanFD = "motion", pathology_ad = "AD pathology"),
                                    rectangle = TRUE)


apoe_abneg <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, 
                                                                  diagnosis %in% c("Normal", "SCD"), 
                                                                 # abnorm_ab==0
                                                                  ), 
                                    gradient_data = grad_df |> filter(study== "biofinder"), 
                                    gradients = c(1, 3),
                                    gradient_colors = gradient_cols,
                                    list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                    empty_row_height = -0.15,
                                    brain_title_size = 6,
                                    axis_text_size = 5,
                                    axis_title_size = 6,
                                    scatter_title_vjust = 0,
                                    base_size_ = 7,
                                    vect = TRUE,
                                    rasterize = TRUE,
                                    ggrastr_dpi = 300,
                                    add_shade = TRUE, 
                                    shade_alpha = 0.1,
                                    shade_size = 0.1,
                                    r_spin_size = 0.75,
                                    point_size = 0.05,
                                    point_alpha = 0.3,
                                    plot_net_legend = FALSE,
                                    net_legend_x = 0.2,
                                    net_legend_y = 0.01,
                                    mod_formula = formula(paste0("~ age + tau_pathology + apoe4 + sex + rsqa__MeanFD")),
                                    covariates = c("sex", "rsqa__MeanFD"),
                                    plt_title = "Cognitively Unimpaired (Ab-)",
                                    brain_names = c("Age", "Tau Pathology", "e4 carrier"),
                                    plt_subtitle = TRUE,
                                    subtit_lookup = c(apoe4 = "e4 carrier", rsqa__MeanFD = "motion", pathology_ad = "AD pathology"),
                                    rectangle = TRUE)




l_marg = 10
apoe_p_and_tauab <- ggdraw() +
  draw_plot(tau_ab_cu_adj$plot + plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=l_marg)),
                                                               plot.subtitle = element_text(margin = margin(l=l_marg)))),
            x = 0, y = 0.505, width = 0.555, height = 0.495) +
  draw_plot(tau_ab_cu$plot + plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=l_marg)),
                                                           plot.subtitle = element_text(margin = margin(l=l_marg))))
            , x = 0.565, y = 0.505, width = 1 - 0.565, height = 0.495) +
  draw_plot(apoe$plot+ plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=l_marg)),
                                                     plot.subtitle = element_text(margin = margin(l=l_marg))))
            , x = 0, y = 0, width = 0.495, height = 0.495) +
  draw_plot(apoe_abneg$plot+ plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=l_marg)),
                                                           plot.subtitle = element_text(margin = margin(l=l_marg))))
            , x = 0.505, y = 0, width = 0.495, height = 0.495) +
  draw_plot_label("A", size = 9) +
  draw_plot_label("B", x = 0.570, size = 9) +
  draw_plot_label("C", x = 0.0, y = 0.495, size = 9) +
  draw_plot_label("D", x = 0.505, y = 0.495, size = 9) 


scaling_factor <- 1
p_name <- "apoe_plot_v1_w_tau_ab.pdf"
ggsave(file.path(figure_path, p_name), apoe_p_and_tauab,
       width = img_width*scaling_factor, height = img_width*0.8*scaling_factor, units = "mm", dpi = 300, device = "pdf", bg = "white")


tau_ab_cu_adj$tmaps |> mutate(panel = "A") |> 
  bind_rows(tau_ab_cu$tmaps |> mutate(panel = "B")) |> 
  bind_rows(apoe$tmaps |> mutate(panel = "C")) |> 
  bind_rows(apoe_abneg$tmaps |> mutate(panel = "D")) |> 
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient"))) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))



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
                                         empty_row_height = -0.2,
                                         brain_title_size = 6,
                                         axis_text_size = 5,
                                         axis_title_size = 6,
                                         scatter_title_vjust = 0,
                                         base_size_ = 7,
                                         vect = TRUE,
                                         rasterize = TRUE,
                                         ggrastr_dpi = 300,
                                         add_shade = TRUE, 
                                         shade_alpha = 0.1,
                                         shade_size = 0.1,
                                         r_spin_size = 0.75,
                                         point_size = 0.05,
                                         point_alpha = 0.3,
                                         plot_net_legend = TRUE,
                                         net_legend_x = 0.2,
                                         net_legend_y = 0.01,
                                         mod_formula = formula(paste0("~ age + sex + rsqa__MeanFD")),
                                         covariates = c("sex", "rsqa__MeanFD"),
                                         plt_title = "Parcelwise tau on parcelwise FCS regression",
                                         plt_subtitle = TRUE,
                                         rectangle = TRUE)


scaling_factor <- 1
p_name <- "nodal_tau.pdf"
ggsave(file.path(figure_path, p_name), nodal_tau$plot,
       width = img_width*scaling_factor, height = img_width*0.9*scaling_factor, units = "mm", dpi = 300, device = "pdf", bg = "white")


source_data <- nodal_tau$tmaps |> 
  select(-n, -model_formula) |> 
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient")))
write_csv(source_data, file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))



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
                                                empty_row_height = -0.15,
                                                brain_title_size = 6,
                                                axis_text_size = 5,
                                                axis_title_size = 6,
                                                scatter_title_vjust = 0,
                                                base_size_ = 7,
                                                rasterize = TRUE,
                                                ggrastr_dpi = 300,
                                                add_shade = TRUE, 
                                                shade_alpha = 0.1,
                                                shade_size = 0.1,
                                                r_spin_size = 0.75,
                                                point_size = 0.05,
                                                point_alpha = 0.2,
                                                plot_net_legend = TRUE,
                                                net_legend_x = 0.08,
                                                net_legend_y = 0.01,
                                                brain_title_vjust = -3,
                                                brain_subtitle_vjust = -4,
                                                l_width = 0.025,
                                                vect = TRUE,
                                                mod_formula = formula(paste0(" ~ age + pathology_ad + ct_tot + wmhs_tot + lacunes + microbleeds + a_syn + sex + rsqa__MeanFD")),
                                                covariates = c("sex", "rsqa__MeanFD"),
                                                plt_title = "Effects of age and AD pathology on FCS are not explained by other age-related factors",
                                                brain_names = c("Age", "AD Pathology", "-CT total", "WMHS tot", "Lacunes+", "Microbleeds+", "aSyn+"),
                                                #brain_title_vjust = -5,
                                                #brain_subtitle_vjust = -7.5,
                                                plt_subtitle = TRUE,
                                                subtit_lookup = c(pathology_ad = "AD pathology", ct_tot = "-ct_tot", a_syn = "a_syn_pos", rsqa__MeanFD = "motion"),
                                                group_n_title = TRUE,
                                                rectangle = TRUE)


img_width <-  180 
p_name <- "atrophy_and_cerebrovasc.pdf"
ggsave(file.path(figure_path, p_name), cerebro_vascular$plot,
       width = img_width*scaling_factor, height = img_width*0.45*scaling_factor, units = "mm", dpi = 300, device = "pdf", bg = "white")

cerebro_vascular$tmaps |>
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient"))) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))




###########################
# Cognition (imputation)
###########################

library(psych)
library(mice)



###############
# Linear stuff
###############

nic <- read_delim("stuff_for_revisions/nicola__20240326_143046.csv") 


cog <- nic |> select(1:5, mPACC_v1:fcsrt_immediate,
                     cognitive_test_date
)



cog_bl <- cog |> 
  filter(!is.na(cognitive_test_date)) |> 
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

df <- cog_bl %>%
  mutate(cognitive_test_date = as.Date(cognitive_test_date))

baseline_dates <- df %>%
  group_by(sid) %>%
  summarise(
    baseline_date = min(cognitive_test_date, na.rm = TRUE),
    .groups = "drop"
  )

df <- df %>%
  left_join(baseline_dates, by = "sid") %>%
  mutate(
    days_from_baseline = as.numeric(cognitive_test_date - baseline_date)
  )


get_baseline_within_window <- function(date_diff, values, window_days = 365) {
  
  # index of baseline row
  i0 <- which(date_diff == 0)
  
  # if baseline is observed, use it
  if (length(i0) == 1 && !is.na(values[i0])) {
    return(values[i0])
  }
  
  # otherwise look within window
  candidates <- which(
    !is.na(values) &
      abs(date_diff) <= window_days
  )
  
  if (length(candidates) == 0) {
    return(NA_real_)
  }
  
  # take closest in time
  closest <- candidates[which.min(abs(date_diff[candidates]))]
  values[closest]
}


cog_vars <- c(
  "adas_delayed_word_recall",
  "adas_immediate_word_recall_average",
  "symbol_digit",
  "trailmaking_a",
  "trailmaking_b",
  "letter_s",
  "animal_fluency",
  "vosp_cube",
  "bnt_15_2",
  "aqt_color_form"
)


baseline_df <- df %>%
  group_by(sid) %>%
  summarise(
    baseline_date = first(baseline_date),
    
    across(
      all_of(cog_vars),
      ~ get_baseline_within_window(
        date_diff = days_from_baseline,
        values    = .x,
        window_days = 365
      ),
      .names = "{.col}"
    ),
    .groups = "drop"
  )

biofinder_cog <- baseline_df |> 
  inner_join(biofinder_df |> filter(fmri_bl), by = "sid") 

imp_df <- biofinder_cog %>%
  select(
    sid,
    age,
    sex,
    pathology_ad,
    education,
    diagnosis,
    mPACC_v1,
    any_of(cog_vars)
  ) |> 
  mutate(n_missing = rowSums(is.na(across(all_of(cog_vars))))) |> 
  filter(n_missing < 5) |> 
  select(-n_missing)


id_vars   <- "sid"

library(mice)

# --- initialise -------------------------------------------------------
ini  <- mice(imp_df, maxit = 0, print = FALSE)
meth <- ini$method
pred <- ini$predictorMatrix

# --- ID: never imputed, never used -----------------------------------
meth[id_vars] <- ""
pred[id_vars, ] <- 0
pred[, id_vars] <- 0

# --- variables that should NOT be imputed (rows = 0) ------------------
no_impute_vars <- c(
  "sex",
  "diagnosis",
  "age",
  "education",
  "mPACC_v1",
  "pathology_ad"
)

pred[no_impute_vars, ] <- 0

# --- predictors for cognitive variables -------------------------------
cog_predictors <- c(
  "age",
  "sex",
  "education",
  "pathology_ad",
  "mPACC_v1"
)

pred[cog_vars, cog_predictors] <- 1
pred[cog_vars, cog_vars]       <- 1

# --- no self-prediction ----------------------------------------------
diag(pred) <- 0


imp <- mice(
  imp_df,
  m      = 20,
  maxit  = 30,
  method = meth,
  predictorMatrix = pred,
  seed   = 42,
  print  = FALSE
)
imp_list <- complete(imp, action = "all")

ave_imp <- imp_list[[1]] |> column_to_rownames("sid") |> select(all_of(cog_vars))
ave_imp[, ] <- 0
for(i in seq_along(imp_list)) {
  ave_imp <- ave_imp + imp_list[[i]] |> column_to_rownames("sid") |> select(all_of(cog_vars))
}
ave_imp <- ave_imp/30
ave_imp <- ave_imp |> rownames_to_column("sid")


#####################
# Cognition modelling
#####################

cog_mat_cu <- ave_imp |> 
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
  ) |> inner_join(biofinder_cog |> select(-all_of(cog_vars))) |> 
  drop_na(pathology_ad) |>
  filter(diagnosis %in% c("Normal"), age <70 ) |> 
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
  ) |> 
  column_to_rownames("sid")


cog_fac <- psych::omega(cog_mat_cu, nfactors = 4)


om_factors <- cog_fac[["scores"]]
refs <- cog_mat_cu |> select(adas_delayed_word_recall, vosp_cube, animal_fluency, trailmaking_a)


rename_factors_by_corr <- function(omega_factors, ref_vars) {
  global_cog <- omega_factors[, 1]
  omega_factors <- omega_factors[, -1]
  label_map <- c(
    adas_delayed_word_recall = "memory",
    trailmaking_a           = "executive",
    animal_fluency           = "language",
    vosp_cube                = "visuospatial"
  )
  cor_mat <- cor(ref_vars, omega_factors, use = "pairwise.complete.obs")
  strongest_ref <- apply(abs(cor_mat), 2, function(x) names(x)[which.max(x)])
  new_names <- label_map[strongest_ref]
  colnames(omega_factors) <- new_names
  omega_factors <- cbind(global_cog, omega_factors)
  return(omega_factors)
}



cog_mat <- ave_imp |> 
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
  ) |> inner_join(biofinder_cog |> select(-all_of(cog_vars))) |> 
  drop_na(pathology_ad) |>
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
  ) |> 
  column_to_rownames("sid")


domain_scores <- psych::factor.scores(cog_mat, f = cog_fac$schmid$sl, method = "regression") 
domain_scores <- domain_scores$scores[, 1:5]
colnames(domain_scores) <- c("global_cog", "executive", "memory", "language", "visuospatial")


biofinder_cog_comp <- biofinder_cog |> 
  inner_join(domain_scores |> as_tibble(rownames = NA) |> rownames_to_column("sid")) 


cog_diff_plot <- function(b_size = 7, subt_size = rel(0.75)) {
  fit <- biofinder_cog_comp |> 
    mutate(
      # memory = resid(lm(memory ~ global_cog)),
      # executive = resid(lm(executive ~ global_cog))
    ) |> 
    lm(memory ~ scale(age) + scale(pathology_ad) + sex , data = _) 
  
  # Extract coefficients
  coef_df <- broom::tidy(fit) |> 
    filter(term %in% c("scale(age)", "scale(pathology_ad)")) |> 
    select(term, estimate, std.error) |> 
    mutate(ci_low = estimate - 1.96*std.error,
           ci_high = estimate + 1.96*std.error,
           outcome = "memory",
           term = recode(term,
                         "scale(age)" = "Age",
                         "scale(pathology_ad)" = "Pathology"))
  
  # Plot
  mem_plot <- ggplot(coef_df, aes(x = term, y = estimate)) +
    geom_point(size = 1) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    labs(y = "Memory", x = "") +
    theme_bw(base_size = b_size) 
  
  
  
  fit_exec <- biofinder_cog_comp |> 
    mutate(
      # memory = resid(lm(memory ~ global_cog)),
      # executive = resid(lm(executive ~ global_cog))
    ) |> 
    lm(executive ~ scale(age) + scale(pathology_ad) + sex, data = _) 
  
  # Extract coefficients
  coef_df <- broom::tidy(fit_exec) |> 
    filter(term %in% c("scale(age)", "scale(pathology_ad)")) |> 
    select(term, estimate, std.error) |> 
    mutate(ci_low = estimate - 1.96*std.error,
           ci_high = estimate + 1.96*std.error,
           outcome = "executive",
           term = recode(term,
                         "scale(age)" = "Age",
                         "scale(pathology_ad)" = "Pathology"))
  
  source_data <- bind_rows(
    mem_plot$data,
    coef_df
  )
  
  # Plot
  exec_plot <- ggplot(coef_df, aes(x = term, y = estimate)) +
    geom_point(size = 1) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
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
          plot.subtitle = element_text(family = "mono", size = subt_size, hjust = 0),
          plot.title = element_text(size = rel(0.95), hjust = 0))
  list(
    plot = cog_eff_plot,
    data = source_data
  )
}

cog_diff_res <- cog_diff_plot(subt_size = rel(0.75))
cog_diff <- cog_diff_res$plot



cognitive_composites <- plot_gradient_relationships(biofinder_cog_comp |> 
                                                      filter(fmri_bl), 
                                                    gradient_data = grad_df |> filter(study== "biofinder"), 
                                                    gradients = c(1, 3),
                                                    gradient_colors = gradient_cols,
                                                    list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                                    brain_title_size = 6,
                                                    axis_text_size = 5,
                                                    axis_title_size = 6,
                                                    scatter_title_vjust = 0,
                                                    rasterize = TRUE,
                                                    ggrastr_dpi = 300,
                                                    add_shade = TRUE, 
                                                    shade_alpha = 0.1,
                                                    shade_size = 0.1,
                                                    point_size = 0.05,
                                                    point_alpha = 0.2,
                                                    l_width = 0.05,
                                                    empty_row_height = -0.225,
                                                    base_size_ = 7,
                                                    brain_title_vjust = 1,
                                                    r_spin_size = 0.75,
                                                    vect = TRUE,
                                                    mod_formula = formula(paste0(" ~ age + pathology_ad + global_cog + executive + memory + sex + rsqa__MeanFD")),
                                                    covariates = c("age", "pathology_ad", "sex", "rsqa__MeanFD"),
                                                    plt_title = "Higher executive or memory function show positive gradient alignment",
                                                    brain_names = c("Global Cognition", "Executive", "Memory"),
                                                    plt_subtitle = TRUE,
                                                    group_n_title = FALSE,
                                                    rectangle = TRUE,
                                                    plot_net_legend = TRUE,
                                                    net_legend_x = 0.15,
                                                    net_legend_y = 0.01)




source("src/make_omega_plot.R")
omega_plot_res <- plot_cognition_factor_model()
omega_plot <- omega_plot_res$plot


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

pw_mediation_bar_res <- plot_mediation_barplot(subject_data = df |> 
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
                                               plt_title = "Parcelwise mediation",
                                               terms_nice = c("Age", "Age", "ADpath", "ADpath"),
                                               point_size = 0.75,
                                               plot_spacing = 0.3,
                                               empty_row_height = -0.15,
                                               padding_med_str = 0.5,
                                               l_width = 0.05,
                                               med_lab_size = 1.5,
                                               base_size_ = 7,
                                               n_lab_size = 4.5,
                                               add_shade = TRUE,
                                               shade_alpha = 0.1,
                                               shade_size = 0.1,
                                               rasterize = TRUE,
                                               rectangle = TRUE)
pw_mediation_bar <- pw_mediation_bar_res$plot



cognition_full_new <- 
  ggdraw() +
  draw_plot(omega_plot, x = 0, y = 0.505, width = 0.40, height = 0.495) +
  draw_plot(cognitive_composites$plot, x = 0.41, y = 0.505, width = 0.59, height = 0.495) +
  draw_plot(cog_diff #+ plot_annotation(theme = theme(plot.margin = margin(10, 5, 0, 5)))
            , x = 0.0, y = 0.0, width = 0.35, height = 0.495) +
  draw_plot(pw_mediation_bar, x = 0.36, y = 0, width = 0.64, height = 0.495) +
  draw_plot_label("A", x = 0.38, y = 1, size = 8, color = "#4e4e4e") +
  draw_plot_label("B", x = 0.98, y = 1, size = 8, color = "#4e4e4e") +
  draw_plot_label("C", x = 0.33, y = 0.495, size = 8, color = "#4e4e4e") +
  draw_plot_label("D", x = 0.98, y = 0.495, size = 8, color = "#4e4e4e") +
  #draw_text("Parcelwise mediation", x = 0.365, y = 0.495, hjust = 0, vjust = 1.5, size = 20) +
  #draw_text("Covariates: \nAge/ADpath \nSex \nMotion", x = 0.3675, y = 0.485, hjust = 0, vjust = 1.5, size = 12)
  draw_text("Covariates: \nAge/ADpath \nSex \nMotion", x = 0.3675, y = 0.46, hjust = 0, vjust = 1.5, size = 5)

img_width <-  180 
figure_path <- "paper/revision_figures/imputed_cog"
p_name <- "cognition_full.svg"
ggsave(file.path(figure_path, p_name), cognition_full_new,
       width = img_width*scaling_factor, height = img_width*0.8*scaling_factor, 
       units = "mm", dpi = 300, 
       device = "svg", 
       bg = "transparent")


omega_plot_res$data |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_panelA.csv")))

cognitive_composites$tmaps |> 
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient"))) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_panelB.csv")))

cog_diff_res$data |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_panelC.csv")))

pw_mediation_bar_res$data |> select(-c(method, affinity, sim_method, threshold)) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_panelD.csv")))



###############
# Interaction 
###############

int_plot <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl), 
                                        gradient_data = grad_df |> filter(study== "biofinder"), 
                                        gradients = c(1, 3),
                                        gradient_colors = gradient_cols,
                                        list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                        brain_title_size = 6,
                                        axis_text_size = 5,
                                        axis_title_size = 6,
                                        scatter_title_vjust = 0,
                                        rasterize = TRUE,
                                        ggrastr_dpi = 300,
                                        add_shade = TRUE, 
                                        shade_alpha = 0.1,
                                        shade_size = 0.1,
                                        point_size = 0.05,
                                        point_alpha = 0.2,
                                        l_width = 0.05,
                                        empty_row_height = -0.15,
                                        base_size_ = 7,
                                        brain_title_vjust = -3,
                                        r_spin_size = 0.775,
                                        vect = TRUE,
                                        plt_subtitle = TRUE,
                                        group_n_title = FALSE,
                                        rectangle = TRUE,
                                        plot_net_legend = TRUE,
                                        net_legend_x = 0.15,
                                        net_legend_y = 0.01,
                                        mod_formula = formula(" ~ scale(age) * pathology_ad + sex + rsqa__MeanFD"),
                                        covariates = c("sex", "rsqa__MeanFD"),
                                        brain_names = c("Age (scaled)", "AD Pathology", "Age × ADpath"),
                                        plt_title = "Full sample")



####################
# Cognitive subtypes
####################



domains <- c("executive", "memory", "language", "visuospatial")
biofinder_cog_comp_res <- biofinder_cog_comp


ctrl_stats <- biofinder_cog_comp_res |>
  filter(diagnosis %in% c("Normal")) |>
  summarise(
    across(
      all_of(domains),
      list(
        mean = \(x) mean(x, na.rm = TRUE),
        sd   = \(x) sd(x, na.rm = TRUE)
      )
    )
  )

df_cog <- biofinder_cog_comp_res |>
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
    min_domain = which.min(c_across(all_of(domains))),
    cog_subtype = case_when(
      imp_memory ~ "amnestic",
      !imp_memory & n_impairments > 0 ~ "non-amnestic",
      TRUE ~ NA_character_
    )
  ) |>
  ungroup() |>
  mutate(cog_subtype = ifelse(diagnosis %in% c("MCI", "AD"), cog_subtype, NA),
         min_domain = ifelse(diagnosis %in% c("MCI", "AD"), min_domain, NA))


source("src/util.R")


plot_gam_gradient_effects <- function(
    gam_nonamnestic,
    gam_amnestic,
    gam_load,
    gam_eoad,
    grad_df,
    study = "biofinder",
    gradient = "gradient1",
    n_nonamnestic = NULL,
    n_amnestic = NULL,
    n_load = NULL,
    n_eoad = NULL
) {
  
  extract_gam_cors <- function(gam_obj, sample_name, n = NULL) {
    
    point_effects <- do.call(rbind, gam_obj$pat_derivs)
    
    pathology <- point_effects[nrow(point_effects), ]
    point_effects <- point_effects[-nrow(point_effects), , drop = FALSE]
    
    grad_vec <- grad_df |>
      dplyr::filter(study == !!study) |>
      dplyr::pull(.data[[gradient]])
    
    cors <- apply(point_effects, 2, function(x) cor(x, grad_vec, use = "pairwise.complete.obs"))
    
    label <- if (!is.null(n)) {
      paste0(sample_name, " (n = ", n, ")")
    } else {
      sample_name
    }
    
    tibble::tibble(
      pathology = as.numeric(pathology),
      gradient_correlation = as.numeric(cors),
      sample = label,
      sample_raw = sample_name,
      n = n,
      gradient = gradient
    )
  }
  
  source_data_age <- dplyr::bind_rows(
    extract_gam_cors(gam_load, "MCI/AD >65", n_load),
    extract_gam_cors(gam_eoad, "MCI/AD <65", n_eoad)
  ) |>
    dplyr::mutate(
      panel = "Full control sample, subsets of patients",
      sample = factor(
        sample,
        levels = unique(sample)
      )
    )
  
  source_data_cog <- dplyr::bind_rows(
    extract_gam_cors(gam_amnestic, "Amnestic", n_amnestic),
    extract_gam_cors(gam_nonamnestic, "Non-amnestic", n_nonamnestic)
  ) |>
    dplyr::mutate(
      panel = "Cognitive subtype",
      sample = factor(
        sample,
        levels = unique(sample)
      )
    )
  
  source_data <- dplyr::bind_rows(
    source_data_age,
    source_data_cog
  )
  
  age_plot <- source_data_age |>
    ggplot2::ggplot(ggplot2::aes(x = pathology, y = gradient_correlation)) +
    ggplot2::geom_line(show.legend = FALSE, linewidth = 0.5) +
    ggplot2::facet_wrap(~sample) +
    ggplot2::theme_bw(base_size = 7) +
    ggplot2::labs(
      y = "Gradient Correlation (r)",
      x = "AD Pathology"
    ) +
    ggplot2::ggtitle("Full control sample, subsets of patients") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = ggplot2::rel(1))
    )
  
  cog_plot <- source_data_cog |>
    ggplot2::ggplot(ggplot2::aes(x = pathology, y = gradient_correlation)) +
    ggplot2::geom_line(show.legend = FALSE, linewidth = 0.5) +
    ggplot2::facet_wrap(~sample) +
    ggplot2::theme_bw(base_size = 7) +
    ggplot2::labs(
      y = "Gradient Correlation (r)",
      x = "AD Pathology"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = ggplot2::rel(1))
    )
  
  plot <- patchwork::wrap_plots(
    list(age_plot, cog_plot),
    ncol = 1,
    guides = "collect",
    axis_titles = "collect",
    axes = "collect"
  ) +
    patchwork::plot_annotation(
      theme = ggplot2::theme(
        plot.background = ggplot2::element_rect(
          color = "black",
          linewidth = 0.5
        )
      )
    )
  
  list(
    plot = plot,
    source_data = source_data
  )
}

gam_preds_wo_amnestic <- gam_pred_nodes(df_cog |> filter(is.na(cog_subtype) | cog_subtype == "non-amnestic" ),
                                        fc_matrix = fc_measures_bf$affinity,
                                        roi_names = rois,
                                        id_var = "image_file",
                                        print_ticker = TRUE,
                                        model_formula = formula("FC ~ s(age) + s(pathology_ad, k = 5) + sex + rsqa__MeanFD"))

gam_preds_only_amnestic <- gam_pred_nodes(df_cog |>  filter(is.na(cog_subtype) | cog_subtype == "amnestic" ),
                                          fc_matrix = fc_measures_bf$affinity,
                                          roi_names = rois,
                                          id_var = "image_file",
                                          print_ticker = TRUE,
                                          model_formula = formula("FC ~ s(age) + s(pathology_ad, k = 5) + sex + rsqa__MeanFD"))

eoad_gam <- gam_pred_nodes(biofinder_df %>% filter(fmri_bl) |> filter(!(diagnosis %in% c("MCI", "AD") & age>=65)),
                           fc_matrix = fc_measures_bf$affinity,
                           roi_names = rois, 
                           id_var = "image_file",
                           print_ticker = TRUE,
                           model_formula = formula("FC ~ s(age) + s(pathology_ad, k = 5) + sex + rsqa__MeanFD"))


load_gam  <- gam_pred_nodes(biofinder_df %>% filter(fmri_bl) |> filter(!(diagnosis %in% c("MCI", "AD") & age<=65)),
                            fc_matrix = fc_measures_bf$affinity,
                            roi_names = rois, 
                            id_var = "image_file",
                            print_ticker = TRUE,
                            model_formula = formula("FC ~ s(age) + s(pathology_ad, k = 5) + sex + rsqa__MeanFD"))

# write_rds(list(gam_preds_wo_amnestic = gam_preds_wo_amnestic, 
#                gam_preds_only_amnestic = gam_preds_only_amnestic,
#                eoad_gam = eoad_gam,
#                load_gam = load_gam), "data/modelled_data/gams_supplementary.rds")


n_nonamn <- df_cog |> filter(is.na(cog_subtype) | cog_subtype == "non-amnestic" ) |> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()
n_amn <- df_cog |>  filter(is.na(cog_subtype) | cog_subtype == "amnestic" )|> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()

n_eoad <- biofinder_df %>% filter(fmri_bl) |> filter(!(diagnosis %in% c("MCI", "AD") & age>=65)) |> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()
n_load <- biofinder_df %>% filter(fmri_bl) |> filter(!(diagnosis %in% c("MCI", "AD") & age<=65)) |> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()


plot_gam_gradient_effects_res <- plot_gam_gradient_effects(
  gam_preds_wo_amnestic,
  gam_preds_only_amnestic,
  load_gam,
  eoad_gam,
  grad_df = grad_df,
  gradient = "gradient1",
  n_nonamnestic = n_nonamn$n,
  n_amnestic = n_amn$n,
  n_load = n_load$n,
  n_eoad = n_eoad$n
)

gam_plots <- plot_gam_gradient_effects_res$plot

age_plots <- ggdraw()+ 
  draw_plot(int_plot$plot, x = 0, y = 0, width = 0.60,  height = 1) + 
  draw_plot(gam_plots, x = 0.61, y = 0.0, width = 0.39, height = 1) 

scaling_factor <- 1
img_width <- 180
p_name <- "age_plot.pdf"
ggsave(file.path(figure_path, p_name), age_plots,
       width = img_width*scaling_factor, height = img_width*0.5*scaling_factor, 
       units = "mm", dpi = 300, device = "pdf", bg = "white")

plot_gam_gradient_effects_res$source_data |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_LeftPanel.csv")))

int_plot$tmaps |> 
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient"))) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_RightPanel.csv")))


##################################
# Causality comment
##################################


orig_analysis <- nodal_regression_fits(
  biofinder_df %>% filter(fmri_bl), 
  fc_measures_bf$affinity,
  roi_names = rois,
  vectorised = TRUE,
  model_formula = formula(paste0(" ~ age + pathology_ad + sex  + rsqa__MeanFD"))
)

orig_analysis <- get_nodal_ests(orig_analysis)


subj_corr_g1 <- apply(fc_measures_bf$affinity[biofinder_df %>% filter(fmri_bl) |> pull(image_file), ], 1 , function(x) cor(x, g1))
subj_corr_g3 <- apply(fc_measures_bf$affinity[biofinder_df %>% filter(fmri_bl) |> pull(image_file), ], 1 , function(x) cor(x, g3))

tibble(`SA corr` = subj_corr_g1, `RE corr` = subj_corr_g3) |> 
  pivot_longer(everything(), names_to = "Gradient", values_to = "Correlation") |> 
  ggplot(aes(x = Correlation)) +
  geom_histogram(color = "black", linewidth = 0.1) +
  theme_bw(base_size = 7) +
  labs(title = "Individual correlation with gradients" , y = "") +
  facet_wrap(~Gradient) +
  theme(plot.title = element_text(size = rel(0.85)),
        plot.background = element_blank()) -> subject_corr
  
sexed <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl), 
                                     gradient_data = grad_df |> filter(study== "biofinder"), 
                                     gradients = c(1, 3),
                                     gradient_colors = gradient_cols,
                                     list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                     brain_title_size = 6,
                                     axis_text_size = 5,
                                     axis_title_size = 6,
                                     scatter_title_vjust = 0,
                                     rasterize = TRUE,
                                     ggrastr_dpi = 300,
                                     add_shade = TRUE, 
                                     shade_alpha = 0.1,
                                     shade_size = 0.1,
                                     point_size = 0.05,
                                     point_alpha = 0.2,
                                     l_width = 0.05,
                                     empty_row_height = -0.15,
                                     base_size_ = 7,
                                     brain_title_vjust = -3,
                                     r_spin_size = 0.775,
                                     vect = TRUE,
                                     plt_subtitle = TRUE,
                                     group_n_title = FALSE,
                                     rectangle = TRUE,
                                     plot_net_legend = TRUE,
                                     net_legend_x = 0.15,
                                     net_legend_y = 0.01,
                                     mod_formula = formula(" ~ age + pathology_ad + education + sex + rsqa__MeanFD"),
                                     covariates = c("rsqa__MeanFD"),
                                     plt_title = "Full Sample",
                                     title_lmargin = 10,
                                     brain_names = c("Age", "AD Pathology", "Education", "Sex")
                                     )


mean_aff <- fc_measures_bf$affinity[biofinder_df %>% filter(fmri_bl) |> pull(image_file), ] |> colMeans()
t_path <- orig_analysis |> filter(term == "pathology_ad") |> pull(statistic)
t_age <- orig_analysis |> filter(term == "age") |> pull(statistic)

scale_fac <- 7/20

cor_path <- ppcor::pcor(data.frame(SA_axis = g1, tmap_path = t_path, Affinity = mean_aff))$estimate |> 
  as.data.frame() |> 
  rownames_to_column("Var1") |> 
  pivot_longer(-Var1, names_to = "Var2", values_to = "correlation") |> 
  mutate(
    Var1 = str_replace_all(Var1, "_", " "),
    Var2 = str_replace_all(Var2, "_", " "),
    Var1 = factor(Var1, levels = unique(Var1)),
    Var2 = factor(Var2, levels = rev(unique(Var1)))  # reverse order
  ) |> 
  ggplot(aes(x = Var1, y = Var2, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = muted("blue"), high = muted("red"), midpoint = 0,   limits = c(-1, 1)) +
  scale_x_discrete(position = "top") +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 6*scale_fac, color = "darkgray") +
  coord_fixed() +
  labs(title = "Partial correlation pathology",
       x = NULL, y = NULL, fill = "Correlation") +
  theme_minimal(base_size = 7) +
  theme(
    legend.position = "",
    axis.text = element_text(size = rel(0.65)),
    axis.text.x = element_text(angle = 45, vjust = -1, hjust = 0),
    panel.grid = element_blank(),
    plot.title = element_text(size = rel(0.85)))

cor_age <- ppcor::pcor(data.frame(RE_axis = g3, tmap_age = t_age, Affinity = mean_aff))$estimate |> 
  as.data.frame() |> 
  rownames_to_column("Var1") |> 
  pivot_longer(-Var1, names_to = "Var2", values_to = "correlation") |> 
  mutate(
    Var1 = str_replace_all(Var1, "_", " "),
    Var2 = str_replace_all(Var2, "_", " "),
    Var1 = factor(Var1, levels = unique(Var1)),
    Var2 = factor(Var2, levels = rev(unique(Var1)))  # reverse order
  ) |> 
  ggplot(aes(x = Var1, y = Var2, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = muted("blue"), high = muted("red"), midpoint = 0,   limits = c(-1, 1)) +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 6*scale_fac, color = "darkgray") +
  scale_x_discrete(position = "top") +
  coord_fixed() +
  labs(title = "Partial correlation age",
       x = NULL, y = NULL, fill = "Correlation") +
  theme_minimal(base_size = 7) +
  theme(
    legend.position = "",
    axis.text = element_text(size = rel(0.65)),
    axis.text.x = element_text(angle = 45, vjust = -1, hjust = 0),
    panel.grid = element_blank(),
    plot.title = element_text(size = rel(0.85))
  )


part_cor <- (cor_path + cor_age) / subject_corr
caus_plot <-
  ggdraw() +
  draw_plot(sexed$plot, x = 0, y = 0, width = 0.6, height = 1) +
  draw_plot_label("A", size = 7) +
  draw_plot(part_cor, x = 0.575, y = 0, width = 0.425, height = 1) +
  draw_plot_label("B", x = 0.6, size = 7) +
  draw_plot_label("C", x = 0.6, y = 0.55, size = 7) 

img_width <- 180
figure_path <- "paper/revision_figures"
scaling_factor <- 1
p_name <- "caus_plot.pdf"
ggsave(file.path(figure_path, p_name), caus_plot,
       width = img_width*scaling_factor, height = img_width*0.45*scaling_factor, 
       units = "mm", dpi = 300, device = "pdf", bg = "white")

sexed$tmaps |> 
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient"))) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_panelA.csv")))

ppcor::pcor(data.frame(SA_axis = g1, tmap_path = t_path, Affinity = mean_aff))$estimate |> 
  as.data.frame() |> 
  rownames_to_column("Var1") |> 
  pivot_longer(-Var1, names_to = "Var2", values_to = "correlation") |> 
  mutate(panel = "Partial correlation pathology") |> 
  bind_rows(ppcor::pcor(data.frame(RE_axis = g3, tmap_age = t_age, Affinity = mean_aff))$estimate |> 
              as.data.frame() |> 
              rownames_to_column("Var1") |> 
              pivot_longer(-Var1, names_to = "Var2", values_to = "correlation") |> 
              mutate(panel = "Partial correlation age")) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_panelB.csv")))
  
tibble(`SA corr` = subj_corr_g1, `RE corr` = subj_corr_g3) |> 
  pivot_longer(everything(), names_to = "Gradient", values_to = "Correlation") |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_panelC.csv")))


###########################
# Gradients from AD
###########################


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
                                        brain_title_size = 6,
                                        axis_text_size = 5,
                                        axis_title_size = 5,
                                        scatter_title_vjust = 0,
                                        rasterize = TRUE,
                                        ggrastr_dpi = 300,
                                        add_shade = TRUE, 
                                        shade_alpha = 0.1,
                                        shade_size = 0.1,
                                        point_size = 0.05,
                                        point_alpha = 0.2,
                                        l_width = 0.05,
                                        empty_row_height = -0.25,
                                        base_size_ = 7,
                                        brain_title_vjust = -1,
                                        r_spin_size = 0.775,
                                        vect = TRUE,
                                        plt_subtitle = TRUE,
                                        group_n_title = FALSE,
                                        rectangle = TRUE,
                                        plot_net_legend = TRUE,
                                        net_legend_x = 0.2,
                                        net_legend_y = 0.015,
                                        mod_formula = formula(" ~ age + pathology_ad + `-mPACC` + sex + rsqa__MeanFD"),
                                        covariates = c("sex", "rsqa__MeanFD"),
                                        plt_title = "Diagnosed MCI/AD, gradients derived from full sample ",
                                        brain_names = c("Age", "AD Pathology", "-mPACC"))

img_width = 88
p_name <- "all_grad_plot.pdf"
ggsave(file.path(figure_path, p_name), all_grad$plot,
       width = img_width*scaling_factor,
       height = img_width*0.8*scaling_factor, 
       units = "mm", dpi = 300, device = "pdf", bg = "white")

all_grad$tmaps |> 
  left_join(grad_all |> filter(study %in% c("everyone")) |> 
              select(region, study, starts_with("gradient"))) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))

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
                                        brain_title_size = 6,
                                        axis_text_size = 5,
                                        axis_title_size = 5,
                                        scatter_title_vjust = 0,
                                        rasterize = TRUE,
                                        ggrastr_dpi = 300,
                                        add_shade = TRUE, 
                                        shade_alpha = 0.1,
                                        shade_size = 0.1,
                                        point_size = 0.05,
                                        point_alpha = 0.2,
                                        l_width = 0.05,
                                        empty_row_height = -0.25,
                                        base_size_ = 7,
                                        brain_title_vjust = -1,
                                        r_spin_size = 0.75,
                                        vect = TRUE,
                                        plt_subtitle = TRUE,
                                        group_n_title = FALSE,
                                        rectangle = TRUE,
                                        mod_formula = formula(" ~ age + pathology_ad + sex + rsqa__MeanFD"),
                                        covariates = c("sex", "rsqa__MeanFD"),
                                        plt_title = "Gradients derived from ADNI in BioFINDER sample",
                                        brain_names = c("Age", "AD Pathology"))

adni_bf_grad <- plot_gradient_relationships(adni_df %>% filter(fmri_bl), 
                                            gradient_data = grad_df |> filter(study== "biofinder"), 
                                            gradients = c(1, 3),
                                            gradient_colors = gradient_cols,
                                            list_of_parcel_data = list(nodal_affinity = fc_measures_adni$affinity),
                                            brain_title_size = 6,
                                            axis_text_size = 5,
                                            axis_title_size = 5,
                                            scatter_title_vjust = 0,
                                            rasterize = TRUE,
                                            ggrastr_dpi = 300,
                                            add_shade = TRUE, 
                                            shade_alpha = 0.1,
                                            shade_size = 0.1,
                                            point_size = 0.05,
                                            point_alpha = 0.2,
                                            l_width = 0.05,
                                            empty_row_height = -0.25,
                                            base_size_ = 7,
                                            brain_title_vjust = -1,
                                            r_spin_size = 0.75,
                                            vect = TRUE,
                                            plt_subtitle = TRUE,
                                            group_n_title = FALSE,
                                            rectangle = TRUE,
                                            id_var = "file_func",
                                            mod_formula = formula(" ~ age + pathology_ad + sex + rsqa__MeanFD"),
                                            covariates = c("sex", "rsqa__MeanFD"),
                                            plt_title = "Gradients derived from BioFINDER in ADNI sample",
                                            brain_names = c("Age", "AD Pathology"))



bf_path_all <- plot_gradient_relationships(biofinder_df_patscore %>% filter(fmri_bl), 
                                           gradient_data = grad_df |> filter(study== "biofinder"), 
                                           gradients = c(1, 3),
                                           gradient_colors = gradient_cols,
                                           list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                           brain_title_size = 6,
                                           axis_text_size = 5,
                                           axis_title_size = 5,
                                           scatter_title_vjust = 0,
                                           rasterize = TRUE,
                                           ggrastr_dpi = 300,
                                           add_shade = TRUE, 
                                           shade_alpha = 0.1,
                                           shade_size = 0.1,
                                           point_size = 0.05,
                                           point_alpha = 0.2,
                                           l_width = 0.05,
                                           empty_row_height = -0.25,
                                           base_size_ = 7,
                                           brain_title_vjust = -1,
                                           r_spin_size = 0.75,
                                           vect = TRUE,
                                           plt_subtitle = TRUE,
                                           group_n_title = FALSE,
                                           rectangle = TRUE,
                                           mod_formula = formula(" ~ age + path_all + sex + rsqa__MeanFD"),
                                           covariates = c("sex", "rsqa__MeanFD"),
                                           plt_title = "BioFINDER, pathology score derived from both cohorts",
                                           brain_names = c("Age", "AD Pathology (all)"))

adni_path_all <- plot_gradient_relationships(adni_df_patscore %>% filter(fmri_bl), 
                                            gradient_data = grad_df |> filter(study== "adni"), 
                                            gradients = c(1, 3),
                                            gradient_colors = gradient_cols,
                                            list_of_parcel_data = list(nodal_affinity = fc_measures_adni$affinity),
                                            brain_title_size = 6,
                                            axis_text_size = 5,
                                            axis_title_size = 5,
                                            scatter_title_vjust = 0,
                                            rasterize = TRUE,
                                            ggrastr_dpi = 300,
                                            add_shade = TRUE, 
                                            shade_alpha = 0.1,
                                            shade_size = 0.1,
                                            point_size = 0.05,
                                            point_alpha = 0.2,
                                            l_width = 0.05,
                                            empty_row_height = -0.25,
                                            base_size_ = 7,
                                            brain_title_vjust = -1,
                                            r_spin_size = 0.75,
                                            vect = TRUE,
                                            plt_subtitle = TRUE,
                                            group_n_title = FALSE,
                                            rectangle = TRUE,
                                            id_var = "file_func",
                                            mod_formula = formula(" ~ age + path_all + sex + rsqa__MeanFD"),
                                            covariates = c("sex", "rsqa__MeanFD"),
                                            plt_title = "ADNI, pathology score derived from both cohorts",
                                            brain_names = c("Age", "AD Pathology (all)"))


switched_grads_all_path <- ggdraw() +
  draw_plot(bf_adni_grad$plot + plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=10)),
                                                               plot.subtitle = element_text(margin = margin(l=10)))),
            x = 0, y = 0.505, width = 0.495, height = 0.495) +
  draw_plot(adni_bf_grad$plot + plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=10)),
                                                           plot.subtitle = element_text(margin = margin(l=10))))
            , x = 0.505, y = 0.505, width = 0.495, height = 0.495) +
  draw_plot(bf_path_all$plot+ plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=10)),
                                                     plot.subtitle = element_text(margin = margin(l=10))))
            , x = 0, y = 0, width = 0.495, height = 0.495) +
  draw_plot(adni_path_all$plot+ plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=10)),
                                                           plot.subtitle = element_text(margin = margin(l=10))))
            , x = 0.505, y = 0, width = 0.495, height = 0.495) +
  draw_plot_label("A", size = 9) +
  draw_plot_label("B", x = 0.505, size = 9) +
  draw_plot_label("C", x = 0.0, y = 0.495, size = 9) +
  draw_plot_label("D", x = 0.505, y = 0.495, size = 9) 

scaling_factor <- 1
img_width <- 180
p_name <- "switched_grad_all_path.pdf"
ggsave(file.path(figure_path, p_name), switched_grads_all_path,
       width = img_width*scaling_factor, 
       height = img_width*0.85*scaling_factor, 
       units = "mm", dpi = 300, device = "pdf", bg = "white")


bf_adni_grad$tmaps |> mutate(panel = "A") |> 
  bind_rows(adni_bf_grad$tmaps |> mutate(panel = "B")) |> 
  bind_rows(bf_path_all$tmaps |> mutate(panel = "C")) |> 
  bind_rows(adni_path_all$tmaps |> mutate(panel = "D")) |> 
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient"))) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))
  


###########################
# Age and pathology collinearity
###########################

pal <- ggsci::pal_futurama()(3)


age_path <- biofinder_df |> 
  dplyr::filter(fmri_bl, !is.na(abnorm_ab)) |> 
  dplyr::mutate(
    Ab_pos = as.logical(abnorm_ab),
    Diagnosis = dplyr::case_when(
      diagnosis %in% c("Normal", "SCD") ~ "CU",
      TRUE ~ diagnosis
    ) |> factor(levels = c("CU", "MCI", "AD"))
  ) |> 
  ggplot(aes(age, pathology_ad, color = Diagnosis)) +
  geom_point(alpha = 0.2, size = 0.7) +
  labs(x = "Age", y = "Pathology Score") +
  scale_color_manual(
    values = c(pal[3], pal[1], pal[2]),
    labels = parse(text = c(
      "CU",
      "MCI~(A*beta^'+')",
      "AD~(A*beta^'+')"
    ))
  ) +
  theme_bw(base_size = 7) +
  theme(legend.position = "top")

p_name <- "age_path_dist.pdf"
ggsave(file.path(figure_path, p_name), age_path,
       width = 88, height = 60, units = "mm", 
       dpi = 300, device = "pdf", bg = "white")

biofinder_df |> 
  dplyr::filter(fmri_bl, !is.na(abnorm_ab)) |> 
  dplyr::mutate(
    Ab_pos = as.logical(abnorm_ab),
    Diagnosis = dplyr::case_when(
      diagnosis %in% c("Normal", "SCD") ~ "CU",
      TRUE ~ diagnosis
    ) |> factor(levels = c("CU", "MCI", "AD"))
  ) |> 
  select(age, pathology_ad, Diagnosis) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))

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


