###########################
# Cognition
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

