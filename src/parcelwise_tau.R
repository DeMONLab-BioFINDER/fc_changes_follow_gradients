tau <- readRDS("/media/strg/repos/fmri_project1/data/tau_pet/tau_pet_extracted/tau_pet_schaefer_1000_r.rds")


health_cog_sep <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, diagnosis %in% c("Normal", "SCD")), 
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
                                              base_size_ = 24,
                                              vect = TRUE,
                                              mod_formula = formula(paste0("~ age + tau_pathology + abnorm_ab + sex + rsqa__MeanFD")),
                                              covariates = c("sex", "rsqa__MeanFD"),
                                              r2_size = rel(6),
                                              filter_criteria = quo(),
                                              show_networks = FALSE,
                                              tag_prefix = "",
                                              tag_sep = "",
                                              layout_construction = "horizontal",
                                              include_gradient_plots = TRUE,
                                              right_term_side = FALSE,
                                              plt_title = "",
                                              plt_subtitle = TRUE,
                                              cache_runs = FALSE)




nodal_regression_fits_roiwise_pred <- function(subject_data,
                                               fc_matrix, 
                                               roi_names,
                                               tau_matrix = NULL,
                                               logistic = FALSE,
                                               scale_fc = FALSE,
                                               vectorised = FALSE,
                                               id_var = "image_file",
                                               id_tau = "tau_file",
                                               dep_var = "FC",
                                               model_formula = formula(" ~ age + tau_pathology + sex + rsqa__MeanFD")) {
  
  # --- Vectorized mode (fast matrix regression) ---
  if (vectorised) {
    bio_filtered <- subject_data %>% 
      drop_na(all_of(all.vars(model_formula)), cho12)
    X_base <- model.matrix(model_formula, data = bio_filtered)
    Y <- as.data.frame(fc_matrix[bio_filtered[[id_var]], ])
    
    if (!is.null(tau_matrix)) {
      tau_df <- as.data.frame(tau_matrix[bio_filtered[[id_tau]], ])
    }
    
    similarity_fits <- list()
    for (roi in roi_names) {
      if (!is.null(tau_matrix)) {
        X <- cbind(X_base, tau_df[[roi]])
        colnames(X)[ncol(X)] <- "tau_roi"
      } else {
        X <- X_base
      }
      fit <- lm.fit(X, Y[[roi]])
      fit$formula <- update(model_formula, paste0("~ tau_roi + ", as.character(model_formula)[2]))
      similarity_fits[[roi]] <- fit
    }
    
    return(similarity_fits)
  }
  
  # --- Non-vectorized mode (loop with tidy df) ---
  inter_sub_long <- as.data.frame(fc_matrix) %>%
    rownames_to_column(id_var) %>%
    inner_join(subject_data, by = id_var) %>%
    pivot_longer(starts_with('7Networks'), names_to = 'region', values_to = dep_var)
  
  if (!is.null(tau_matrix)) {
    tau_long <- as.data.frame(tau_matrix) %>%
      rownames_to_column(id_var) %>%
      pivot_longer(starts_with('7Networks'), names_to = 'region', values_to = 'tau_val')
    
    inter_sub_long <- inter_sub_long %>%
      left_join(tau_long, by = c(id_var, "region"))
  }
  
  similarity_fits <- list()
  ticker <- 0
  pb <- txtProgressBar(min = 0, max = length(roi_names), initial = 0) 
  
  for (roi in roi_names) {
    reg_dat <- inter_sub_long %>% filter(region == roi)
    if (scale_fc) reg_dat <- reg_dat %>% mutate(dep_var = scale(.data[[dep_var]]))
    
    if (!is.null(tau_matrix)) {
      # Add tau term dynamically
      formula_with_tau <- update(model_formula, paste("~ tau_val +", as.character(model_formula)[2]))
    } else {
      formula_with_tau <- model_formula
    }
    
    if (!logistic) {
      fit <- lm(formula_with_tau, data = reg_dat)
    } else {
      fit <- glm(formula_with_tau, data = reg_dat, family = binomial(link = "logit"))
    }
    
    similarity_fits[[roi]] <- fit
    ticker <- ticker + 1
    setTxtProgressBar(pb, ticker)
  }
  
  close(pb)
  return(similarity_fits)
}

test <- nodal_regression_fits_roiwise_pred(biofinder_df |> filter(fmri_bl),
                                           fc_matrix = fc_measures_bf$affinity,
                                           roi_names = rois,
                                           #tau_matrix = tau, 
                                           vectorised = TRUE,
                                           model_formula = formula(" ~ age + cho56 + sex + rsqa__MeanFD"))




ests <- get_nodal_ests(test)



ests <- list_of_ests[[1]] %>% rename_with(~ names(list_of_ests[1]), statistic)

if (length(list_of_ests) > 1) {
  for(i in 2:length(list_of_ests)) {
    ests <- inner_join(ests, list_of_ests[[i]] %>% rename_with(~ names(list_of_ests[i]), statistic), by = c("term", "region"))
  }
}

plot_brain_ests <- function(ests, tag = "a", covariates = c("sex", "rsqa__MeanFD"), atlas_geometry = readRDS("data/atlas_data/schaef_ggseg2.rds")) {
  parcel_line_size = 0.1
  #bsize = 16
  
  terms_of_interest <- unique(ests$term)[!(unique(ests$term) %in% c(covariates, "(Intercept)"))]
  plots_of_terms <- list()

  i = 1
  for (term_of_i in terms_of_interest) {
    
    p <-  ests %>%
      filter(term == term_of_i) %>%
      right_join(atlas_geometry$atlas, by = "region") %>%
      ggplot() +
      geom_sf(aes(
        fill = statistic,
        geometry = geometry), linewidth= 0.1,
        show.legend = FALSE) +
      (if (add_shade)
        geom_sf(data = atlas_geometry$shade, size = shade_size, alpha = shade_alpha)
       else NULL) +
      theme_void(base_size = base_size_)+
      labs(fill = 't', title = str_to_title(str_replace(term_of_i, "_", " "))
      ) +
      theme(#legend.position = "",
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        plot.title = element_text(color = "black", hjust = 0.5)
      ) +
      scale_fill_gradient2(
        low = muted("blue"),
        mid = "white",
        high = muted("red") 
      )
    
    i = i + 1
    

    plots_of_terms[[term_of_i]]<- p
  }
  
  return(plots_of_terms)
}

b_e <- plot_brain_ests(ests)

wrap_plots(b_e)




make_gradient_relationship_plots <- function(
    ests,
    gradient_data,
    grad_char,
    terms_of_interest,
    net_names,
    gradient_colors,
    layout_construction = "horizontal",
    right_term_side = FALSE,
    gray_out = FALSE,
    spintest = TRUE,
    perms = NULL,
    r2_size = rel(4.6),
    side_color_bar = TRUE,
    base_size_ = 11
) {
  require(tidyverse)
  require(ggpmisc)
  require(ggtext)
  require(ggside)
  require(ggnewscale)
  
  plots <- list()
  count <- 1
  
  for (g in grad_char) {
      for (term_ in terms_of_interest) {
        
        plot_data <- gradient_data %>%
          inner_join(ests %>% filter(term == term_), by = "region")
        
        i <- t
        j <- g
        
        # Swap axes if layout is "vertical"
        if (layout_construction == "vertical") {
          i <- g
          j <- t
        }
        
        lab_grad <- "Gradient score"
        lab_term <- paste(str_to_title(str_replace(t, "_", " ")), "t-value")
        
        if (layout_construction == "horizontal") {
          x_lab <- lab_term
          y_lab <- lab_grad
        } else {
          y_lab <- lab_term
          x_lab <- lab_grad
        }
        
        x_min <- min(plot_data[[i]], na.rm = TRUE)
        x_max <- max(plot_data[[i]], na.rm = TRUE)
        y_min <- min(plot_data[[j]], na.rm = TRUE)
        y_max <- max(plot_data[[j]], na.rm = TRUE)
        
        p <- ggplot(plot_data, aes(x = .data[[i]], y = .data[[j]], color = if (!gray_out) name else "gray")) +
          geom_point(alpha = ifelse(gray_out, 0.05, 0.2)) +
          stat_poly_line(se = FALSE, color = ifelse(gray_out, "#808080", "#323232")) +
          labs(x = x_lab, y = y_lab, color = "Network") +
          xlim(x_min, x_max) +
          theme_bw(base_size = base_size_) +
          theme(
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA)
          ) +
          scale_color_manual(values = net_names %>% select(name, col) %>% deframe()) +
          scale_y_continuous(limits = c(y_min, y_max), position = ifelse(right_term_side, "right", "left"))
        
        # ---------------- Spin test or standard correlation ----------------
        if (spintest) {
          r <- cor(plot_data[[i]], plot_data[[j]], method = "pearson", use = "complete.obs")
          p_val <- if (!is.null(perms)) perm_sphere_p(plot_data[[i]], plot_data[[j]], perm.id = perms, corr.type = "pearson") else NA
          p_lab <- scales::label_pvalue(
            accuracy = 0.001,
            prefix = c("<i>p<sub>spin</sub></i> &lt; ", "<i>p<sub>spin</sub></i> = ", "<i>p<sub>spin</sub></i> &gt; ")
          )(p_val)
          label_corr <- paste0("<i>r</i> = ", round(r, 2), " ,&#8203;&nbsp;", p_lab)
          
          p <- p + labs(subtitle = label_corr) +
            theme(plot.subtitle = element_markdown(
              size = rel(0.75),
              hjust = 0.95,
              margin = margin(b = 3),
              color = ifelse(gray_out, scales::alpha("black", 0.2), "black")
            ))
        } else {
          p <- p + stat_poly_eq(
            aes(label = paste(after_stat(rr.label),
                              str_remove(after_stat(rr.confint.label), "95% CI "),
                              sep = "*\" \"*")),
            parse = TRUE, color = "#323232", label.x = "left", label.y = "top", size = r2_size
          )
        }
        
        # ---------------- Gray-out overlay ----------------
        if (gray_out) {
          p <- p + ggpp::annotate(
            "text_npc",
            npcx = "middle",
            npcy = "middle",
            label = "?",
            size = 28
          )
        }
        
        # ---------------- Side color bars ----------------
        if (side_color_bar) {
          x_range_df <- data.frame(x = seq(x_min, x_max, length.out = 100))
          colnames(x_range_df)[1] <- i
          
          y_range_df <- data.frame(y = seq(y_min, y_max, length.out = 100))
          colnames(y_range_df)[1] <- j
          
          if (layout_construction == "horizontal") {
            alpha_switch_xside <- ifelse(g == tail(grad_char, 1), 1, 0)
            alpha_switch_yside <- ifelse(term_ == terms_of_interest[1], 1, 0)
            x_axis_colorbar <- c(muted("blue"), muted("red"))
            y_axis_colorbar <- gradient_colors[[g]]
          } else { # vertical
            alpha_switch_xside <- ifelse(term_ == tail(terms_of_interest, 1), 1, 0)
            alpha_switch_yside <- ifelse(g == grad_char[1], 1, 0)
            x_axis_colorbar <- gradient_colors[[g]]
            y_axis_colorbar <- c(muted("blue"), muted("red"))
          }
          
          if (right_term_side && layout_construction == "vertical") {
            alpha_switch_yside <- ifelse(g == tail(grad_char, 1), 1, 0)
          }
          
          p <- p +
            geom_xsidetile(
              data = x_range_df,
              aes(x = .data[[i]], y = 0, fill = .data[[i]]),
              alpha = alpha_switch_xside,
              show.legend = FALSE,
              inherit.aes = FALSE
            ) +
            scale_fill_gradient2(low = x_axis_colorbar[1], mid = "white", high = x_axis_colorbar[2]) +
            ggnewscale::new_scale_fill() +
            geom_ysidetile(
              data = y_range_df,
              aes(y = .data[[j]], x = 0, fill = .data[[j]]),
              alpha = alpha_switch_yside,
              show.legend = FALSE,
              inherit.aes = FALSE
            ) +
            scale_fill_gradient2(low = y_axis_colorbar[1], mid = "white", high = y_axis_colorbar[2]) +
            theme_ggside_void() +
            ggside(x.pos = "bottom", y.pos = ifelse(right_term_side, "right", "left")) +
            theme(ggside.panel.scale = 0.04)
        }
        
        plots[[count]] <- p
        count <- count + 1
      }
    
  }
  
  return(plots)
}

unique_terms <- ests$term |> 
  unique() |> 
  (\(x) x[!(x %in% c("(Intercept)", "rsqa__MeanFD", "tau_roi"))])()

scatters <- make_gradient_relationship_plots(
  ests,
  grad_df |> filter(study == "biofinder"),
  grad_char = c("gradient1", "gradient3"),
  terms_of_interest = unique_terms,
  net_names = net_names,
  gradient_colors = gradient_cols
)



plot_gradient_relationships <- function(subject_data,
                                        gradient_data,
                                        list_of_parcel_data, 
                                        atlas_geometry = readRDS("data/atlas_data/schaef_ggseg2.rds"),
                                        add_shade = FALSE,
                                        shade_alpha = 0.01,
                                        shade_size = 0.01,
                                        vect = FALSE,
                                        base_size_ = 11,
                                        grad_name = TRUE,
                                        gradients = c(1, 3),
                                        gradient_colors = NULL,
                                        empty_row_height = 0,
                                        padding = 0,
                                        r2_size = rel(4.6),
                                        spintest = TRUE,
                                        perms = readRDS("data/atlas_data/permutations_1000_hungarian.rds"),
                                        id_var = "image_file",
                                        mod_formula = formula(paste0("FC ~ age")),
                                        covariates = c("sex", "rsqa__MeanFD"),
                                        filter_criteria = quo(),
                                        plt_title = NULL,
                                        plt_subtitle = FALSE,
                                        tag_sep = "",
                                        tag_prefix = "",
                                        layout_construction = "horizontal",
                                        right_term_side = FALSE,
                                        include_gradient_plots = TRUE,
                                        gray_out = FALSE,
                                        side_color_bar = TRUE,
                                        plot_spacing = 0.3,
                                        show_networks = FALSE,
                                        network_geometry = NULL,
                                        cache_runs = FALSE,
                                        longitudinal = FALSE,
                                        logistic_fit = FALSE,
                                        scale_fc = FALSE,
                                        sub_id = "sid",
                                        longitudinal_formula = formula(paste0("FC ~ age + (1 | ", sub_id, ")"))
){
  
  source("src/util.R")
  require(tidyverse)
  require(scales)
  require(patchwork)
  require(ggpmisc)
  require(sf)
  require(ggtext)
  
  
  old <- theme_set(theme_bw(base_size = base_size_))
  theme_update(panel.background = element_rect(fill = "transparent", colour = NA),
               plot.background = element_rect(fill = "transparent", colour = NA),
               legend.background = element_rect(fill = "transparent", colour = NA),
               legend.box.background = element_rect(fill = "transparent", colour = NA))
  
  net_names <- data.frame(name = c('Vis', 'SomMot', 'DorsAttn','SalVentAttn','Limbic', 'Cont', 'Default'),
                          col = c("#781286", "#4682B4", "#00760E", "#C43AFA", "#c7cc7a", "#E69422", "#CD3E4E"), #"#DCF8A4"
                          label = c(1:7))
  
  grad_char <- paste0("gradient", gradients)

  
  
  


  ## Gradient plots 
  gradient_plots <- list()
  tag_labs <- paste0(tag_prefix, tag_sep, "a",tag_sep, tolower(as.roman(1:length(gradients))))
  i = 1
  for (stud in unique(gradient_data$study)) {
    for (grad in grad_char) {
      
      if (include_gradient_plots) {
        gradient_plots[[paste0(stud,"_", grad)]] <- gradient_data %>% filter(study==stud) %>% 
          right_join(atlas_geometry$atlas, by = "region") %>%
          ggplot() +
          geom_sf(aes(
            fill = .data[[grad]],
            geometry = geometry), linewidth= 0.1,
            show.legend = FALSE)+
          (if (add_shade)
            geom_sf(data = atlas_geometry$shade, size = shade_size, alpha = shade_alpha)
           else NULL) +
          theme_void(base_size = base_size_)+
          labs(fill = "", title = str_to_title(str_replace(paste0("_", grad), "_", " ")),
               #tag = tag_labs[i]
          ) + 
          (if (grad_name)
            labs(title = case_when(
              grad == "gradient1" ~ "Sens-Assoc",
              grad == "gradient2" ~ "Vis-Mot",
              grad == "gradient3" ~ "Social-Exec",
              TRUE ~ grad
            ))  else NULL) +
          theme(legend.position = "",
                panel.background = element_rect(fill = "transparent", colour = NA),
                plot.background = element_rect(fill = "transparent", colour = NA),
                legend.background = element_rect(fill = "transparent", colour = NA),
                legend.box.background = element_rect(fill = "transparent", colour = NA),
                plot.title = element_text(color = "black", hjust = 0.5)
          ) +
          scale_fill_gradient2(
            low = gradient_colors[[grad]][1],
            mid = "white",
            high = gradient_colors[[grad]][2] 
          ) 
        
        i = i +1
      } else {
        gradient_plots[[paste0(stud,"_", grad)]] <- gradient_data %>% filter(study==stud) %>% 
          inner_join(atlas_geometry$atlas, by = "region") %>%
          ggplot() +
          geom_sf(aes(
            fill = .data[[grad]],
            geometry = geometry), alpha = 0, linewidth= 0.1, color = NA,
            show.legend = FALSE)+
          theme_void(base_size = base_size_)+
          labs(fill = "", title = str_to_title(str_replace(paste0("_", grad), "_", " ")),
               #tag = tag_labs[i]
          ) +
          theme(legend.position = "",
                panel.background = element_rect(fill = "transparent", colour = NA),
                plot.background = element_rect(fill = "transparent", colour = NA),
                legend.background = element_rect(fill = "transparent", colour = NA),
                legend.box.background = element_rect(fill = "transparent", colour = NA),
                plot.title = element_text(color = NA, hjust = 0.5)
          ) +
          scale_fill_gradient2(
            low = gradient_colors[[grad]][1],
            mid = "white",
            high = gradient_colors[[grad]][2] 
          ) 
        i = i +1
      }
      
      
    }
  }
  
  
  
  vars_of_interest <- analysis_name
  n_analysis <- length(analysis_name)
  terms_of_interest <- unique(ests$term)[!(unique(ests$term) %in% c(covariates, "(Intercept)"))]
  n_terms <- length(terms_of_interest)
  n_plts_row <- length(terms_of_interest)*length(vars_of_interest)
  tag_labs <- c(length(analysis_name))
  
  n_gradients <- length(grad_char)
  
  labels <- letters[2:(n_analysis+1)]  
  tag_labs  <- c()  
  
  for (g in 0:(n_gradients - 1)) {
    num_seq <- (g * n_terms + 1):((g + 1) * n_terms) 
    for (label in labels) {
      temp <- paste0(tag_prefix, tag_sep, label, tag_sep, num_seq)  
      tag_labs  <- c(tag_labs , temp) 
    }
  }
  
  plots <- list()
  count <- 1
  for (g in grad_char) {
    for (t in vars_of_interest) {
      for (term_ in terms_of_interest){
        plot_data  <- gradient_data %>% inner_join(ests %>% filter(term == term_), by = "region") 
        
        i <- t
        j <- g
        
        # Swap axes if layout is "vertical"
        if (layout_construction == "vertical") {
          i <- g
          j <- t
        }
        
        lab_grad <- "Gradient score"
        lab_term <- paste(str_to_title(str_replace(t, "_", " ")), "t-value")
        
        if (layout_construction == "horizontal") {
          x_lab <- lab_term
          y_lab <- lab_grad
        } else {
          y_lab <- lab_term
          x_lab <- lab_grad
        }
        
        
        
        x_min <- min(plot_data[i])
        x_max <- max(plot_data[i])
        
        y_min <- min(plot_data[j])
        y_max <- max(plot_data[j])
        
        
        
        p <- plot_data %>% 
          ggplot(aes(x = .data[[i]], y =.data[[j]],
                     color = if (!gray_out) name else "gray")) +
          geom_point(alpha = ifelse(gray_out, 0.05, 0.2)) +
          stat_poly_line(se = FALSE, 
                         color = ifelse(gray_out, "#808080", "#323232")) +
          labs(#title = term_,
            x = x_lab,
            y = y_lab,
            #tag = tag_labs[count],
            color = "Network") +
          xlim(x_min, x_max)+
          #ylim(y_min, y_max)+
          theme_bw(base_size = base_size_) +
          theme(
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            #axis.title = element_text(size = rel(1.2))
          ) +
          scale_color_manual(values = net_names %>% select(name, col) %>% deframe()) +
          scale_y_continuous(limits = c(y_min, y_max), position = ifelse(right_term_side, "right", "left"))
        
        if (spintest) {
          r <- cor(plot_data[[i]], plot_data[[j]], method = "pearson")
          p_val <- perm_sphere_p(plot_data[[i]], plot_data[[j]], perm.id = perms, corr.type='pearson')
          
          #p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("italic(p[spin]) < ", "italic(p[spin]) == ", "italic(p[spin]) > "))(p_val)
          #label_corr <- paste0("italic(r) == ", r |> round(2),"*','~", p_lab)
          p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("<i>p<sub>spin</sub></i> &lt; ", "<i>p<sub>spin</sub></i> = ", "<i>p<sub>spin</sub></i> &gt; "))(p_val)
          label_corr <- paste0("<i>r</i> = ", r |> round(2), " ,&#8203;&nbsp;" , p_lab)
          
          
          p <- p + labs(subtitle = label_corr) +
            theme(plot.subtitle = element_markdown(size = rel(0.75), hjust = 0.95, margin = margin(b = 3),
                                                   color = ifelse(gray_out, scales::alpha("black", 0.2), "black")))
          
          #   subtitle = label = label_corr, alpha = ifelse(gray_out, 0.2, 1)
          # p <- p + ggpp::annotate(geom = "text_npc", label = label_corr,
          #                         alpha = ifelse(gray_out, 0.2, 1),
          #                           size = r2_size,
          #                           npcx = "left", npcy = "top",
          #                           family = "sans",
          #                           parse = TRUE)
        } else {
          p <- p + stat_poly_eq(aes(label = paste(after_stat(rr.label),
                                                  str_remove(after_stat(rr.confint.label), "95% CI "),
                                                  sep = "*\" \"*")),
                                parse = TRUE, color = "#323232", label.x = "left", label.y = "top", size = r2_size)
        }
        
        if (gray_out) {
          p <- p + ggpp::annotate("text_npc",
                                  npcx = "middle", 
                                  npcy = "middle",
                                  label = "?",
                                  size = 28
          )
        }
        
        
        if (side_color_bar) {
          require(ggside)
          x_range_df <- data.frame(x = seq(x_min, x_max, length.out = 100))
          colnames(x_range_df)[1] <- i
          
          y_range_df <- data.frame(y = seq(y_min, y_max, length.out = 100))
          colnames(y_range_df)[1] <- j
          
          if(layout_construction == "horizontal") {
            alpha_switch_xside <- ifelse(g == tail(grad_char, 1), 1, 0)
            alpha_switch_yside <- ifelse(term_ == terms_of_interest[1], 1, 0)
            x_axis_colorbar <- c(muted("blue"), muted("red"))
            y_axis_colorbar <- gradient_colors[[g]]
          }
          if(layout_construction == "vertical") {
            alpha_switch_xside <- ifelse(term_ == tail(terms_of_interest, 1), 1, 0)
            alpha_switch_yside <- ifelse(g == grad_char[1], 1, 0)
            x_axis_colorbar <- gradient_colors[[g]]
            y_axis_colorbar <- c(muted("blue"), muted("red"))
          }
          
          if(right_term_side & layout_construction == "vertical") alpha_switch_yside <- ifelse(g == tail(grad_char, 1), 1, 0)
          
          
          
          p <- p +
            geom_xsidetile(data = x_range_df, aes(x = .data[[i]], y = 0, fill = .data[[i]]), 
                           alpha = alpha_switch_xside,
                           show.legend = FALSE, 
                           inherit.aes = FALSE) +
            scale_fill_gradient2(
              low = x_axis_colorbar[1],
              mid = "white",
              high = x_axis_colorbar[2] 
            ) +
            ggnewscale::new_scale_fill() +
            geom_ysidetile(data = y_range_df, aes(y = .data[[j]], x = 0, fill = .data[[j]]), 
                           alpha = alpha_switch_yside,
                           show.legend = FALSE,
                           inherit.aes = FALSE) +
            scale_fill_gradient2(
              low = y_axis_colorbar[1],
              mid = "white",
              high = y_axis_colorbar[2]
            ) +
            theme_ggside_void() +
            ggside(x.pos = "bottom", y.pos = ifelse(right_term_side, "right", "left")) +
            theme(ggside.panel.scale = 0.04)
        }
        
        plots[[count]] <- p
        count <- count + 1
      }
    }
  }
  
  n_terms = length(unique(ests$term)[!(unique(ests$term) %in% c(covariates, "(Intercept)"))])
  n_analysis <- length(list_of_parcel_data)
  
  list_of_brain_plots_ests <- list()
  letter_tag <- letters[2:(n_analysis+1)]
  i = 1
  for (analysis in names(list_of_ests)){
    list_of_brain_plots_ests[[analysis]] <- plot_brain_ests(list_of_ests[[analysis]], tag = letter_tag[i])
    i = i + 1
  }
  
  
  if (layout_construction == "horizontal") {
    n_plot_cols <- length(1:((n_terms*n_analysis) + (n_analysis - 1) + 1))
    col_widths <- c(1, rep(c(rep(1, n_terms), plot_spacing), n_analysis))[-(n_plot_cols+1)]
    empty_cols <- seq(1, n_plot_cols, by = n_terms+1)[-1]
    
    n_plot_rows <- n_gradients + 1 + (n_gradients)  
    row_heights <- rep(1, n_plot_rows)  
    empty_rows <- seq(2, n_plot_rows, by = 2)  
    
    row_heights[empty_rows] <- empty_row_height
    
    layout <- c(
      area(3, 1)  # First gradient
    )
    
    if (n_gradients>1) {
      for (g in seq(5, n_plot_rows, by =2)){
        layout <- c(layout, area(g, 1))
      }
    }
    
    for (col in 2:n_plot_cols){
      if (col %in% empty_cols) {
      } else {
        layout <- c(layout, area(1, col))
      }
    }
    
    for(i in seq(3, (n_gradients*2-1)+2, by = 2)){
      for (j in 2:((n_terms*n_analysis) + (n_analysis - 1) + 1)) {
        if ((j %in% empty_cols)) {
        } else {
          layout <- c(layout, area(i, j))
        }
      }
    }
    
    for(empt in empty_cols){
      layout <- c(layout, area(1, empt , b = n_gradients + 1))
    }
    
    if (padding != 0) {
      for(pad in 1:padding){
        layout <- c(layout, area(1, n_plot_cols + pad, b = n_gradients + 1))
      }
    }
    
    
  }
  
  
  if (layout_construction == "vertical") {
    n_plot_rows <- (n_terms * n_analysis) + (n_analysis - 1) + 1
    row_heights <- c(1, rep(c(rep(1, n_terms), plot_spacing), n_analysis))[-(n_plot_rows+1)]
    col_widths <- NULL
    empty_rows <- seq(1, n_plot_rows, by = n_terms + 1)[-1]
    
    if (right_term_side) {
      layout <- c(
        area(1, 1)
      )
      
      for (g in 2:(length(gradients))) {
        layout <- c(layout, area(1, g))
      }
      
      for (row in 2:n_plot_rows) {
        if (!(row %in% empty_rows)) {
          layout <- c(layout, area(row, (length(gradients) + 1)))
        }
      }
      
      for (j in 1:(length(gradients))) {
        for (i in 2:n_plot_rows) {
          if (!(i %in% empty_rows)) {
            layout <- c(layout, area(i, j))
          }
        }
      }
      
    } else {
      layout <- c(
        area(1, 2)
      )
      
      for (g in 3:(length(gradients) + 1)) {
        layout <- c(layout, area(1, g))
      }
      
      for (row in 2:n_plot_rows) {
        if (!(row %in% empty_rows)) {
          layout <- c(layout, area(row, 1))
        }
      }
      
      for (j in 2:(length(gradients) + 1)) {
        for (i in 2:n_plot_rows) {
          if (!(i %in% empty_rows)) {
            layout <- c(layout, area(i, j))
          }
        }
      }
    }
    
  }
  
  
  filt_char = as_label(filter_criteria )
  
  
  brain_plots <- unlist(list_of_brain_plots_ests, recursive = FALSE)
  
  plots_to_include <- c(gradient_plots, brain_plots, plots)
  
  if (plt_subtitle) {
    f <- mod_formula
    rhs <- paste0(f)[2]
    rhs <- gsub("\\*", "×", rhs)
    expr_str <- paste0("FC[parcel] ~ '~' ~ ", rhs)
    #subtit_expr <- parse(text = expr_str)[[1]]
    subtit_expr <- bquote(FC[parcel] ~ "~" ~ .(rhs))
  } else {
    subtit_expr <- ""
  }
  
  if (!is.null(plt_title)) {
    n = ests %>% pull(n) %>% unique()
    plt_title <- paste0(plt_title, "(N = ", n, ")")
  } 
  
  p <- Reduce(`+`, plots_to_include) +
    plot_annotation(title = plt_title, subtitle = subtit_expr,
                    #tag_levels = c('A', '1'),
                    theme = theme(
                      plot.title = element_text(#size = 28,
                        hjust = 0),
                      plot.subtitle = element_text(size = rel(0.9),
                                                   hjust = 0,
                                                   vjust = -0.05,
                                                   family = "mono",
                                                   face = "italic")
                    )) +
    plot_layout(design = layout, axis_titles = "collect", axes = "collect", guides = "collect",
                widths = col_widths,
                heights = row_heights
    ) & 
    theme(legend.position = "",
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          plot.tag.position  = if (right_term_side) c(0.1, 0.96) else c(0.9, 0.96),
          #plot.tag = element_text(size = 8, hjust = 0, vjust = 0)
    )
  list(plot = p, n = ests %>% pull(n) %>% unique(), model_formula = ests %>% pull(model_formula) %>% unique(), tmaps = ests)
  
}