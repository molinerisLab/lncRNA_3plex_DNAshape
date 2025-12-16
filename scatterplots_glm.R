library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)
library(viridis)
generate_all_lnc_scatterplots <- function(stability_col, df_path){
  set.seed(42) 
  
  # --- 1. LOAD ALL DATA ---
  full_df <- read.table(gzfile(df_path), header = TRUE)
  full_df$pos_neg <- factor(full_df$pos_neg, levels = c("neg", "pos"))
  
  # --- 2. DEFINE TRAINING SET (ROBUST ONLY) ---
  training_targets <- c("NEAT1", "TERC", "MEG3", "lncSmad7", "HOTAIR", "CDKN2B-AS1")
  
  # Create the Training Dataframe (subset)
  df_training_subset <- full_df[full_df$lncRNA %in% training_targets, ]
  
  # Balance the Training Set (Downsampling)
  min_rows <- min(table(df_training_subset$lncRNA))
  df_training_subset %>% group_by(lncRNA) %>% sample_n(min_rows) -> df_balanced_train
  
  # --- 3. TRAIN MODELS ---
  formula_full <- as.formula(paste("pos_neg ~", stability_col, "+ HelT + MGW + ProT + Roll"))
  formula_stab <- as.formula(paste("pos_neg ~", stability_col))
  
  model_full <- glm(formula_full, family = binomial(link = "logit"), data = df_balanced_train, na.action = na.exclude)
  model_stab <- glm(formula_stab, family = binomial(link = "logit"), data = df_balanced_train)
  
  # --- 4. PREDICT ON *ALL* lncRNAs ---
  df_all_positives <- subset(full_df, pos_neg == "pos")
  
  df_all_positives$score_full_shape <- predict(model_full, newdata = df_all_positives, type = "response")
  df_all_positives$score_stability  <- predict(model_stab, newdata = df_all_positives, type = "response")
  
  # --- 5. PLOTTING ---
  
  # Combined Plot (All points together)
  p_combined <- ggplot(df_all_positives, aes(x = score_stability, y = score_full_shape)) +
    geom_point(alpha = 0.3, color = "dodgerblue4", size = 1.2) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + 
    theme_bw() +
    labs(
      title = "Impact of DNA Shape (All lncRNAs)",
      subtitle = paste("Trained on Robust Set (", length(training_targets), "lncRNAs) -> Predicted on All"),
      x = "Prediction Score: Stability Only",
      y = "Prediction Score: Stability + Shape"
    ) +
    coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1))

  # Faceted Plot (One panel per lncRNA found in the file)
  p_faceted <- p_combined + 
    facet_wrap(~lncRNA, ncol = 5) + 
    labs(title = "Impact of DNA Shape by lncRNA (All Available)")

  return(list(combined = p_combined, faceted = p_faceted))
}

stability_metric = "Stability_best"
file_path = "best_param_3plex/ALL_shape.3plex_stability.matrix.gz"

plots <- generate_all_lnc_scatterplots(stability_metric, file_path)

# Save Combined
ggsave("Scatter_Positives_ALL_Combined.pdf", plot = plots$combined, width = 8, height = 8)

ggsave("Scatter_Positives_ALL_Faceted.pdf", plot = plots$faceted, width = 15, height = 10)



####################### Plot Positives and Negatives toghether ###########################


generate_comparison_scatterplots <- function(stability_col, df_path){
  set.seed(42) 
  
  # --- 1. LOAD ALL DATA ---
  full_df <- read.table(gzfile(df_path), header = TRUE)
  
  # Set Factor Levels: "pos" is 1st, "neg" is 2nd.
full_df$pos_neg <- factor(full_df$pos_neg, levels = c("neg", "pos"))
  
  # --- 2. DEFINE TRAINING SET (ROBUST ONLY) ---
  training_targets <- c("NEAT1", "TERC", "MEG3", "lncSmad7", "HOTAIR", "CDKN2B-AS1")
  
  # Create Training Subset
  df_training_subset <- full_df[full_df$lncRNA %in% training_targets, ]
  
  # Balance the Training Set (Downsampling)
  min_rows <- min(table(df_training_subset$lncRNA))
  df_balanced_train <- df_training_subset %>% 
    group_by(lncRNA) %>% 
    sample_n(min_rows) %>%
    ungroup()
  
  # --- 3. TRAIN MODELS ---
  formula_full <- as.formula(paste("pos_neg ~", stability_col, "+ HelT + MGW + ProT + Roll"))
  formula_stab <- as.formula(paste("pos_neg ~", stability_col))
  
  model_full <- glm(formula_full, family = binomial(link = "logit"), data = df_balanced_train, na.action = na.exclude)
  model_stab <- glm(formula_stab, family = binomial(link = "logit"), data = df_balanced_train)
  
  # --- 4. PREDICT ON EVERYTHING ---
  full_df$score_full_shape <- predict(model_full, newdata = full_df, type = "response")
  full_df$score_stability  <- predict(model_stab, newdata = full_df, type = "response")
  
  # --- 5. PREPARE PLOTTING DATA (SORTING) ---
  # We arrange by 'pos_neg'. Since levels are c("pos", "neg"), 
  full_df_sorted <- full_df %>% arrange(pos_neg)
  
  # Create a second subset for the "Training Only" plots
  train_only_sorted <- full_df_sorted %>% filter(lncRNA %in% training_targets)
  
    # --- 6. PLOTTING FUNCTION ---
  make_plot <- function(data, title_suffix) {
    ggplot(data, aes(x = score_stability, y = score_full_shape, color = pos_neg)) +
      # Alpha and size
      geom_point(alpha = 0.4, size = 1.2) + 
      # Add diagonal line
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + 
      theme_bw() +
      scale_color_manual(values = c("pos" = "#0571b0", "neg" = "#f4a582")) +
      labs(
        title = paste("DNA Shape Impact:", title_suffix),
        subtitle = "Negatives (Orange) plotted on top of Positives (Blue)",
        x = "Prediction Score: Stability Only",
        y = "Prediction Score: Stability + Shape",
        color = "Class"
      ) +
      coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1))
  }
  
  # --- 7. GENERATE PLOTS ---
  
  # Set A: ALL lncRNAs
  p_combined_all <- make_plot(full_df_sorted, "All lncRNAs (Combined)")
  p_faceted_all  <- make_plot(full_df_sorted, "All lncRNAs (Faceted)") + 
    facet_wrap(~lncRNA, ncol = 5)
  
  # Set A-NEW: Faceted by lncRNA AND pos_neg
  p_faceted_all_posneg <- make_plot(full_df_sorted, "All lncRNAs (Faceted by Class)") +
    facet_grid(lncRNA ~ pos_neg) +
    labs(subtitle = "Rows: lncRNA, Columns: Class (neg/pos)") +
    theme(
      panel.spacing.x = unit(2, "lines"),  # horizontal spacing between neg/pos columns
      panel.spacing.y = unit(1.5, "lines")   # vertical spacing between lncRNA rows
    )

  # Set B: TRAINING lncRNAs Only
  p_combined_train <- make_plot(train_only_sorted, "Training Set Only (Combined)")
  p_faceted_train  <- make_plot(train_only_sorted, "Training Set Only (Faceted)") + 
    facet_wrap(~lncRNA, ncol = 3)
  
  # Set B-NEW: Faceted by lncRNA AND pos_neg for training set
  p_faceted_train_posneg <- make_plot(train_only_sorted, "Training Set (Faceted by Class)") +
    facet_grid(pos_neg ~ lncRNA) +
    labs(subtitle = "Rows: Class (neg/pos), Columns: lncRNA")

  p_faceted_meg3 <- make_plot(
    full_df_sorted %>% filter(lncRNA == "MEG3"), 
    "MEG3 Only (Faceted by Class)"
  ) +
    facet_grid(lncRNA ~ pos_neg) +
    labs(subtitle = "Rows: Class (neg/pos), Columns: lncRNA")
  
  p_faceted_train_by_class <- make_plot(train_only_sorted, "Training Set (Faceted by Class Only)") +
    facet_wrap(~ pos_neg, ncol = 2) +
    labs(subtitle = "Training lncRNAs grouped by Class")

  p_faceted_all_by_class <- make_plot(full_df_sorted, "All lncRNAs (Faceted by Class Only)") +
    facet_wrap(~ pos_neg, ncol = 2) +
    labs(subtitle = "All lncRNAs grouped by Class")

   p_faceted_all_by_class <- make_plot(full_df_sorted, "All lncRNAs (Faceted by Class Only)") +
    facet_wrap(~ pos_neg, ncol = 2) +
    labs(subtitle = "All lncRNAs grouped by Class")

  p_faceted_all_by_class_density <- ggplot(full_df_sorted, aes(x = score_stability, y = score_full_shape)) +
    geom_pointdensity(adjust = 4, size = 0.5, alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
        scale_color_viridis_c(
      option = "D",
      name = "Local Density"
    ) +
    facet_wrap(~ pos_neg, ncol = 2) +
    theme_bw() +
    labs(
      title = "DNA Shape Impact: All lncRNAs (Density by Class)",
      subtitle = "Points colored by local density",
      x = "Prediction Score: Stability Only",
      y = "Prediction Score: Stability + Shape"
    ) +
    coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
    theme(
      strip.background = element_rect(fill = "#f0f0f0"),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11)
    )

  return(list(
    all_combined = p_combined_all,
    all_faceted = p_faceted_all,
    all_faceted_posneg = p_faceted_all_posneg,
    all_faceted_by_class = p_faceted_all_by_class,  
    all_faceted_by_class_density = p_faceted_all_by_class_density,  
    train_combined = p_combined_train,
    train_faceted = p_faceted_train,
    train_faceted_posneg = p_faceted_train_posneg,
    train_faceted_by_class = p_faceted_train_by_class, 
    meg3_faceted = p_faceted_meg3
  ))
}

# --- EXECUTION ---
stability_metric = "Stability_best"
file_path = "best_param_3plex/ALL_shape.3plex_stability.matrix.gz"

plots <- generate_comparison_scatterplots(stability_metric, file_path)

# --- SAVE PLOTS ---

# 1. All lncRNAs
ggsave("Scatter_ALL_Combined_PosNeg.pdf", plot = plots$all_combined, width = 8, height = 8)
ggsave("Scatter_ALL_Faceted_PosNeg.pdf", plot = plots$all_faceted, width = 15, height = 12)
ggsave("Scatter_ALL_Faceted_ByClass.pdf", plot = plots$all_faceted_posneg, width = 10, height = 40)
ggsave("Scatter_ALL_Faceted_ClassOnly.pdf", plot = plots$all_faceted_by_class, width = 10, height = 5)  
ggsave("Scatter_ALL_Faceted_ClassOnly_Density.pdf", plot = plots$all_faceted_by_class_density, width = 10, height = 5)  # NEW
ggsave("Scatter_MEG3_Faceted_ByClass.pdf", plot = plots$meg3_faceted, width = 6, height = 6)

# 2. Training Only
ggsave("Scatter_Training_Combined_PosNeg.pdf", plot = plots$train_combined, width = 8, height = 8)
ggsave("Scatter_Training_Faceted_PosNeg.pdf", plot = plots$train_faceted, width = 10, height = 8)
ggsave("Scatter_Training_Faceted_ByClass.pdf", plot = plots$train_faceted_posneg, width = 12, height = 6)
ggsave("Scatter_Training_Faceted_ClassOnly.pdf", plot = plots$train_faceted_by_class, width = 10, height = 5)