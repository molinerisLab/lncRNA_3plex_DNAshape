library(dplyr)
library(ggplot2)

generate_both_scatterplots <- function(stability_col, df_path){
  # 1. Reproducibility
  set.seed(42) 
  
  # --- LOAD DATA ---
  full_df <- read.table(gzfile(df_path), header = TRUE)
  
  # Filter for target lncRNAs
  target_lncRNAs <- c("NEAT1", "TERC", "MEG3", "lncSmad7", "HOTAIR", "CDKN2B-AS1")
  df <- full_df[full_df$lncRNA %in% target_lncRNAs, ]
  df$pos_neg <- factor(df$pos_neg, levels = c("pos", "neg"))
  
  # --- TRAINING STEP ---
  min_rows <- min(table(df$lncRNA))
  df %>% group_by(lncRNA) %>% sample_n(min_rows) -> df_training
  
  # Fit models
  formula_full <- as.formula(paste("pos_neg ~", stability_col, "+ HelT + MGW + ProT + Roll"))
  formula_stab <- as.formula(paste("pos_neg ~", stability_col))
  
  model_full <- glm(formula_full, family = binomial(link = "logit"), data = df_training, na.action = na.exclude)
  model_stab <- glm(formula_stab, family = binomial(link = "logit"), data = df_training)
  
  # --- PREDICTION STEP ---
  df_positives <- subset(df, pos_neg == "pos")
  df_positives$score_full_shape <- predict(model_full, newdata = df_positives, type = "response")
  df_positives$score_stability  <- predict(model_stab, newdata = df_positives, type = "response")
  
  # --- PLOTTING ---
  
  # 1. Create the Base Plot (Combined)
  p_combined <- ggplot(df_positives, aes(x = score_stability, y = score_full_shape)) +
    geom_point(alpha = 0.4, color = "dodgerblue4", size = 1.5) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + 
    theme_bw() +
    labs(
      title = "Impact of DNA Shape on Prediction Confidence (All lncRNAs)",
      subtitle = "Positive Binding Sites Only (Model trained on Balanced Data)",
      x = "Prediction Score: Stability Only (Baseline)",
      y = "Prediction Score: Stability + Shape (Full Model)"
    ) +
    coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1))

  # 2. Create the Faceted Plot (by adding facet_wrap to the base plot)
  p_faceted <- p_combined + 
    facet_wrap(~lncRNA, ncol = 3) +
    labs(title = "Impact of DNA Shape by lncRNA") # Update title for facet version

  # Return both plots in a list
  return(list(combined = p_combined, faceted = p_faceted))
}

# --- EXECUTION ---
stability_metric = "Stability_best"
file_path = "best_param_3plex/ALL_shape.3plex_stability.matrix.gz"

# Run the function once
plots <- generate_both_scatterplots(stability_metric, file_path)

# Save the Combined plot
print(plots$combined)
ggsave("Scatter_Positives_Combined.pdf", plot = plots$combined, width = 7, height = 7)

# Save the Faceted plot
print(plots$faceted)
ggsave("Scatter_Positives_Faceted.pdf", plot = plots$faceted, width = 10, height = 7)