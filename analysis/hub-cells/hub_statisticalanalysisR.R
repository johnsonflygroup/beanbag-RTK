---
title: "hub analysis final"
author: "######"
date: "2025-10-23"
output:
   pdf_document:
    fig_caption: yes
    latex_engine: xelatex
always_allow_html: true
header-includes:
  - \usepackage{float}
  - \floatplacement{figure}{H}
  - \usepackage{booktabs}
  - \usepackage{caption}
  - \captionsetup[table]{position=bottom}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)

library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(emmeans)
library(multcompView)
library(broom)
library(kableExtra)
library(stringr)
library(tidyr)
library(readxl)
library(tibble)



# Significance stars
signif_stars <- function(p) {
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01, "**",
  ifelse(p < 0.05, "*",
  ifelse(p < 0.1, ".", ""))))
}


# Pretty y-axis labels with correct units
y_labels <- c(
  nuclei_in_hub = "Nuclei in Hub",
  hub_volume_um3 = "Hub Volume ($\\mu$m³)",
  hub_sphericity = "Hub Sphericity",
  hub_solidity = "Hub Solidity",
  hub_surface_area = "Surface Area ($\\mu$m²)",
  hub_equivalent_diameter = "Equivalent Diameter ($\\mu$m)",
  hub_major_axis_length = "Major Axis Length ($\\mu$m)",
  hub_minor_axis_length = "Minor Axis Length ($\\mu$m)",
  hub_extent = "Hub Extent",
  log_volume = "Log Hub Volume ($\\mu$m³)"
)

# Analysis + plotting function
analyze_and_plot_report <- function(df, response_var) {
  formula <- as.formula(paste(response_var, "~ genotype"))
  model <- aov(formula, data = df)
  y_label <- y_labels[[response_var]]
  if (is.null(y_label)) y_label <- response_var


  anova_tidy <- broom::tidy(model)
  
  em <- emmeans(model, ~ genotype)
  tukey <- contrast(em, method = "tukey") %>% summary(infer = TRUE)
  tukey$Significance <- signif_stars(tukey$p.value)
  
  cld_res <- multcomp::cld(em, Letters = letters, adjust = "tukey")
  
  max_vals <- df %>%
    group_by(genotype) %>%
    summarise(max_value = max(.data[[response_var]], na.rm = TRUE), .groups = "drop")
  
  offset <- 0.1 * max(df[[response_var]], na.rm = TRUE)
  
  plot_data <- cld_res %>%
    left_join(max_vals, by = c("genotype")) %>%
    mutate(y_pos = max_value + offset)


  
  p <- ggplot(df, aes(x = genotype, y = .data[[response_var]], fill = genotype)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    theme_bw(base_size = 14) +
    labs(x = "Genotype", y = y_label) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    geom_text(data = plot_data,
              aes(x = genotype, y = y_pos, label = .group),
              inherit.aes = FALSE,
              size = 4.5,
              fontface = "bold")
  
  list(model = model, anova_tidy = anova_tidy, tukey = tukey, plot = p)
}

plot_metric_if_valid <- function(df, response_var, label) {
  formula <- as.formula(paste(response_var, "~ genotype"))
  model <- aov(formula, data = df)
  sw <- shapiro.test(residuals(model))
  lev <- car::leveneTest(formula, data = df)
  lev_p <- lev[1, "Pr(>F)"]

  # Label
  y_label <- y_labels[[response_var]]
  if (is.null(y_label)) y_label <- response_var

  # If ANOVA assumptions pass
  if (sw$p.value > 0.05 & lev_p > 0.05) {
    cat("\n\n## Valid ANOVA:", response_var, "\n\n")
    result <- analyze_and_plot_report(df, response_var)
    print(result$plot)

  } else {
    # Run Kruskal-Wallis
    kw <- kruskal.test(formula, data = df)

    if (kw$p.value < 0.05) {
      cat("\n\n## Valid Kruskal-Wallis:", response_var, "\n\n")

      # Pairwise Wilcoxon
      pw <- pairwise.wilcox.test(df[[response_var]], df$genotype, p.adjust.method = "BH")
      pw_matrix <- pw$p.value

      # Compact letter display from emmeans (fallback to ANOVA model)
      em <- emmeans(model, ~ genotype)
      cld <- multcomp::cld(em, Letters = letters, adjust = "tukey")

      max_vals <- df %>% group_by(genotype) %>% summarise(max_value = max(.data[[response_var]], na.rm = TRUE), .groups = "drop")
      offset <- 0.1 * max(df[[response_var]], na.rm = TRUE)

      plot_data <- cld %>% left_join(max_vals, by = "genotype") %>% mutate(y_pos = max_value + offset)

      p <- ggplot(df, aes(x = genotype, y = .data[[response_var]], fill = genotype)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
        theme_bw(base_size = 14) +
        labs(x = "Genotype", y = y_label) +
        theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
        geom_text(data = plot_data, aes(x = genotype, y = y_pos, label = .group), inherit.aes = FALSE, size = 4.5, fontface = "bold")

      print(p)

    } else {
      cat("\n\n## Skipped (No valid test):", response_var, "\n\n")
    }
  }
}

analyze_metric <- function(df, response_var, label) {
  checks <- check_anova_assumptions(df, response_var)
  y_label <- y_labels[[response_var]] %||% response_var
  formula <- as.formula(paste(response_var, "~ genotype"))

  if (checks$valid_anova) {
    cat("\n\n## Valid ANOVA:", response_var, "\n\n")
    res <- analyze_and_plot_report(df, response_var)
    print(res$plot)
    if (!dir.exists("figures")) dir.create("figures")
    ggsave(filename = paste0("figures/", response_var, "_plot.pdf"), plot = res$plot, width = 8, height = 6)
    
    df %>%
      group_by(genotype) %>%
      summarise(
        Mean = mean(.data[[response_var]], na.rm = TRUE),
        SD = sd(.data[[response_var]], na.rm = TRUE),
        N = sum(!is.na(.data[[response_var]]))
      ) %>%
      kable(format = "latex", booktabs = TRUE, caption = paste("Genotype Summary for", label)) %>%
      kable_styling(latex_options = "hold_position") %>%
      print()
    
    print_anova_tables(res$anova_tidy, res$tukey, y_label)
    
  } else {
    kw <- kruskal.test(formula, data = df)
    
    if (kw$p.value < 0.05) {
      cat("\n\n## Valid Kruskal-Wallis:", response_var, "\n\n")
      plot_metric_if_valid(df, response_var, y_labels)
      
      df %>%
        group_by(genotype) %>%
        summarise(
          Mean = mean(.data[[response_var]], na.rm = TRUE),
          SD = sd(.data[[response_var]], na.rm = TRUE),
          N = sum(!is.na(.data[[response_var]]))
        ) %>%
        kable(format = "latex", booktabs = TRUE, caption = paste("Genotype Summary for", label)) %>%
        kable_styling(latex_options = "hold_position") %>%
        print()
      
      print_kw_table(kw, response_var, label)
      print_pw_table(df, response_var, label)
    } else {
      cat("\n\n## Skipped (No valid test):", response_var, "\n\n")
    }
  }
}


# Check ANOVA assumptions
check_anova_assumptions <- function(df, response_var) {
  formula <- as.formula(paste(response_var, "~ genotype"))
  model <- aov(formula, data = df)
  sw <- shapiro.test(residuals(model))
  lev <- car::leveneTest(formula, data = df)
  lev_p <- lev[1, "Pr(>F)"]
  list(
    model = model,
    normality = sw$p.value > 0.05,
    variance = lev_p > 0.05,
    valid_anova = sw$p.value > 0.05 & lev_p > 0.05
  )
}

# Print ANOVA tables
print_anova_tables <- function(anova_tidy, tukey, label) {
  anova_tidy %>%
    select(term, df, statistic, p.value) %>%
    rename(`F value` = statistic) %>%
    mutate(p.value = ifelse(p.value < 0.001, "<0.001", signif(p.value, 3))) %>%
    kable(format = "latex", booktabs = TRUE, caption = paste("ANOVA Results for", label)) %>%
    kable_styling(latex_options = "hold_position") %>%
    print()

  tukey %>%
    select(contrast, estimate, lower.CL, upper.CL, p.value, Significance) %>%
    rename(
      Comparison = contrast,
      Difference = estimate,
      `Lower CI` = lower.CL,
      `Upper CI` = upper.CL,
      `Adjusted p-value` = p.value
    ) %>%
    mutate(`Adjusted p-value` = signif(`Adjusted p-value`, 3)) %>%
    kable(format = "latex", booktabs = TRUE, caption = paste("Tukey Post Hoc Comparisons for", label)) %>%
    kable_styling(latex_options = "hold_position") %>%
    print()
}

# Print Kruskal-Wallis and Wilcoxon tables
print_kw_table <- function(kw, response_var, label) {
  tibble(
    Metric = label,
    `Chi-squared` = round(kw$statistic, 3),
    df = kw$parameter,
    `p-value` = signif(kw$p.value, 3),
    Significance = signif_stars(kw$p.value)
  ) %>%
    kable(format = "latex", booktabs = TRUE, caption = paste("Kruskal-Wallis Test for", label)) %>%
    kable_styling(latex_options = "hold_position") %>%
    print()
}

print_pw_table <- function(df, response_var, label) {
  pw <- pairwise.wilcox.test(df[[response_var]], df$genotype, p.adjust.method = "BH")
  pw_table <- as.data.frame(pw$p.value) %>%
    rownames_to_column("Group1") %>%
    pivot_longer(-Group1, names_to = "Group2", values_to = "Adjusted p-value") %>%
    filter(!is.na(`Adjusted p-value`)) %>%
    mutate(Significance = signif_stars(`Adjusted p-value`)) %>%
    mutate(`Adjusted p-value` = signif(`Adjusted p-value`, 3)) %>%
    arrange(Group1, Group2)

  kable(pw_table, caption = paste("Pairwise Wilcoxon Comparisons for", label)) %>%
    kable_styling(latex_options = "hold_position") %>%
    print()
}
```

```{r load-data}

data <- read.csv("D:/UNI/HUB_CELLS/EXPERIMENTAL/Hub cells/Analysisfinal_all/hub_analysis_all.csv")
data$genotype <- factor(  str_to_upper(data$genotype),
  levels = c("WT", "T2A.WT", "DF.WT", "T2A.DF")
)
voxel_volume_um3 <- 0.0661
data$hub_volume_um3 <- data$hub_volume_voxels * voxel_volume_um3

```

---


```{r metrics-analysis-loop, results='asis', fig.width=10, fig.height=6}
metrics <- c("nuclei_in_hub", "hub_volume_um3", "hub_sphericity", "hub_solidity", "hub_extent","hub_surface_area","hub_equivalent_diameter","hub_major_axis_length",
  "hub_minor_axis_length" )
for (metric in metrics) {
  label <- y_labels[[metric]] %||% metric
  analyze_metric(data, metric, label)
}


```
```{r test-summary-table, results='asis'}
# Initialize tracking vectors
test_summary <- tibble(
  Metric = character(),
  Test_Type = character(),
  Normality = character(),
  Equal_Variance = character(),
  ANOVA_Valid = character()
)

for (metric in metrics) {
  label <- y_labels[[metric]] %||% metric
  formula <- as.formula(paste(metric, "~ genotype"))
  model <- aov(formula, data = data)
  sw <- shapiro.test(residuals(model))
  lev <- car::leveneTest(formula, data = data)
  lev_p <- lev[1, "Pr(>F)"]
  kw <- kruskal.test(formula, data = data)

  normality <- ifelse(sw$p.value > 0.05, "Passed", "Failed")
  variance <- ifelse(lev_p > 0.05, "Passed", "Failed")
  anova_valid <- sw$p.value > 0.05 & lev_p > 0.05

  test_type <- if (anova_valid) {
    "ANOVA"
  } else if (kw$p.value < 0.05) {
    "Kruskal-Wallis"
  } else {
    "Skipped"
  }

  test_summary <- add_row(test_summary,
    Metric = label,
    Test_Type = test_type,
    Normality = normality,
    Equal_Variance = variance,
    ANOVA_Valid = ifelse(anova_valid, "Yes", "No")
  )
}

kable(test_summary, caption = "Statistical Test Summary by Metric") %>%
  kable_styling(latex_options = "hold_position")
```

```{r log-volume-anova-tables, results='asis'}
library(broom)

# Log-transform and run ANOVA
data$log_volume <- log(data$hub_volume_um3 + 1)
model_log_volume <- aov(log_volume ~ genotype, data = data)

# Format ANOVA table
broom::tidy(model_log_volume) %>%
  select(term, df, statistic, p.value) %>%
  rename(`F value` = statistic) %>%
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", formatC(p.value, format = "e", digits = 2))) %>%
  kable(format = "latex", booktabs = TRUE, caption = "ANOVA Results for Log-Transformed Hub Volume") %>%
  kable_styling(latex_options = "hold_position") %>%
  print()

# Tukey post hoc
library(emmeans)
library(multcompView)

emmeans(model_log_volume, ~ genotype) %>%
  contrast(method = "tukey") %>%
  summary(infer = TRUE) %>%
  mutate(Significance = ifelse(`p.value` < 0.001, "***",
                        ifelse(`p.value` < 0.01, "**",
                        ifelse(`p.value` < 0.05, "*",
                        ifelse(`p.value` < 0.1, ".", "")))) ) %>%
  select(contrast, estimate, lower.CL, upper.CL, p.value, Significance) %>%
  rename(
    Comparison = contrast,
    Difference = estimate,
    `Lower CI` = lower.CL,
    `Upper CI` = upper.CL,
    `Adjusted p-value` = p.value
  ) %>%
  mutate(`Adjusted p-value` = signif(`Adjusted p-value`, 3)) %>%
  kable(format = "latex", booktabs = TRUE, caption = "Tukey Post Hoc Comparisons for Log-Transformed Hub Volume") %>%
  kable_styling(latex_options = "hold_position") %>%
  print()
```

## Assumption Checks and ANOVA for Transformed Data

```{r log-volume-assumption-checks, results='asis', fig.width=7, fig.height=6}
library(car)

cat("\n\n## Assumption Checks for: log_volume\n\n")

model_log <- aov(log_volume ~ genotype, data = data)

# Q-Q plot
qqnorm(residuals(model_log), main = "Q-Q Plot: log_volume")
qqline(residuals(model_log), col = "blue", lwd = 2)

# Shapiro-Wilk test
sw_log <- shapiro.test(residuals(model_log))
cat("\n**Shapiro-Wilk test for normality**\n\n")
cat(paste0("W = ", round(sw_log$statistic, 4), ", p-value = ", signif(sw_log$p.value, 4), "\n"))

# Levene’s test
lev_log <- leveneTest(log_volume ~ genotype, data = data)
cat("\n**Levene’s test for equal variances**\n\n")
print(knitr::kable(lev_log, format = "markdown"))

# Residuals vs Fitted plot
plot(model_log, which = 1, main = "Residuals vs Fitted: log_volume")
```
## Assumption Summary (Transformed Data)

```{r automated-assumption-summary-transformed, results='asis'}
# Only include transformed metrics here
metrics_transformed <- c("log_volume", "hub_sphericity", "hub_solidity")
normality_results <- c()
variance_results <- c()
anova_validity <- c()

for (metric in metrics_transformed) {
  formula <- as.formula(paste(metric, "~ genotype"))
  model <- aov(formula, data = data)
  
  sw <- shapiro.test(residuals(model))
  normality_results <- c(normality_results, ifelse(sw$p.value > 0.05, "Passed", "Failed"))
  
  lev <- leveneTest(formula, data = data)
  lev_p <- lev[1, "Pr(>F)"]
  variance_results <- c(variance_results, ifelse(lev_p > 0.05, "Passed", "Failed"))
  
  anova_validity <- c(anova_validity, ifelse(sw$p.value > 0.05 & lev_p > 0.05, "Yes", "No"))
}

assumption_summary_transformed <- tibble(
  Metric = c("Hub Volume (log)", "Sphericity", "Solidity"),
  Normality = normality_results,
  Equal_Variance = variance_results,
  ANOVA_Valid = anova_validity
)

kable(assumption_summary_transformed, caption = "Assumption Check Summary (Transformed Data)") %>%
  kable_styling(latex_options = "hold_position")
```
```{r plot-log-volume-if-valid, results='asis', fig.width=10, fig.height=6}
# Only plot if assumptions passed
if (sw_log$p.value > 0.05 & lev_log[1, "Pr(>F)"] > 0.05) {
  res_log_volume <- analyze_and_plot_report(data, "log_volume")
  print(res_log_volume$plot)
}
```

```{r correlation-heatmap-pooled, results='asis', fig.width=10, fig.height=8}
library(Hmisc)
library(dplyr)
library(tidyr)
library(ggplot2)

corr_metrics <- c("nuclei_in_hub", "hub_volume_um3", "hub_sphericity", "hub_solidity",
                  "hub_surface_area", "hub_equivalent_diameter",
                  "hub_major_axis_length", "hub_minor_axis_length", "hub_extent")

# Drop rows with missing values
complete_data <- data[, corr_metrics] %>% na.omit()

# Run rcorr
corr_res <- rcorr(as.matrix(complete_data), type = "spearman")

# Extract matrices
r <- corr_res$r
p <- corr_res$P

# Tidy format
cor_df <- as.data.frame(as.table(r)) %>%
  rename(Metric_1 = Var1, Metric_2 = Var2, Spearman_rho = Freq) %>%
  mutate(p_value = as.vector(p)) %>%
  filter(Metric_1 != Metric_2)

# Filter for significant correlations
cor_df_sig <- cor_df %>% filter(p_value < 0.05)

# Rename metrics for display
pretty_names <- c(
  nuclei_in_hub = "Nuclei in Hub",
  hub_volume_um3 = "Hub Volume",
  hub_sphericity = "Hub Sphericity",
  hub_solidity = "Hub Solidity",
  hub_surface_area = "Surface Area",
  hub_equivalent_diameter = "Equivalent Diameter",
  hub_major_axis_length = "Major Axis Length",
  hub_minor_axis_length = "Minor Axis Length",
  hub_extent = "Hub Extent"
)

cor_df_sig <- cor_df_sig %>%
  mutate(
    Metric_1 = pretty_names[Metric_1],
    Metric_2 = pretty_names[Metric_2]
  )


# Plot
ggplot(cor_df_sig, aes(x = Metric_1, y = Metric_2, fill = Spearman_rho)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#00BFC4",   # sky blue 
    mid = "white",
    high = "#F8766D",  # coral pink 
    midpoint = 0,
    limits = c(-1, 1),
    name = expression(rho)
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) +
  labs(title = "Significant Spearman Correlations Across All Genotypes",
       x = "", y = "")
```