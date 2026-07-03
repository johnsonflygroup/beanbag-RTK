---
title: "Statistical Analysis gut"
author: "####"
date: "2025-08-10"
output:
  pdf_document:
    fig_caption: yes
header-includes:
  - \usepackage{float}
  - \floatplacement{figure}{H}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)

library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(emmeans)
library(multcompView)
library(broom)
library(kableExtra)
library(car)
library(stringr)
library(effectsize)



# Significance stars function
signif_stars <- function(p) {
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01, "**",
  ifelse(p < 0.05, "*",
  ifelse(p < 0.1, ".", ""))))
}

# Pretty y-axis labels
y_labels <- c(
  n_nuclei = "Number of Nuclei",
  total_integrated_intensity = "Integrated DAPI Intensity (a.u.)",
  mean_intensity = "Mean DAPI Intensity per Nucleus (a.u.)"
)

```


```{r load-data}
data <- read_excel("data_for_r.xlsx")
data$genotype <- factor(data$genotype)
data$region <- factor(data$region)
data <- data %>%
  mutate(mean_intensity = total_integrated_intensity / n_nuclei)
data$genotype <- factor(
str_to_upper(data$genotype),
  levels = c("WT", "T2A.WT", "DF.WT", "T2A.DF")
)
data$region <- factor(
  data$region, levels = c("R1", "R2", "R3", "R4", "R5", "R6")
  )

```

```{r define-analysis-function}
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
  max_vals <- df %>% group_by(genotype) %>% summarise(max_value = max(.data[[response_var]], na.rm = TRUE), .groups = "drop")
  offset <- 0.1 * max(df[[response_var]], na.rm = TRUE)
  plot_data <- cld_res %>% left_join(max_vals, by = "genotype") %>% mutate(y_pos = max_value + offset)

  p <- ggplot(df, aes(x = genotype, y = .data[[response_var]], fill = genotype)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    theme_bw(base_size = 14) +
    labs(x = "Genotype", y = y_label) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    geom_text(data = plot_data, aes(x = genotype, y = y_pos, label = .group), inherit.aes = FALSE, size = 4.5, fontface = "bold")

  list(model = model, anova_tidy = anova_tidy, tukey = tukey, plot = p)
}

analyze_two_way <- function(df, response_var) {
  formula <- as.formula(paste(response_var, "~ genotype * region"))
  model <- aov(formula, data = df)
  y_label <- y_labels[[response_var]]
  if (is.null(y_label)) y_label <- response_var

  # Assumption checks
  sw <- shapiro.test(residuals(model))
  lev <- car::leveneTest(formula, data = df)
  lev_p <- lev[1, "Pr(>F)"]
  assumptions <- tibble(
    Metric = y_label,
    Normality = ifelse(sw$p.value > 0.05, "Passed", "Failed"),
    Equal_Variance = ifelse(lev_p > 0.05, "Passed", "Failed"),
    ANOVA_Valid = ifelse(sw$p.value > 0.05 & lev_p > 0.05, "Yes", "No")
  )

  # ANOVA summary
  anova_tidy <- broom::tidy(model)

  # Plot
  p <- ggplot(df, aes(x = genotype, y = .data[[response_var]], fill = genotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  facet_wrap(~ region, ncol = 3) +
  theme_bw(base_size = 14) +
  labs(x = "Genotype", y = y_label) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


  list(model = model, assumptions = assumptions, anova_tidy = anova_tidy, plot = p)
}
```


```{r nuclei-analysis, results='asis', fig.width=10, fig.height=6}
res_nuclei <- analyze_and_plot_report(data, "n_nuclei")

res_nuclei$anova_tidy %>%
  select(term, df, statistic, p.value) %>%
  rename(`F value` = statistic) %>%
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", formatC(p.value, format = "e", digits = 2))) %>%
  kable(format = "latex", booktabs = TRUE, caption = "ANOVA Results for Number of Nuclei") %>%
  kable_styling(latex_options = "hold_position") %>%
  print()

eta_nuclei <- eta_squared(res_nuclei$model)
print(knitr::kable(eta_nuclei, caption = "Effect Size (η²) for Genotype on Number of Nuclei"))

res_nuclei$tukey %>%
  select(contrast, estimate, lower.CL, upper.CL, p.value, Significance) %>%
  rename(
    Comparison = contrast,
    Difference = estimate,
    `Lower CI` = lower.CL,
    `Upper CI` = upper.CL,
    `Adjusted p-value` = p.value
  ) %>%
  mutate(`Adjusted p-value` = signif(`Adjusted p-value`, 3)) %>%
  kable(format = "latex", booktabs = TRUE, caption = "Tukey Post Hoc Comparisons for Number of Nuclei") %>%
  kable_styling(latex_options = "hold_position") %>%
  print()

print(res_nuclei$plot)
```


```{r intensity-analysis, results='asis', fig.width=10, fig.height=6}
res_intensity <- analyze_and_plot_report(data, "mean_intensity")

res_intensity$anova_tidy %>%
  select(term, df, statistic, p.value) %>%
  rename(`F value` = statistic) %>%
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", formatC(p.value, format = "e", digits = 2)))%>%
  kable(format = "latex", booktabs = TRUE, caption = "ANOVA Results for Mean DAPI Intensity per Nucleus") %>%
  kable_styling(latex_options = "hold_position") %>%
  print()

eta_intensity <- eta_squared(res_intensity$model)
print(knitr::kable(eta_intensity, caption = "Effect Size (η²) for Genotype on Mean DAPI Intensity"))


res_intensity$tukey %>%
  select(contrast, estimate, lower.CL, upper.CL, p.value, Significance) %>%
  rename(
    Comparison = contrast,
    Difference = estimate,
    `Lower CI` = lower.CL,
    `Upper CI` = upper.CL,
    `Adjusted p-value` = p.value
  ) %>%
  mutate(`Adjusted p-value` = signif(`Adjusted p-value`, 3)) %>%
  kable(format = "latex", booktabs = TRUE, caption = "Tukey Post Hoc Comparisons for Mean DAPI Intensity per Nucleus") %>%
  kable_styling(latex_options = "hold_position") %>%
  print()
```
\clearpage

```{r intensity-plot, fig.width=10, fig.height=7}
print(res_intensity$plot)
```

```{r two-way-mean-intensity, results='asis', fig.width=10, fig.height=6}
res_mean <- analyze_two_way(data, "mean_intensity")

print(res_mean$assumptions)

res_mean$anova_tidy %>%
  select(term, df, statistic, p.value) %>%
  rename(`F value` = statistic) %>%
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", signif(p.value, 3))) %>%
  kable(format = "latex", booktabs = TRUE, caption = "Two-Way ANOVA for Mean DAPI Intensity") %>%
  kable_styling(latex_options = "hold_position") %>%
  print()

eta_two_way <- eta_squared(res_mean$model)
print(knitr::kable(eta_two_way, caption = "Effect Sizes (η²) for Genotype, Region, and Interaction"))


print(res_mean$plot)
```

```{r two-way-nuclei, results='asis', fig.width=10, fig.height=6}
res_mean <- analyze_two_way(data, "n_nuclei")

print(res_mean$assumptions)

res_mean$anova_tidy %>%
  select(term, df, statistic, p.value) %>%
  rename(`F value` = statistic) %>%
  mutate(p.value = ifelse(p.value < 0.001, "<0.001", signif(p.value, 3))) %>%
  kable(format = "latex", booktabs = TRUE, caption = "Two-Way ANOVA for Number of Nuclei") %>%
  kable_styling(latex_options = "hold_position") %>%
  print()

eta_two_way <- eta_squared(res_mean$model)
print(knitr::kable(eta_two_way, caption = "Effect Sizes (η²) for Genotype, Region, and Interaction"))


print(res_mean$plot)
```

```{r assumption-checks, results='asis'}
metrics <- c("n_nuclei", "mean_intensity")
normality_results <- c()
variance_results <- c()
anova_validity <- c()

for (metric in metrics) {
  formula <- as.formula(paste(metric, "~ genotype"))
  model <- aov(formula, data = data)
  sw <- shapiro.test(residuals(model))
  lev <- leveneTest(formula, data = data)
  lev_p <- lev[1, "Pr(>F)"]

  normality_results <- c(normality_results, ifelse(sw$p.value > 0.05, "Passed", "Failed"))
  variance_results <- c(variance_results, ifelse(lev_p > 0.05, "Passed", "Failed"))
  anova_validity <- c(anova_validity, ifelse(sw$p.value > 0.05 & lev_p > 0.05, "Yes", "No"))
}

tibble(
  Metric = y_labels[metrics],
  Normality = normality_results,
  Equal_Variance = variance_results,
  ANOVA_Valid = anova_validity
) %>%
  kable(caption = "Assumption Check Summary") %>%
  kable_styling(latex_options = "hold_position")
```
```{r kruskal-tests, results='asis'}
metrics <- c("n_nuclei", "mean_intensity")

kw_results <- lapply(metrics, function(metric) {
  test <- kruskal.test(as.formula(paste(metric, "~ genotype")), data = data)
  tibble(
    Metric = y_labels[[metric]],
    `Chi-squared` = round(test$statistic, 3),
    df = test$parameter,
    `p-value` = signif(test$p.value, 4),
    Significance = signif_stars(test$p.value)
  )
})

kw_table <- do.call(rbind, kw_results)

kable(kw_table, caption = "Kruskal-Wallis Test Results") %>%
  kable_styling(latex_options = "hold_position")
```

```{r nuclei-intensity-correlation, fig.width=7, fig.height=6, fig.align='center', results='asis'}
# Calculate correlation by genotype and region
library(ggplot2)
library(dplyr)

cor_results <- data %>%
  group_by(genotype) %>%
  summarise(
    cor = cor(n_nuclei, total_integrated_intensity, method = "pearson"),
    p = cor.test(n_nuclei, total_integrated_intensity)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      p < 0.1 ~ ".",
      TRUE ~ ""
    )
  )

cor_results_mean <- data %>%
  group_by(genotype) %>%
  summarise(
    cor = cor(n_nuclei, mean_intensity, method = "pearson"),
    p = cor.test(n_nuclei, mean_intensity)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      p < 0.1 ~ ".",
      TRUE ~ ""
    )
  )

print(knitr::kable(cor_results, caption = "Pearson Correlation between Number of Nuclei and Total DAPI Intensity"))
print(knitr::kable(cor_results_mean, caption = "Pearson Correlation between Number of Nuclei and Mean DAPI Intensity"))


# Plot with regression lines by group
ggplot(data, aes(x = n_nuclei, y = total_integrated_intensity, color = genotype)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw(base_size = 14) +
  labs(
    x = "Number of Nuclei",
    y = "Integrated Fluorescent Intensity",
    color = "Genotype",
    title = "Correlation of Nuclei and Intensity by Genotype"
  )

ggplot(data, aes(x = n_nuclei, y = mean_intensity, color = genotype)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw(base_size = 14) +
  labs(
    x = "Number of Nuclei",
    y = "Mean DAPI Intensity per Nucleus",
    color = "Genotype",
    title = "Correlation of Nuclei and Mean Intensity by Genotype"
  )

```
