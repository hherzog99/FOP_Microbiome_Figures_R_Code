# -------------------------
# HO Volume Post-GMA 
# -------------------------

# ---- Install missing packages ----
required_packages <- c("here", "readxl", "ggplot2", "dplyr", "ggsignif")

installed <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# ---- Load packages ----
library(here)
library(readxl)
library(ggplot2)
library(dplyr)
library(ggsignif)

# ---- Create output folder ----
if (!dir.exists(here("figures"))) {
  dir.create(here("figures"), recursive = TRUE)
}

# ---- Read data ----
data <- read_excel(here("data", "ho_quantification.xlsx"))
colnames(data) <- gsub("\\s+", "_", colnames(data))

plot_data <- data %>%
  filter(group %in% c("FOP No ABX", "FOP ABX")) %>%
  select(id, group, HO_Volume) %>%
  mutate(
    group = factor(group, levels = c("FOP No ABX", "FOP ABX")),
    group = recode(group,
                   "FOP No ABX" = "FOP",
                   "FOP ABX"    = "FOP GMA")
  )

# -------------------------
# Excluding sample with likely genotype error
# -------------------------
likely_genotype_error_ids <- c(789)

plot_data_filtered <- plot_data %>%
  filter(!id %in% likely_genotype_error_ids)

excluded_samples <- plot_data %>%
  filter(id %in% likely_genotype_error_ids)

# -------------------------
# Extract groups
# -------------------------
group1 <- plot_data_filtered %>%
  filter(group == "FOP") %>%
  pull(HO_Volume)

group2 <- plot_data_filtered %>%
  filter(group == "FOP GMA") %>%
  pull(HO_Volume)

# -------------------------
# Assumption testing
# -------------------------
shapiro1 <- shapiro.test(group1)
shapiro2 <- shapiro.test(group2)
var_test <- var.test(group1, group2)

# -------------------------
# Choose statistical test
# -------------------------
if (shapiro1$p.value > 0.05 && shapiro2$p.value > 0.05) {
  
  if (var_test$p.value > 0.05) {
    stat_test <- t.test(
      HO_Volume ~ group,
      data = plot_data_filtered,
      var.equal = TRUE,
      alternative = "two.sided"
    )
    test_used <- "Unpaired Student's t-test"
    
  } else {
    stat_test <- t.test(
      HO_Volume ~ group,
      data = plot_data_filtered,
      var.equal = FALSE,
      alternative = "two.sided"
    )
    test_used <- "Welch's two-sample t-test"
  }
  
} else {
  
  stat_test <- wilcox.test(
    HO_Volume ~ group,
    data = plot_data_filtered,
    alternative = "two.sided",
    exact = FALSE
  )
  
  test_used <- "Wilcoxon rank-sum test"
}

# -------------------------
# Summary stats
# -------------------------
summary_stats <- plot_data_filtered %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean = mean(HO_Volume, na.rm = TRUE),
    sd   = sd(HO_Volume, na.rm = TRUE),
    median = median(HO_Volume, na.rm = TRUE),
    min = min(HO_Volume, na.rm = TRUE),
    max = max(HO_Volume, na.rm = TRUE),
    .groups = "drop"
  )

mean_fop <- summary_stats %>% filter(group == "FOP") %>% pull(mean)
mean_gma <- summary_stats %>% filter(group == "FOP GMA") %>% pull(mean)

n_fop <- summary_stats %>% filter(group == "FOP") %>% pull(n)
n_gma <- summary_stats %>% filter(group == "FOP GMA") %>% pull(n)

mean_difference <- mean_fop - mean_gma
percent_reduction <- 100 * mean_difference / mean_fop

# p-value text for figure
p_value <- signif(stat_test$p.value, 3)

# CI if t-test
ci_text <- NA
if (grepl("t-test", test_used)) {
  ci_low <- stat_test$conf.int[1]
  ci_high <- stat_test$conf.int[2]
  ci_text <- paste0(
    "(",
    round(ci_low, 2),
    ", ",
    round(ci_high, 2),
    ")"
  )
}

# -------------------------
# Print report
# -------------------------
cat("\n================ HO VOLUME ANALYSIS REPORT ================\n")
cat("Outcome: Heterotopic Ossification Volume (mm^3)\n")
cat("Groups compared: FOP vs FOP GMA\n\n")

cat("Exclusions:\n")
if (nrow(excluded_samples) > 0) {
  print(excluded_samples)
} else {
  cat("None\n")
}
cat("\n")

cat("Sample sizes:\n")
cat("FOP: n =", n_fop, "\n")
cat("FOP GMA: n =", n_gma, "\n\n")

cat("Assumption testing:\n")
cat("Shapiro-Wilk FOP p =", signif(shapiro1$p.value, 4), "\n")
cat("Shapiro-Wilk FOP GMA p =", signif(shapiro2$p.value, 4), "\n")
cat("Variance F-test p =", signif(var_test$p.value, 4), "\n\n")

cat("Primary comparison:\n")
cat("Test used:", test_used, "\n")
cat("Sidedness: Two-sided\n")
cat("Exact p-value:", signif(stat_test$p.value, 4), "\n\n")

cat("Descriptive statistics:\n")
cat("FOP:", round(mean_fop,2), "±", round(summary_stats$sd[1],2), "\n")
cat("FOP GMA:", round(mean_gma,2), "±", round(summary_stats$sd[2],2), "\n\n")

cat("Effect size:\n")
cat("Mean difference:", round(mean_difference,2), "mm^3\n")
cat("Percent reduction:", round(percent_reduction,1), "%\n")

if (!is.na(ci_text)) {
  cat("95% CI of difference:", ci_text, "\n")
}

cat("===========================================================\n\n")

# -------------------------
# Plot
# -------------------------
y_max <- max(plot_data_filtered$HO_Volume, na.rm = TRUE)
y_limit <- ceiling(y_max * 1.35)

p <- ggplot(plot_data_filtered, aes(x = group, y = HO_Volume)) +
  geom_boxplot(
    coef = Inf,
    fill = "white",
    outlier.shape = NA,
    size = 1.1
  ) +
  geom_jitter(
    position = position_jitter(0.15),
    size = 4,
    shape = 21,
    color = "black",
    fill = "#00008B",
    stroke = 0.9
  ) +
  geom_signif(
    comparisons = list(c("FOP", "FOP GMA")),
    annotations = paste0("p = ", p_value),
    y_position = y_limit * 1.05,
    textsize = 6,
    tip_length = 0.02,
    size = 1.4
  ) +
  labs(
    x = NULL,
    y = expression("HO Volume (mm"^3 * ")")
  ) +
  theme_classic(base_family = "Helvetica") +
  ylim(0, y_limit * 1.15) +
  annotate(
    "text",
    x = 2,
    y = y_limit * 0.50,
    label = paste0("↓ ", round(percent_reduction,1), "% reduction"),
    size = 5,
    fontface = "bold",
    family = "Helvetica"
  )

print(p)

# -------------------------
# Save
# -------------------------
ggsave(
  here("figures", "HO_Volume_PostGMA.tiff"),
  p,
  width = 12,
  height = 14,
  units = "cm",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  here("figures", "HO_Volume_PostGMA.pdf"),
  p,
  width = 12,
  height = 14,
  units = "cm"
)
