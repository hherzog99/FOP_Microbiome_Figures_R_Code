# -------------------------
# Trauma-induced Anti-IL-1 HO boxplot
# -------------------------

# -------------------------
# Install packages (run once)
# -------------------------
install.packages("readxl")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("ggsignif")
install.packages("writexl")
install.packages("here")

# -------------------------
# Load libraries
# -------------------------
library(readxl)
library(ggplot2)
library(dplyr)
library(ggsignif)
library(writexl)
library(here)

# -------------------------
# Helpers
# -------------------------
round_half_up <- function(x, digits = 2) {
  posneg <- sign(x)
  z <- abs(x) * 10^digits
  z <- z + 0.5 + sqrt(.Machine$double.eps)
  z <- trunc(z)
  posneg * z / 10^digits
}

choose_two_group_test <- function(df, group1, group2, value_col = "Traumatic HO Volume",
                                  group_col = "treatment_given", alpha = 0.05) {
  sub <- df %>%
    filter(.data[[group_col]] %in% c(group1, group2)) %>%
    mutate(group_tmp = factor(as.character(.data[[group_col]]),
                              levels = c(group1, group2))) %>%
    filter(!is.na(.data[[value_col]]))
  
  x <- sub %>% filter(group_tmp == group1) %>% pull(.data[[value_col]])
  y <- sub %>% filter(group_tmp == group2) %>% pull(.data[[value_col]])
  
  n1 <- length(x)
  n2 <- length(y)
  
  if (n1 < 2 || n2 < 2) {
    return(data.frame(
      comparison = paste(group1, "vs", group2),
      group1 = group1,
      group2 = group2,
      n1 = n1,
      n2 = n2,
      shapiro_p_group1 = NA_real_,
      shapiro_p_group2 = NA_real_,
      variance_test_p = NA_real_,
      chosen_test = NA_character_,
      sidedness = NA_character_,
      raw_p_value = NA_real_,
      display_p_value = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  shapiro_p1 <- if (n1 >= 3 && n1 <= 5000) shapiro.test(x)$p.value else NA_real_
  shapiro_p2 <- if (n2 >= 3 && n2 <= 5000) shapiro.test(y)$p.value else NA_real_
  
  both_normal <- !is.na(shapiro_p1) && !is.na(shapiro_p2) &&
    shapiro_p1 >= alpha && shapiro_p2 >= alpha
  
  if (both_normal) {
    var_p <- var.test(x, y)$p.value
    
    if (!is.na(var_p) && var_p >= alpha) {
      test_obj <- t.test(x, y, var.equal = TRUE, alternative = "two.sided")
      chosen_test <- "Student t-test"
    } else {
      test_obj <- t.test(x, y, var.equal = FALSE, alternative = "two.sided")
      chosen_test <- "Welch t-test"
    }
  } else {
    var_p <- NA_real_
    test_obj <- wilcox.test(x, y, exact = FALSE, alternative = "two.sided")
    chosen_test <- "Wilcoxon rank-sum test"
  }
  
  raw_p <- test_obj$p.value
  display_p <- round_half_up(raw_p, 2)
  
  data.frame(
    comparison = paste(group1, "vs", group2),
    group1 = group1,
    group2 = group2,
    n1 = n1,
    n2 = n2,
    shapiro_p_group1 = shapiro_p1,
    shapiro_p_group2 = shapiro_p2,
    variance_test_p = var_p,
    chosen_test = chosen_test,
    sidedness = "Two-sided",
    raw_p_value = raw_p,
    display_p_value = display_p,
    stringsAsFactors = FALSE
  )
}

# -------------------------
# Paths
# -------------------------
xlsx_path <- here("data", "anti_il_1_ho_quants.xlsx")
fig_dir   <- here("figures")

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

tiff_file  <- file.path(fig_dir, "anti_il_1_ho_quants.tiff")
pdf_file   <- file.path(fig_dir, "anti_il_1_ho_quants.pdf")
excel_file <- file.path(fig_dir, "anti_il_1_ho_quants_p_value.xlsx")

# -------------------------
# Read data
# -------------------------
data <- read_excel(xlsx_path)

# -------------------------
# Filter + format
# -------------------------
excluded_ids <- c(602, 561, 222, 258, 261)

plot_data <- data %>%
  filter(
    treatment_given %in% c("BSUR01", "IgG2A", "VEH"),
    !mouse_id %in% excluded_ids,
    !is.na(`Traumatic HO Volume`),
    !is.na(`sex`)
  ) %>%
  select(mouse_id, `sex`, treatment_given, `Traumatic HO Volume`) %>%
  mutate(
    treatment_given = factor(treatment_given, levels = c("VEH", "IgG2A", "BSUR01")),
    sex = factor(`sex`, levels = c("M", "F")),
    sex_label = recode(sex, "M" = "Male", "F" = "Female")
  )

excluded_samples <- data %>%
  filter(
    treatment_given %in% c("BSUR01", "IgG2A", "VEH"),
    mouse_id %in% excluded_ids
  )

# -------------------------
# Pairwise comparisons (Wilcoxon for all)
# -------------------------
run_wilcox <- function(df, group1, group2) {
  sub <- df %>%
    filter(treatment_given %in% c(group1, group2)) %>%
    mutate(treatment_given = factor(as.character(treatment_given),
                                    levels = c(group1, group2)))
  
  x <- sub %>% filter(treatment_given == group1) %>% pull(`Traumatic HO Volume`)
  y <- sub %>% filter(treatment_given == group2) %>% pull(`Traumatic HO Volume`)
  
  shapiro_p1 <- if (length(x) >= 3 && length(x) <= 5000) shapiro.test(x)$p.value else NA_real_
  shapiro_p2 <- if (length(y) >= 3 && length(y) <= 5000) shapiro.test(y)$p.value else NA_real_
  
  test <- wilcox.test(`Traumatic HO Volume` ~ treatment_given,
                      data = sub,
                      exact = FALSE,
                      alternative = "two.sided")
  
  data.frame(
    comparison = paste(group1, "vs", group2),
    group1 = group1,
    group2 = group2,
    n1 = sum(sub$treatment_given == group1),
    n2 = sum(sub$treatment_given == group2),
    shapiro_p_group1 = shapiro_p1,
    shapiro_p_group2 = shapiro_p2,
    chosen_test = "Wilcoxon rank-sum test",
    sidedness = "Two-sided",
    raw_p_value = test$p.value,
    display_p_value = round_half_up(test$p.value, 2),
    stringsAsFactors = FALSE
  )
}

stats_veh_bsur <- run_wilcox(plot_data, "VEH", "BSUR01")
stats_igg_bsur <- run_wilcox(plot_data, "IgG2A", "BSUR01")
stats_veh_igg  <- run_wilcox(plot_data, "VEH", "IgG2A")

pairwise_stats <- bind_rows(
  stats_veh_bsur,
  stats_igg_bsur,
  stats_veh_igg
)

print(pairwise_stats)

# make p labels for plot annotations
p_veh_bsur <- stats_veh_bsur$display_p_value
p_igg_bsur <- stats_igg_bsur$display_p_value
p_veh_igg  <- stats_veh_igg$display_p_value

# -------------------------
# Summary stats
# -------------------------
summary_stats <- plot_data %>%
  group_by(treatment_given) %>%
  summarise(
    mean = mean(`Traumatic HO Volume`, na.rm = TRUE),
    sd   = sd(`Traumatic HO Volume`, na.rm = TRUE),
    median = median(`Traumatic HO Volume`, na.rm = TRUE),
    min = min(`Traumatic HO Volume`, na.rm = TRUE),
    max = max(`Traumatic HO Volume`, na.rm = TRUE),
    n    = n(),
    male = sum(sex == "M", na.rm = TRUE),
    female = sum(sex == "F", na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

mean_veh  <- summary_stats %>% filter(treatment_given == "VEH") %>% pull(mean)
mean_igg  <- summary_stats %>% filter(treatment_given == "IgG2A") %>% pull(mean)
mean_bsur <- summary_stats %>% filter(treatment_given == "BSUR01") %>% pull(mean)

pct_label_veh <- paste0(
  "↓ ", round(100 * (mean_veh - mean_bsur) / mean_veh, 1), "% vs VEHICLE"
)
pct_label_igg <- paste0(
  "↓ ", round(100 * (mean_igg - mean_bsur) / mean_igg, 1), "% vs IgG2A"
)

# -------------------------
# Console Stats Reporting
# -------------------------
cat("\n================ ANTI-IL-1 HO ANALYSIS REPORT ================\n")
cat("Outcome: Trauma-Induced Gastrocnemius HO Volume (mm^3)\n")
cat("Groups compared: VEH, IgG2A, BSUR01\n\n")

cat("Excluded mouse IDs:\n")
if (nrow(excluded_samples) > 0) {
  print(excluded_samples)
} else {
  cat("None\n")
}
cat("\n")

cat("Group summary:\n")
print(summary_stats)
cat("\n")

cat("Pairwise comparisons:\n")
print(pairwise_stats)
cat("\n")

cat("Primary inferential tests used for figure annotations:\n")
cat("- VEH vs BSUR01: Wilcoxon rank-sum test, two-sided, p =", signif(stats_veh_bsur$raw_p_value, 4), "\n")
cat("- IgG2A vs BSUR01: Wilcoxon rank-sum test, two-sided, p =", signif(stats_igg_bsur$raw_p_value, 4), "\n")
cat("- VEH vs IgG2A: Wilcoxon rank-sum test, two-sided, p =", signif(stats_veh_igg$raw_p_value, 4), "\n")
cat("==============================================================\n\n")

# -------------------------
# Save stats + p-values to Excel
# -------------------------
write_xlsx(
  list(
    summary_stats = summary_stats,
    pairwise_stats = pairwise_stats,
    plot_data_used = plot_data,
    excluded_samples = excluded_samples
  ),
  path = excel_file
)

# -------------------------
# Plot parameters
# -------------------------
axis_text_size  <- 20
axis_title_size <- 24
title_size      <- 24
linesize <- 1.5
dpi <- 300

y_limit <- ceiling(max(plot_data$`Traumatic HO Volume`, na.rm = TRUE) * 1.45)

sex_fills <- c("Male" = "#00008B", "Female" = "#00008B")

annotation_df <- data.frame(
  xmin = c("VEH", "IgG2A", "VEH"),
  xmax = c("BSUR01", "BSUR01", "IgG2A"),
  annotations = c(
    paste0("p = ", p_veh_bsur),
    paste0("p = ", p_igg_bsur),
    paste0("p = ", p_veh_igg)
  ),
  y_position = c(y_limit * 1.02, y_limit * 1.10, y_limit * 1.18),
  stringsAsFactors = FALSE
)

# -------------------------
# Plot
# -------------------------
p <- ggplot(plot_data, aes(x = treatment_given, y = `Traumatic HO Volume`)) +
  geom_boxplot(
    coef = Inf,
    fill = "white",
    outlier.shape = NA,
    size = 1.1
  ) +
  geom_jitter(
    aes(fill = sex_label),
    position = position_jitter(0.15),
    size = 4,
    shape = 21,
    color = "black",
    stroke = 0.9,
    show.legend = FALSE
  ) +
  geom_signif(
    data = annotation_df,
    aes(xmin = xmin, xmax = xmax, annotations = annotations, y_position = y_position),
    manual = TRUE,
    textsize = 6,
    tip_length = 0.02,
    size = 1.4
  ) +
  labs(
    x = NULL,
    y = expression("Trauma-Induced Gastrocnemius HO Volume (mm"^3*")"),
    title = ""
  ) +
  theme_classic() +
  ylim(0, y_limit * 1.25) +
  scale_fill_manual(values = sex_fills, guide = "none") +
  scale_x_discrete(labels = c(
    "VEH" = "VEHICLE",
    "IgG2A" = "IgG2A",
    "BSUR01" = "BSUR01"
  )) +
  theme(
    axis.line   = element_line(size = linesize),
    axis.text.y = element_text(size = axis_text_size),
    axis.text.x = element_text(size = axis_text_size),
    axis.title.y = element_text(size = axis_title_size),
    plot.title  = element_text(size = title_size, hjust = 0.5)
  ) +
  annotate(
    "text",
    x = 3,
    y = y_limit * 0.17,
    label = pct_label_veh,
    size = 6,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 3,
    y = y_limit * 0.09,
    label = pct_label_igg,
    size = 6,
    fontface = "bold"
  )

print(p)

# -------------------------
# Save
# -------------------------
ggsave(
  filename = tiff_file,
  plot = p,
  width = 18,
  height = 20,
  units = "cm",
  dpi = dpi
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 18,
  height = 20,
  units = "cm",
  dpi = dpi
)
