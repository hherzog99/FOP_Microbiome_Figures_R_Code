# =========================
# GMA Survival KM curve
# =========================

# ===== Install missing packages =====
required_packages <- c(
  "here",
  "survival",
  "survminer",
  "ggplot2",
  "readxl",
  "dplyr",
  "tidyr",
  "stringr"
)

installed <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# ===== Libraries =====
library(here)
library(survival)
library(survminer)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

# ===== Create figures folder if needed =====
if (!dir.exists(here("figures"))) {
  dir.create(here("figures"), recursive = TRUE)
}

# ===== parameters =====
yaxis      <- 26
ytext      <- 22
xaxis      <- 22
linesize   <- 1.8
annot_size <- 7

# MMP9 grey for numbers + legend
mmp9_grey <- "grey40"

# KM colors
km_colors <- c("#0072B2", "#D55E00")

# ===== Read data from data folder =====
metadata_raw <- read_excel(here("data", "FOP_ABX_KM_Curve.xlsx"))

metadata <- metadata_raw %>%
  transmute(
    Sex     = toupper(str_trim(as.character(Sex))),
    Days_Alive = suppressWarnings(as.numeric(Days_Alive)),
    Condition  = str_squish(as.character(Condition))
  )

# Assume all are events unless Event column exists
if (!"Event" %in% names(metadata_raw)) {
  metadata$Event <- 1L
} else {
  metadata$Event <- suppressWarnings(as.integer(metadata_raw$Event))
}

# ===== Normalize groups =====
metadata <- metadata %>%
  mutate(
    cond_norm = toupper(Condition),
    is_GMA    = str_detect(cond_norm, "GMA"),
    Group     = if_else(is_GMA, "FOP GMA", "FOP")
  ) %>%
  filter(!is.na(Days_Alive), !is.na(Event)) %>%
  mutate(
    Event = as.integer(Event),
    Group = factor(Group, levels = c("FOP", "FOP GMA"))
  )

# ===== Survival fit =====
surv_obj <- with(metadata, Surv(time = Days_Alive, event = Event))
fit <- survfit(surv_obj ~ Group, data = metadata)

# ===== Legend labels with n =====
n_per_group <- metadata %>%
  count(Group) %>%
  complete(Group = levels(metadata$Group), fill = list(n = 0)) %>%
  arrange(Group)

legend_labels <- c(
  sprintf("FOP (n=%d)",     n_per_group$n[n_per_group$Group == "FOP"]),
  sprintf("FOP GMA (n=%d)", n_per_group$n[n_per_group$Group == "FOP GMA"])
)

# ===== Log-rank test =====
logrank <- survdiff(Surv(Days_Alive, Event) ~ Group, data = metadata)
logrank_p <- 1 - pchisq(logrank$chisq, df = length(logrank$n) - 1)
p_text <- paste0("Log-rank p = ", formatC(logrank_p, format = "f", digits = 3))

# ===== Mean survival by group =====
mean_survival <- metadata %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    mean_survival_days = mean(Days_Alive, na.rm = TRUE),
    sd_survival_days   = sd(Days_Alive, na.rm = TRUE),
    min_survival_days  = min(Days_Alive, na.rm = TRUE),
    max_survival_days  = max(Days_Alive, na.rm = TRUE),
    events = sum(Event == 1, na.rm = TRUE),
    censored = sum(Event == 0, na.rm = TRUE),
    .groups = "drop"
  )

# ===== Print stats report =====
cat("\n================ SURVIVAL ANALYSIS REPORT ================\n")
cat("Outcome: Days Alive\n")
cat("Groups compared: FOP vs FOP GMA\n\n")

cat("Sample sizes:\n")
for (i in seq_len(nrow(mean_survival))) {
  cat(" ", as.character(mean_survival$Group[i]),
      ": n =", mean_survival$n[i], "\n")
}
cat("\n")

cat("Event counts:\n")
for (i in seq_len(nrow(mean_survival))) {
  cat(" ", as.character(mean_survival$Group[i]),
      ": events =", mean_survival$events[i],
      ", censored =", mean_survival$censored[i], "\n")
}
cat("\n")

cat("Primary survival comparison:\n")
cat("  Test: Log-rank test\n")
cat("  Sidedness: Two-sided\n")
cat("  Test statistic (chi-square):", round(logrank$chisq, 4), "\n")
cat("  Degrees of freedom:", length(logrank$n) - 1, "\n")
cat("  Exact p-value:", signif(logrank_p, 4), "\n\n")

cat("Mean survival by group:\n")
for (i in seq_len(nrow(mean_survival))) {
  cat(" ", as.character(mean_survival$Group[i]),
      ": mean =", round(mean_survival$mean_survival_days[i], 2), "days",
      ", SD =", round(mean_survival$sd_survival_days[i], 2),
      ", range =", paste0("(", round(mean_survival$min_survival_days[i], 2),
                          "-", round(mean_survival$max_survival_days[i], 2), ")"),
      "\n")
}
cat("=========================================================\n\n")

print(mean_survival)

# ===== Annotation placement (DOWN + LEFT) =====
max_x <- max(metadata$Days_Alive, na.rm = TRUE)
left_x <- max_x * 0.015
y_text_annot <- 0.06

# ===== Plot =====
ggsurv <- ggsurvplot(
  fit,
  data = metadata,
  conf.int   = FALSE,
  risk.table = FALSE,
  pval       = FALSE,
  palette    = km_colors,
  legend.labs = legend_labels,
  legend.title = "",
  xlab = "Days Alive",
  ylab = "Survival Probability",
  title = "FOP Survival: FOP vs FOP GMA",
  ggtheme = theme_classic(base_size = 16)
)

ggsurv$plot <- ggsurv$plot +
  theme(
    axis.line    = element_line(linewidth = linesize, color = "black"),
    axis.title.x = element_text(size = yaxis, color = "black"),
    axis.title.y = element_text(size = yaxis, color = "black"),
    plot.title   = element_text(size = yaxis + 2, hjust = 0.5, color = "black"),
    axis.text.x  = element_text(size = xaxis, color = mmp9_grey),
    axis.text.y  = element_text(size = ytext, color = mmp9_grey),
    legend.text  = element_text(size = xaxis, color = mmp9_grey)
  ) +
  annotate(
    "text",
    x = left_x,
    y = y_text_annot,
    label = p_text,
    hjust = 0,
    color = "black",
    size = annot_size
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.6)))

print(ggsurv)

# ===== Save square outputs to figures folder =====
tiff_file <- here("figures", "KM_Curve_GMA.tiff")
pdf_file  <- here("figures", "KM_Curve_GMA.pdf")

# TIFF
tiff(
  filename = tiff_file,
  width = 7.5,
  height = 7.5,
  units = "in",
  res = 600,
  compression = "lzw"
)
print(ggsurv$plot)
dev.off()

# PDF
pdf(
  file = pdf_file,
  width = 7.5,
  height = 7.5,
  useDingbats = FALSE
)
print(ggsurv$plot)
dev.off()
