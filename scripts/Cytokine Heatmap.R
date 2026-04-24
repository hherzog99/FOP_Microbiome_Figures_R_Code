# ==========================================
# Serum cytokines heatmap (log2 fold-change) 
# ==========================================

# --------------------------
# Install packages (run once)
# --------------------------
install.packages("readxl")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ComplexHeatmap")
install.packages("circlize")
install.packages("writexl")
install.packages("here")

# If ComplexHeatmap does not install from CRAN, use:
# install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")

# --------------------------
# Load libraries
# --------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(writexl)
library(here)

# --------------------------
# Config
# --------------------------
xlsx_path <- here("data", "Cytokine.xlsx")
sheet        <- NULL
eps          <- 1e-6
clean_labels <- TRUE
lim          <- 2
alpha_star   <- 0.059

# Output directory
fig_dir <- here("figures")

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

# Output files
png_file   <- file.path(fig_dir, "Serum_Cytokines_log2FC_heatmap.png")
pdf_file   <- file.path(fig_dir, "Serum_Cytokines_log2FC_heatmap.pdf")
excel_file <- file.path(fig_dir, "Cytokine_pvalues_Wilcoxon.xlsx")

# --------------------------
# Helpers
# --------------------------
fixlab <- function(x) {
  x <- gsub("β","beta", x, fixed = TRUE)
  x <- gsub("α","alpha", x, fixed = TRUE)
  x
}

wrap_cytokine_labels <- function(x, width = 18) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

round_half_up <- function(x, digits = 2) {
  posneg <- sign(x)
  z <- abs(x) * 10^digits
  z <- z + 0.5 + sqrt(.Machine$double.eps)
  z <- trunc(z)
  posneg * z / 10^digits
}

p_to_symbol <- function(p){
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01,  "**",
                       ifelse(p < alpha_star, "*",
                              ifelse(p < 0.10, "·", "")))))
}

# --------------------------
# Load data
# --------------------------
df <- if (is.null(sheet)) read_excel(xlsx_path) else read_excel(xlsx_path, sheet = sheet)
stopifnot(all(c("Group","Sex") %in% names(df)))

df <- df %>%
  mutate(
    Group = factor(trimws(as.character(Group)), levels = c("WT","WT Abx","FOP","FOP Abx")),
    Sex   = factor(trimws(as.character(Sex)),   levels = c("F","M"))
  )

cytokine_cols <- setdiff(names(df), c("Group","Sex"))
df[cytokine_cols] <- lapply(df[cytokine_cols], function(x) suppressWarnings(as.numeric(x)))

# --------------------------
# Comparisons
# --------------------------
comp_specs <- tibble::tribble(
  ~g1,        ~g2,       ~comp_key,             ~col_label,
  "FOP",      "WT",      "FOP_vs_WT",           "FOP vs WT",
  "WT Abx",   "WT",      "WTGMA_vs_WT",         "WT GMA vs WT",
  "FOP Abx",  "WT",      "FOPGMA_vs_WT",        "FOP GMA vs WT",
  "FOP Abx",  "WT Abx",  "FOPGMA_vs_WTGMA",     "FOP GMA vs WT GMA",
  "FOP Abx",  "FOP",     "FOPGMA_vs_FOP",       "FOP GMA vs FOP"
)

# --------------------------
# Wilcoxon p-values
# --------------------------
test_one_comparison <- function(df, cytokines, g1, g2, comp_key){
  bind_rows(lapply(cytokines, function(cyt){
    x <- df %>% filter(Group == g1) %>% pull(.data[[cyt]])
    y <- df %>% filter(Group == g2) %>% pull(.data[[cyt]])
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    
    p <- if (length(x) < 2 || length(y) < 2) NA_real_
    else suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
    
    data.frame(
      Cytokine = cyt,
      Comparison = comp_key,
      p_value = p,
      stringsAsFactors = FALSE
    )
  }))
}

stats_all <- bind_rows(lapply(seq_len(nrow(comp_specs)), function(i){
  cs <- comp_specs[i, ]
  test_one_comparison(df, cytokine_cols, cs$g1, cs$g2, cs$comp_key)
}))

if (clean_labels) stats_all$Cytokine <- fixlab(stats_all$Cytokine)

# Round p-values to 2 decimals using half-up rounding
stats_all <- stats_all %>%
  mutate(p_value = round_half_up(p_value, 2))

# --------------------------
# Save Excel
# --------------------------
write_xlsx(
  list(`All_Comparisons` = stats_all),
  path = excel_file
)

# --------------------------
# log2FC matrix
# --------------------------
group_means <- df %>%
  group_by(Group) %>%
  summarise(across(all_of(cytokine_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

M_grp <- t(as.matrix(group_means %>% select(-Group)))
colnames(M_grp) <- as.character(group_means$Group)
rownames(M_grp) <- cytokine_cols

fc_mat <- do.call(cbind, lapply(seq_len(nrow(comp_specs)), function(i){
  cs <- comp_specs[i, ]
  log2((M_grp[, cs$g1] + eps) / (M_grp[, cs$g2] + eps))
}))

colnames(fc_mat) <- comp_specs$col_label
rownames(fc_mat) <- rownames(M_grp)

if (clean_labels) rownames(fc_mat) <- fixlab(rownames(fc_mat))

# Order rows
ord <- order(fc_mat[, "FOP vs WT"], decreasing = TRUE, na.last = TRUE)
fc_mat <- fc_mat[ord, , drop = FALSE]

row_labels <- wrap_cytokine_labels(rownames(fc_mat), width = 18)

# --------------------------
# Significance symbols
# --------------------------
p_wide <- stats_all %>%
  pivot_wider(names_from = Comparison, values_from = p_value)

p_wide2 <- p_wide[match(rownames(fc_mat), p_wide$Cytokine), , drop = FALSE]

label_to_key <- setNames(comp_specs$comp_key, comp_specs$col_label)
symbols_mat <- sapply(colnames(fc_mat), function(lbl){
  p_to_symbol(p_wide2[[label_to_key[[lbl]]]])
})

symbols_mat <- as.matrix(symbols_mat)
rownames(symbols_mat) <- rownames(fc_mat)

# --------------------------
# Color function
# --------------------------
col_fun <- colorRamp2(
  c(-lim, 0, lim),
  c("#2166ac", "white", "#b2182b")
)

# --------------------------
# Heatmap
# --------------------------
ht <- Heatmap(
  fc_mat,
  name = "log2FC",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = NULL,
  
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "Helvetica"),
  column_names_centered = TRUE,
  column_names_max_height = unit(22, "mm"),
  
  row_labels = row_labels,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 14, fontfamily = "Helvetica"),
  row_names_max_width = unit(80, "mm"),
  
  width = unit(115, "mm"),
  show_heatmap_legend = FALSE,
  
  rect_gp = gpar(col = "white", lwd = 0.7),
  
  cell_fun = function(j, i, x, y, w, h, fill){
    lab <- symbols_mat[i, j]
    if (lab != "") {
      grid.text(
        lab, x, y,
        gp = gpar(fontsize = 12, fontface = "bold",
                  col = "black", fontfamily = "Helvetica")
      )
    }
  }
)

# --------------------------
# Legends (RIGHT, vertical)
# --------------------------
lgd_log2fc <- Legend(
  col_fun = col_fun,
  title = "log2FC",
  at = seq(-lim, lim, by = 1),
  direction = "vertical",
  legend_height = unit(45, "mm"),
  title_gp  = gpar(fontsize = 14, fontface = "bold", fontfamily = "Helvetica"),
  labels_gp = gpar(fontsize = 12, fontfamily = "Helvetica")
)

sig_lgd <- Legend(
  title = "Significance",
  labels = c("***  p < 0.001",
             "**   p < 0.01",
             "*    p < 0.05",
             "·    p < 0.1"),
  legend_gp = gpar(col = "black"),
  title_gp  = gpar(fontsize = 14, fontface = "bold", fontfamily = "Helvetica"),
  labels_gp = gpar(fontsize = 12, fontfamily = "Helvetica")
)

lgd_right <- packLegend(
  lgd_log2fc,
  sig_lgd,
  direction = "vertical",
  gap = unit(4, "mm")
)

# --------------------------
# Draw & save
# --------------------------
png(png_file, width = 3200, height = 2400, res = 300)
draw(
  ht,
  heatmap_legend_list = list(lgd_right),
  heatmap_legend_side = "right",
  padding = unit(c(6, 14, 6, 6), "mm")
)
dev.off()

pdf(pdf_file, width = 11.5, height = 9, family = "Helvetica")
draw(
  ht,
  heatmap_legend_list = list(lgd_right),
  heatmap_legend_side = "right",
  padding = unit(c(6, 14, 6, 6), "mm")
)
dev.off()
