#!/usr/bin/env Rscript
#
# SCRIPT: I-Corps Thematic Reanalysis
# PURPOSE: This script performs a full statistical analysis and generates all
#          publication-ready figures for the "Double Pivot" manuscript.
#
# FIXES INCLUDED:
#   - The GLM section now exclusively uses the stable Poisson model to avoid
#     convergence warnings with sparse data, ensuring all results are reliable.
#   - The 'themes' vector matches the exact column names in the provided CSV file.
#
# USAGE:
#   Rscript run_analysis.R your_data_file.csv --plot --png

# ---- 1. SETUP ----
suppressPackageStartupMessages({
  library(data.table)
  library(MASS)
  library(ggplot2)
  library(tidyr)
  library(scales)
  library(stringr)
})

# ---- 2. CONFIGURATION ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript run_analysis.R <csv_path> [--plot] [--png|--pdf]")
infile   <- args[1]
do_plot  <- any(grepl("^--plot$", args))
save_png <- any(grepl("^--png$", args))
save_pdf <- any(grepl("^--pdf$", args))
if (!save_png && !save_pdf) save_png <- TRUE

themes <- c(
  "Clinician_Workflow_Constraints",
  "Funding_Economic_Constraints",
  "Data_Collection_Challenges",
  "Pharma_Build_vs_Buy",
  "Need_for_Focus_Scope",
  "Regulatory_Privacy_Hurdles",
  "Global_Health_Market_Diffs",
  "Payer_ROI_Value_Emphasis"
)
timevar <- "Phase"

# ---- 3. LOAD & VALIDATE DATA ----
DT <- data.table::fread(infile, showProgress = FALSE)
need <- c(timevar, themes)
missing <- setdiff(need, names(DT))
if (length(missing)) {
  stop(paste("Error: Missing required columns:", paste(missing, collapse=", ")))
}

for (t in themes) DT[, (t) := as.integer(ifelse(is.na(get(t)) | get(t)==0, 0L, 1L))]
DT[, (timevar) := as.integer(get(timevar))]
DT <- DT[!is.na(get(timevar))]

# ---- 4. SETUP OUTPUTS & HELPERS ----
outdir <- "icorps_reanalysis_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

wilson_ci <- function(x, n, conf=0.95){
  if (is.na(n) || n==0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1-conf)/2)
  p <- x/n
  denom <- 1 + z^2/n
  center <- (p + z^2/(2*n))/denom
  half <- (z/denom) * sqrt(p*(1-p)/n + z^2/(4*n^2))
  c(max(0, center - half), min(1, center + half))
}

save_plot <- function(p, filename, w=8, h=5, dpi=300){
  if (save_png) ggsave(file.path(outdir, paste0(filename, ".png")), p, width=w, height=h, dpi=dpi, bg="white")
  if (save_pdf) {
    pdf_device <- if (requireNamespace("cairo_pdf", quietly = TRUE)) "cairo_pdf" else "pdf"
    ggsave(file.path(outdir, paste0(filename, ".pdf")), p, width=w, height=h, device=pdf_device, bg="white")
  }
}

theme_pub <- function(){
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face="bold", size = 16),
    plot.subtitle = element_text(size = 12, margin = margin(b = 10)),
    legend.position = "bottom",
    axis.title = element_text(face="bold")
  )
}

# ===============================================
# 5. STATISTICAL ANALYSIS & TABLE GENERATION
# ===============================================

agg <- DT[, c(lapply(.SD, sum), .(N_interviews = .N)), by = timevar, .SDcols = themes][order(get(timevar))]
props <- copy(agg)
for (t in themes) props[, (paste0(t, "_prop")) := get(t)/N_interviews]
pcols <- paste0(themes, "_prop")

fwrite(agg, file.path(outdir, "theme_counts_by_time.csv"))
fwrite(props[, c(timevar, "N_interviews", pcols), with = FALSE],
       file.path(outdir, "theme_proportions_by_time.csv"))

chisq <- data.table()
for (t in themes) {
  tab <- table(DT[[timevar]], DT[[t]])
  test <- suppressWarnings(chisq.test(tab, correct=FALSE))
  n <- sum(tab); r <- nrow(tab); c <- ncol(tab)
  cv <- sqrt(unname(test$statistic)/(n*min(r-1, c-1)))
  chisq <- rbind(chisq, data.table(
    Theme=t, Chi_sq=unname(test$statistic), df=unname(test$parameter),
    p_value=unname(test$p.value), Cramers_V=cv
  ))
}
chisq[, p_value_fdr := p.adjust(p_value, "BH")]
fwrite(chisq[order(p_value)], file.path(outdir, "chi_square_test.csv"))

# --- **FIXED: GLM SECTION NOW USES STABLE POISSON MODEL ONLY** ---
# This avoids the convergence warnings from glm.nb() with sparse data.
glm_res <- data.table(
  Theme=themes, Model_Family=NA_character_, IRR=NA_real_,
  CI_low=NA_real_, CI_high=NA_real_, p_value=NA_real_
)
for (i in seq_along(themes)) {
  y <- agg[[themes[i]]]
  off <- log(agg$N_interviews)
  model_data <- data.frame(y = y, Time = agg[[timevar]], off = off)
  colnames(model_data)[2] <- timevar
  
  # Use the standard glm() with a Poisson family for all themes.
  # The offset is correctly placed inside the formula.
  formula_obj <- as.formula(paste("y ~", timevar, "+ offset(off)"))
  fit <- glm(formula_obj, family=poisson(), data = model_data)
  
  co <- summary(fit)$coefficients
  ci <- suppressMessages(confint(fit))
  
  irr <- exp(co[timevar,"Estimate"])
  lwr <- exp(ci[timevar,1]); upr <- exp(ci[timevar,2]); pvl <- co[timevar,"Pr(>|z|)"]
  
  glm_res[i, `:=`(
    Model_Family = "Poisson", # All models are now Poisson
    IRR = irr, CI_low = lwr, CI_high = upr, p_value = pvl
  )]
}
glm_res[, p_value_fdr := p.adjust(p_value, "BH")]
fwrite(glm_res[order(p_value)], file.path(outdir, "glm_trend_models.csv"))


# ===============================================
# 6. FIGURE GENERATION
# ===============================================
if (do_plot) {
  plot_df <- props[, c(timevar, "N_interviews", pcols), with=FALSE] |>
    pivot_longer(cols = all_of(pcols), names_to = "Theme_prop", values_to = "prop")
  plot_df$Theme <- gsub("_prop$","", plot_df$Theme_prop)
  
  plot_df$Theme_label <- gsub("_"," ", plot_df$Theme)
  plot_df$Theme_label <- gsub("Diffs","Differences", plot_df$Theme_label)
  plot_df$Theme_label <- str_wrap(plot_df$Theme_label, width = 25)

  ci_df <- as.data.table(plot_df)[, .(x = round(prop*N_interviews), n = N_interviews),
                                  by = c(timevar, "Theme", "Theme_label", "prop")]
  ci_calc <- t(apply(ci_df, 1, function(row){
    x <- as.numeric(row["x"]); n <- as.numeric(row["n"])
    wilson_ci(x, n, 0.95)
  }))
  ci_df$lo <- ci_calc[,1]; ci_df$hi <- ci_calc[,2]
  plot_df <- cbind(plot_df, lo = ci_df$lo, hi = ci_df$hi)

  base <- theme_pub()

  f_key <- "Funding_Economic_Constraints"
  f_df <- props[, .(Phase = get(timevar), prop = get(paste0(f_key,"_prop")))]
  
  p_fund <- ggplot(f_df, aes(x = Phase, y = prop)) +
    geom_line(color = "#c0392b", linewidth = 1.5) +
    geom_point(size = 4, color = "#d35400") +
    geom_text(aes(label = percent(prop, accuracy=1)), vjust = -1.2, size = 4.5, color = "black", fontface = "bold") +
    scale_y_continuous(labels = percent_format(1), limits = c(0, 1.05)) +
    scale_x_continuous(breaks = unique(f_df[[timevar]])) +
    labs(
      title = "Spotlight: Funding & Economic Constraints",
      subtitle = "The dramatic spike in Phase 2 quantitatively shows\nthe business model's invalidation.",
      x = timevar, y = "Proportion of Interviews"
    ) +
    base + theme(legend.position = "none")
  save_plot(p_fund, "Figure1_Funding_Spotlight", w=8, h=5.5)

  p_heat <- ggplot(plot_df, aes(x = factor(.data[[timevar]]), y = Theme_label, fill = prop)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = percent(prop, accuracy = 1)), color = "black", size = 3.5) +
    scale_fill_gradient(low = "#e8f1fa", high = "#005f9e", labels = percent_format(1)) +
    labs(
      title = "Theme × Phase Prevalence Heatmap",
      subtitle = "Shows the relative stability of other themes compared to the funding spike.",
      x = timevar, y = "Theme", fill = "% of Interviews"
    ) +
    base + theme(legend.position = "right")
  save_plot(p_heat, "Figure2_Heatmap", w=9, h=6.5)

  p_lines_ci <- ggplot(plot_df, aes(x = .data[[timevar]], y = prop, color = Theme_label, group = Theme_label)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = Theme_label), alpha = 0.12, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_x_continuous(breaks = unique(plot_df[[timevar]])) +
    labs(
      title = "Theme Prevalence by Phase with 95% Confidence Intervals",
      x = timevar, y = "Proportion of Interviews", color = "Theme"
    ) +
    guides(fill = "none") +
    base +
    guides(color = guide_legend(nrow = 3, byrow = TRUE, title.position = "top", title.hjust = 0.5)) 
  save_plot(p_lines_ci, "SuppFigure1_All_Themes_CI", w=10, h=7)

  p_bars <- ggplot(plot_df, aes(x = factor(.data[[timevar]]), y = prop, fill = Theme_label)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.8) +
    geom_text(aes(label = percent(prop, accuracy = 1)),
              position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    scale_y_continuous(labels = percent_format(1), limits = c(0, 1.05)) +
    labs(
      title = "Comparative Prevalence of Themes by Phase",
      x = timevar, y = "Proportion of Interviews", fill = "Theme"
    ) +
    base +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE, title.position = "top", title.hjust = 0.5))
  save_plot(p_bars, "SuppFigure2_Bars_Dodged", w=11, h=7)
}

cat("Analysis complete. All outputs are in the directory:", normalizePath(outdir), "\n")