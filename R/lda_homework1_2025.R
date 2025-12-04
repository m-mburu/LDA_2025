## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE,
                        message = FALSE,
                        warning = FALSE,
                        fig.width = 7,
                        fig.height = 4.3)

library(kableExtra)

use_latex_packages()


#list.files(lda_data_path)

#' Tidy model with 95% CIs and a formatted estimate column
#'
#' Works for lm/glm via broom and lme4/glmmTMB via broom.mixed.
#' @param model a fitted model (lm, glm, lmerMod, glmerMod, glmmTMB, etc.)
#' @param effects which effects to return
#' for mixed models ("fixed", "ran_pars", "ran_vals")
#' @param conf_level CI level
#' @param exponentiate exponentiate estimates/CIs (useful for logit/Poisson)
#' @param digits rounding for formatted strings
#' @return tibble with columns term, estimate (formatted),
#' std.error, statistic, p.value, conf.low, conf.high
tidy_with_ci <- function(model,
                         effects = c("fixed"),
                         conf_level = 0.95,
                         exponentiate = FALSE,
                         digits = 2) {
  mm_class <- class(model)

  is_mixed <-
    any(c(
      "lmerMod", "lmerModLmerTest",
      "glmerMod", "nlmerMod"
    ) %in% mm_class) ||
      inherits(model, "glmmTMB")

  if (is_mixed) {
    if (!requireNamespace("broom.mixed", quietly = TRUE)) {
      stop("broom.mixed is required for mixed models.\n
        Install with install.packages('broom.mixed').")
    }
    out <- broom.mixed::tidy(
      model,
      effects = effects,
      conf.int = TRUE,
      conf.level = conf_level,
      exponentiate = exponentiate
    )
    # keep only common columns
    out <- out[, intersect(
      c(
        "effect", "group", "term", "estimate", "std.error",
        "statistic", "p.value", "conf.low", "conf.high"
      ),
      names(out)
    )]
  } else if (inherits(model, "gls")) {
    # nlme::gls objects are not supported by broom::tidy; extract manually
    s <- summary(model)
    if (!is.null(s$tTable)) {
      tt <- as.data.frame(s$tTable)
      # tTable columns: Value, Std.Error, t-value, p-value
      out <- data.frame(
        term = rownames(tt),
        estimate = tt[["Value"]],
        std.error = tt[["Std.Error"]],
        statistic = tt[["t-value"]],
        p.value = tt[["p-value"]],
        stringsAsFactors = FALSE
      )
      # Wald (normal) CIs
      z <- stats::qnorm(1 - (1 - conf_level) / 2)
      out$conf.low <- out$estimate - z * out$std.error
      out$conf.high <- out$estimate + z * out$std.error
      # ensure column names match later expectations
      out <- out[, c(
        "term", "estimate", "std.error",
        "statistic", "p.value", "conf.low", "conf.high"
      )]
    } else {
      stop("Unable to extract coefficient table from gls summary output.")
    }
  } else {
    if (!requireNamespace("broom", quietly = TRUE)) {
      stop("broom is required. Install with install.packages('broom').")
    }
    out <- broom::tidy(
      model,
      conf.int = TRUE,
      conf.level = conf_level,
      exponentiate = exponentiate
    )
    out <- out[, intersect(
      c(
        "term", "estimate", "std.error",
        "statistic", "p.value", "conf.low", "conf.high"
      ),
      names(out)
    )]
  }

  # format numeric columns
  num_cols <- intersect(c("estimate", "std.error", "statistic", "
  p.value", "conf.low", "conf.high"), names(out))
  out[num_cols] <- lapply(out[num_cols], function(x) {
    if (is.numeric(x)) round(x, digits) else x
  })

  # formatted estimate with CI
  if (all(c("estimate", "conf.low", "conf.high") %in%
    names(out))) {
    # Use sprintf to format numbers to a consistent width for alignment
    # make 3f to be digitsf
    out$estimate_str <- sprintf(paste0("%7.", digits, "f"), out$estimate)
    out$conf.low_str <- sprintf(paste0("%7.", digits, "f"), out$conf.low)
    out$conf.high_str <- sprintf(paste0("%7.", digits, "f"), out$conf.high)

    est_ci <- paste0(
      out$estimate_str, " [",
      out$conf.low_str, ", ", out$conf.high_str, "]"
    )
    out$estimate <- est_ci
  }

  # optional: pretty p-values (e.g., <0.001)
  if ("p.value" %in% names(out)) {
    out$p.value <- ifelse(is.na(out$p.value), NA_character_,
      ifelse(out$p.value < 10^(-digits),
        paste0("<", format(10^(-digits), scientific = FALSE)),
        format(round(out$p.value, digits), nsmall = digits)
      )
    )
  }

  # reorder columns: put term first, then estimate string
  keep_order <- c(
    "effect", "group", "term", "estimate", "std.error",
    "statistic", "p.value"
  )
  out <- out[, intersect(keep_order, names(out))]

  tibble::as_tibble(out)
}


# Tidy adjacent model comparisons from base anova() output
# Works for anova(m1, m2, ...): lm/glm, lmer/glmer, gls/lme, etc.
tidy_anova_compare <- function(anova_df, digits = 2) {
  x <- anova_df
  is_num <- vapply(x, is.numeric, logical(1))
  x[is_num] <- lapply(x[is_num], function(v) {
    ifelse(is.na(v), v, round(v, digits))
  })
  x$call <- NULL
  x
}



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(here)
library(tidyverse)
library(data.table)
library(haven)
library(janitor)
library(knitr)
library(lme4)   
library(lmerTest)
library(pbkrtest)
library(gridExtra)

lda_path <- here()
lda_data_path <- here(lda_path, "data")
#list.files(lda_data_path)



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load the data
alzheimer_data <- read_sas(here(
    lda_data_path,
    "alzheimer25.sas7bdat"
)) |> clean_names()




## ----factor-conversion-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
setDT(alzheimer_data)
# log bmi
#alzheimer_data[, bmi := log(bmi)]
factor_info <- list(
  sex = list(
    levels = c(0, 1),
    labels = c("Male", "Female")
  ),
  edu = list(
    levels = c(1, 2, 3, 4),
    labels = c("Primary", "Lower secondary", "Upper secondary", "Higher")
  ),
  job = list(
    levels = c(0, 1),
    labels = c("No Paid Job", "Paid Job")
  ),
  wzc = list(
    levels = c(0, 1),
    labels = c("Lives at Home", "Nursing Home")
  )
)


for (col_name in names(factor_info)) {
  info <- factor_info[[col_name]]
  alzheimer_data[, (col_name) := factor(get(col_name),
    levels = info$levels,
    labels = info$labels
  )]
}

# Verify the changes by checking the structure of the new columns
#kable(alzheimer_data[1:5, c("sex", "edu", "job", "wzc")])


## ----reshape-data----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

id_vars <- c("patid", "trial", "sex", "age", 
             "edu", "bmi", "inkomen", "job", "adl", "wzc",
             "cdrsb0",  "abpet0", "taupet0")



alzheimer_long <- melt(alzheimer_data,
  id.vars = id_vars,
  measure.vars = patterns(
    cdrsb = "^cdrsb[1-6]",
    bprs = "^bprs[0-6]",
    abpet = "^abpet[1-6]",
    taupet = "^taupet[1-6]"
  ),
  variable.name = "time",
  value.name = c("cdrsb", "bprs", "abpet", "taupet"),
  variable.factor = FALSE
)

alzheimer_long[, time := as.integer(time)]

#alzheimer_long[, time := time - 1]


# Display the first few rows of the new long-format data
#kable(head(alzheimer_long))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(table1)

# Use the original wide-format data for a true baseline summary
baseline_data <- copy(alzheimer_long) 

# Define a named list of labels for clarity
var_labels <- list(
  age    = "Age (years)",
  sex    = "Sex",
  edu    = "Education Level",
  bmi    = "Body Mass Index (kg/m²)",
  job    = "Job Status",
  adl    = "Activities of Daily Living",
  wzc    = "Living Situation",
  bprs  = "Baseline BPRS Score",
  cdrsb0 = "Baseline CDR-SB Score",
  abpet0 = "Baseline Amyloid PET",
  taupet0= "Baseline Tau PET"
)

# Use a loop to apply the labels to the copied data
for (var in names(var_labels)) {
  if (var %in% names(baseline_data)) {
    Hmisc::label(baseline_data[[var]]) <- var_labels[[var]]
  }
}


tab1_summary <- table1(
  ~ bprs + age + sex + edu + bmi + job + adl +
    wzc  + cdrsb0 + abpet0 + taupet0 | time,
  data = baseline_data[!is.na(bprs)],
  caption = "Baseline Patient Characteristics",
  footnote = "Continuous variables are presented as Mean (SD);
        Categorical variables as n (%)."
) |> as.data.frame(row.names = "Variable")

tab1_summary$Overall <- NULL  # Remove overall column if present
# Print the table
# kable(tab1_summary,
#   format = "latex",
#   booktabs = TRUE,
#   #longtable = TRUE,
#   align = "lccccccc"
# ) |>
#   # scale down to fit page width
#   kableExtra::kable_styling(
#     latex_options = c(
      
#       "scale_down"
#     ),
#     font_size = 9
#   )



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-missingness-analysis
#| fig-cap: "Missingness Analysis: Baseline BPRS by Number of Missed Visits (left) and Percentage Missing Over Time (right)"
#| fig.width: 8.5
#| fig.height: 3.5

bprs_cols <- paste0("bprs", 0:6)
alzheimer_data[, n_missing := rowSums(is.na(.SD)), .SDcols = bprs_cols]

# Group by the number of missing visits
# and calculate summary stats for baseline BPRS
dropout_summary <- alzheimer_data[!is.na(bprs0), .(
  mean_bprs0 = mean(bprs0),
  median_bprs0 = median(bprs0),
  visit_n_patients = .N
), by = n_missing][order(n_missing)]



missing_by_time <- alzheimer_long[, .(
  n_patients = uniqueN(patid),
  n_obs = sum(!is.na(bprs))
), by = time][order(time)]

missing_by_time[, pct_missing := round(100 *
                                         (n_patients - n_obs)/
                                         n_patients, 1)]



boxplot_meanbprso0 <- ggplot(alzheimer_data[!is.na(bprs0)], 
                             aes(x = factor(n_missing),
                                 y = bprs0)) +
  geom_boxplot(aes(fill = factor(n_missing)),
               alpha = 0.7) +
  labs(
    title = "Baseline BPRS Score by
    Number of Missed Visits",
    x = "Number of Missed 
    BPRS Measurements",
    y = "Baseline BPRS 
    Score (bprs0)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

missing_by_time_plot <- ggplot(missing_by_time, 
                               aes(x = time, y = pct_missing)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 2) +
  labs(
    title = "Percentage of Missing BPRS \n Measurements Over Time",
    x = "Time (Years)",
    y = "Percentage of Missing BPRS Measurements (%)"
  ) +
  theme_minimal()

grid.arrange(boxplot_meanbprso0,
             missing_by_time_plot,
             ncol = 2
             )


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: loess-smoothing
#  LOESS model

loess_fit <- loess(bprs ~ time,
  data = alzheimer_long[!is.na(bprs)],
  na.action = na.omit
)

# Score the data to get predictions (equivalent to `score data=...`)

bprs_smoothed <- predict(loess_fit, newdata = alzheimer_long)

# Create and sort the output data

loess_results <- alzheimer_long[!is.na(bprs), .(patid, time, bprs, sex)]
loess_results[, predicted_bprs := predict(loess_fit, newdata = .SD)]

# Sort the results by time
setorder(loess_results, time)


# Visualize 
simple_avg_plot <- ggplot(alzheimer_long, aes(x = time, y = bprs)) +
  geom_point(alpha = 0.1, color = "gray") + # Plot raw data points transparently
  geom_line(
    data = loess_results,
    aes(y = predicted_bprs), color = "blue", size = 1
  ) +
  labs(
    title = "LOESS Smoothed Trend of BPRS Over Time",
    x = "Time (Years)",
    y = "BPRS Score"
  ) +
  theme_minimal()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-loess-by-categorical-vars
#| fig-width: 10.3
#| fig-height: 8
#| fig-cap: "LOESS Smoothed BPRS Over Time by Subgroups. Overall average (top),  sex, education, job status, living situation, and clinical center."


library(gridExtra)

subgroup_vars <- c("sex", 
                   "edu", 
                   "job",
                   "wzc", 
                   "trial")

# factor trial

alzheimer_long[, trial := factor(as.character(trial))]

mycolors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A",
  "#984EA3", "#FF7F00", "#FFFF33",
  "#A65628", "#1B9E77", "#D95F02",
  "#7570B3", "#E7298A", "#66A61E",
  "#E6AB02", "#A6761D", "#7FC97F",
  "#BEAED4", "#FDC086", "#FFFF99",
  "#386CB0", "#F0027F", "#BF5B17",
  "#FBB4AE", "#B3CDE3", "#CCEBC5",
  "#DECBE4", "#FED9A6", "#FFFFCC", 
  "#E5D8BD"
)
# Titles for each plot
plot_titles <- c(
  "sex" = "LOESS Smoothed BPRS by Sex",
  "edu" = "LOESS Smoothed BPRS by Education",
  "job" = "LOESS Smoothed BPRS by Job Status",
  "wzc" = "LOESS Smoothed BPRS by Living Situation",
  "trial" = "LOESS Smoothed BPRS by Clinical Center",
  "taupet0" = "LOESS Smoothed BPRS by Baseline Tau PET"
)

#alzheimer_long[, taupet0 := factor(as.character(taupet0))]
# List to store plots
plot_list <- list()

# Loop to create plots
for (var in subgroup_vars) {
  p <- ggplot(
    alzheimer_long[!is.na(bprs)],
    aes_string(x = "time", y = "bprs", color = var)
  ) +
    geom_smooth(method = "loess", se = FALSE) +
    labs(
      title = plot_titles[var],
      x = "Time (Years)",
      y = "BPRS Score"
    ) +
    theme_minimal() +
    scale_color_manual(
      name = "",
      values =mycolors
    ) +
    theme(legend.position = "bottom")


  if (var == "trial") {
    p <- p + guides(color = guide_legend(ncol = 10, byrow = TRUE))
  }
  plot_list[[var]] <- p
}

names(plot_list) <- subgroup_vars

trial_centers_plot <- plot_list[["trial"]]
plot_list[["trial"]] <- NULL
# Arrange the plots in a grid
plot_list[["simple_avg_plot"]] <- simple_avg_plot


# create ordered grob list with overall plot first
glist <- plot_list[c(
  "simple_avg_plot",
  setdiff(
    names(plot_list),
    "simple_avg_plot"
  )
)]
glist <- unname(glist)


ncol <- 2
n_other <- length(glist) - 1
nrows <- 1 + ceiling(n_other / ncol)
layout_mat <- matrix(NA_integer_, nrow = nrows, ncol = ncol)
layout_mat[1, ] <- 1L
idx <- 2L
for (r in 2:nrows) {
  for (c in 1:ncol) {
    if (idx <= length(glist)) {
      layout_mat[r, c] <- idx
      idx <- idx + 1L
    }
  }
}

## add figure 2 caption
bottom_caption <- grid::textGrob(
  "Figure 2: LOESS Smoothed BPRS Over Time by Subgroups.
  Overall average (top),  sex, education, job status, 
  living situation, and clinical center.",
  gp = grid::gpar(fontsize = 10),
  x = 0.5,
  hjust = 0.5
)

grid.arrange(
  grobs = glist,
  layout_matrix = layout_mat
  #bottom = bottom_caption
)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-loess-by-time-cross-section
#| fig-width: 10.3
#| fig-height: 5
#| fig-cap: "LOESS-smoothed BPRS vs continuous predictors by visit year."

library(data.table)
library(ggplot2)

# Select and melt via data.table
continuous_vars <- c("age", "bmi", "cdrsb0", "abpet0")

dt_long_cs <- melt(
  alzheimer_long[!is.na(bprs)],
  id.vars = c("patid", "time", "bprs"),
  measure.vars = continuous_vars,
  variable.name = "predictor",
  value.name = "xval",
  variable.factor = TRUE
)[!is.na(xval)]

# Nice facet labels
label_map <- c(
  age    = "Age (years)",
  bmi    = "Body Mass Index (kg/m²)",
  cdrsb0 = "Baseline CDR-SB",
  abpet0 = "Baseline Amyloid PET"
)
dt_long_cs[, predictor_label :=
  factor(as.character(predictor),
         levels = names(label_map),
         labels = unname(label_map))]

ggplot(dt_long_cs, aes(x = xval, y = bprs, color = factor(time))) +
  geom_smooth(method = "loess", se = FALSE, span = 0.8) +
  # geom_point(alpha = 0.15, size = 0.6) +   # uncomment to show points
  labs(
    title = "BPRS vs Continuous Predictors by Visit Year",
    x = NULL,
    y = "BPRS Score",
    color = "Time (Year)"
  ) +
  scale_color_manual(values = mycolors, guide = guide_legend(ncol = 7)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~ predictor_label, scales = "free_x", ncol = 2)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-trial-centersplot
#| fig-width: 7
#| fig-height: 3.0
#| fig-cap: "LOESS Smoothed BPRS Over Time by Clinical Center."

# use grid extra to display trial centers plot alone
trial_centers_plot <- trial_centers_plot +
    theme(
        legend.position = "none"
    )
grid.arrange(
    trial_centers_plot
)



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: variance-structure

alzheimer_long_full <- alzheimer_long[!is.na(bprs)]
mean_model <- lm(bprs ~ factor(time), data = alzheimer_long_full)


alzheimer_long_full[, resid_sq := residuals(mean_model)^2]


var_function_plot <- ggplot(
  alzheimer_long_full,
  aes(x = time, y = resid_sq)
) +
  geom_point(alpha = 0.3) +
  geom_smooth(
    method = "loess", se = FALSE,
    color = "black", size = 1.2
  ) +
  labs(
    title = "Smoothed Variance Function (BPRS)",
    x = "Time (Years)",
    y = "Squared Residual"
  ) +
  theme_minimal()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-correlation-structure
#| fig-width: 8
#| fig-height: 3.0
#| fig-cap: "Smoothed Variance Function (left) and Semi-variogram (right) of BPRS Residuals."

library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

# Residuals from a linear average trend: bprs ~ time  (no covariates)
dt <- alzheimer_long %>%
  select(patid, time, bprs) %>%
  filter(!is.na(bprs), !is.na(time))

m_lin <- lm(bprs ~ time, data = dt)
dt_res <- dt %>% mutate(r = residuals(m_lin))

#  WITHIN-subject pairs: u = |t_k - t_j|, v = 0.5 * (r_k - r_j)^2
pairs_within <- dt_res %>%
  inner_join(dt_res, by = "patid", suffix = c(".j", ".k")) %>%
  filter(time.k > time.j) %>%
  transmute(
    u = time.k - time.j,               # time lag
    v = 0.5 * (r.k - r.j)^2            # semi-variance point
  )

# Discrete-by-lag estimate and LOESS smooth
vario_by_u <- pairs_within %>%
  group_by(u) %>%
  summarise(v = mean(v, na.rm = TRUE), n = n(), .groups = "drop") %>%
  arrange(u)

lo <- loess(v ~ u, data = pairs_within, span = 0.7)
grid <- tibble(u = seq(min(pairs_within$u),
  max(pairs_within$u),
  length.out = 200
)) %>%
  mutate(vsmooth = as.numeric(predict(lo, newdata = cur_data())))

# 3) BETWEEN-subject pairs for "total variance" reference (≈ σ^2 + τ^2 + ν^2)
#    (SAS used all cross-subject half-squared differences and took the mean)
pairs_between_mean <- {
  # light-weight approx: sample some cross pairs for speed if data is large
  set.seed(1)
  # sample up to 20k cross pairs; increase if dataset is small
  target_n <- 20000

  # draw indices
  idx_j <- sample.int(nrow(dt_res), min(target_n, nrow(dt_res)), replace = TRUE)
  idx_k <- sample.int(nrow(dt_res), min(target_n, nrow(dt_res)), replace = TRUE)

  cross_df <- tibble(
    patid_j = dt_res$patid[idx_j],
    r_j     = dt_res$r[idx_j],
    patid_k = dt_res$patid[idx_k],
    r_k     = dt_res$r[idx_k]
  ) %>% filter(patid_j != patid_k)

  mean(0.5 * (cross_df$r_k - cross_df$r_j)^2, na.rm = TRUE)
}
total_hat <- pairs_between_mean 


# Plot: raw points by lag, LOESS smooth, and total-variance reference line
coor_plot <- ggplot() +
  geom_point(data = vario_by_u, aes(u, v, size = n), alpha = 0.5) +
  geom_line(data = grid, aes(u, vsmooth), linewidth = 1) +
  geom_hline(yintercept = total_hat, linetype = 2) +
  labs(
    title = "Semi-variogram of BPRS residuals, linear trend",
    x = "Time lag u (years)",
    y = "Semivariogram v(u)"
  ) +
  guides(size = "none")+
  theme_minimal()

# use grid arrange for variance and correlation plots
grid.arrange(var_function_plot,
  coor_plot,
  ncol = 2
)



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: individual-profiles
uniqu_ppid <- alzheimer_data[, unique(patid)]
set.seed(200)
random_ppid <- sample(uniqu_ppid, 50)

# Filter long data for these patients
plot_data <- alzheimer_long[patid %in% random_ppid]

# Make patid a factor for cleaner colouring
plot_data[, patid := factor(patid)]

# Expand our 28-colour palette to 50 colours
mycolors50 <- colorRampPalette(mycolors)(50)

spaggheti_plot <- ggplot(
  plot_data,
  aes(
    x = time, y = bprs,
    group = patid, colour = patid
  )
) +
  geom_line(size = 0.5, alpha = 0.9) +
  geom_point(size = 1, alpha = 0.9) +
  scale_colour_manual(values = mycolors50) +
  labs(
    title = "Individual BPRS Profiles \n(50 Random Patients)",
    x = "Time (Years)",
    y = "BPRS"
  ) +
  theme_minimal() +
  theme(legend.position = "none")



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-per-patient-slope
#| fig-width: 8
#| fig-height: 3.5
#| fig-cap: "Individual BPRS Profiles (left) and Distribution of Per-Patient Slopes (right)."

# Calculate per-patient slopes
slopes <- alzheimer_long[!is.na(bprs), .(
    slope = if (.N > 1) {
        coef(lm(bprs ~ time))[2]
    } else {
        NA_real_
    }
), by = patid]
# Plot histogram of slopes
per_patient_slope <- ggplot(slopes, aes(x = slope)) +
    geom_histogram(
        binwidth = 0.5,
        color = "black",
        alpha = 0.7
    ) +
    labs(
        title = "Distribution of Per-Patient BPRS Slopes",
        x = "Slope of BPRS over Time",
        y = "Number of Patients"
    ) +
    theme_minimal()

grid.arrange(
    spaggheti_plot,
    per_patient_slope,
    ncol = 2
)




## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-increment-analysis
#| fig-cap: "Year-to-Year Change in BPRS. Mean increments with standard deviation bars."
#| fig.width: 5.3
#| fig.height: 3.5

# Summary statistics of increments
alz_incr <- alzheimer_long %>%
  arrange(patid, time) %>%
  group_by(patid) %>%
  mutate(increment = bprs - lag(bprs)) %>%
  ungroup()

incr_summary <- alz_incr %>%
  filter(!is.na(increment)) %>%
  group_by(time) %>%
  summarise(
    n = n(),
    mean_increment = mean(increment, na.rm = TRUE),
    sd_increment   = sd(increment, na.rm = TRUE)
  )

#incr_summary

ggplot(incr_summary, aes(time, mean_increment)) +
  geom_line( colour = "#1b9e77") +
  geom_point(size = 1, colour = "#d95f02") +
  geom_errorbar(aes(ymin = mean_increment - sd_increment,
                    ymax = mean_increment + sd_increment),
                width = 0.15, colour = "grey30") +
  theme_minimal() +
  labs(
    title = "Year-to-Year Change in BPRS",
    subtitle = "Mean increments with standard deviation bars",
    x = "Time",
    y = "Mean increment (BPRS)"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40")
  )




## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-separate-analyses
#| tbl-cap: "Separate Linear Models for Each Follow-up Year"

baseline_covariates_str <- "sex + age + bmi + job + adl + wzc +
cdrsb0 + abpet0 + taupet0"

# Split the string into a vector of predictor names
all_predictors <- trimws(strsplit(baseline_covariates_str, "\\+")[[1]])

list_models <- list()

label_map <- c(
  "(Intercept)"         = "Intercept",
  "age"                 = "Age (years)",
  "bmi"                 = "Body Mass Index",
  "adl"                 = "Activities of Daily Living Score",
  "cdrsb0"              = "Baseline CDR-SB Score",
  "abpet0"              = "Baseline Amyloid PET",
  "taupet0"             = "Baseline Tau PET",
  "sexFemale"           = "Sex: Female vs Male",
  "jobPaid Job"         = "Paid Job vs No Paid Job",
  "wzcNursing Home"     = "Lives in Nursing Home vs Home"
)



# Loop through each follow-up year from 0 to 6
for (year in 0:6) {
  bprs_outcome <- paste0("bprs", year)

  # 1. Create the data subset for the current year
  current_data <- alzheimer_data[!is.na(get(bprs_outcome)), ]

  # 2. Identify valid predictors by checking for variance
  valid_predictors <- c()
  for (pred in all_predictors) {
    col_data <- current_data[[pred]]
    is_valid <- FALSE

    if (is.numeric(col_data)) {
      # For numeric: check if standard deviation is greater than 0
      is_valid <- sd(col_data, na.rm = TRUE) > 0
    } else if (is.factor(col_data)) {
      # For factors: check if there is more than one unique level
      is_valid <- uniqueN(col_data, na.rm = TRUE) > 1
    }

    # Add to list if valid, otherwise print a message
    if (is_valid & !is.na(is_valid)) {
      valid_predictors <- c(valid_predictors, pred)
    }
  }

  # 3. Create the new formula string with only the valid predictors
  formula_str <- paste(
    bprs_outcome, "~",
    paste(valid_predictors, collapse = " + ")
  )

  # 4. Fit the model using the cleaned formula
  model <- lm(as.formula(formula_str), data = current_data)

  list_models[[bprs_outcome]] <- model
}



format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
}

years <- 0:6

dt_long <- rbindlist(lapply(years, function(y) {
  outcome <- paste0("bprs", y)
  mod <- list_models[[outcome]]

  sm <- summary(mod)$coef
  ci <- confint(mod)

  dt <- data.table(
    covariate = rownames(sm),
    estimate  = sm[, "Estimate"],
    lower     = ci[, 1],
    upper     = ci[, 2],
    #p         = sm[, "Pr(>|t|)"],
    year      = paste0("Year", y)
  )

  dt[,
     .(covariate,
       year,
       cell = sprintf(
         "%.2f [%.2f, %.2f]",
         estimate, lower, upper
       ))]
}))


# Wide table: one row per covariate, columns = Year1..Year6
summary_dt <- dcast(
  dt_long,
  covariate ~ year,
  value.var = "cell"
)[order(covariate)]

# Map the variable names
summary_dt$covariate_clean <- label_map[ summary_dt$covariate ]

summary_dt$covariate_clean <- ifelse(
  is.na(summary_dt$covariate_clean),
  summary_dt$covariate,   # fallback: keep original if no match in map
  summary_dt$covariate_clean
)

# Replace the 'covariate column' with the clean version
summary_dt <- summary_dt %>%
  mutate(
    covariate = ifelse(
      is.na(covariate_clean),
      covariate,
      covariate_clean
    )
  ) %>%
  select(-covariate_clean)


# Optional nicer name for clinicians
#setnames(summary_dt, "covariate", "Baseline covariate")
kable(summary_dt,
  format = "latex",
  booktabs = TRUE,
  align = "lrrrrrrr", # l=left, r=right alignment for columns
  linesep = ""
) |>
  kable_styling(latex_options = c("hold_position", "scale_down"))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(dplyr)
library(nlme)
library(emmeans)


fm_adj_full <- bprs ~ time + edu + sex + 
  age + bmi + job + adl + wzc +  cdrsb0 +  abpet0 +
  time:sex + time:age + time:bmi + time:job + time:adl + time:wzc +
  time:cdrsb0 + time:abpet0 

gls_full_mean_structure <- gls(fm_adj_full,
  data = alzheimer_long,
  correlation = corSymm(form = ~ time | patid),
  weights = varIdent(form = ~ 1 | time),
  method = "ML", na.action = na.exclude
)

#summary(gls_full_mean_structure)

# tidy_with_ci(gls_full_mean_structure) |>
#   kable(
#     format = "latex",
#     booktabs = TRUE,
#     #longtable = TRUE,
#     align = "lrrrr"
#   ) |>
#   column_spec(2, width = "5cm") |>
#   # add hold position to keep table close to text
#   kable_styling(font_size = 9, position = "center")   |>
#   kable_styling(latex_options = c("hold_position", "scale_down"))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-question-3-model
#| tbl-cap: "Simplfied mean structiure Unstructured Covariance Results"
#| tbl-pos: "H"

fm_adj <- bprs ~ time +  sex + 
  age + bmi + job + adl + wzc +  cdrsb0 + 
   + time:cdrsb0 

#taupet0 + abpet0 
gls_reduced_mean_structure <- gls(fm_adj,
  data = alzheimer_long,
  correlation = corSymm(form = ~ time | patid),
  weights = varIdent(form = ~ 1 | time),
  method = "ML", na.action = na.exclude
)

#summary(gls_reduced_mean_structure )

tidy_with_ci(gls_reduced_mean_structure) |>
  kable(
    format = "latex",
    booktabs = TRUE,
    #longtable = TRUE,
    align = "lrrrr"
  ) |>
  column_spec(2, width = "5cm") |>
  kable_styling(font_size = 9, position = "center")   |>
  kable_styling(latex_options = c("hold_position", "scale_down"))



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-lr-test
#| tbl-cap: "Likelihood Ratio Test between Full and Reduced Mean Structure Models"

lr_test <- anova(
  gls_full_mean_structure,
  gls_reduced_mean_structure
) |> tidy_anova_compare()

# lr_test |>
#   kable(
#     format = "latex",
#     booktabs = TRUE,
#     #longtable = TRUE,
#     align = "lrrrr"
#   ) |>
#   column_spec(2, width = "5cm") |>
#   kable_styling(font_size = 9, position = "left")   |>
#   kable_styling(latex_options = c("hold_position", "scale_down"))



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-reduce-covariance-structure
#| tbl-cap: "Marginal Model with AR(1) Covariance Results"

gls_reduced_cov_structure <- gls(fm_adj,
  data = alzheimer_long,
  correlation = corAR1(form = ~ time | patid),
  weights = varIdent(form = ~ 1 | time),
  method = "ML", na.action = na.exclude
)

# #summary(gls_reduced_cov_structure )

# tidy_with_ci(gls_reduced_cov_structure) |>
#   kable(
#     format = "latex",
#     booktabs = TRUE,
#     #longtable = TRUE,
#     align = "lrrrr"
#   ) |>
#   kable_styling(font_size = 9, position = "left")   |>
#   kable_styling(latex_options = c("hold_position", "scale_down"))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-lr-test-covariance
#| tbl-cap: "Likelihood Ratio Test between Reduced Mean Structure and Reduced Covariance Structure Models"

lr_test_cov <- anova(
  gls_reduced_mean_structure,
  gls_reduced_cov_structure
) |> tidy_anova_compare() 


# lr_test_cov |>
#   kable(
#     format = "latex",
#     booktabs = TRUE,
#     #longtable = TRUE,
#     align = "lrrrr"
#   ) |>
#   column_spec(2, width = "5cm") |>
#   kable_styling(font_size = 9, position = "left")   |>
#   kable_styling(latex_options = c("hold_position", "scale_down"))



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-gls-reduced-residual-diagnostics
#| fig-cap: "GLS Reduced Mean Structure Diagnostics: QQ plot, Residuals vs Fitted, Residuals vs Time."
#| fig-width: 9
#| fig-height: 4.8

library(ggplot2)
library(gridExtra)

# Use Pearson (standardized) residuals
res_df <- data.frame(
  fitted    = fitted(gls_reduced_mean_structure),
  resid_raw = residuals(gls_reduced_mean_structure, type = "response"),
  resid_pearson = residuals(gls_reduced_mean_structure, type = "pearson"),
  time      = alzheimer_long$time
)

# QQ plot of Pearson residuals
p_qq <- ggplot(res_df, aes(sample = resid_pearson)) +
  stat_qq(alpha = 0.6) +
  stat_qq_line(color = "red") +
  labs(title = "QQ Plot",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_minimal()

# Residuals vs Fitted
p_rvf <- ggplot(res_df, aes(x = fitted, y = resid_pearson)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red") +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(title = "Residuals vs Fitted",
       x = "Fitted Values",
       y = "Residuals") +
  theme_minimal()

# Residuals vs Time (scatter + smooth)
p_rvt <- ggplot(res_df, aes(x = time, y = resid_pearson)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red") +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(title = "Residuals vs Time",
       x = "Time (Years)",
       y = "Residuals") +
  theme_minimal()

# Optional: boxplot by time
p_box <- ggplot(res_df, aes(x = factor(time), y = resid_pearson)) +
  geom_boxplot(outlier.alpha = 0.4) +
  geom_hline(yintercept = 0, color = "red") +
  labs(title = "Residuals by Time",
       x = "Time",
       y = "Residuals") +
  theme_minimal()

# Arrange (choose 3 or 4 panels)
grid.arrange(p_qq, p_rvf, p_rvt, p_box, ncol = 2)



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-two-stage-stage2
#| fig-cap: "Individual BPRS Profiles (left) and Distribution of Individual Intercepts and Slopes (right)."
#| fig.width: 8
#| fig.height: 3
#| fig.pos: "H"

# We fit patient with more than one observation
stage1_results <- alzheimer_long[!is.na(bprs),
  if (.N > 1) {
    model <- lm(bprs ~ time)
    .(intercept = coef(model)[1], slope = coef(model)[2])
  },
  by = patid
]


# melt for plotting

stage1_results_melt <- melt(stage1_results,
  id.vars = "patid",
  variable.name = "parameter",
  value.name = "estimate"
)

# distribution of intercepts and slopes
two_stage <- ggplot(stage1_results_melt, aes(x = estimate)) +
  geom_histogram(
    binwidth = 1, color = "black",
    fill = "lightblue", alpha = 0.7
  ) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(
    title = "Distribution of Individual Intercepts and Slopes",
    x = "Estimate",
    y = "Number of Patients"
  ) +
  theme_minimal()

two_stage

## save the plot in  results/
ggsave(here("results/fig-two-stage-stage2.png"),
  plot = two_stage,
  width = 8,
  height = 3
)



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-two-stage-stage2
#| tbl-cap: "Two-Stage Intercept Model Results"
#| tbl-pos: "H"

# Merge the Stage 1 results with the baseline data
stage2_data <- merge(stage1_results,
  alzheimer_data,
  by = "patid"
)

# baseline covariates to test
baseline_covariates <- " sex + age + bmi + job + adl +
wzc + cdrsb0 + abpet0 + taupet0"

# What predicts a patient's starting BPRS level?
intercept_formula <- as.formula(paste(
  "intercept ~",
  baseline_covariates
))
intercept_model <- lm(intercept_formula,
  data = stage2_data
)

tidy_with_ci(intercept_model) |>
  kable(
    format = "latex",
    booktabs = TRUE,
    align = "lrrrr",
    # l=left, r=right alignment for columns
    linesep = ""
  ) |>
  column_spec(2, width = "5cm") |> # Set the width of the 2nd column (estimate)
  kable_styling(latex_options = c("hold_position", "scale_down"))



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-slope-model
#| tbl-cap: "Two-Stage Slope Model Results"
#| tbl-pos: "H"

slope_formula <- as.formula(paste(
  "slope ~",
  baseline_covariates
))
slope_model <- lm(slope_formula,
  data = stage2_data
)

tidy_with_ci(slope_model) |>
  kable(
    format = "latex",
    booktabs = TRUE,
    align = "lrrrr", # l=left, r=right alignment for columns
    linesep = ""
  ) |>
  column_spec(2, width = "5cm") |> # Set the width of the 2nd column (estimate)
  kable_styling(latex_options = c("hold_position", "scale_down"))



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-lmer-reduced-mean-structure
#| tbl-cap: "Mixed Model with Random Intercept and Slope Results"
#| tbl-pos: "H"


library(lme4 )
fm_adj <- bprs ~ time + sex + 
  age + bmi + job + adl + wzc +  cdrsb0 + 
   time:cdrsb0 +
   1 + (time | patid) + (1 | trial)
lmer_reduced_mean_structure <- lmer(fm_adj,
  data = alzheimer_long,
  na.action = na.exclude
)

#summary(lmer_reduced_mean_structure)
lmer_tidy <-tidy_with_ci(lmer_reduced_mean_structure) 

lmer_tidy |>
 kable(
   format = "latex",
   caption = "Mixed Model with Random Intercept
    and Slope Results",
   booktabs = TRUE,
    longtable = TRUE,
   align = "lrrrr"
 ) |>
 column_spec(2, width = "5cm") |>
 kable_styling(font_size = 9, position = "center") |>
  kable_styling(latex_options = c("hold_position", "scale_down"))

# save lmer_tidy as csv in results with row.names()

write.csv(lmer_tidy,
  file = here("results/lmer_reduced_mean_structure_results.csv"),
  row.names = FALSE
)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-residual-diagnostics
#| fig-cap: "Residual Diagnostics: Residuals vs Fitted (left) and Normal Q-Q Plot (right)."
#| fig-pos: "H"
#| fig.width: 9
#| fig.height: 3



resid_re <- resid(lmer_reduced_mean_structure)
fitted_re <- fitted(lmer_reduced_mean_structure)

# ---- ggplot version of Residuals vs Fitted ----
p1_diag_q5 <- ggplot(
    data.frame(fitted = fitted_re, resid = resid_re),
    aes(x = fitted, y = resid)
) +
    geom_point() +
    geom_hline(yintercept = 0, colour = "red") +
    labs(title = "Residuals vs Fitted Values") +
    labs(x = "Fitted values", y = "Residuals") +
    theme_minimal()

# ---- ggplot version of QQ plot ----
p2_diag_q5 <- ggplot(data.frame(resid = resid_re), aes(sample = resid)) +
    stat_qq() +
    stat_qq_line(colour = "red") +
    labs(title = "Normal Q-Q Plot") +
    theme_minimal()

# ---- Arrange with caption ----
grid.arrange(
    p1_diag_q5, p2_diag_q5,
    ncol = 2
)

# save grid
#
png(here("results/fig-residual-diagnostics.png"),
    width = 8, height = 6,
    units = "in", res = 100
)

grid.arrange(
    p1_diag_q5, p2_diag_q5,
    ncol = 2
)

dev.off()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-random-effects-diagnostics
#| fig-cap: "Diagnostics of Random Effects: Q-Q plots for random intercepts (left) and random slopes (middle), and their joint distribution (right)."
#| fig-pos: "H"
#| fig.width: 9
#| fig.height: 4.5


library(ggplot2)
library(gridExtra)
library(grid)

re <- ranef(lmer_reduced_mean_structure)$patid
re_df <- data.frame(
  intercept = re[,"(Intercept)"],
  slope = re[,"time"]
)

p1 <- ggplot(re_df, aes(sample = intercept)) +
  stat_qq() +
  stat_qq_line(colour = "red") +
  labs(title = "Q-Q Plot Random Intercepts",
       x = "Theoretical Quantiles", 
       y = "Sample Quantiles") + theme_minimal()

p2 <- ggplot(re_df, aes(sample = slope)) +
  stat_qq() +
  stat_qq_line(colour = "red") +
  labs(title = "Q-Q Plot Random Slopes",
       x = "Theoretical Quantiles", 
       y = "Sample Quantiles") + theme_minimal()

p3 <- ggplot(re_df, aes(x = intercept, y = slope)) +
  geom_point() +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_vline(xintercept = 0, colour = "grey") +
  labs(title = "Random Intercepts vs Random Slopes",
       x = "Random Intercept",
       y = "Random Slope") + theme_minimal()

grid.arrange(
  p1, p2, p3,
  ncol = 2
)



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-two-stage-vs-mixed
#| fig-cap: "Comparison of Patient-Specific Intercepts and Slopes from the Two-Stage and Mixed-Model methods."
#| fig-subcap: c("Intercepts", "Slopes")
#| fig-pos: "H"
#| fig.width: 8
#| fig.height: 3.5


# Build per-patient newdata at time 0 and 1 (baseline covariates)
base_covs <- alzheimer_data[, .(
  patid, trial, sex,
  age, bmi, job, adl, wzc, cdrsb0, abpet0
)]

nd0 <- copy(base_covs)
nd0[, time := 0L]
nd1 <- copy(base_covs)
nd1[, time := 1L]

# Subject-specific predictions (uses BLUPs since patid/trial exist in model)
pred0 <- predict(lmer_reduced_mean_structure,
  newdata = nd0,
  re.form = NULL,
  allow.new.levels = FALSE
)
pred1 <- predict(lmer_reduced_mean_structure,
  newdata = nd1,
  re.form = NULL,
  allow.new.levels = FALSE
)

# Extract intercept (at time 0) and slope (per year ≈ pred at 1 - pred at 0)
subj_spec <- data.table(
  patid = nd0$patid,
  intercept_mixed = pred0,
  slope_mixed = pred1 - pred0
)

# also show random-effects-only deviations (BLUPs)
re_pat <- data.table(
  patid = rownames(ranef(lmer_reduced_mean_structure)$patid),
  ranef(lmer_reduced_mean_structure)$patid
)

setnames(
  re_pat, c("patid", "(Intercept)", "time"),
  c("patid", "b0i", "b1i")
)

# Compare to two-stage estimates (Question 4)
comp <- merge(stage1_results,
  subj_spec,
  by = "patid",
  all.x = TRUE
)

# Correlations
corr_int <- cor(comp$intercept,
  comp$intercept_mixed,
  use = "complete.obs"
)
corr_slope <- cor(comp$slope,
  comp$slope_mixed,
  use = "complete.obs"
)

# Plots
p_int <- ggplot(comp, aes(intercept, intercept_mixed)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(
    title = sprintf("Intercepts: Two-stage vs Mixed (r = %.2f)", corr_int),
    x = "Two-stage intercept", y = "Mixed-model intercept"
  ) +
  theme_minimal()

p_slope <- ggplot(comp, aes(slope, slope_mixed)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(
    title = sprintf("Slopes: Two-stage vs Mixed (r = %.2f)", corr_slope),
    x = "Two-stage slope", y = "Mixed-model slope"
  ) +
  theme_minimal()

gridExtra::grid.arrange(
  p_int,
  p_slope,
  ncol = 2
)




## ----eval=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# summary(lmer_reduced_mean_structure)


## ----eval=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# library(lme4)
# library(broom.mixed)
# library(knitr)
# library(kableExtra)
# 
# #tidy_with_ci(lmer_reduced_mean_structure) |>
# #  kable(
#  #   align = "lrrrr"
# #  ) |>
#  # column_spec(2, width = "5cm") |>
#  # kable_styling(font_size = 10, position = "center") |>
#  # kable_styling(latex_options = c("hold_position", "scale_down"))
# model_table <- broom.mixed::tidy(lmer_reduced_mean_structure, effects = "fixed")
# 
# kable(model_table, digits = 3, caption = "Table 4: Fixed Effects Estimates from Mixed Model") %>%
#   kable_styling(full_width = FALSE)
# 
# 
# # Extract variance–covariance components
# vc <- VarCorr(lmer_reduced_mean_structure)
# 
# # Convert to a clean data frame
# re_table <- as.data.frame(vc)
# 
# # Clean up and rename columns
# re_table_clean <- re_table[, c("grp", "var1",
#                                "var2", "vcov", "sdcor")]
# names(re_table_clean) <- c("Group", "Term1", "Term2", "Value", "SD_or_Corr")
# 
# kable(
#   re_table_clean,
#   digits = 3,
#   caption = "Table 5: Random-Effects Variance–Covariance Components"
# ) %>%
#   kable_styling(full_width = FALSE)
# 
# 


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-random-effects-variance-components
#| tbl-cap: "Random Effects: Variance Components from Mixed Model"
#| tbl-pos: "H"

# Helper to tidy random-effects variance/correlation like lmer summary
tidy_re_varcorr <- function(model, digits = 3) {
  vc <- lme4::VarCorr(model)
  rows <- list()

  for (grp in names(vc)) {
    m <- vc[[grp]]
    sd <- attr(m, "stddev")
    trm <- names(sd)
    var <- sd^2

    corMat <- tryCatch(attr(m, "correlation"),
      error = function(e) NULL
    )

    for (k in seq_along(trm)) {
      term <- trm[k]
      corr_val <- NA_real_
      if (!is.null(corMat) && nrow(corMat) > 1 &&
        term != "(Intercept)" &&
        "(Intercept)" %in% colnames(corMat)) {
        corr_val <- as.numeric(corMat[term, "(Intercept)"])
      }
      rows[[length(rows) + 1L]] <- data.frame(
        Groups = grp,
        Name = term,
        Variance = var[k],
        Std.Dev. = sd[k],
        Corr = corr_val,
        check.names = FALSE
      )
    }
  }

  out <- do.call(rbind, rows)

  # Add Residual variance
  resid_sd <- sigma(model)
  resid_var <- resid_sd^2
  resid_row <- data.frame(
    Groups = "Residual",
    Name = "",
    Variance = resid_var,
    Std.Dev. = resid_sd,
    Corr = NA_real_,
    check.names = FALSE
  )
  out <- rbind(out, resid_row)

  num_cols <- c("Variance", "Std.Dev.", "Corr")
  out[num_cols] <- lapply(
    out[num_cols],
    function(x) {
      ifelse(is.na(x),
        NA, round(x, digits)
      )
    }
  )
  out
}

display_random_effects <- function(
    model,
    caption = "Random effects: variance components",
    digits = 3) {
  df <- tidy_re_varcorr(model, digits = digits)
  knitr::kable(
    df,
    format = "latex",
    booktabs = TRUE,
    #caption = caption,
    align = c("l", "l", "r", "r", "r")
  ) |>
    kableExtra::kable_styling(latex_options = "hold_position")
}

# Display variance components for the fitted mixed model
display_random_effects(lmer_reduced_mean_structure,
  caption = "Random effects:
                        variance components (lmer)",
  digits = 3
)

var_comp <- tidy_re_varcorr(lmer_reduced_mean_structure,
 digits = 3)

 ## save csv 
write.csv(var_comp,
  file = here("results/lmer_random_effects_variance_components.csv"),
  row.names = FALSE)


## ----echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-question-3
#| tbl-cap: "Gls Full Model Results"
tidy_with_ci(gls_full_mean_structure) |>
  kable(
    format = "latex",
    booktabs = TRUE,
    #longtable = TRUE,
    align = "lrrrr"
  ) |>
  column_spec(2, width = "5cm") |>
  # add hold position to keep table close to text
  kable_styling(font_size = 9, position = "center")   |>
  kable_styling(latex_options = c("hold_position", "scale_down"))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-q1-baseline-summary-tables
#| tbl-cap: "Baseline Patient Characteristics"
tab1_summary$Overall <- NULL  
# Print the table
kable(tab1_summary,
  format = "latex",
  booktabs = TRUE,
  #longtable = TRUE,
  align = "lccccccc"
) |>
  # scale down to fit page width
  kableExtra::kable_styling(
    latex_options = c(
      
      "scale_down"
    ),
    font_size = 9
  )


