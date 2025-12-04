## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE,
                        message = FALSE,
                        warning = FALSE,
                        fig.width = 7,
                        fig.height = 4.3)

library(kableExtra)

use_latex_packages()


tidy_with_ci <- function(model,
                         effects = c("fixed", "ran_pars"),
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

  # formatted estimate with CI
  if (all(c("estimate", "conf.low", "conf.high") %in% names(out))) {
    # Identify rows that have valid CIs (i.e., not ran_pars)
    has_ci <- !is.na(out$conf.low)

    # Create the formatted estimate string, applying CI only where valid
    out$estimate <- ifelse(
      has_ci,
      paste0(
        sprintf(paste0("%.", digits, "f"), out$estimate), " [",
        sprintf(paste0("%.", digits, "f"), out$conf.low), ", ",
        sprintf(paste0("%.", digits, "f"), out$conf.high), "]"
      ),
      sprintf(paste0("%.", digits, "f"), out$estimate) # Just format the estimate for ran_pars
    )
  }

  # format numeric columns that are not the main estimate
  num_cols <- intersect(c("std.error", "statistic"), names(out))
  out[num_cols] <- lapply(out[num_cols], function(x) {
    if (is.numeric(x)) round(x, digits) else x
  })

  # optional: pretty p-values (e.g., <0.001)
  if ("p.value" %in% names(out)) {
    # Format p-values, leaving NAs as they are
    out$p.value <- ifelse(is.na(out$p.value), NA_character_,
      ifelse(out$p.value < 10^(-digits),
        paste0("<0.", paste(rep("0", digits - 1), collapse = ""), "1"),
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
             "edu", "bmi", "inkomen", "job", "adl", "wzc")




alzheimer_long <- melt(alzheimer_data,
  id.vars = id_vars,
  measure.vars = patterns(
    cdrsb = "^cdrsb[0-6]",
    bprs = "^bprs[0-6]",
    abpet = "^abpet[0-6]",
    taupet = "^taupet[0-6]"
  ),
  variable.name = "time",
  value.name = c("cdrsb", "bprs", "abpet", "taupet"),
  variable.factor = FALSE
)

alzheimer_long[, time := as.integer(time)]

alzheimer_long[, CDRSBbin := ifelse(cdrsb > 10, 1, 0)]
#alzheimer_long[, time := time - 1]


# Display the first few rows of the new long-format data
#kable(head(alzheimer_long))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-cdrsb-by-time
#| fig.height: 3.5
#| fig.width: 10
#| fig.cap: "Distribution of CDRSB Scores Over Time and by sex"

library(ggplot2)
library(scales)

# 1. Stacked Bar Plot: CDRSB by Time
# Calculate counts and percentages within each time point
plot_data_time <- alzheimer_long[!is.na(cdrsb), .N, by = .(time, CDRSBbin)]
plot_data_time[, pct := N / sum(N), by = time]

cdrsb_time_plot <- ggplot(plot_data_time, aes(
    x = factor(time),
    y = pct, fill = factor(CDRSBbin)
)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    #scale_fill_brewer(palette = "Dark2") +
    labs(
        title = "Distribution of CDRSB Scores Over Time",
        x = "Time (Visit)",
        y = "Percentage",
        fill = "CDRSB Score"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

# 2. Stacked Bar Plot: CDRSB by Site (trial)
# Calculate counts and percentages within each trial
plot_data_site <- alzheimer_long[!is.na(cdrsb), .N,
by = .(trial, CDRSBbin)]
plot_data_site[, pct := N / sum(N), by = trial]

cdrsb_site_plot <- ggplot(plot_data_site, aes(
    x = factor(trial),
    y = pct, fill = factor(CDRSBbin)
)) +
    geom_bar(stat = "identity", position = "fill", width = .5) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_brewer(palette = "Dark2") +
    labs(
        title = "Distribution of CDRSB Scores by Site",
        x = "Site (Trial)",
        y = "Percentage",
        fill = "CDRSB Score"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90))




outcome_by_sex <- alzheimer_long %>%
    group_by(sex, time) %>%
    summarise(p_high = mean(CDRSBbin, na.rm = TRUE))

sex_time <- ggplot(outcome_by_sex, aes(x = time, y = p_high, color = sex)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
        title = "Proportion of Patients with CDRSB > 10 \n by Sex Over Time",
        x = "Time (Visit)",
        y = "Proportion with CDRSB > 10",
        color = "Sex"
    ) +
    guides(color = guide_legend(title = "Sex")) +
    theme_minimal()+
    theme(legend.position = "bottom")

grid.arrange(cdrsb_time_plot,
    sex_time,
    ncol = 2
)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Identify baseline CDRSBbin for each patient at time == 0
baseline_df <- alzheimer_long[time == 1, .(patid, CDRSBbin0 = CDRSBbin)]

# Merge the baseline group back into the long dataset
long_with_base <- merge(
  alzheimer_long,
  baseline_df,
  by = "patid",
  all.x = TRUE
)

# Create baseline severity groups based on CDRSBbin0
long_with_base[, baseline_group := fcase(
  CDRSBbin0 == 1, "High Severity at Baseline (CDRSB > 10)",
  CDRSBbin0 == 0, "Low Severity at Baseline (CDRSB <= 10)",
  default = NA_character_
)]

# Count total N in each group at baseline
baseline_counts <- long_with_base[!is.na(baseline_group),
                                  .(n_baseline = uniqueN(patid)),
                                  by = baseline_group]

# Count retention by group and time (number of unique patients observed)
retention_by_group <- long_with_base[
  !is.na(cdrsb) & !is.na(baseline_group),
  .(n_remaining = uniqueN(patid)),
  by = .(baseline_group, time)
][order(baseline_group, time)]

# Merge and compute retention percentage
retention_by_group <- merge(retention_by_group, baseline_counts, by = "baseline_group")
retention_by_group[, pct_remaining := n_remaining / n_baseline]

# Plot
retention_plot <- ggplot(retention_by_group,
       aes(x = time, y = pct_remaining,
           color = baseline_group, group = baseline_group)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Retention Over Time by Baseline Severity",
    subtitle = "Patients with high baseline severity faster until year visit 6",
    x = "Time (Visit)",
    y = "Percent of Initial Cohort Remaining",
    color = "Baseline Severity"
  ) +
  guides(color = guide_legend(title = NULL)) +
  theme_minimal() +
  theme(legend.position = "bottom")




## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-age-bmi-by-cdrsb
#| fig.cap: "Distribution of Age and BMI by CDRSBbin"
#| fig.height: 4.5
#| fig.width: 10
age_by_crdrsb <- ggplot(
    subset(alzheimer_long, !is.na(CDRSBbin) & !is.na(age)),
    aes(fill = factor(CDRSBbin), y = age, x = factor(time))
) +
    geom_boxplot() +
    labs(
        x = "Time",
        y = "Age",
        title = "Distribution of Age by CDRSBbin"
    ) +
    guides(fill = guide_legend(title = "CDRSBbin")) +
    theme_minimal()+
    theme(legend.position = "bottom")




### bmi by cdrsbbin
bmi_by_crdrsb <- ggplot(
    subset(alzheimer_long, !is.na(CDRSBbin) & !is.na(bmi)),
    aes(fill = factor(CDRSBbin), y = bmi, x = factor(time))
) +
    geom_boxplot() +
    labs(
        x = "Time",
        y = "BMI",
        title = "Distribution of BMI by CDRSBbin"
    ) +
    guides(fill = guide_legend(title = "CDRSBbin")) +
    theme_minimal()+
    theme(legend.position = "bottom")

grid.arrange(age_by_crdrsb,
    bmi_by_crdrsb,
    ncol = 2
)



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-semi-variogram
#| fig.cap: "Empirical Semivariogram of CDRSBbin Over Time and Patient drop out over time by CDRSB Baseline Severity"
#| fig.width: 10
#| fig.height: 3.5

max_time <- max(alzheimer_long$time, na.rm = TRUE)
lags <- 1:max_time

semivar_df <- lapply(lags, function(h) {
  tmp <- alzheimer_long %>%
    group_by(patid) %>%
    arrange(time, .by_group = TRUE) %>%
    mutate(
      y_t  = CDRSBbin,
      y_th = dplyr::lead(CDRSBbin, n = h),
      diff2 = (y_t - y_th)^2
    ) %>%
    filter(!is.na(diff2))  # keep only valid pairs at lag h

  data.frame(
    lag = h,
    semivariance = 0.5 * mean(tmp$diff2)
  )
}) %>%
  bind_rows()

# Plot semivariogram
semi_variogram <- ggplot(semivar_df, aes(x = lag, y = semivariance)) +
  geom_point() +
  geom_line() +
  labs(
    title = "Empirical Semivariogram of CDRSBbin Over Time",
    x = "Time lag (years)",
    y = "Semivariance"
  ) +
  theme_minimal()

grid.arrange(semi_variogram, retention_plot, ncol = 2)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: gee-model
#| tbl-cap: "GEE Model Results (Odds Ratios)"
#| tbl-pos: "H"

library(geepack)
# Ensure data is sorted for GEE
setorder(alzheimer_long, patid, time)

gee_model <- geeglm(
    CDRSBbin ~ time + sex + age + bmi +
     adl + abpet + taupet +
     time:abpet + time:taupet,
    id = patid,
    data = alzheimer_long,
    family = binomial(link = "logit"),
    corstr = "ar1",
    std.err = "san.se" # Robust standard errors
)

#summary(gee_model)

# Tidy output with Odds Ratios
tidy_with_ci(gee_model, exponentiate = TRUE) %>%
    kable(
        format = "latex",
        booktabs = TRUE,
        # longtable = TRUE,
        align = "lcccccc",
    ) |>
    # scale down to fit page width
    kableExtra::kable_styling(
        latex_options = c(
            "scale_down"
        ),
        font_size = 9
    )



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
# center time

glmm_model_orig <- glmer(
  CDRSBbin ~ time + sex + age + bmi + abpet + taupet+
   adl  + time:taupet + 1 + (time | patid),
  data = alzheimer_long,
  family = binomial(link = "logit"),
  nAGQ = 1,
  control = ctrl
)





## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
# center time

glmm_model_ppid <- glmer(
  CDRSBbin ~ time + sex + age + bmi + abpet + taupet+
   adl  + time:taupet +  (1 | patid),
  data = alzheimer_long,
  family = binomial(link = "logit"),
  nAGQ = 1,
  control = ctrl
)





## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-glmm-model-center
#| tbl-cap: "GLMM Model Results (Odds Ratios)"
#| tbl-pos: "H"


alzheimer_long <- alzheimer_long |>
  dplyr::mutate(
    time_c   = time   - mean(time,   na.rm = TRUE),
    age_c    = age    - mean(age,    na.rm = TRUE),
    bmi_c    = bmi    - mean(bmi,    na.rm = TRUE),
    adl_c    = adl    - mean(adl,    na.rm = TRUE),
    abpet_c  = abpet  - mean(abpet,  na.rm = TRUE),
    taupet_c = taupet - mean(taupet, na.rm = TRUE)
  )

ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
glmm_center <- glmer(
  CDRSBbin ~ time_c + sex + age_c + bmi_c + adl_c +
    abpet_c + taupet_c + time_c:taupet_c +
    (time_c | patid),
  data   = alzheimer_long,
  family = binomial(link = "logit"),
  nAGQ   = 1,
  control = ctrl
)
#summary(glmm_center)

tidy_with_ci(glmm_center, exponentiate = TRUE) %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    # longtable = TRUE,
    align = "lcccccc"
  ) |>
  # scale down to fit page width
  kableExtra::kable_styling(
    latex_options = c(
      "scale_down"
    ),
    font_size = 9
  )


## ----eval=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# # 1) Fit a plain logistic regression with the same fixed effects
# glm_raw <- glm(
#   CDRSBbin ~ time + sex + age + bmi + adl +
#     abpet + taupet + time:taupet,
#   data   = alzheimer_long,
#   family = binomial
# )
#


## ----eval=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# X <- model.matrix(glm_raw)[, -1]  # remove intercept
# cor(X)
#


## ----eval=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# library(car)
#
# vif(glm_raw)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-glmm-ranef
#| tbl-cap: "GLMM Random Effects Variance Components"
#| tbl-pos: "H"

tidy_re_varcorr(glmm_center, digits = 4) %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    # longtable = TRUE,
    align = c("l", "l", "r", "r", "r")
  ) |>
  # scale down to fit page width
  kableExtra::kable_styling(
    latex_options = c(
      "scale_down"
    ),
    font_size = 9
  )


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-eb-hist-intercept-slope
#| fig-cap: "Empirical Bayes random effects: intercepts (left) and slopes (right)"
#| fig.width: 9
#| fig.height: 4

library(lme4)
library(data.table)
library(ggplot2)
library(gridExtra)

# Ensure model with (time_c | patid) is fitted above
mm <- glmm_center

re_pat <- ranef(mm)$patid
re_dt  <- as.data.table(re_pat, keep.rownames = "patid")
setnames(re_dt, c("(Intercept)", "time_c"), c("random_intercept", "random_slope"))

p_int <- ggplot(re_dt, aes(x = random_intercept)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(title = "Patient random intercepts", x = "BLUP (log-odds)", y = "Count") +
  theme_minimal()

p_slope <- ggplot(re_dt, aes(x = random_slope)) +
  geom_histogram(bins = 30, fill = "darkorange", color = "white") +
  labs(title = "Patient random slopes", x = "BLUP (per unit time)", y = "Count") +
  theme_minimal()

grid.arrange(p_int, p_slope, ncol = 2)


## ----eval=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: extreme-subjects
#| tbl-cap: "Subjects with Highest Random Intercepts (High Risk)"
#| tbl-pos: "H"

# top_risk <- re_dt %>% arrange(desc(random_intercept))
#
# head(top_risk, 10) %>%
#   kable(
#     format = "latex",
#     booktabs = TRUE,
#     # longtable = TRUE,
#     align = "lcccccc"
#   ) |>
#   # scale down to fit page width
#   kableExtra::kable_styling(
#     latex_options = c(
#       "scale_down"
#     ),
#     font_size = 9
#   )
#
#


## ----eval=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: extreme-subjects-low
#| tbl-cap: "Subjects with Lowest Random Intercepts (Low Risk)"
#| tbl-pos: "H"
# low_risk <- re_dt %>% arrange(random_intercept)
# head(low_risk, 10) %>%
#   kable(
#     format = "latex",
#     booktabs = TRUE,
#     # longtable = TRUE,
#     align = "lcccccc"
#   ) |>
#   # scale down to fit page width
#   kableExtra::kable_styling(
#     latex_options = c(
#       "scale_down"
#     ),
#     font_size = 9
#   )
#


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-glmm-model-orig
#| tbl-cap: "GLMM Model Results on Original Scale (Odds Ratios)"
#| tbl-pos: "H"

tidy_with_ci(glmm_model_orig, exponentiate = TRUE) %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    # longtable = TRUE,
    align = "lcccccc"
  ) |>
  # scale down to fit page width
  kableExtra::kable_styling(
    latex_options = c(
      "scale_down"
    ),
    font_size = 9
  )


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-glmm-model-ppid
#| tbl-cap: "GLMM Model Results on Original Scale (Odds Ratios)"
#| tbl-pos: "H"

tidy_with_ci(glmm_model_ppid, exponentiate = TRUE) %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    # longtable = TRUE,
    align = "lcccccc"
  ) |>
  # scale down to fit page width
  kableExtra::kable_styling(
    latex_options = c(
      "scale_down"
    ),
    font_size = 9
  )


save.image(here("R","LDA2.RData"))

