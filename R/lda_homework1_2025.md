LDA Homework 1 (2025-2026)
================
  Student A (student number)  
Student B (student number)  
Student C (student number)  
Student D (student number)   

  
% replaces ‘color’









  
% optional

  
% optional



  
% for \[H\] (hold_position)

  
% if customizing header/footer

### Question 1

#### Check if Data is Balanced

From the description of the experiment we note that there is dropout
during follow up. Before using any summary methods to check the mean,
variance and correlation structure we first check if the data is
balanced. We see that there is missingness in the data as shown in
figure 1 below. Patients with higher baseline BPRS scores tend to miss
more follow-up visits. There is also an increase of of the percentage of
missingness over time, reaching above 40% by year 6. Because simple by
time averages overweight those who remain and under represent those who
dropped out (patients with higher baseline BPRS scores), they can paint
a biased trend if simple bprs by time averages are used. This pattern
suggests data are not missing completely at random, so we use smoothing
methods that are better suited for unbalanced data

![](lda_homework1_2025_files/figure-gfm/missingness-analysis-1.png)<!-- -->

#### Zero-variance check for baseline Tau PET (taupet0)

From the table, taupet0 still varies through Year 2, but from Year 3
onwards it has zero variance (only one unique value among those with
BPRS observed). That means taupet0 cannot explain any between-patient
differences after Year 2 and coefficients for those years are not
estimable. We should not use taupet0 from year 3 as this variable would
not provide information on growth trajectories from here onwards.

#### Average loess smoothing

#### Categorical Variables

To explore potential subgroup differences and to see how various
variables are associated with the evolution of BPRS, we create LOESS
smoothed plots of BPRS over time stratified by categorical variables:
sex, education, job status, living situation, and clinical center.
Across these variables, BPRS increases over time. Most subgroup curves
are roughly parallel, showing main effects without interactions for
variables like sex, job status, and living situation. The four education
categories overlap almost completely and there is no difference in
trend. This means there is unlikely a significant effect of education on
BPRS trajectory.

![](lda_homework1_2025_files/figure-gfm/loess-by-subgroups-1.png)<!-- -->

#### Continuous Variables

To explore continuous predictors, we plotted LOESS-smoothed BPRS against
age, BMI and baseline CDR-SB (cdrsb0) and abept0 using separate curves
for each follow-up year. Higher age is associated with higher BPRS at
every visit. BMI shows a curved pattern that peaks in the mid-20s,
suggesting that a non-linear term would describe it better than a
straight line. Baseline CDR-SB shows only a weak and slightly curved
association with BPRS. The LOESS plot of BPRS versus ABPET0 by visit
year shows an unstable curved pattern, with peak at year 3 that is not
repeated at other visits, which might be because of high missingness of
this variable in follow up years

![](lda_homework1_2025_files/figure-gfm/loess-by-time-cross-section-1.png)<!-- -->

#### Trial Centers Plot

To check whether clinical center (trial) affects the mean structure, we
plot(figure 4) LOESS-smoothed BPRS over time separately by trial centre.
We see clear separate curves between centers, indicating that trial
should be included as a random effect in the model to account for
between center variability.

![](lda_homework1_2025_files/figure-gfm/trial-centersplot-1.png)<!-- -->

#### Variance & Correlation Structure

The smoothed squared residuals shows a flat residual curve, we dont
observe a heteroscedastic pattern, meaning the spread of patients BPRS’s
readings trajectory around there expected readings doesn’t change a lot.
The semi-variogram also rises with the time gap between measurements,
showing that observations from the same patient are more similar when
close together and become less similar as the gap widens. These pattern
indicates a within-patient correlation that weakens with distance in
time therefore, the models we fit need to allow correlation structure
that decays over time. The plots are shown in figure 5 below.

![](lda_homework1_2025_files/figure-gfm/correlation-structure-1.png)<!-- -->

#### Question 1: Conclusions and modeling implications

- Mean structure: BPRS rises steadily over time. Subgroup curves are
  mostly parallel with level shifts by sex, job, living situation, and
  clear vertical shifts by clinical center; education shows little
  separation.
- For continuous variables age, showed a straight line pattern while BMI
  showed non linear relationship with BPRS.
- Variance structure: The smoothed squared-residual function is roughly
  flat across visits (no visual heteroscedasticity observed). We will
  fit a model with un assumed heterogeneous variance and compare this
  with a constant variance model using AIC.
- Correlation structure: The semivariogram increases with time lag,
  indicating strong within-patient correlation at short lags that decays
  as visits are farther apart.
- Missingness: Dropout increases over time and is more common for
  patients with higher baseline BPRS, so simple time-point averages
  would be biased
- Modeling implications: For population-average inference, we will fit
  repeated measures model with unstructured covariance structure and
  compare with AR(1) correlation
- For subject-specific inference, we will fit a mixed model with patient
  random intercept and slope and include a random intercept for clinical
  center, because there is a clear clinical site differences
- we don’t recommend using taupet 0 as a covariate in longitudinal
  models as it has zero variance from year 3 onwards, for the
  participants who’s BPRS measurement were taken at year 3 and beyond.
  Hence it will not show how taupet0 affects BPRS over the years

### Question 2

#### Individual profiles & per-patient slopes

To visualise heterogeneity in longitudinal evolution, BPRS trajectories
of 40 random patients were plotted. These individual profiles show
variation in baseline symptom severity, together with increasing trends
over time. The profiles also illustrate the unbalanced nature of the
data, as subjects contribute different numbers of observations. These
observations help justify the choice of summary statistics for
longitudinal data.

![](lda_homework1_2025_files/figure-gfm/per-patient-slope-1.png)<!-- -->

#### Separate analyses for each follow-up year

For each follow-up year $t \in \{1,\dots,6\}$, we fit an ordinary linear
model relating the BPRS score at year $t$ to baseline covariates:

$$
\begin{aligned}
  ext{BPRS}_{it} = &\; \alpha_{0t} + \alpha_{1t} \cdot \text{sex}_i + \alpha_{2t} \cdot \text{age}_i + \alpha_{3t} \cdot \text{bmi}_i + \alpha_{4t} \cdot \text{job}_i + \alpha_{5t} \cdot \text{adl}_i \\
& + \alpha_{6t} \cdot \text{wzc}_i + \alpha_{7t} \cdot \text{cdrsb0}_i + \alpha_{8t} \cdot \text{abpet0}_i + \alpha_{9t} \cdot \text{taupet0}_i + \varepsilon_{it} .
\end{aligned}
$$

where, for each fixed year $t$:

- $\alpha_{0t}$ is the expected BPRS at year $t$ for the reference
  patient (Male, No Paid Job, Lives at Home), with continuous covariates
  at 0.
- $\alpha_{1t}$ is the difference at year $t$ for Female vs Male.
- $\alpha_{2t}$ is the change at year $t$ per 1-year increase in age.
- $\alpha_{3t}$ is the change at year $t$ per 1-unit increase in BMI.
- $\alpha_{4t}$ is the difference at year $t$ for Paid Job vs No Paid
  Job.
- $\alpha_{5t}$ is the change at year $t$ per 1-unit increase in ADL.
- $\alpha_{6t}$ is the difference at year $t$ for Nursing Home vs Lives
  at Home.
- $\alpha_{7t}$ is the change at year $t$ per 1-unit increase in
  baseline CDR-SB (cdrsb0).
- $\alpha_{8t}$ is the change at year $t$ per 1-unit increase in
  baseline Amyloid PET (abpet0).
- $\alpha_{9t}$ is the change at year $t$ per 1-unit increase in
  baseline Tau PET (taupet0).
- $\varepsilon_{it}$ is the error term at year $t$.

### Question 3

#### Motivation Repeated-Measures Model

In Question 2, In Question 2, we analysed each year separately and
presented yearly results (Table 2). While this approach is simple, it
has limitations. Running separate models for each year means repeating
the same hypothesis test multiple times, which increases the risk of
Type I error. It also reduces statistical power, because each model uses
only one year of data instead of leveraging all observations together.
Additionally, separate analyses do not address the main research
question: how BPRS scores changes over time and how baseline
characteristics influence that change. They also ignore the fact that
repeated measurements from the same patient are correlated, which can
lead to misleading conclusions.

To overcome these limitations when estimating the population effects, we
use a repeated-measures multivariate model. This approach treats each
patient’s yearly BPRS readings dependent. A normal multivariate linear
model assumes all observations are independent, which is not realistic
here because the same patient BPRS measurements are taken yearly up to
year 6. The repeated measures model correctly handles the correlation
within a person, the uneven follow-up due to dropout, and the changing
variability across visits.This model uses all available data and gives
more reliable averages and standard errors. We start from a saturated
repeated measures model where time, all baseline characteristics, and
their interactions with time are included. This gives us a way to assess
how the variable affect both the starting BPRS level and the change over
time, so that we do not miss important effects when we later simplify
the model. For the repeated measurements within a patient, we allow a
very flexible pattern for how BPRS scores at different visits are
correlated to each other, without forcing a specific correlation pattern
in advance. We do not use this saturated model to report final clinical
results, we use it as a screening model and then simplify the mean
structure.

To estimate population effects, we use a repeated-measures multivariate
model. This approach treats each patient’s yearly BPRS readings
dependent. A normal multivariate linear model assumes all observations
are independent, which is not realistic here because the same patient
BPRS measurements are taken yearly up to year 6. The repeated measures
model correctly handles the correlation within a person, the uneven
follow-up due to dropout, and the changing variability across
visits.This model uses all available data and gives more reliable
averages and standard errors. We start from a saturated repeated
measures model where time, all baseline characteristics, and their
interactions with time are included. This gives us a way to assess how
the variable affect both the starting BPRS level and the change over
time, so that we do not miss important effects when we later simplify
the model. For the repeated measurements within a patient, we allow a
very flexible pattern for how BPRS scores at different visits are
correlated to each other, without forcing a specific correlation pattern
in advance. We do not use this saturated model to report final clinical
results, we use it as a screening model and then simplify the mean
structure.

In question 2 we have yearly results(Table 2) but analysing each year
separately would mean running several models for the same
question/hypothesis. This increases the chance of getting a significant
result just by chance (Type 1 error), because we are not using all
observations together it also leads to lower power. Separate analyses
also do not answer the main question of the study, which is how BPRS
changes over time and how baseline factors relate to that change. They
ignore that measurements from the same patient are correlated, which can
give misleading results. The repeated-measures model avoids these
problems by using all time points at once and by correctly accounting
for the fact that each patient contributes several, correlated
measurements.

$$
\begin{aligned}
E[Y_{ij} | \boldsymbol{X}_i] = & \beta_0 + \beta_1 \cdot \text{time}_{j} + \beta_2 \cdot \text{edu}_i + \beta_3 \cdot \text{sex}_i + \beta_4 \cdot \text{age}_i + \beta_5 \cdot \text{bmi}_i + \beta_6 \cdot \text{job}_i \\
& + \beta_7 \cdot \text{adl}_i + \beta_8 \cdot \text{wzc}_i + \beta_9 \cdot \text{cdrsb0}_i + \beta_{10} \cdot \text{abpet0}_i \\
& + \beta_{11} \cdot (\text{time}_{j} \times \text{sex}_i) + \beta_{12} \cdot (\text{time}_{j} \times \text{age}_i) + \beta_{13} \cdot (\text{time}_{j} \times \text{bmi}_i) \\
& + \beta_{14} \cdot (\text{time}_{j} \times \text{job}_i) + \beta_{15} \cdot (\text{time}_{j} \times \text{adl}_i) + \beta_{16} \cdot (\text{time}_{j} \times \text{wzc}_i) \\
& + \beta_{17} \cdot (\text{time}_{j} \times \text{cdrsb0}_i) + \beta_{18} \cdot (\text{time}_{j} \times \text{abpet0}_i)
\end{aligned}
$$

where:

- $Y_{ij}$ is the BPRS score for patient $i$ at time $j$.
- $\beta_0$ is the intercept: expected BPRS at time 0 for the reference
  patient (Male, Primary education, No Paid Job, Lives at Home), with
  continuous covariates at 0.
- $\beta_1$ is the main effect of time: average yearly change in BPRS
  for the reference patient.
- $\beta_2$ is the (vector of) main effect(s) of education: contrasts
  for Lower secondary, Upper secondary, and Higher vs Primary (written
  compactly as $\beta_2\,\text{edu}_i$).
- $\beta_3$ is the main effect of sex: Female vs Male at time 0.
- $\beta_4$ is the main effect of age: per 1-year increase at baseline.
- $\beta_5$ is the main effect of BMI: per 1-unit increase at baseline.
- $\beta_6$ is the main effect of job: Paid Job vs No Paid Job at time
  0.
- $\beta_7$ is the main effect of ADL: per 1-unit increase at baseline.
- $\beta_8$ is the main effect of living situation (wzc): Nursing Home
  vs Lives at Home at time 0.
- $\beta_9$ is the main effect of baseline CDR-SB (cdrsb0): per 1-unit
  increase at baseline.
- $\beta_{10}$ is the main effect of baseline Amyloid PET (abpet0): per
  1-unit increase at baseline.
- $\beta_{11}$ is the time-by-sex interaction: additional yearly change
  for Female vs Male.
- $\beta_{12}$ is the time-by-age interaction: modification of yearly
  change per 1-year increase in age.
- $\beta_{13}$ is the time-by-BMI interaction: modification of yearly
  change per 1-unit increase in BMI.
- $\beta_{14}$ is the time-by-job interaction: additional yearly change
  for Paid Job vs No Paid Job.
- $\beta_{15}$ is the time-by-ADL interaction: modification of yearly
  change per 1-unit increase in ADL.
- $\beta_{16}$ is the time-by-living situation (wzc) interaction:
  additional yearly change for Nursing Home vs Lives at Home.
- $\beta_{17}$ is the time-by-cdrsb0 interaction: modification of yearly
  change per 1-unit increase in baseline CDR-SB.
- $\beta_{18}$ is the time-by-abpet0 interaction: modification of yearly
  change per 1-unit increase in baseline Amyloid PET.

Notes: Categorical predictors (sex, edu, job, wzc) are encoded as
factors in R. Their coefficients represent contrasts relative to the
reference levels shown above; for education, the single symbol $\beta_2$
denotes a vector of three contrasts (levels 2–4 versus Primary).

The within-patient errors
$\boldsymbol{e}_i = (e_{i1}, \dots, e_{in_i})^T$ are assumed to follow a
multivariate normal distribution,
$\boldsymbol{e}_i \sim N(\boldsymbol{0}, \boldsymbol{\Sigma}_i)$, where
$\boldsymbol{\Sigma}_i$ is a covariance matrix that accounts for
correlation and heteroscedasticity. In this model, we use an
unstructured correlation matrix (`corSymm`) and allow for different
variances at each time point (`varIdent`).

Note: taupet0 terms (main effect and time interaction) were
intentionally excluded from this model due to near-zero variability at
later follow-up waves (see earlier diagnostic section), which prevents
stable estimation of their effects in a marginal mean model.

The complex repeated-measures model that includes time, all baseline
characteristics, and all their interactions with time, while allowing
for unstructured covariance structure. This complex model confirms that
BPRS increases over time and that age, BMI, ADL, job status, and living
in a nursing home are all associated to BPRS levels. However, most of
interactions between time and other covariates are not significant,
which tells us we can safely move to a simpler mean structure without
losing important clinical information.

##### Simplified Mean Structure

For this step we fit a repeated-measures model with a simplified mean
structure but the same unstructured covariance as before. In this model,
time and all main baseline characteristics except education remain in
the mean structure, and only one time interaction (time with baseline
CDR-SB) is kept. The estimated correlation matrix shows a positive
correlation between repeated BPRS scores within a patient, especially
for visits that are closer in time, and the visit-specific standard
deviations increase from baseline to later years, contradicting with our
earlier variance plots showing homoscedastic variance. We will use this
reduced-mean, unstructured covariance model as the starting point when
we compare with simpler covariance structures for instance AR(1)

$$
\begin{aligned}
E[Y_{ij} | \boldsymbol{X}_i] = & \beta_0 + \beta_1 \cdot \text{time}_{j} + \beta_2 \cdot \text{sex}_i + \beta_3 \cdot \text{age}_i + \beta_4 \cdot \text{bmi}_i + \beta_5 \cdot \text{job}_i \\
& + \beta_6 \cdot \text{adl}_i + \beta_7 \cdot \text{wzc}_i + \beta_8 \cdot \text{cdrsb0}_i + \beta_{9} \cdot \text{taupet0}_i + \beta_{10} \cdot \text{abpet0}_i \\
& + \beta_{11} \cdot (\text{time}_{j} \times \text{cdrsb0}_i)
\end{aligned}
$$

where:

- $Y_{ij}$ is the BPRS score for patient $i$ at time $j$.
- $\beta_0$ is the intercept: expected BPRS at time 0 for the reference
  patient (Male, No Paid Job, Lives at Home), with continuous covariates
  at 0.
- $\beta_1$ is the main effect of time: average yearly change in BPRS
  for the reference patient.
- $\beta_2$ is the main effect of sex: Female vs Male at time 0.
- $\beta_3$ is the main effect of age: per 1-year increase at baseline.
- $\beta_4$ is the main effect of BMI: per 1-unit increase at baseline.
- $\beta_5$ is the main effect of job: Paid Job vs No Paid Job at time
  0.
- $\beta_6$ is the main effect of ADL: per 1-unit increase at baseline.
- $\beta_7$ is the main effect of living situation (wzc): Nursing Home
  vs Lives at Home at time 0.
- $\beta_8$ is the main effect of baseline CDR-SB (cdrsb0): per 1-unit
  increase at baseline.
- $\beta_9$ is the main effect of baseline Tau PET (taupet0): per 1-unit
  increase at baseline.
- $\beta_{10}$ is the main effect of baseline Amyloid PET (abpet0): per
  1-unit increase at baseline.
- $\beta_{11}$ is the time-by-cdrsb0 interaction: modification of yearly
  change per 1-unit increase in baseline CDR-SB.

The covariance structure remains the same as in the full model.

##### LR test between full and reduced mean structure models

The likelihood ratio test compares the complex mean structure model
(with all time–covariate interactions) to the reduced model (only time
$\times$ CDR-SB). The test statistic is 13.78 with a p-value of 0.25,
which means there is no evidence that the extra interaction terms
improve the fit. This shows the simpler mean structure explains the data
just as well as the saturated one, while using fewer parameters, so we
will keep the reduced model as our working mean structure.

#### Reduced Covariance Structure

For this step, we keep the same simplified mean structure and only
change how we model the correlation over time. The aim is to check
whether a simpler covariance structure is sufficient. We fit a
repeated-measures model with an AR(1) correlation structure, where each
visit is highly correlated with the previous one and the correlation
weakens as the time gap increases. We still allow variability to differ
by visit. This means BPRS measurements taken close in time tend to be
similar, and later visits show more spread.

When we compare this AR(1) model to the more flexible model with a
unstructured covariance structure, the likelihood ratio test shows
(L.Ratio $\approx 1233$, $p < 0.0001$) and AIC/BIC worsen for AR(1).
This indicates that the data pattern is too complex to be captured by
this AR(1) correlation structure, so for population-average inference we
keep the unstructured covariance from the previous model.

##### Reduced mean structure vs reduced covariance structure

#### Discussion of Question 3 Results

For Question 3, our final repeated-measures model uses a reduced mean
structure with a unstructured covariance. On the mean side, BPRS
increases on average by about 7 points per year. After adjustment,
higher age, higher BMI, worse ADL, job status, and living in a nursing
home are all associated with higher BPRS levels, while women have
slightly lower BPRS than men. The interaction of time and CDR-SB term is
is significant meaning patients who start with better cognition tend to
show a faster increase in BPRS, whereas those who already have worse
cognition at baseline worsen more slowly. This directly answers our main
question at the population level: BPRS worsens over time, and the speed
of worsening is associated with baseline cognitive status.

For the biomarkers, we focus on ABPET0 and CDR-SB. ABPET0 is retained in
the marginal model, but its coefficient is small and the confidence
interval is wide, so its association with BPRS is uncertain. The LOESS
plot of BPRS versus ABPET0 by visit year shows an unstable curved
pattern, with a peak at year 3 that is not repeated at other visits.
This could be due to data missingness. TAUPET0, on the other hand, was
dropped because from year 3 onwards it has no variability, so it cannot
explain differences in BPRS.

The covariance part of the model lets each visit have its own
variability and allows all time points within a patient to be correlated
without imposing a correlation structure. This matches what we saw in
the exploratory plots repeated BPRS measurements from the same patient
are strongly correlated, especially when visits are close in time, and
the spread increases over follow-up. Using this unstructured structure
gives more realistic standard errors for the population averages.

Finally, this marginal model has some limitations compared with a
mixed-effects model. It provides population-average effects but not
patient-specific trajectories, and it does not explicitly model
differences between clinical centres. The LOESS plot by centre shows
centres with roughly parallel but clearly non overlapping curves for
majority of centres, suggesting centre level differences in BPRS levels
over time. In the current model, this extra variation is absorbed into
the residual covariance rather than being separated out. In Question 5
we therefore move to a mixed-effects model with random intercepts and
slopes for patients, and a random intercept for trial, to better
describe individual and centre-specific disease trajectories while still
addressing the main objective of how BPRS changes over time, and how
those changes relate to baseline clinical characteristics, cognition,
and biomarkers.

### Question 4

In the first stage of a two-stage analysis, we fit a separate linear
regression model for each patient to estimate their individual intercept
and slope over time. For each patient $i$, the model is:

$$
Y_{ij} = \beta_{0i} + \beta_{1i} \cdot \text{time}_j + e_{ij}
$$

where:

- $Y_{ij}$ is the BPRS score for patient $i$ at time $j$.
- $\beta_{0i}$ is the patient-specific intercept (estimated baseline
  BPRS).
- $\beta_{1i}$ is the patient-specific slope (estimated rate of change
  in BPRS).
- $e_{ij}$ are the individual error terms, assumed to be independent
  with mean zero.

![](lda_homework1_2025_files/figure-gfm/two-stage-stage1-1.png)<!-- -->

#### Stage 2: Relate Intercepts and Slopes to Covariates

##### Model for the Intercepts

In the second stage, we model the estimated patient-specific intercepts
($\hat{\beta}_{0i}$) as a function of baseline covariates to understand
which factors predict the initial BPRS score. The model is:

$$
\begin{aligned}
\hat{\beta}_{0i} = & \gamma_{00} + \gamma_{01} \cdot \text{sex}_i + \gamma_{02} \cdot \text{age}_i + \gamma_{03} \cdot \text{bmi}_i + \gamma_{04} \cdot \text{job}_i + \gamma_{05} \cdot \text{adl}_i \\
& + \gamma_{06} \cdot \text{wzc}_i + \gamma_{07} \cdot \text{cdrsb0}_i + \gamma_{08} \cdot \text{abpet0}_i + \gamma_{09} \cdot \text{taupet0}_i + \epsilon_{0i}
\end{aligned}
$$

where:

- $\hat{\beta}_{0i}$ is the estimated intercept for patient $i$ from
  Stage 1.
- $\gamma_{00}$ is the average intercept for the reference individual
  (Male, No Paid Job, Lives at Home), with continuous covariates at 0.
- $\gamma_{01}$ is the intercept difference for Female vs Male.
- $\gamma_{02}$ is the intercept change per 1-year increase in age.
- $\gamma_{03}$ is the intercept change per 1-unit increase in BMI.
- $\gamma_{04}$ is the intercept difference for Paid Job vs No Paid Job.
- $\gamma_{05}$ is the intercept change per 1-unit increase in ADL.
- $\gamma_{06}$ is the intercept difference for Nursing Home vs Lives at
  Home.
- $\gamma_{07}$ is the intercept change per 1-unit increase in baseline
  CDR-SB (cdrsb0).
- $\gamma_{08}$ is the intercept change per 1-unit increase in baseline
  Amyloid PET (abpet0).
- $\gamma_{09}$ is the intercept change per 1-unit increase in baseline
  Tau PET (taupet0).
- $\epsilon_{0i}$ is the error term.

#### Model for the Slopes

Similarly, we model the estimated patient-specific slopes
($\hat{\beta}_{1i}$) as a function of baseline covariates to see what
predicts the rate of change in BPRS. The model is:

$$
\begin{aligned}
\hat{\beta}_{1i} = & \gamma_{10} + \gamma_{11} \cdot \text{sex}_i + \gamma_{12} \cdot \text{age}_i + \gamma_{13} \cdot \text{bmi}_i + \gamma_{14} \cdot \text{job}_i + \gamma_{15} \cdot \text{adl}_i \\
& + \gamma_{16} \cdot \text{wzc}_i + \gamma_{17} \cdot \text{cdrsb0}_i + \gamma_{18} \cdot \text{abpet0}_i + \gamma_{19} \cdot \text{taupet0}_i + \epsilon_{1i}
\end{aligned}
$$

where:

- $\hat{\beta}_{1i}$ is the estimated slope for patient $i$ from Stage
  1.
- $\gamma_{10}$ is the average slope for the reference individual (Male,
  No Paid Job, Lives at Home), with continuous covariates at 0.
- $\gamma_{11}$ is the slope difference for Female vs Male.
- $\gamma_{12}$ is the slope change per 1-year increase in age.
- $\gamma_{13}$ is the slope change per 1-unit increase in BMI.
- $\gamma_{14}$ is the slope difference for Paid Job vs No Paid Job.
- $\gamma_{15}$ is the slope change per 1-unit increase in ADL.
- $\gamma_{16}$ is the slope difference for Nursing Home vs Lives at
  Home.
- $\gamma_{17}$ is the slope change per 1-unit increase in baseline
  CDR-SB (cdrsb0).
- $\gamma_{18}$ is the slope change per 1-unit increase in baseline
  Amyloid PET (abpet0).
- $\gamma_{19}$ is the slope change per 1-unit increase in baseline Tau
  PET (taupet0).
- $\epsilon_{1i}$ is the error term.

### Question 5

#### Plausible random effects model

To understand how psychiatric symptoms evolve over time in Alzheimer’s
disease, we developed a statistical model that accounts for repeated
measurements in each patient. This model allows each patient to have
their own starting level of symptoms and their own rate of change over
time. The outcome of interest was the Brief Psychiatric Rating Scale
(BPRS), measured annually over six years.

We model the longitudinal BPRS scores using a linear mixed-effects model
with random effects at both the patient and trial level:

$$
\begin{aligned}
\text{BPRS}_{it} = (\beta_0 + b_{0i} + u_{0j}) + (\beta_1 + b_{1i})\,\text{time}_{it} + \beta_2 \text{sex}_i
+ \beta_3 \text{age}_i + \beta_4 \text{bmi}_i \\
+ \beta_5 \text{job}_i + \beta_6 \text{adl}_i + \beta_7 \text{wzc}_i + \beta_8 \text{cdrsb0}_i + \beta_9 (\text{time}_{it} \cdot \text{cdrsb0}_i) + \epsilon_{it}
\end{aligned}
$$

where:

- $\text{BPRS}_{it}$ is the BPRS score for patient $i$ in trial $j$ at
  time $t$  
- $\beta_0, \dots, \beta_9$ are **fixed (population-level) effects**  
- $b_{0i}$ is the patient-specific random intercept  
- $b_{1i}$ is the **patient-specific random slope** for time  
- $u_{0j}$ is the **trial-specific random intercept**  
- $\epsilon_{it}$ is the **residual error term**

The patient-level random effects are assumed to follow a bivariate
normal distribution:

$$
\begin{pmatrix}
b_{0i} \\
b_{1i}
\end{pmatrix}
\sim
N\!\left(
\begin{pmatrix}
0 \\
0
\end{pmatrix},
\begin{pmatrix}
\sigma_{b0}^{2} & \sigma_{b0,b1} \\
\sigma_{b0,b1} & \sigma_{b1}^{2}
\end{pmatrix}
\right),
$$

and the trial-level random intercept follows:

$$
u_{0j} \sim N(0, \sigma_{\text{trial}}^{2}).
$$

#### Model Motivation and Justification

Here, We chose to fit linear mixed effects model because our dataset
contains repeated BPRS measurements for the same patients over 6 years
and across multiple clinical trials, with clear variability both in
baseline symptom severity and rate of change between patients. An
unstructured covariance matrix was used for the random effects, allowing
the intercept and slope to be correlated. The model included:

1)  Random intercepts and random slopes for each patient, reflecting
    that individuals differ both in their baseline psychiatric symptom
    severity and their rate of change over time.
2)  A random intercept for trial was added to account for systematic
    differences across the 25 clinical trials.
3)  Fixed effects included all clinically relevant baseline covariates
    (age, sex, BMI, job status, living situation, CDR-SB), ensuring
    adjustment for major patient characteristics.
4)  An interaction term between time and CDRSB was included to assess
    whether cognitive impairment at baseline affects the rate of
    psychiatric symptom progression.

#### Estimation Method

The model was fitted using Restricted Maximum Likelihood (REML), which
provides more reliable and less biased estimates of variance components,
especially because our dataset had many clusters (patients,
trials).Because our model includes random intercepts for both patients
and trials, REML is appropriate to ensure stable estimates of these
random effect variances.

#### Assumptions made

1)  Linearity: BPRS changes linearly over the follow-up period.
2)  Normality: Random effects and residual errors are approximately
    normally distributed.
3)  Independence: After accounting for random effects, residuals are
    independent.
4)  Homoscedasticity: The residual errors have constant variance across
    time and subjects.
5)  Missing data (dropout) are at least missing at random (MAR)
    conditional on the covariates and random effects.
6)  Correct Covariance Structure: The specified random-effects structure
    (random intercepts and slopes) appropriately captures within-subject
    and within-trial correlation.

#### Limitations

1)  More complex to interpret for non-statisticians.The presence of
    fixed effects, random effects, and correlation structures makes
    interpretation harder compared to ordinary regression
2)  Results depend on correct specification of random-effects structure.
    Misspecifying random intercepts/slopes or the covariance structure
    can lead to biased estimates or incorrect standard errors.
3)  Convergence issues can arise with complex models. Models with
    multiple random effects, unstructured covariance matrices, or small
    sample sizes may fail to converge or produce unreliable estimates.
4)  Assumptions such as normality of random effects and linearity may
    not fully hold in real data. Although LME models are robust,
    violation of assumptions especially the normality of random effects
    can affect estimates, particularly with small sample sizes.
5)  MAR assumption cannot be fully verified.The model assumes missing
    data are Missing At Random, but this assumption is untestable and
    can lead to biased results if dropout relates to unobserved
    outcomes.

#### Diagnostics: Evaluating Model Appropriateness

To ensure that the linear mixed-effects model is statistically reliable,
several diagnostic checks were performed. These checks help determine
whether the model’s assumptions are met, and whether the results are
reliable.

**Normality of Residuals (Q–Q Plot)** Figure 8 (right plot) presents the
Normal Q–Q plot which compares the distribution of the model residuals
(errors) to what would be expected if they followed a normal
distribution, which is an underlying assumption of linear mixed-effects
models.The residuals aligned closely with the expected diagonal line,
showing two deviations at each tail.The assumption of normally
distributed residuals is largely satisfied. This supports the
reliability of the estimated effects.

![](lda_homework1_2025_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

**Residuals vs. Fitted Values Plot** Figure 8(right) plot assesses
whether residuals are randomly distributed around zero, which reflects
the assumption of constant variability (“homoscedasticity”).Residuals
are spread evenly around zero with no funnel shape or curvature. This
indicates Variance is stable across time, the model does not
systematically under or over predict symptoms and the linearity
assumption is appropriate. The model fits the data well and the
variability is appropriately captured.

**Distribution of Patient-Specific Intercepts and Slopes** To evaluate
whether the random-effects model appropriately captured differences
between patients, we examined the distribution of the estimated
subject-specific intercepts (patients’ baseline BPRS scores) and
subject-specific slopes (their rates of symptom change over time)(Figure
9). The random intercepts(Figure 9: left Q-Q plot) follow the
theoretical normal line very closely, with only minor deviations at the
extreme upper and lower tails. This suggests that the normality
assumption for the random intercepts is well satisfied. The distribution
is symmetric and shows no major outliers, indicating that patient-level
baseline differences in BPRS are reasonably captured by a normal random
intercept.The random slopes(Figure 9:Middle Q-Q plot) also align closely
with the diagonal reference line, again with only slight deviations in
the two tails. This indicates that the random slopes are approximately
normally distributed, supporting the assumption that individual
differences in rate of change over time follow a normal pattern. No
strong violations or heavy-tailed behavior are present. The scatterplot
of random intercepts versus random slopes (Figure 9: Scatterplot) shows
cloud of points with no strong linear relationship. This implies that
patients’ baseline BPRS levels are not strongly associated with their
rate of change over time. In other words, starting severity does not
systematically predict whether a patient improves faster or slower. The
spread is balanced around zero for both intercepts and slopes,
confirming that the random-effects structure is behaving as expected.

![](lda_homework1_2025_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

#### Comparison With the Two-Stage Analysis

To better understand how each patient’s symptoms (BPRS scores) changed
over time, we estimated individual starting points (intercepts) and
rates of change (slopes) using two different approaches. First, in
Question 4, we applied a two-stage method, fitting a separate linear
regression for each patient to obtain individual intercepts and slopes.
Second, we used a random-effects mixed model, which estimates all
patient-specific intercepts and slopes simultaneously. We then compared
the two methods to see whether the advanced mixed model provides
meaningful improvement.

To compare the two-stage approach with the mixed-effects model, we
examined patient-specific intercepts and slopes estimated by each
method. The intercepts showed very strong agreement between the two
approaches (Figure 10, left), with a correlation of r = 0.94 and points
closely following the 45-degree line. This indicates that baseline BPRS
levels are estimated consistently in both methods. In contrast, the
agreement for slopes was moderate (Figure 10, right; r = 0.77), and the
two-stage estimates displayed greater variability and several extreme
values. This pattern reflects the inherent instability of estimating
individual change rates from limited per-patient data. The mixed-effects
model, by using all availabledata simultaneously, shrinks extreme slopes
toward the population mean and produces more stable and reliable
estimates. Overall, the mixed model provides a clearer and less noisy
picture of individual symptom trajectories, particularly for the rate of
change over time.

Overall, the two methods give broadly similar patient-specific results,
but the mixed model provides a clearer and more stable representation of
individual symptom trajectories, particularly for the rate of change
over time

**Subject-specific intercepts and slopes from mixed model vs two-stage
estimates**
![](lda_homework1_2025_files/figure-gfm/subject-specific-int-slope-1.png)<!-- -->

**Why are the mixed-model estimates still preferred?**

Compared to the two-stage approach, mixed-effects models are preferred
because they correctly account for the correlation among repeated
measurements within each patient, whereas the two-stage method treats
each patient’s data independently. Mixed models also allow a flexible
covariance structure for intercepts and slopes, leading to more accurate
uncertainty estimates. In contrast, the two-stage approach produces
unstable slopes when patients have few observations and cannot model
within-patient correlation. Overall, the mixed model provides more
reliable and statistically valid patient-specific estimates.

#### Results

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: fm_adj
    ##    Data: alzheimer_long
    ## 
    ## REML criterion at convergence: 27121.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.6390 -0.5982 -0.0086  0.5904  3.3222 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev. Corr 
    ##  patid    (Intercept) 0.8067   0.8982        
    ##           time        0.6863   0.8284   -0.61
    ##  trial    (Intercept) 8.0037   2.8291        
    ##  Residual             2.6980   1.6426        
    ## Number of obs: 6220, groups:  patid, 1253; trial, 25
    ## 
    ## Fixed effects:
    ##                   Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)     -6.029e+01  4.169e+00  9.534e+02 -14.461  < 2e-16 ***
    ## time             7.220e+00  4.026e-02  9.957e+02 179.318  < 2e-16 ***
    ## sexFemale       -1.609e-01  8.760e-02  1.169e+03  -1.837   0.0665 .  
    ## age              1.697e+00  4.006e-02  1.044e+03  42.364  < 2e-16 ***
    ## bmi              1.784e-01  2.729e-02  1.204e+03   6.537 9.26e-11 ***
    ## jobPaid Job     -5.136e+00  8.042e-01  1.059e+03  -6.386 2.54e-10 ***
    ## adl              6.535e-02  1.178e-01  1.054e+03   0.555   0.5791    
    ## wzcNursing Home  1.915e+00  9.790e-02  1.258e+03  19.565  < 2e-16 ***
    ## cdrsb0           2.002e-03  7.603e-03  1.083e+03   0.263   0.7924    
    ## time:cdrsb0     -2.178e-02  4.149e-03  9.707e+02  -5.248 1.88e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) time   sexFml age    bmi    jbPdJb adl    wzcNrH cdrsb0
    ## time        -0.009                                                        
    ## sexFemale    0.026 -0.006                                                 
    ## age         -0.981 -0.003 -0.054                                          
    ## bmi         -0.773 -0.005 -0.002  0.693                                   
    ## jobPaid Job  0.955 -0.001  0.025 -0.960 -0.697                            
    ## adl         -0.975  0.002 -0.001  0.980  0.708 -0.983                     
    ## wzcNursngHm  0.052 -0.016  0.020 -0.086 -0.003  0.037 -0.009              
    ## cdrsb0      -0.017  0.444  0.044  0.006 -0.014 -0.013  0.011  0.016       
    ## time:cdrsb0  0.006 -0.677  0.005  0.002  0.003  0.000  0.001  0.000 -0.651

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Table 4: Fixed Effects Estimates from Mixed Model
</caption>
<thead>
<tr>
<th style="text-align:left;">
effect
</th>
<th style="text-align:left;">
term
</th>
<th style="text-align:right;">
estimate
</th>
<th style="text-align:right;">
std.error
</th>
<th style="text-align:right;">
statistic
</th>
<th style="text-align:right;">
df
</th>
<th style="text-align:right;">
p.value
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
fixed
</td>
<td style="text-align:left;">
(Intercept)
</td>
<td style="text-align:right;">
-60.293
</td>
<td style="text-align:right;">
4.169
</td>
<td style="text-align:right;">
-14.461
</td>
<td style="text-align:right;">
953.390
</td>
<td style="text-align:right;">
0.000
</td>
</tr>
<tr>
<td style="text-align:left;">
fixed
</td>
<td style="text-align:left;">
time
</td>
<td style="text-align:right;">
7.220
</td>
<td style="text-align:right;">
0.040
</td>
<td style="text-align:right;">
179.318
</td>
<td style="text-align:right;">
995.729
</td>
<td style="text-align:right;">
0.000
</td>
</tr>
<tr>
<td style="text-align:left;">
fixed
</td>
<td style="text-align:left;">
sexFemale
</td>
<td style="text-align:right;">
-0.161
</td>
<td style="text-align:right;">
0.088
</td>
<td style="text-align:right;">
-1.837
</td>
<td style="text-align:right;">
1169.059
</td>
<td style="text-align:right;">
0.067
</td>
</tr>
<tr>
<td style="text-align:left;">
fixed
</td>
<td style="text-align:left;">
age
</td>
<td style="text-align:right;">
1.697
</td>
<td style="text-align:right;">
0.040
</td>
<td style="text-align:right;">
42.364
</td>
<td style="text-align:right;">
1043.954
</td>
<td style="text-align:right;">
0.000
</td>
</tr>
<tr>
<td style="text-align:left;">
fixed
</td>
<td style="text-align:left;">
bmi
</td>
<td style="text-align:right;">
0.178
</td>
<td style="text-align:right;">
0.027
</td>
<td style="text-align:right;">
6.537
</td>
<td style="text-align:right;">
1204.048
</td>
<td style="text-align:right;">
0.000
</td>
</tr>
<tr>
<td style="text-align:left;">
fixed
</td>
<td style="text-align:left;">
jobPaid Job
</td>
<td style="text-align:right;">
-5.136
</td>
<td style="text-align:right;">
0.804
</td>
<td style="text-align:right;">
-6.386
</td>
<td style="text-align:right;">
1058.656
</td>
<td style="text-align:right;">
0.000
</td>
</tr>
<tr>
<td style="text-align:left;">
fixed
</td>
<td style="text-align:left;">
adl
</td>
<td style="text-align:right;">
0.065
</td>
<td style="text-align:right;">
0.118
</td>
<td style="text-align:right;">
0.555
</td>
<td style="text-align:right;">
1053.977
</td>
<td style="text-align:right;">
0.579
</td>
</tr>
<tr>
<td style="text-align:left;">
fixed
</td>
<td style="text-align:left;">
wzcNursing Home
</td>
<td style="text-align:right;">
1.915
</td>
<td style="text-align:right;">
0.098
</td>
<td style="text-align:right;">
19.565
</td>
<td style="text-align:right;">
1258.449
</td>
<td style="text-align:right;">
0.000
</td>
</tr>
<tr>
<td style="text-align:left;">
fixed
</td>
<td style="text-align:left;">
cdrsb0
</td>
<td style="text-align:right;">
0.002
</td>
<td style="text-align:right;">
0.008
</td>
<td style="text-align:right;">
0.263
</td>
<td style="text-align:right;">
1082.728
</td>
<td style="text-align:right;">
0.792
</td>
</tr>
<tr>
<td style="text-align:left;">
fixed
</td>
<td style="text-align:left;">
time:cdrsb0
</td>
<td style="text-align:right;">
-0.022
</td>
<td style="text-align:right;">
0.004
</td>
<td style="text-align:right;">
-5.248
</td>
<td style="text-align:right;">
970.693
</td>
<td style="text-align:right;">
0.000
</td>
</tr>
</tbody>
</table>
<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Table 5: Random-Effects Variance–Covariance Components
</caption>
<thead>
<tr>
<th style="text-align:left;">
Group
</th>
<th style="text-align:left;">
Term1
</th>
<th style="text-align:left;">
Term2
</th>
<th style="text-align:right;">
Value
</th>
<th style="text-align:right;">
SD_or_Corr
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
patid
</td>
<td style="text-align:left;">
(Intercept)
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.807
</td>
<td style="text-align:right;">
0.898
</td>
</tr>
<tr>
<td style="text-align:left;">
patid
</td>
<td style="text-align:left;">
time
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.686
</td>
<td style="text-align:right;">
0.828
</td>
</tr>
<tr>
<td style="text-align:left;">
patid
</td>
<td style="text-align:left;">
(Intercept)
</td>
<td style="text-align:left;">
time
</td>
<td style="text-align:right;">
-0.455
</td>
<td style="text-align:right;">
-0.611
</td>
</tr>
<tr>
<td style="text-align:left;">
trial
</td>
<td style="text-align:left;">
(Intercept)
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
8.004
</td>
<td style="text-align:right;">
2.829
</td>
</tr>
<tr>
<td style="text-align:left;">
Residual
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
2.698
</td>
<td style="text-align:right;">
1.643
</td>
</tr>
</tbody>
</table>

The linear mixed-effects model was fitted to 6,220 observations from
1,253 patients across 25 trials. The REML criterion at convergence was
27,121.3, and scaled residuals ranged from −3.64 to 3.32.The
random-effects structure included patient-level random intercepts and
random slopes for time. These patient-specific intercepts and slopes
were correlated. The model also included trial-level random intercepts
to account for differences across trials, in addition to a residual
variance component. The full set of estimated variances, standard
deviations, and correlations for all random-effect terms is presented in
Table 5.

Several fixed effects were statistically significant, including the
intercept, time, age, BMI, job status, WZC category, and the
time-by-cdrsb0 interaction, all with (p\<0.001). The effects of sex,
ADL, and the main effect of cdrsb0 were not statistically significant.
Full estimates, standard errors, and significance levels for all fixed
effects are presented in Table 4.

##### Discussion

Using a mixed-effects model, we examined how psychiatric symptoms change
over time in patients with Alzheimer’s disease. This modelling approach
accounts for individual differences in both baseline symptom severity
and rates of progression.

The analysis showed that psychiatric symptoms worsen steadily over time,
with patients on average increasing by about 7 BPRS points per year.
However, not all patients start at the same level or decline at the same
speed. Our model captured this variation and revealed important patient
characteristics that help explain these differences.

Patients who were older, had higher BMI, poorer daily functioning (ADL),
lived in a nursing home, or did not have a paid job tended to have more
severe psychiatric symptoms at baseline. In addition, patients with more
impaired cognition at diagnosis (higher CDR-SB) showed faster worsening
of psychiatric symptoms over time.

Overall, this approach provides a detailed picture of how psychiatric
symptoms evolve in Alzheimer’s disease. It shows which baseline factors
are linked with faster decline, and it gives individualized estimates of
symptom progression. These findings can help to better understand
patient variability, identify patients who may be at higher risk for
rapid deterioration, and guide monitoring and care strategies over the
course of the disease.

#### Reflection on parameterization

The chosen parameterization had two main components:

Mean structure (fixed effects): We specified a linear effect of time,
adjusted for sex, age, BMI, job status, ADL, nursing home residence
(WZC), baseline CDR-SB, amyloid PET, and the interaction between time
and CDR-SB. This reflects a clinical hypothesis that:

1)  symptom severity depends on demographic, functional, and disease
    severity markers at baseline, and

2)  cognitive status modifies the trajectory of psychiatric symptoms
    (time × CDRSB0).

Although the assumption of linear change simplifies interpretation,
non-linear patterns could be explored in future work if justified..

For the covariance structure, we selected patient-level random
intercepts and random slopes for time with an unstructured covariance
matrix. This parameterization captures individual differences in
baseline severity and rates of change while allowing them to be
correlated. Trial-level random intercepts accounted for between-trial
variation. Overall, this structure provides a good balance between
flexibility and interpretability. Alternative parameterizations, such as
non-linear time or more complex random effects, could be considered, but
the chosen model is appropriate and coherent for the current analysis.

#### Recommendations for future similar experiments

From a modeling and design perspective, these are recommendations that
emerged:

1)  systematic recording of dropout reasons: Since dropout can bias
    longitudinal analyses, future studies should routinely document why
    patients leave the study (death, institutionalization, refusal,
    health decline). This would allow sensitivity analyses under
    different missing-data assumptions (MAR vs MNAR).

2)  Richer covariate information on treatment and care: Including
    medication use (e.g., antipsychotics, antidepressants), psychosocial
    interventions, caregiver burden, and comorbidities might help
    explain additional between-patient variation and could identify
    modifiable risk factors for rapid deterioration.

3)  Parallel assessment of multiple outcome domains Since cognition,
    function, and psychiatric symptoms are closely related in
    Alzheimer’s disease, future experiments could be designed to jointly
    model these domains rather than analyzing each in isolation. This
    would provide a more integrated picture of disease progression and
    may better inform clinical decision-making.
