# LDA Homework 3 (2025-2026)
Anita Kerubo Ogero (2469491), Okolie Tess Eucharia (2468507), Moses
Mburu (2469245), Dorothy Chepkoech (2469284)

## Methodology

### Data structure and notation

Let $i = 1,\dots,N$ index patients and $j = 0,\dots,J$ index scheduled
follow-up visits (years since baseline). Two outcomes were analysed:

- Continuous outcome: $Y_{ij}$ (BPRS).
- Binary outcome: $D_{ij}\in\{0,1\}$ (event/condition indicator defined
  in the data dictionary).

Baseline covariates are denoted by $\mathbf{x}_i$ and include sex, age,
BMI, job status, ADL, living situation (WZC), baseline cognition
(CDR-SB), and additional baseline biomarkers when informative.
Trial/centre is indexed by $\ell(i)\in\{1,\dots,L\}$.

All analyses were performed in long format with one row per $(i,j)$
observation.

### 1) Missingness exploration and classification of patterns

#### Missingness indicators and patterns

For each outcome, we defined an observation indicator

$$
R_{ij}=\begin{cases}
1,& \text{if }Y_{ij}\text{ (or }D_{ij}\text{) is observed} \\
0,& \text{if missing.}
\end{cases}
$$

We classified missingness patterns as:

- **Monotone missingness (dropout):**
  $$R_{ij}=0 \Rightarrow R_{ik}=0 \quad \text{for all } k>j.$$

- **Non-monotone (intermittent) missingness:** there exist (j\<k) such
  that (R\_{ij}=0) and (R\_{ik}=1).

For each patient, we defined a dropout time

$$
T_i=\max{j: R_{ij}=1},
$$

and summarised (i) the distribution of (T_i), (ii) the proportion of
patients with monotone dropout, and (iii) transitions (0) (direct
evidence of non-monotone missingness).

#### Descriptive missingness diagnostics

We used complementary diagnostics to characterise missingness:

1.  Missingness by visit: $\Pr(R_{ij}=0)$ across $j$.
2.  Baseline severity vs missingness: distributions of baseline outcome
    values stratified by number of missed visits and/or dropout time
    $T_i$.
3.  Graphical displays (heatmaps/spaghetti with missingness) to
    distinguish dropout from intermittent gaps.

#### Modelling the missingness process (to support MAR arguments)

To assess whether missingness plausibly depends on observed history, we
modelled the probability of being observed at visit $j$ using a
discrete-time model conditional on observed information up to $j-1$. A
generic form is:

$$
\text{logit}\{\Pr(R_{ij}=1 \mid \mathcal{H}_{i,j-1})\}
= \alpha_0 + \alpha_1 j + \boldsymbol{\alpha}_2^\top \mathbf{x}_i + \alpha_3 Y_{i,j-1} + \alpha_4 \Delta Y_{i,j-1},
$$

where $\mathcal{H}_{i,j-1}$ denotes observed history up to $j-1$, and
$\Delta Y_{i,j-1}=Y_{i,j-1}-Y_{i,j-2}$ (when available). This was used
as an empirical check that observation probability is explainable by
observed covariates and past outcomes (consistent with MAR, not a
proof).

### 2) Continuous outcome model: selected linear mixed-effects model

#### Model specification

For BPRS ($Y_{ij}$), we fitted a linear mixed-effects model with
patient-specific random intercepts and slopes, plus a trial-level random
intercept. The fixed-effects mean structure included time, baseline
covariates, and a time-by-CDR-SB interaction:

$$
\begin{aligned}
Y_{ij} &= (\beta_0 + b_{0i} + u_{0\ell(i)}) + (\beta_1 + b_{1i})t_{ij} \\
&+ \beta_2\text{sex}_i + \beta_3\text{age}_i + \beta_4\text{bmi}_i \\
&+ \beta_5\text{job}_i + \beta_6\text{adl}_i + \beta_7\text{wzc}_i \\
&+ \beta_8\text{cdrsb0}_i + \beta_9(t_{ij}\cdot \text{cdrsb0}_i) \\
&+ \varepsilon_{ij}.
\end{aligned}
$$

Random effects were assumed jointly normal:

$$
\begin{pmatrix} b_{0i} \\ b_{1i} \end{pmatrix}
\sim N\left(
\begin{pmatrix} 0 \\ 0 \end{pmatrix},
\begin{pmatrix}
\sigma_{b0}^2 & \sigma_{b0b1} \\
\sigma_{b0b1} & \sigma_{b1}^2
\end{pmatrix}
\right),
\qquad
u_{0\ell} \sim N(0,\sigma_{\text{trial}}^2),
\qquad
\varepsilon_{ij}\sim N(0,\sigma^2).
$$

**Motivation for this parameterisation.** The model targets (i) the
average change in BPRS over time, (ii) heterogeneity in baseline level
and progression rate between patients via ((b\_{0i}, b\_{1i})), and
(iii) systematic level differences between trials via (u\_{0}). The
interaction (t\_{ij}\_i) directly addresses effect modification of
progression by baseline cognition.

#### Estimation and inference

The model was estimated using restricted maximum likelihood (REML) for
variance components. Fixed effects were interpreted as mean differences
in BPRS (main effects) and mean differences in yearly change (time
interaction). Random-effects variance components were reported to
quantify between-patient and between-trial variability.

### 3) Binary outcome models: random-effects model and marginal model

Let $D_{ij}$ denote the binary outcome.

#### 3a) Random-effects model (GLMM)

We fitted a logistic mixed model (subject-specific interpretation) with
a patient-level random intercept and a trial-level random intercept:

$$
\begin{aligned}
\text{logit}\{\Pr(D_{ij}=1 \mid b_{0i}, u_{0\ell(i)})\} &= \eta_{ij} \\
&= (\gamma_0 + b_{0i} + u_{0\ell(i)}) + \gamma_1 t_{ij} \\
&+ \boldsymbol{\gamma}_2^\top \mathbf{x}_i \\
&+ \gamma_3 (t_{ij}\cdot \text{cdrsb0}_i),
\end{aligned}
$$

with

$$
b_{0i}\sim N(0,\tau_b^2),\qquad
u_{0\ell}\sim N(0,\tau_{\text{trial}}^2).
$$

**Motivation.** This model quantifies how covariates affect the
*subject-specific* odds of the event while accounting for clustering
within patient and systematic differences across trials.

#### 3b) Marginal model (GEE)

For a population-average interpretation, we fitted a generalized
estimating equation model with logit link:

$$
\text{logit}\{\Pr(D_{ij}=1)\}
= \delta_0 + \delta_1 t_{ij} + \boldsymbol{\delta}_2^\top \mathbf{x}_i + \delta_3 (t_{ij}\cdot \text{cdrsb0}_i).
$$

Let $\boldsymbol{\mu}_i = E(\mathbf{D}_i)$ and $\mathbf{A}_i$ be the
diagonal matrix of marginal variances $\mu_{ij}(1-\mu_{ij})$. A working
correlation matrix $\mathbf{R}(\alpha)$ was specified, giving

$$
\mathbf{V}_i = \mathbf{A}_i^{1/2},\mathbf{R}(\alpha),\mathbf{A}_i^{1/2}.
$$

Regression parameters $\boldsymbol{\delta}$ were obtained by solving the
GEE:

$$
\sum_{i=1}^N \mathbf{D}_i^\top \mathbf{V}_i^{-1}(\mathbf{D}_i-\boldsymbol{\mu}_i)=\mathbf{0},
$$

with robust (sandwich) standard errors for valid inference under
possible misspecification of $\mathbf{R}(\alpha)$.

### 4) One MAR-valid version for each of the three models and justification

Let MAR mean:

$$
\Pr(R_{ij}=1 \mid \text{full data}) = \Pr(R_{ij}=1 \mid \text{observed data}),
$$

i.e., missingness depends only on observed outcomes/covariates, not on
the unobserved value at the same visit after conditioning.

#### Linear mixed model (continuous outcome) under MAR

The REML/ML likelihood for the linear mixed model uses all observed
measurements and yields valid inference under MAR when the mean and
random-effects structure are correctly specified. No ad-hoc imputation
was required for MAR validity.

**Chosen MAR version:** likelihood-based LMM fitted to all available
(Y\_{ij}) using REML (and ML for any likelihood-based comparisons if
needed).

#### Logistic GLMM (binary outcome) under MAR

Similarly, maximum likelihood estimation of the GLMM yields valid
inference under MAR, again assuming the model for
$\Pr(D_{ij}=1\mid \cdot)$ and random-effects distribution is correctly
specified.

**Chosen MAR version:** likelihood-based logistic GLMM fitted to all
available $D_{ij}$ using maximum likelihood.

#### GEE (binary outcome) under MAR: inverse-probability weighted GEE (IPW-GEE)

Standard GEE with missing outcomes is generally valid under MCAR; to
obtain a version valid under MAR, we used inverse-probability weights
derived from the observation model for $R_{ij}$. Let

$$
\pi_{ij} = \Pr(R_{ij}=1 \mid \mathcal{H}_{i,j-1}),
$$

estimated from the discrete-time model for observation. Define
stabilized weights

$$
w_{ij} = \frac{\Pr(R_{ij}=1 \mid j)}{\widehat{\pi}_{ij}}.
$$

The weighted GEE solves:

$$
\sum_{i=1}^N \mathbf{D}_i^\top \mathbf{W}_i,\mathbf{V}_i^{-1},(\mathbf{D}_i-\boldsymbol{\mu}_i)=\mathbf{0},
$$

where $\mathbf{W}_i$ is diagonal with entries $w_{ij}$ for observed
outcomes. Robust sandwich standard errors were used.

**Chosen MAR version:** IPW-GEE with stabilized weights because it
directly targets population-average effects while adjusting for
dropout/intermittent missingness driven by observed history.

### 5) Agreement between exploratory findings and fitted models

To evaluate consistency between empirical features and model
implications, we compared:

#### Empirical vs model-implied variance

From exploratory work, the empirical variance function was estimated by
smoothing squared residuals from a simple mean trend model:

$$
\widehat{v}(t) = \text{LOESS}\left(r_{ij}^2 \text{ on } t_{ij}\right),
\quad r_{ij}=Y_{ij}-\widehat{E}(Y_{ij}\mid t_{ij}).
$$

For the selected linear mixed model, the implied marginal variance at
time (t) is

$$
\text{Var}(Y_{ij}\mid t)=
\text{Var}(b_{0i}+b_{1i}t) + \sigma^2
=
\sigma_{b0}^2 + 2t,\sigma_{b0b1} + t^2\sigma_{b1}^2 + \sigma^2.
$$

We examined whether $\widehat{v}(t)$ is compatible with the quadratic
form implied by random intercept/slope plus constant residual variance.
Discrepancies were discussed in terms of (i) sampling variability at
late visits, (ii) the role of random slopes in generating time-dependent
variance even with constant $\sigma^2$, and (iii) whether additional
residual heteroscedasticity would be needed.

#### Empirical vs model-implied correlation

We estimated an empirical semi-variogram based on within-patient
residual differences:

$$
\widehat{\gamma}(u) = \frac{1}{2},\widehat{E}\left[(r_{ik}-r_{ij})^2 ,\middle|, |t_k-t_j|=u\right].
$$

For the mixed model, correlation is induced by shared random effects;
for the GEE model, correlation is represented through the working
structure $\mathbf{R}(\alpha)$. We assessed whether the empirical
correlation decay with lag is broadly consistent with the fitted
random-effects correlation (or working correlation), noting that robust
inference in GEE does not require perfect correlation specification.

### 6) Sensitivity analysis for the continuous outcome

Because MAR is not empirically testable, we performed a brief
sensitivity analysis to departures from MAR using a pattern-mixture
style $\delta$-adjustment. After defining dropout time $T_i$, we
considered that unobserved post-dropout outcomes might be systematically
higher (or lower) than predicted under MAR by a constant shift $\delta$.
Conceptually:

$$
Y_{ij}^{(\delta)}=
\begin{cases}
Y_{ij}, & R_{ij}=1, \\
Y_{ij}^{\text{MAR}}+\delta, & R_{ij}=0 \text{ and } j>T_i,
\end{cases}
$$

where $Y_{ij}^{\text{MAR}}$ represents the model-based conditional
expectation under MAR. We refitted the selected linear mixed model
across a small, pre-specified grid of $\delta$ values (e.g.,
$\delta \in \{-5,-2,0,2,5\}$) and assessed robustness of key estimands
(time effect and $t\times \text{cdrsb0}$ interaction). This provides a
transparent check of how conclusions would change if missing outcomes
after dropout were systematically worse or better than MAR predicts.

### Software implementation

Analyses were implemented in R using likelihood-based mixed modelling
for the continuous outcome (linear mixed model), logistic mixed
modelling for the binary outcome (GLMM), and marginal modelling via GEE
with robust standard errors. The MAR-adjusted marginal binary analysis
used inverse-probability weighted GEE with weights estimated from an
explicit observation model based on observed history.
