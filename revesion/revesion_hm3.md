Revision Homework 3
================

Let $Y_{ij}$ denote the outcome for subject $i$ at visit $j$, and define
the missingness indicator as:

$$
R_{ij} =
\begin{cases}
1, & \text{if } Y_{ij} \text{ is observed} \\
0, & \text{if } Y_{ij} \text{ is missing}
\end{cases}
$$

We also decompose the vector of outcomes for subject $i$ as:

$$
Y_i = (Y_i^{\text{obs}},\ Y_i^{\text{mis}})
$$

where $Y_i^{\text{obs}}$ and $Y_i^{\text{mis}}$ are the observed and
missing components, respectively.

### 1. Missingness Pattern

- **Complete**: No missing data for subject $i$, i.e., $R_{ij} = 1$ for
  all $j$.

- **Monotone**: Once missing, always missing thereafter. There exists
  $t_i$ such that:

  $$
  R_{ij} = 1 \text{ for } j \le t_i, \quad R_{ij} = 0 \text{ for } j > t_i.
  $$

- **Non-monotone**: Intermittent missingness, where missing visits are
  followed by observed ones. There exist $j < k$ such that:

  $$
  R_{ij} = 0 \quad \text{and} \quad R_{ik} = 1.
  $$

### 2. Dropout Pattern

- **Complete**: All scheduled visits are attended.
- **Dropout**: Observed until a certain point, then missing thereafter
  (monotone).
- **Intermittent**: Misses at least one visit but returns later
  (non-monotone).

### 3. Model Framework

We classify joint models for $(Y_i, R_i)$ into the following frameworks:

#### a) SEM – Selection Models

$$
p(Y_i, R_i \mid X_i) = p(Y_i \mid X_i; \theta)\ p(R_i \mid Y_i, X_i; \psi)
$$

Interpretation: Model outcome process first, then dropout as a function
of observed and unobserved outcomes. Often used in MNAR sensitivity
analysis.

#### b) PMM – Pattern Mixture Models

$$
p(Y_i, R_i \mid X_i) = p(R_i \mid X_i; \psi)\ p(Y_i \mid R_i, X_i; \theta)
$$

Interpretation: Stratify the analysis by dropout pattern (e.g., dropout
time), and model outcomes within pattern strata.

#### c) SPM – Shared Parameter Models

$$
p(Y_i, R_i \mid X_i) = \int p(Y_i \mid b_i, X_i)\ p(R_i \mid b_i, X_i)\ p(b_i)\ db_i
$$

Interpretation: Missingness and outcome share latent random effects
$b_i$, e.g., frailty.

### 4. Missingness Mechanism

This characterizes the assumptions on how the dropout process depends on
outcomes and covariates:

- **MCAR (Missing Completely at Random)**:

  $$
  p(R_i \mid Y_i, X_i) = p(R_i)
  $$

  Missingness is independent of all data. Rarely plausible.

- **MAR (Missing at Random)**: $$
  p(R_i \mid Y_i, X_i) = p(R_i \mid Y_i^{\text{obs}}, X_i)
  $$

  Missingness depends only on observed outcomes and covariates. This
  assumption allows likelihood-based and Bayesian inference without
  explicitly modeling $R_i$.

- **MNAR (Missing Not at Random)**: $$
  p(R_i \mid Y_i, X_i) \text{ depends on } Y_i^{\text{mis}}
  $$ Missingness depends on unobserved values, requiring joint modeling
  or sensitivity analysis.

### 5. Ignorability

This concept tells us whether we can ignore the missing data mechanism
in inference:

- **Ignorable**: Holds when the mechanism is MAR *and* the parameter
  sets $\theta$ (for $Y$) and $\psi$ (for $R$) are distinct. In this
  case, valid inference for $\theta$ can proceed using the observed-data
  likelihood:

  $$
  L(\theta) \propto p(Y^{\text{obs}} \mid X; \theta)
  $$

- **Non-ignorable**: Either the mechanism is MNAR or parameters are not
  distinct. Then, inference must model $R$ jointly with $Y$, or include
  sensitivity analysis.

### 6. Inference Paradigm

- **Frequentist**: Parameters are fixed; uncertainty quantified via
  sampling distributions. Estimation uses ML, REML, or estimating
  equations. Example: LMMs, GEE.

- **Likelihood-based**: Emphasizes the likelihood function:

  $$
  L(\theta; Y^{\text{obs}})
  $$ Inference is based on likelihood profiles, likelihood ratio tests,
  and AIC/BIC.

- **Bayesian**: Parameters are random variables

  $$
  p(\theta \mid Y^{\text{obs}}) \propto L(\theta; Y^{\text{obs}})\ p(\theta)
  $$ Posterior inference naturally incorporates missing data via
  imputation or marginalization.

### Summary

This taxonomy enables us to:

- Describe the nature of dropout/missingness in real-world data.
- Choose appropriate models (e.g., LMM, GLMM, GEE) under the MAR
  assumption.
- Perform sensitivity analyses (SEM, PMM, SPM) under MNAR.
- Select a valid inferential strategy (frequentist, likelihood, or
  Bayesian).

For the Alzheimer’s dataset, classification of subjects into
*completers*, *dropouts*, and *intermittent missers* will guide our
modeling approach. Under MAR, likelihood-based methods (linear mixed
models, GLMM, GEE) remain valid and are adopted for primary analyses.
MNAR concerns will be addressed using shared parameter or
pattern-mixture approaches in sensitivity analyses.

### 30.15 Modeling frameworks and missing-data mechanisms (notes linked to the Alzheimer cohort)

#### Notation (what the slide is talking about)

For patient $i=1,\dots,N$, let the longitudinal outcome vector be

$$
y_i = (y_{i0}, y_{i1}, \dots, y_{iJ}),
$$

where $j=0$ is baseline and $j=1,\dots,J$ are follow-up visits (e.g.,
yearly). Let $X_i$ be baseline covariates (age, sex, education, baseline
severity markers, etc.). Define the missingness indicator vector

$$
r_i = (r_{i0}, r_{i1}, \dots, r_{iJ}), \qquad
r_{ij} =
\begin{cases}
1, & \text{if } y_{ij}\ \text{is observed},\
0, & \text{if } y_{ij}\ \text{is missing}.
\end{cases}
$$

Split outcomes into observed and missing parts:
$y_i = (y_i^{o}, y_i^{m})$.

The “big object” on the slide is the joint distribution:

$$
f(y_i, r_i \mid X_i, \theta, \psi),
$$

where $\theta$ indexes the outcome model parameters and $\psi$ indexes
the missingness model parameters.

In the Alzheimer data, $y_i$ could be the continuous BPRS trajectory,
while $r_i$ captures the fact that patients increasingly miss later
visits (often because of dropout, illness progression,
institutionalisation, or death). This is exactly why the joint object
$f(y_i,r_i\mid X_i)$ matters: later years are a selected subset of
patients.

### 1) Selection models (SEM): “model $y$ first, then model $r$ given $y$”

Factorisation:

$$
f(y_i, r_i \mid X_i, \theta, \psi) =
f(y_i \mid X_i, \theta),
f(r_i \mid X_i, y_i^{o}, y_i^{m}, \psi).
$$

Interpretation in Alzheimer setting:

- $f(y_i\mid X_i,\theta)$ describes the disease course (e.g., mean BPRS
  change over time, subject heterogeneity).
- $f(r_i\mid \cdot)$ describes the observation/dropout process,
  potentially depending on the patient’s underlying symptom severity
  trajectory.

### Missing-data mechanisms expressed inside the selection model

This is the cleanest way to state MCAR/MAR/MNAR.

- **MCAR (as on the slide, “covariate-dependent MCAR”)**

$$
f(r_i \mid X_i, y_i^{o}, y_i^{m}, \psi) = f(r_i \mid X_i, \psi).
$$

- **MAR**

$$
f(r_i \mid X_i, y_i^{o}, y_i^{m}, \psi) = f(r_i \mid X_i, y_i^{o}, \psi).
$$

A practical MAR-style dropout model (discrete-time hazard):

$$
\Pr(r_{i,j+1}=1 \mid \mathcal{H}_{ij}) =
\operatorname{logit}^{-1}\Big(
\alpha_0 + \alpha_1 y_{ij} + \alpha_2^\top X_i + \alpha_3 j
\Big),
$$

- **MNAR**

$$
f(r_i \mid X_i, y_i^{o}, y_i^{m}, \psi)\ \text{depends on } y_i^{m}.
$$

Example MNAR-flavoured model:

$$
\Pr(r_{i,j+1}=1 \mid \mathcal{H}_{ij}, y_{i,j+1}) =
\operatorname{logit}^{-1}\Big(
\alpha_0 + \alpha_1 y_{i,j+1} + \alpha_2^\top X_i
\Big).
$$

### 2) Pattern-mixture models (PMM): “model $y$ within missingness patterns, then mix”

Factorisation:

$$
f(y_i, r_i \mid X_i, \theta, \psi) =
f(y_i \mid X_i, r_i, \theta),
f(r_i \mid X_i, \psi).
$$

Define dropout time:

$$
t_i = \max{j: r_{ij}=1}.
$$

Let $g_i=t_i$ (last observed visit). Then PMM models:

$$
f(y_i \mid X_i, r_i, \theta) \equiv f(y_i \mid X_i, g_i, \theta).
$$

Delta-adjustment sensitivity (for BPRS):

$$
y_{ij}^{m}(\delta) = \widehat{y}^{MAR}_{ij} + \delta, \qquad j>t_i.
$$

### 3) Shared-parameter models (SPM): “connect $y$ and $r$ through latent patient health”

Latent subject effects $b_i$ drive both processes.

$$
f(y_i, r_i \mid X_i, \theta, \psi) =
\int
f(y_i \mid X_i, b_i, \theta),
f(r_i \mid X_i, b_i, \psi),
f(b_i)
db_i.
$$

Outcome model:

$$
y_{ij} = x_{ij}^\top \beta + z_{ij}^\top b_i + \varepsilon_{ij},
\qquad \varepsilon_{ij}\sim N(0,\sigma^2).
$$

Dropout model:

$$
\Pr(r_{i,j+1}=1 \mid X_i, b_i) =
\operatorname{logit}^{-1}\Big(
\gamma_0 + \gamma_1^\top b_i + \gamma_2^\top X_i + \gamma_3 j
\Big).
$$

#### Practical notes for Alzheimer modeling

- MAR justifies standard mixed models for BPRS.
- Sensitivity via PMM-$\delta$ approach is clinically meaningful.
- SPM offers a biologically plausible missingness narrative in
  Alzheimer’s follow-up.

To edit the provided RMarkdown file for proper MathJax rendering, you
need to ensure that all LaTeX equations are properly formatted for
MathJax. This involves using `$$ ... $$` for block equations and
`\(...\)` for inline equations. Below is the revised version of your
RMarkdown file:

Let $Y_{ij}$ denote the outcome for subject $i$ at visit $j$, and define
the missingness indicator as:

$$
R_{ij} =
\begin{cases}
1, & \text{if } Y_{ij} \text{ is observed}, \\
0, & \text{if } Y_{ij} \text{ is missing}.
\end{cases}
$$

We also decompose the vector of outcomes for subject $i$ as:

$$
Y_i = (Y_i^{\text{obs}},\ Y_i^{\text{mis}})
$$

where $Y_i^{\text{obs}}$ and $Y_i^{\text{mis}}$ are the observed and
missing components, respectively.

### 1. Missingness Pattern

- **Complete**: No missing data for subject $i$, i.e., $R_{ij} = 1$ for
  all $j$.

- **Monotone**: Once missing, always missing thereafter. There exists
  $t_i$ such that:

  $$
  R_{ij} = 1 \text{ for } j \le t_i, \quad R_{ij} = 0 \text{ for } j > t_i.
  $$

- **Non-monotone**: Intermittent missingness, where missing visits are
  followed by observed ones. There exist $j < k$ such that:

  $$
  R_{ij} = 0 \quad \text{and} \quad R_{ik} = 1.
  $$

### 2. Dropout Pattern

- **Complete**: All scheduled visits are attended.
- **Dropout**: Observed until a certain point, then missing thereafter
  (monotone).
- **Intermittent**: Misses at least one visit but returns later
  (non-monotone).

### 3. Model Framework

We classify joint models for $(Y_i, R_i)$ into the following frameworks:

#### a) SEM – Selection Models

$$
p(Y_i, R_i \mid X_i) = p(Y_i \mid X_i; \theta)\ p(R_i \mid Y_i, X_i; \psi)
$$

Interpretation: Model outcome process first, then dropout as a function
of observed and unobserved outcomes. Often used in MNAR sensitivity
analysis.

#### b) PMM – Pattern Mixture Models

$$
p(Y_i, R_i \mid X_i) = p(R_i \mid X_i; \psi)\ p(Y_i \mid R_i, X_i; \theta)
$$

Interpretation: Stratify the analysis by dropout pattern (e.g., dropout
time), and model outcomes within pattern strata.

#### c) SPM – Shared Parameter Models

$$
p(Y_i, R_i \mid X_i) = \int p(Y_i \mid b_i, X_i)\ p(R_i \mid b_i, X_i)\ p(b_i)\ db_i
$$

Interpretation: Missingness and outcome share latent random effects
$b_i$, e.g., frailty.

### 4. Missingness Mechanism

This characterizes the assumptions on how the dropout process depends on
outcomes and covariates:

- **MCAR (Missing Completely at Random)**:

  $$
  p(R_i \mid Y_i, X_i) = p(R_i)
  $$

  Missingness is independent of all data. Rarely plausible.

- **MAR (Missing at Random)**: $$
  p(R_i \mid Y_i, X_i) = p(R_i \mid Y_i^{\text{obs}}, X_i)
  $$

  Missingness depends only on observed outcomes and covariates. This
  assumption allows likelihood-based and Bayesian inference without
  explicitly modeling $R_i$.

- **MNAR (Missing Not at Random)**: $$
  p(R_i \mid Y_i, X_i) \text{ depends on } Y_i^{\text{mis}}
  $$ Missingness depends on unobserved values, requiring joint modeling
  or sensitivity analysis.

### 5. Ignorability

This concept tells us whether we can ignore the missing data mechanism
in inference:

- **Ignorable**: Holds when the mechanism is MAR *and* the parameter
  sets $\theta$ (for $Y$) and $\psi$ (for $R$) are distinct. In this
  case, valid inference for $\theta$ can proceed using the observed-data
  likelihood:

  $$
  L(\theta) \propto p(Y^{\text{obs}} \mid X; \theta)
  $$

- **Non-ignorable**: Either the mechanism is MNAR or parameters are not
  distinct. Then, inference must model $R$ jointly with $Y$, or include
  sensitivity analysis.

### 6. Inference Paradigm

- **Frequentist**: Parameters are fixed; uncertainty quantified via
  sampling distributions. Estimation uses ML, REML, or estimating
  equations. Example: LMMs, GEE.

- **Likelihood-based**: Emphasizes the likelihood function:

  $$
  L(\theta; Y^{\text{obs}})
  $$ Inference is based on likelihood profiles, likelihood ratio tests,
  and AIC/BIC.

- **Bayesian**: Parameters are random variables

  $$
  p(\theta \mid Y^{\text{obs}}) \propto L(\theta; Y^{\text{obs}})\ p(\theta)
  $$ Posterior inference naturally incorporates missing data via
  imputation or marginalization.

$$
R_{ij} =
\begin{cases}
1, & \text{if } Y_{ij} \text{ is observed}, \\
0, & \text{if } Y_{ij} \text{ is missing}.
\end{cases}
$$

We also decompose the vector of outcomes for subject $i$ as:

$$
Y_i = (Y_i^{\text{obs}},\ Y_i^{\text{mis}})
$$

where $Y_i^{\text{obs}}$ and $Y_i^{\text{mis}}$ are the observed and
missing components, respectively.

### 1. Missingness Pattern

- **Complete**: No missing data for subject $i$, i.e., $R_{ij} = 1$ for
  all $j$.

- **Monotone**: Once missing, always missing thereafter. There exists
  $t_i$ such that:

  $$
  R_{ij} = 1 \text{ for } j \le t_i, \quad R_{ij} = 0 \text{ for } j > t_i.
  $$

- **Non-monotone**: Intermittent missingness, where missing visits are
  followed by observed ones. There exist $j < k$ such that:

  $$
  R_{ij} = 0 \quad \text{and} \quad R_{ik} = 1.
  $$

### 2. Dropout Pattern

- **Complete**: All scheduled visits are attended.
- **Dropout**: Observed until a certain point, then missing thereafter
  (monotone).
- **Intermittent**: Misses at least one visit but returns later
  (non-monotone).

### 3. Model Framework

We classify joint models for $(Y_i, R_i)$ into the following frameworks:

#### a) SEM – Selection Models

$$
p(Y_i, R_i \mid X_i) = p(Y_i \mid X_i; \theta)\ p(R_i \mid Y_i, X_i; \psi)
$$

Interpretation: Model outcome process first, then dropout as a function
of observed and unobserved outcomes. Often used in MNAR sensitivity
analysis.

#### b) PMM – Pattern Mixture Models

$$
p(Y_i, R_i \mid X_i) = p(R_i \mid X_i; \psi)\ p(Y_i \mid R_i, X_i; \theta)
$$

Interpretation: Stratify the analysis by dropout pattern (e.g., dropout
time), and model outcomes within pattern strata.

#### c) SPM – Shared Parameter Models

$$
p(Y_i, R_i \mid X_i) = \int p(Y_i \mid b_i, X_i)\ p(R_i \mid b_i, X_i)\ p(b_i)\ db_i
$$

Interpretation: Missingness and outcome share latent random effects
$b_i$, e.g., frailty.

### 4. Missingness Mechanism

This characterizes the assumptions on how the dropout process depends on
outcomes and covariates:

- **MCAR (Missing Completely at Random)**:

  $$
  p(R_i \mid Y_i, X_i) = p(R_i)
  $$

  Missingness is independent of all data. Rarely plausible.

- **MAR (Missing at Random)**:

  $$
  p(R_i \mid Y_i, X_i) = p(R_i \mid Y_i^{\text{obs}}, X_i)
  $$

  Missingness depends only on observed outcomes and covariates. This
  assumption allows likelihood-based and Bayesian inference without
  explicitly modeling $R_i$.

- **MNAR (Missing Not at Random)**:

  $$
  p(R_i \mid Y_i, X_i) \text{ depends on } Y_i^{\text{mis}}
  $$

  Missingness depends on unobserved values, requiring joint modeling or
  sensitivity analysis.

### 5. Ignorability

This concept tells us whether we can ignore the missing data mechanism
in inference:

- **Ignorable**: Holds when the mechanism is MAR *and* the parameter
  sets $\theta$ (for $Y$) and $\psi$ (for $R$) are distinct. In this
  case, valid inference for $\theta$ can proceed using the observed-data
  likelihood:

  $$
  L(\theta) \propto p(Y^{\text{obs}} \mid X; \theta)
  $$

- **Non-ignorable**: Either the mechanism is MNAR or parameters are not
  distinct. Then, inference must model $R$ jointly with $Y$, or include
  sensitivity analysis.

### 6. Inference Paradigm

- **Frequentist**: Parameters are fixed; uncertainty quantified via
  sampling distributions. Estimation uses ML, REML, or estimating
  equations. Example: LMMs, GEE.

- **Likelihood-based**: Emphasizes the likelihood function:

  $$
  L(\theta; Y^{\text{obs}})
  $$

  Inference is based on likelihood profiles, likelihood ratio tests, and
  AIC/BIC.

- **Bayesian**: Parameters are random variables

  $$
  p(\theta \mid Y^{\text{obs}}) \propto L(\theta; Y^{\text{obs}})\ p(\theta)
  $$

  Posterior inference naturally incorporates missing data via imputation
  or marginalization.

You’re going to get hit with some version of: *“You say your exploratory
variance pattern matches the model… how do you justify that in a
defensible way in 30 seconds?”* Here’s a clean defense that sounds like
a statistician, not a student reading slides.

### The core idea to defend

You’re not claiming the two variance functions are *identical*. You’re
claiming they are **compatible at the level that matters for
inference**: same qualitative shape + differences small enough that they
don’t change the fixed‐effect conclusions.

## A tight defense script (20–30 seconds)

“We compared the empirical variance profile from the raw BPRS data to
the variance implied by our fitted mixed model. They’re ‘similar’ in two
senses: first, the variance evolves with time in the same direction and
with the same broad shape; second, the model-based residual diagnostics
did not show systematic time-dependent lack-of-fit. Minor discrepancies
are expected because the empirical variance is noisy and is computed on
changing sample sizes due to dropout. Our inference targets the mean
structure, and we use likelihood-based estimation under MAR with
robustness checks, so small variance-function differences are not
driving our conclusions.”

That is hard to attack in a 5-minute Q&A.

### What to show (if you have 1 slide / 30 seconds)

### 1) One visual comparison

- Plot empirical variance by visit (or time) for BPRS (from observed
  data).
- Overlay the model-implied variance (from your chosen LMM).

You don’t need perfection. You need the audience to see “same
direction + roughly same magnitude”.

### 2) One diagnostic line

Say one of these (pick one you actually did):

- “No clear trend in standardized residual spread across time.”
- “Residual vs fitted plots show roughly constant spread.”
- “Model with heterogeneous residual variance didn’t materially change
  fixed effects.”

### The logic behind “why differences are not fatal” (say it like this)

### Point A: Empirical variance is a *moving target*

Empirical variance at visit (j) is based on whoever is still observed at
(j). With dropout, that subgroup changes.

So even under a correct model, the raw empirical variance can drift
because:

- (n_j) decreases (more noise),
- composition changes (more severe patients drop out, etc.).

That’s the line: *variance estimates are noisier and selection-affected
at later visits*.

### Point B: Mixed models don’t require the empirical variance curve to match perfectly

The LMM is mainly about modeling:

- the mean trajectory (E(Y\_{ij}X_i)),
- and within-subject correlation via random effects.

Residual variance misspecification mostly affects **efficiency** (SEs)
more than **bias** of fixed effects, especially if you did sensible
diagnostics and sensitivity checks.

If they push: you say you checked robustness (below).

------------------------------------------------------------------------

## The “if they push harder” fallback (still short)

### If someone asks: “Would it be a problem if they differ?”

Use this:

“It depends on the kind of difference. If the empirical variance shows
strong heteroscedasticity over time and the model assumes constant
residual variance, SEs could be off and model fit could degrade. That’s
why we checked residual patterns and also considered a model allowing
time-specific residual variance; the fixed-effect conclusions were
stable, so any mismatch doesn’t appear to be driving inference.”

Even if you didn’t fit the heteroscedastic model, you can still say you
**checked residuals** and **the spread looked stable**—but don’t claim a
model you didn’t run.

## A 10-minute structure that makes Q5 easy to defend

**Slide 1 (30s):** Missingness pattern + why MAR primary + MNAR
sensitivity **Slide 2 (2 min):** Continuous outcome model (LMM) **Slide
3 (2 min):** Binary outcome models (GLMM + GEE) **Slide 4 (2 min):**
Missingness diagnostics (MCAR unlikely; MAR plausible) **Slide 5 (2
min):** Sensitivity (()-shift) and robustness **Slide 6 (1 min):** Q5
variance agreement (visual overlay + one diagnostic line)

This way, variance agreement is *not* a random side quest—it’s a quick
validation step.

## Two killer one-liners for Q&A

- **On “similarity”:** “We’re comparing a noisy empirical summary under
  dropout to a parametric model; the criterion is qualitative agreement
  plus no residual evidence of time-varying misfit.”
- **On credibility:** “And we backed it up with an MNAR sensitivity
  analysis, so our conclusions aren’t resting on a fragile variance
  assumption.”

That last line is a confidence grenade—in a good way.

If you tell me what your LMM residual structure is (homoscedastic vs
time-specific vs AR(1) etc.), I can give you the exact 1–2 sentences
that match your model precisely, so you don’t overclaim under pressure.

To clean and format the provided content for RMarkdown with proper
MathJax rendering and multiline equation spacing, here is the revised
version of your file:

Let $Y_{ij}$ denote the outcome for subject $i$ at visit $j$, and define
the missingness indicator as:

$$
R_{ij} =
\begin{cases}
1, & \text{if } Y_{ij} \text{ is observed}, \\
0, & \text{if } Y_{ij} \text{ is missing}.
\end{cases}
$$

We also decompose the vector of outcomes for subject $i$ as:

$$
Y_i = (Y_i^{\text{obs}},\ Y_i^{\text{mis}})
$$

where $Y_i^{\text{obs}}$ and $Y_i^{\text{mis}}$ are the observed and
missing components, respectively.

### 1. Missingness Pattern

- **Complete**: No missing data for subject $i$, i.e., $R_{ij} = 1$ for
  all $j$.

- **Monotone**: Once missing, always missing thereafter. There exists
  $t_i$ such that:

  $$
  R_{ij} = 1 \text{ for } j \le t_i, \quad R_{ij} = 0 \text{ for } j > t_i.
  $$

- **Non-monotone**: Intermittent missingness, where missing visits are
  followed by observed ones. There exist $j < k$ such that:

  $$
  R_{ij} = 0 \quad \text{and} \quad R_{ik} = 1.
  $$

### 2. Dropout Pattern

- **Complete**: All scheduled visits are attended.
- **Dropout**: Observed until a certain point, then missing thereafter
  (monotone).
- **Intermittent**: Misses at least one visit but returns later
  (non-monotone).

### 3. Model Framework

We classify joint models for $(Y_i, R_i)$ into the following frameworks:

#### a) SEM – Selection Models

$$
p(Y_i, R_i \mid X_i) = p(Y_i \mid X_i; \theta)\ p(R_i \mid Y_i, X_i; \psi)
$$

Interpretation: Model outcome process first, then dropout as a function
of observed and unobserved outcomes. Often used in MNAR sensitivity
analysis.

#### b) PMM – Pattern Mixture Models

$$
p(Y_i, R_i \mid X_i) = p(R_i \mid X_i; \psi)\ p(Y_i \mid R_i, X_i; \theta)
$$

Interpretation: Stratify the analysis by dropout pattern (e.g., dropout
time), and model outcomes within pattern strata.

#### c) SPM – Shared Parameter Models

$$
p(Y_i, R_i \mid X_i) = \int p(Y_i \mid b_i, X_i)\ p(R_i \mid b_i, X_i)\ p(b_i)\ db_i
$$

Interpretation: Missingness and outcome share latent random effects
$b_i$, e.g., frailty.

### 4. Missingness Mechanism

This characterizes the assumptions on how the dropout process depends on
outcomes and covariates:

- **MCAR (Missing Completely at Random)**:

  $$
  p(R_i \mid Y_i, X_i) = p(R_i)
  $$

  Missingness is independent of all data. Rarely plausible.

- **MAR (Missing at Random)**:

  $$
  p(R_i \mid Y_i, X_i) = p(R_i \mid Y_i^{\text{obs}}, X_i)
  $$

  Missingness depends only on observed outcomes and covariates. This
  assumption allows likelihood-based and Bayesian inference without
  explicitly modeling $R_i$.

- **MNAR (Missing Not at Random)**:

  $$
  p(R_i \mid Y_i, X_i) \text{ depends on } Y_i^{\text{mis}}
  $$

  Missingness depends on unobserved values, requiring joint modeling or
  sensitivity analysis.

### 5. Ignorability

This concept tells us whether we can ignore the missing data mechanism
in inference:

- **Ignorable**: Holds when the mechanism is MAR *and* the parameter
  sets $\theta$ (for $Y$) and $\psi$ (for $R$) are distinct. In this
  case, valid inference for $\theta$ can proceed using the observed-data
  likelihood:

  $$
  L(\theta) \propto p(Y^{\text{obs}} \mid X; \theta)
  $$

- **Non-ignorable**: Either the mechanism is MNAR or parameters are not
  distinct. Then, inference must model $R$ jointly with $Y$, or include
  sensitivity analysis.

### 6. Inference Paradigm

- **Frequentist**: Parameters are fixed; uncertainty quantified via
  sampling distributions. Estimation uses ML, REML, or estimating
  equations. Example: LMMs, GEE.

- **Likelihood-based**: Emphasizes the likelihood function:

  $$
  L(\theta; Y^{\text{obs}})
  $$

  Inference is based on likelihood profiles, likelihood ratio tests, and
  AIC/BIC.

- **Bayesian**: Parameters are random variables:

  $$
  p(\theta \mid Y^{\text{obs}}) \propto L(\theta; Y^{\text{obs}})\ p(\theta)
  $$

  Posterior inference naturally incorporates missing data via imputation
  or marginalization.

### Summary

This taxonomy enables us to:

- Describe the nature of dropout/missingness in real-world data.
- Choose appropriate models (e.g., LMM, GLMM, GEE) under the MAR
  assumption.
- Perform sensitivity analyses (SEM, PMM, SPM) under MNAR.
- Select a valid inferential strategy (frequentist, likelihood, or
  Bayesian).

For *sensitivity to MNAR* in your Alzheimer/BPRS setting, you want a
model that (i) starts from your MAR analysis, (ii) perturbs the
unobserved part in a controlled way, and (iii) is easy to defend in 5
minutes of Q&A.

The best default is a **pattern–mixture sensitivity analysis with a
$\delta$-adjustment** for the *continuous* BPRS outcome. If you have
time for a second sensitivity “scenario,” a **shared-parameter (joint)
model** is the next most defensible.

### Recommended sensitivity model for BPRS (continuous): Pattern–mixture with $\delta$-shift

You fit your main LMM under MAR (direct likelihood). Then you define
post-dropout counterfactual values by shifting the MAR-based
predictions.

Let $t_i$ be the last observed visit for subject $i$. For $j > t_i$,
define:

$$
y_{ij}^m(\delta) = \widehat{y}_{ij}^{MAR} + \delta .
$$

**Interpretation (in Alzheimer language):** After dropout, the
unobserved BPRS may be systematically worse (or better) than what MAR
predicts.

**How to choose $\delta$ (defensible choices):**

- Use clinically interpretable increments on the BPRS scale, or
- Use standardized shifts like
  $\delta \in \{0, \pm 0.25 s, \pm 0.5 s\}$, where $s$ is the baseline
  SD of BPRS.

You then refit (or re-evaluate) the estimand of interest (time effect,
treatment/group effect) across $\delta$ values and report the “tipping
point” where conclusions change.

A slightly richer version allows the MNAR departure to grow with time
since dropout:

$$
y_{ij}^m(\delta) = \widehat{y}_{ij}^{MAR} + \delta (j - t_i), \qquad j > t_i .
$$

This is often very plausible in Alzheimer follow-up: longer unobserved
time $\Rightarrow$ larger potential departure.

**Why this is easy to defend:**

- It directly targets the *unidentifiable* part: $y_i^m$.
- It’s transparent: $\delta = 0$ is MAR; $\delta \neq 0$ quantifies MNAR
  departure.
- It matches the clinical story: dropout linked to worsening.

### Optional second sensitivity model (if you want one): Shared-parameter (joint) model

If you want a more “mechanistic” MNAR sensitivity, introduce random
effects $b_i$ shared between the BPRS model and dropout.

**Outcome model:**

$$
y_{ij} = x_{ij}^\top \beta + z_{ij}^\top b_i + \varepsilon_{ij}, \qquad \varepsilon_{ij} \sim N(0, \sigma^2).
$$

**Dropout/observation model:**

$$
\Pr(r_{i,j+1} = 1 \mid X_i, b_i) =
\operatorname{logit}^{-1}\Big(\gamma_0 + \gamma_1^\top b_i + \gamma_2^\top X_i + \gamma_3 j\Big).
$$

This encodes “latent severity drives both worse BPRS and higher
dropout.” It’s defensible, but heavier to implement and explain quickly.

## What about the binary outcome sensitivity?

Your assignment’s explicit MNAR sensitivity is for the continuous
outcome (BPRS). For the binary outcome, the clean approach is:

- **Primary MAR analyses:** GLMM (random-effects logistic) + GEE
  (marginal), as requested.
- **If you want a sensitivity-style robustness check for the binary
  outcome:** Do it as a *MAR robustness* check (not MNAR) by comparing
  GLMM vs GEE and/or using weighted GEE (IPW) if you model observation
  probabilities under MAR.

True MNAR sensitivity for the binary outcome
(selection/shared-parameter) is possible, but it’s usually overkill for
a “brief” sensitivity component unless the lecturer explicitly asks.

### Bottom line (what I’d put in your report)

- **Sensitivity for BPRS:** Pattern–mixture $\delta$-adjustment (primary
  recommendation).
- **Optional extra:** Shared-parameter joint model as a second MNAR
  scenario (only if you have bandwidth).
- **Binary outcome:** Keep to GLMM + GEE under MAR; optionally add
  weighted GEE as a robustness check.

If you tell me your BPRS coding (range and whether higher = worse) and
the visit schedule ($j$) (e.g., 0–6 years), I’ll write the exact
sensitivity subsection text in your methodology style, including a
principled $\delta$ grid that doesn’t look arbitrary.

### 30.21 Illustration: why **CC** and **LOCF** are risky, and why we work under **MAR** for the Alzheimer analyses

This slide is doing one job: it shows that *how you handle missing
follow-up values changes both the point estimate and the uncertainty*,
even in a toy example. In Alzheimer follow-up (BPRS, cognition,
function), missingness is common and often related to observed disease
history, so the “MCAR/simple” fixes can mislead.

## Notation (generic, but matches our setting)

Let $Y_i$ be an outcome at a given follow-up (continuous BPRS or a
binary indicator), and let

$$
R_i =
\begin{cases}
1, & \text{if } Y_i \text{ is observed}, \\
0, & \text{if } Y_i \text{ is missing}.
\end{cases}
$$

A simple estimand is the marginal mean at that follow-up,

$$
\theta = \mathbb{E}(Y_i).
$$

For a binary outcome, $\theta$ is a probability; for BPRS it is a mean
score.

### Complete Case (CC): “analyse only those with complete follow-up”

The CC estimator uses only $R_i = 1$ subjects:

$$
\widehat{\theta}_{CC} = \frac{1}{n_{\text{obs}}} \sum_{i:R_i=1} Y_i
= \mathbb{E}(Y_i \mid R_i=1)\ \text{(estimated)}.
$$

So the *target changes* from $\mathbb{E}(Y_i)$ to
$\mathbb{E}(Y_i \mid R_i=1)$. The bias is

$$
\mathbb{E}(\widehat{\theta}_{CC}) - \theta
=
\mathbb{E}(Y_i \mid R_i=1) - \mathbb{E}(Y_i).
$$

- If missingness is **MCAR**, then $R_i \perp Y_i$ and CC can be
  unbiased.
- If missingness depends on observed history (very plausible in
  Alzheimer follow-up), then $R_i$ is associated with $Y_i$ and CC tends
  to be biased.
- Even when bias is modest, CC is **inefficient** because it throws away
  partial trajectories, which is why the slide shows “too wide”
  intervals.

**Alzheimer link:** Patients with worse observed progression often miss
later visits (institutionalisation, death, inability to attend). CC then
over-represents healthier/stabler patients at later visits, typically
pulling estimates toward “less severe” trajectories.

### Last Observation Carried Forward (LOCF): “freeze the trajectory after missingness”

LOCF replaces a missing follow-up value with the last observed value.
Let $t_i$ be the last observed visit for subject $i$, and define the
imputed value

$$
Y_i^{\text{LOCF}} = Y_{i,t_i}.
$$

Then

$$
\widehat{\theta}_{\text{LOCF}} = \frac{1}{n} \sum_{i=1}^n Y_i^{\text{LOCF}}.
$$

Two problems (exactly what the slide labels):

1.  **Bias**: LOCF implicitly assumes “no change after the last observed
    visit.” In progressive disease settings, that is usually
    implausible. For BPRS, it effectively says symptoms remain flat
    after dropout.

2.  **Too narrow uncertainty**: LOCF treats imputed values as if they
    were observed, so it ignores imputation uncertainty. That’s why the
    slide shows a confidence interval that is “too narrow” even when the
    point estimate is off.

**Alzheimer link:** If BPRS tends to worsen (or fluctuate with
worsening) among those who later drop out, LOCF tends to understate
worsening (or distort the time trend), while also overstating precision.

### MAR-based analysis: “use a model that conditions on observed history and integrates over what is missing”

Under MAR, missingness can depend on observed information but not on
unobserved outcomes once the observed history is accounted for:

$$
p(R_i \mid Y_i, X_i) = p(R_i \mid Y_i^{o}, X_i).
$$

For the primary analyses, MAR is operationalised using methods that are
valid under MAR and use all available longitudinal information:

### (i) Direct likelihood (mixed models)

For continuous BPRS (LMM), inference is based on the observed-data
likelihood:

$$
L(\theta) = \prod_{i=1}^N \int f(y_i^{o}, y_i^{m} \mid X_i, \theta)\, dy_i^{m}.
$$

This is exactly what the slide’s “MAR: direct likelihood” line is
pointing at: you don’t fill in values by hand; you fit the model and
integrate over the missing part.

### (ii) For binary outcomes: GLMM and marginal (GEE) approaches

- A **random-effects (GLMM)** uses a subject-specific model and is
  likelihood-based (MAR-valid under correct specification).
- A **marginal model (GEE)** targets population-averaged effects; under
  MAR, a common extension is **weighted GEE** (IPW) when missingness is
  associated with observed history.

### How to defend this in 10 minutes (presentation-friendly)

A clean defence that lands fast:

1.  **What the slide demonstrates:** CC and LOCF are “simple” but they
    change the estimand or bake in unrealistic assumptions, and they can
    distort uncertainty.

2.  **Why it matters here:** In Alzheimer follow-up, missingness is
    rarely plausibly MCAR; it is typically related to observed disease
    history and patient status.

3.  **What we therefore do:** Primary inference is built on MAR-valid
    methods that use observed history efficiently (likelihood-based
    mixed models for BPRS; GLMM and GEE for the binary outcome), and we
    reserve MNAR assumptions for sensitivity analysis rather than
    treating them as truth.

If someone pushes you with “but MAR isn’t testable,” the crisp reply is:
**correct—so MAR is the working assumption for primary inference, and
MNAR is addressed via sensitivity analysis.** That stance is exactly
what the course’s “overview” slide is advocating.

### Quantifying the bias of CC and LOCF (and why we avoid them for the Alzheimer BPRS analysis)

This slide is using a **simple two-visit longitudinal setting**
(baseline and one follow-up) to show why “simple” approaches such as
**complete cases (CC)** and **last observation carried forward (LOCF)**
can mislead when dropout is related to patient evolution — which is very
plausible in an Alzheimer cohort.

#### Set-up and notation (linked to our cohort)

Let $Y_{ij}$ be the continuous outcome for patient $i$ at visit $j$,
with:

- $j=0$: baseline
- $j=1$: follow-up

Let $T_i \in \{0,1\}$ denote a treatment/group indicator (in our
setting: this can be read more generally as a **group indicator** or a
covariate-defined contrast).

Define an indicator for being a **dropout** vs a **completer**:

- **Dropouts**: observed only at baseline ($t_{ij}=0$), with probability
  $p_0$
- **Completers**: observed at baseline and follow-up ($t_{ij}=0,1$),
  with probability $p_1 = 1-p_0$

The key idea is that **dropouts and completers may follow different mean
trajectories**. The slide writes this as two separate mean models.

For **dropouts**:

$$
E(Y_{ij} \mid T_i, \text{dropout}) =
\beta_0 + \beta_1 T_i + \beta_2 t_{ij} + \beta_3 T_i t_{ij}.
$$

For **completers**:

$$
E(Y_{ij} \mid T_i, \text{complete}) =
\gamma_0 + \gamma_1 T_i + \gamma_2 t_{ij} + \gamma_3 T_i t_{ij}.
$$

In an Alzheimer cohort, this is not an artificial scenario: patients who
discontinue follow-up often do so because of **worsening symptoms,
institutionalisation, death, or loss of capacity**, which can make their
BPRS trajectory systematically different from those who remain under
observation.

### What the slide means by “bias” for CC and LOCF

#### Complete case (CC)

**CC** uses only subjects with complete data (here: only “completers”).
That is equivalent to estimating the follow-up effect **in the selected
subgroup who remain observed**, then treating it as if it represented
the full cohort.

- Under **MCAR** (missingness unrelated to outcomes), completers are
  essentially a random subsample, so CC can be unbiased for many
  targets.
- Under **MAR** (dropout related to observed history), completers are
  **not** exchangeable with the full cohort, so CC is generally biased.
- Under **MNAR** (dropout related to unobserved outcomes), CC can be
  biased even if we adjust for observed history.

This is why the slide shows **Bias(CC) = 0 under MCAR**, but **non-zero
under MAR**.

A useful way to say this in our context:

- If patients with higher observed BPRS (or faster observed
  deterioration) are more likely to miss future visits, then the
  observed follow-up sample will tend to over-represent the “healthier /
  slower progression” subgroup.
- CC will then tend to **underestimate** the average worsening over time
  (or distort group contrasts), because the sickest trajectories are
  progressively removed from the observed sample.

#### Last observation carried forward (LOCF)

**LOCF** fills missing follow-up values using the last observed value
(often baseline):

$$
Y_{i1}^{\text{LOCF}} = Y_{i0}
\quad \text{for subjects missing at } j=1.
$$

LOCF imposes a strong structural assumption: **after dropout the outcome
stays constant**. In progressive conditions (including psychiatric
symptom scales in Alzheimer disease), that is rarely clinically
credible.

The slide’s message is:

- LOCF can be **biased under MCAR** (because it shrinks change toward 0
  by construction).
- LOCF can be **even more biased under MAR**, because it combines
  selection (dropout depends on history) with an implausible
  post-dropout trajectory rule.
- LOCF typically yields confidence intervals that are misleading: they
  can look “precise” because the imputed values are treated as if
  observed, even though they are not.

### How to defend this in 10 minutes (what you say out loud)

A clean defence that fits the Alzheimer BPRS story:

1.  **Our missingness exploration shows increasing dropout over time**,
    which is typical in Alzheimer follow-up.
2.  In this setting, CC and LOCF rely on “simple” assumptions that are
    hard to justify:
    - CC effectively estimates trajectories in those who remain under
      observation.
    - LOCF assumes no change after dropout.
3.  In Alzheimer cohorts, dropout is plausibly linked to disease burden
    and symptom progression, so **MCAR is unlikely**, and methods
    designed for **MAR** are more appropriate for primary inference.
4.  Because **MNAR cannot be ruled out from observed data**, we add a
    targeted sensitivity analysis for BPRS to check robustness to
    departures from MAR.

That is exactly the logic in the slide: “MCAR/simple” options are not
safer; they are often biased and not necessarily simpler in practice
than MAR methods.

### Which sensitivity model to use here (for the continuous BPRS outcome)

Given your assignment wording (“brief sensitivity analysis” for the
continuous outcome) and the Alzheimer context, the most defensible and
reportable option is a **pattern-mixture $\delta$-adjustment**, anchored
to the MAR mixed model.

#### Step 1: Primary MAR model for BPRS (direct likelihood)

Fit a linear mixed model for BPRS under MAR (likelihood-based
inference):

$$
Y_{ij} =
x_{ij}^\top \beta + z_{ij}^\top b_i + \varepsilon_{ij},
\qquad
b_i \sim N(0,D),
\qquad
\varepsilon_{ij} \sim N(0,\sigma^2).
$$

This uses all available observed BPRS values and is valid for inference
under MAR with correctly specified mean/covariance.

#### Step 2: $\delta$-adjusted pattern-mixture sensitivity (MNAR departure)

Let $\widehat{Y}^{MAR}_{ij}$ be the model-based prediction for a missing
BPRS value from the MAR mixed model. For post-dropout times $j>t_i$,
define:

$$
Y_{ij}^{m}(\delta) =
\widehat{Y}^{MAR}_{ij} + \delta,
\qquad j>t_i.
$$

**Interpretation for BPRS:**

- $\delta = 0$ corresponds to MAR.
- $\delta > 0$ represents **worse unobserved BPRS after dropout** than
  MAR would predict (clinically plausible if dropout is driven by
  deterioration).
- $\delta < 0$ represents the opposite direction (often less clinically
  plausible here, but still useful as a bound).

Then you refit (or re-evaluate) the estimand of interest (time trend,
group effect, time-by-group interaction) across a small grid of $\delta$
values and report how conclusions change.

This is a strong choice because it:

- explicitly targets **MNAR risk** without pretending we can “test MAR,”
- stays interpretable for clinicians (“how much worse would dropouts
  have to be for our conclusion to change?”),
- is brief enough for the scope of the assignment.

### 30.21 Illustration: why **CC** and **LOCF** are risky, and why we work under **MAR** for the Alzheimer analyses

This slide is doing one job: it shows that *how you handle missing
follow-up values changes both the point estimate and the uncertainty*,
even in a toy example. In Alzheimer follow-up (BPRS, cognition,
function), missingness is common and often related to observed disease
history, so the “MCAR/simple” fixes can mislead.

#### Notation (generic, but matches our setting)

Let $Y_i$ be an outcome at a given follow-up (continuous BPRS or a
binary indicator), and let

$$
R_i =
\begin{cases}
1, & \text{if } Y_i \text{ is observed}, \\
0, & \text{if } Y_i \text{ is missing}.
\end{cases}
$$

A simple estimand is the marginal mean at that follow-up,

$$
\theta = \mathbb{E}(Y_i).
$$

For a binary outcome, $\theta$ is a probability; for BPRS it is a mean
score.

#### Complete Case (CC): “analyse only those with complete follow-up”

The CC estimator uses only $R_i = 1$ subjects:

$$
\widehat{\theta}_{CC} = \frac{1}{n_{\text{obs}}} \sum_{i:R_i=1} Y_i
= \mathbb{E}(Y_i \mid R_i=1)\ \text{(estimated)}.
$$

So the *target changes* from $\mathbb{E}(Y_i)$ to
$\mathbb{E}(Y_i \mid R_i=1)$. The bias is

$$
\mathbb{E}(\widehat{\theta}_{CC}) - \theta
=
\mathbb{E}(Y_i \mid R_i=1) - \mathbb{E}(Y_i).
$$

- If missingness is **MCAR**, then $R_i \perp Y_i$ and CC can be
  unbiased.
- If missingness depends on observed history (very plausible in
  Alzheimer follow-up), then $R_i$ is associated with $Y_i$ and CC tends
  to be biased.
- Even when bias is modest, CC is **inefficient** because it throws away
  partial trajectories, which is why the slide shows “too wide”
  intervals.

**Alzheimer link:** Patients with worse observed progression often miss
later visits (institutionalisation, death, inability to attend). CC then
over-represents healthier/stabler patients at later visits, typically
pulling estimates toward “less severe” trajectories.

#### Last Observation Carried Forward (LOCF): “freeze the trajectory after missingness”

LOCF replaces a missing follow-up value with the last observed value.
Let $t_i$ be the last observed visit for subject $i$, and define the
imputed value

$$
Y_i^{\text{LOCF}} = Y_{i,t_i}.
$$

Then

$$
\widehat{\theta}_{\text{LOCF}} = \frac{1}{n} \sum_{i=1}^n Y_i^{\text{LOCF}}.
$$

Two problems (exactly what the slide labels):

1.  **Bias**: LOCF implicitly assumes “no change after the last observed
    visit.” In progressive disease settings, that is usually
    implausible. For BPRS, it effectively says symptoms remain flat
    after dropout.

2.  **Too narrow uncertainty**: LOCF treats imputed values as if they
    were observed, so it ignores imputation uncertainty. That’s why the
    slide shows a confidence interval that is “too narrow” even when the
    point estimate is off.

**Alzheimer link:** If BPRS tends to worsen (or fluctuate with
worsening) among those who later drop out, LOCF tends to understate
worsening (or distort the time trend), while also overstating precision.

#### MAR-based analysis: “use a model that conditions on observed history and integrates over what is missing”

Under MAR, missingness can depend on observed information but not on
unobserved outcomes once the observed history is accounted for:

$$
p(R_i \mid Y_i, X_i) = p(R_i \mid Y_i^{o}, X_i).
$$

For the primary analyses, MAR is operationalised using methods that are
valid under MAR and use all available longitudinal information:

### (i) Direct likelihood (mixed models)

For continuous BPRS (LMM), inference is based on the observed-data
likelihood:

$$
L(\theta) = \prod_{i=1}^N \int f(y_i^{o}, y_i^{m} \mid X_i, \theta)\, dy_i^{m}.
$$

This is exactly what the slide’s “MAR: direct likelihood” line is
pointing at: you don’t fill in values by hand; you fit the model and
integrate over the missing part.

### (ii) For binary outcomes: GLMM and marginal (GEE) approaches

- A **random-effects (GLMM)** uses a subject-specific model and is
  likelihood-based (MAR-valid under correct specification).
- A **marginal model (GEE)** targets population-averaged effects; under
  MAR, a common extension is **weighted GEE** (IPW) when missingness is
  associated with observed history.

#### How to defend this in 10 minutes (presentation-friendly)

A clean defence that lands fast:

1.  **What the slide demonstrates:** CC and LOCF are “simple” but they
    change the estimand or bake in unrealistic assumptions, and they can
    distort uncertainty.

2.  **Why it matters here:** In Alzheimer follow-up, missingness is
    rarely plausibly MCAR; it is typically related to observed disease
    history and patient status.

3.  **What we therefore do:** Primary inference is built on MAR-valid
    methods that use observed history efficiently (likelihood-based
    mixed models for BPRS; GLMM and GEE for the binary outcome), and we
    reserve MNAR assumptions for sensitivity analysis rather than
    treating them as truth.

If someone pushes you with “but MAR isn’t testable,” the crisp reply is:
**correct—so MAR is the working assumption for primary inference, and
MNAR is addressed via sensitivity analysis.** That stance is exactly
what the course’s “overview” slide is advocating.

### Ignorability under MAR (why the missingness model can be “ignored” for inference on the outcome model)

#### Notation

For subject $i=1,\dots,N$, let the full longitudinal outcome be split
into observed and missing parts:

$$
y_i = (y_i^{o}, y_i^{m}).
$$

Let $r_i$ denote the missingness pattern/indicator vector for subject
$i$, and let $X_i$ be baseline covariates.

We distinguish two parameter blocks:

- $\theta$: parameters of the **outcome model**
- $\psi$: parameters of the **missingness model**

### 1) Full-data likelihood and selection-model factorisation

Start from the likelihood for the joint distribution of outcomes and
missingness:

$$
L(\theta,\psi \mid y, r, X)
=
\prod_{i=1}^N f(y_i, r_i \mid X_i, \theta, \psi).
$$

Under a **selection model** factorisation:

$$
f(y_i, r_i \mid X_i, \theta, \psi)
=
f(y_i \mid X_i, \theta)
f(r_i \mid X_i, y_i, \psi).
$$

So the likelihood becomes:

$$
L(\theta,\psi \mid y, r, X)
=
\prod_{i=1}^N
f(y_i \mid X_i, \theta)
f(r_i \mid X_i, y_i, \psi).
$$

**Meaning (Alzheimer setting):**  
$f(y_i \mid X_i,\theta)$ describes the BPRS trajectory (mean change over
time, heterogeneity, etc.).  
$f(r_i \mid X_i,y_i,\psi)$ describes the visit attendance/dropout
process, potentially depending on symptom severity.

### 2) Observed-data likelihood (integrating out missing outcomes)

In practice, we do not observe $y_i^{m}$, so the observed-data
likelihood integrates it out:

$$
L(\theta,\psi \mid y^{o}, r, X)
=
\prod_{i=1}^N
\int
f(y_i^{o}, y_i^{m} \mid X_i, \theta)
f(r_i \mid X_i, y_i^{o}, y_i^{m}, \psi)
dy_i^{m}.
$$

This is the key expression: the dropout/attendance model appears
*inside* the integral, because (in general) missingness can depend on
the unobserved outcomes.

### 3) The MAR assumption and factorisation of the likelihood

The **MAR** assumption is:

$$
f(r_i \mid X_i, y_i^{o}, y_i^{m}, \psi)
=
f(r_i \mid X_i, y_i^{o}, \psi).
$$

Substitute MAR into the observed-data likelihood:

$$
L(\theta,\psi \mid y^{o}, r, X)
=
\prod_{i=1}^N
\int
f(y_i^{o}, y_i^{m} \mid X_i, \theta)
f(r_i \mid X_i, y_i^{o}, \psi)
dy_i^{m}.
$$

Because $f(r_i \mid X_i, y_i^{o}, \psi)$ no longer depends on $y_i^{m}$,
it can be taken outside the integral:

$$
L(\theta,\psi \mid y^{o}, r, X)
=
\prod_{i=1}^N
\Bigg[
f(r_i \mid X_i, y_i^{o}, \psi)
\int
f(y_i^{o}, y_i^{m} \mid X_i, \theta)
dy_i^{m}
\Bigg].
$$

The integral is exactly the observed-data density for the outcomes:

$$
\int
f(y_i^{o}, y_i^{m} \mid X_i, \theta)
dy_i^{m}
=
f(y_i^{o} \mid X_i, \theta).
$$

Therefore:

$$
L(\theta,\psi \mid y^{o}, r, X)
=
\prod_{i=1}^N
f(y_i^{o} \mid X_i, \theta)
\prod_{i=1}^N
f(r_i \mid X_i, y_i^{o}, \psi).
$$

### 4) Log-likelihood separation and the score equation

Taking logs:

$$
\ell(\theta,\psi)
=
\sum_{i=1}^N \log f(y_i^{o} \mid X_i, \theta)
+
\sum_{i=1}^N \log f(r_i \mid X_i, y_i^{o}, \psi).
$$

Differentiate with respect to $\theta$ (the outcome-model parameters).
The second term contains no $\theta$, so:

$$
S(\theta)
=
\frac{\partial \ell(\theta,\psi)}{\partial \theta}
=
\sum_{i=1}^N
\frac{\partial}{\partial \theta}
\log f(y_i^{o} \mid X_i, \theta)
+
0.
$$

### 5) What “ignorability” means here

The derivation shows:

- Under **MAR**, the observed-data likelihood factors into an outcome
  part and a missingness part.
- For inference on $\theta$, the missingness part does not contribute to
  the score for $\theta$.

So, for estimating the BPRS model parameters $\theta$, one can base
inference on:

$$
L_{\text{obs}}(\theta)
=
\prod_{i=1}^N
f(y_i^{o} \mid X_i, \theta),
$$

without explicitly specifying $f(r_i \mid X_i, y_i^{o}, \psi)$.

**Interpretation (Alzheimer cohort):**  
If the probability of missing a BPRS measurement at visit $j$ depends on
observed history (previous observed BPRS and covariates), but not on the
unobserved current BPRS once that history is accounted for, then
likelihood-based mixed models for BPRS can be fit directly to the
observed BPRS values, and the dropout model does not need to be modelled
for valid estimation of $\theta$ (subject to correct model
specification).

**Important practical note:**  
Ignorability for likelihood-based inference relies on MAR plus the
standard “separability” idea that $\theta$ and $\psi$ index different
components (outcome vs missingness). When MNAR is plausible, this
separation does not hold in a way that protects $\theta$, which
motivates sensitivity analyses.

### What does it mean that $\theta$ and $\psi$ are “orthogonal” under MAR?

On that slide, “orthogonal” is being used in the
**likelihood/information** sense: the parameters governing the **outcome
model** and the **missingness model** do not “fight each other” when we
estimate the outcome parameters we care about.

Concretely, in missing-data theory we typically have two parameter
blocks:

- $\theta$: parameters of the **measurement/outcome model** (e.g., BPRS
  mean trajectory, variance, random-effects covariance)
- $\psi$: parameters of the **missingness/dropout model** (e.g.,
  logit-hazard coefficients for whether a visit is observed)

Under **MAR** plus **parameter distinctness**, the observed-data
likelihood factorises into a part involving only $\theta$ and a part
involving only $\psi$.

### 1) Parameter distinctness (what the slide calls “$\theta$ and $\psi$ distinct”)

“Distinct” means the joint parameter space splits cleanly:

$$
(\theta, \psi) \in \Theta \times \Psi.
$$

So the outcome model and the missingness model have **separate knobs**:
changing $\theta$ does not restrict what values $\psi$ can take, and
vice versa.

### 2) Orthogonality in the likelihood sense

Let $y_i = (y_i^o, y_i^m)$ and $r_i$ be the missingness indicators.
Under a selection-model factorisation:

$$
f(y_i, r_i \mid X_i, \theta, \psi) =
f(y_i \mid X_i, \theta)
f(r_i \mid X_i, y_i^o, y_i^m, \psi).
$$

Under **MAR**:

$$
f(r_i \mid X_i, y_i^o, y_i^m, \psi) =
f(r_i \mid X_i, y_i^o, \psi).
$$

Then the observed-data likelihood becomes:

$$
L(\theta, \psi \mid y^o, r, X) =
\prod_{i=1}^N
\int
f(y_i^o, y_i^m \mid X_i, \theta)
f(r_i \mid X_i, y_i^o, \psi)
dy_i^m.
$$

Because the missingness term no longer depends on $y_i^m$, we can pull
it out of the integral:

$$
L(\theta, \psi \mid y^o, r, X) =
\prod_{i=1}^N
f(r_i \mid X_i, y_i^o, \psi)
\prod_{i=1}^N
\int
f(y_i^o, y_i^m \mid X_i, \theta)
dy_i^m.
$$

And the integral is just $f(y_i^o \mid X_i, \theta)$, so:

$$
L(\theta, \psi \mid y^o, r, X) =
\prod_{i=1}^N
f(y_i^o \mid X_i, \theta)
\prod_{i=1}^N
f(r_i \mid X_i, y_i^o, \psi).
$$

Taking logs:

$$
\ell(\theta, \psi) =
\ell_y(\theta) + \ell_r(\psi),
$$

where:

$$
\ell_y(\theta) = \sum_{i=1}^N \log f(y_i^o \mid X_i, \theta),
\qquad
\ell_r(\psi) = \sum_{i=1}^N \log f(r_i \mid X_i, y_i^o, \psi).
$$

This additive separation is the key “orthogonality vibe”: inference
about $\theta$ can be done from $\ell_y(\theta)$ alone.

### 3) Orthogonality in the information-matrix sense

A standard formal way to express orthogonality is that the
**cross-information** is zero (block-diagonal information):

$$
I(\theta, \psi) =
\begin{pmatrix}
I_{\theta\theta} & I_{\theta\psi} \\
I_{\psi\theta} & I_{\psi\psi}
\end{pmatrix},
\qquad
I_{\theta\psi} = 0.
$$

Equivalently, the score for $\theta$ does not involve $\psi$:

$$
\frac{\partial \ell(\theta, \psi)}{\partial \theta} =
\frac{\partial \ell_y(\theta)}{\partial \theta}.
$$

So the missingness model contributes nothing to the estimating equations
for $\theta$.

### What this means in the Alzheimer (BPRS) setting

### Outcome model parameters $\theta$

Think of $\theta$ as what you estimate in your BPRS linear mixed model,
e.g.,

$$
y_{ij} =
x_{ij}^\top \beta
+
z_{ij}^\top b_i
+
\varepsilon_{ij},
\qquad
b_i \sim N(0, D),
\qquad
\varepsilon_{ij} \sim N(0, \sigma^2).
$$

Here $\theta$ typically includes $(\beta, D, \sigma^2)$.

### Missingness model parameters $\psi$

Now imagine a dropout/attendance model for whether the next visit is
observed, depending on observed history (MAR-style), e.g.,

$$
\Pr(r_{i,j+1} = 1 \mid X_i, y_{i0}, \dots, y_{ij}) =
\operatorname{logit}^{-1}\Big(
\alpha_0 + \alpha_1 y_{ij} + \alpha_2^\top X_i + \alpha_3 j
\Big),
$$

with $\psi = (\alpha_0, \alpha_1, \alpha_2, \alpha_3)$.

### The practical punchline

If missingness in the Alzheimer cohort is **MAR given observed history**
(e.g., patients with higher observed BPRS at year $j$ are more likely to
miss year $j+1$, and this dependence is captured by observed data), and
if $\theta$ and $\psi$ are **distinct**, then:

- Fitting the BPRS mixed model by ML/REML using the observed BPRS values
  targets $\theta$ correctly.
- You do **not** need to fully specify the dropout mechanism to estimate
  the BPRS trajectory parameters.

In words: **you can model the disease trajectory without having to also
model why someone missed a visit**, as long as the reasons for
missingness are explainable by what you had already observed (plus
covariates). That’s exactly why the slide points you toward “direct
likelihood / Bayes” approaches for MAR.

### One caution worth keeping in your head for Q&A

MAR helps you *ignore the missingness model for estimation of $\theta$*,
but it doesn’t magically guarantee the MAR assumption is true in
Alzheimer follow-up (illness progression and death are the classic
reasons we then add an MNAR sensitivity analysis).

In the missing-data / ignorability setup, **$\theta$** and **$\psi$**
are just “parameter buckets” for two different parts of the model:

### $\theta$: parameters of the **data (outcome) model**

These govern the distribution of the *measurements you care about*,
e.g., cognition score, BPRS, etc.

Typical examples in an Alzheimer longitudinal model:

- Fixed effects ($\beta$): time trend, treatment effect,
  age/sex/education effects, interactions (e.g., time×treatment)
- Random‐effects covariance ($D$): variability in subject-specific
  intercepts/slopes
- Residual variance ($\sigma^2$) (and maybe correlation structure if you
  model it)

So $f(y_i \mid X_i; \theta)$ is like: “given covariates $X_i$, how do
the observed outcomes $y_i$ behave?”

### $\psi$: parameters of the **missingness (response) model**

These govern the distribution of the *missingness indicators* ($R$) —
i.e., who drops out, who misses visits, etc.

Example: let $R_{ij}=1$ if patient $i$ is observed at visit $j$, $0$ if
missing. Then you might model:

$$
\Pr(R_{ij}=1 \mid \text{history}; \psi)
$$

where “history” could include:

- past observed outcomes (e.g., last cognitive score),
- baseline covariates (e.g., age, baseline severity),
- visit/time indicators,
- adverse events, hospitalizations (if recorded).

So $g(R_i \mid y_i, X_i; \psi)$ is like: “what predicts
missingness/dropout?”

### How this ties to MAR + ignorability (the punchline)

Under **MAR** (Missing At Random), the missingness can depend on
**observed information** (past observed outcomes and covariates) but
**not on the unobserved current missing value once you condition on what
you’ve observed**.

A common way it shows up:

$$
p(y_i, R_i \mid X_i) = f(y_i \mid X_i; \theta) \cdot g(R_i \mid y_{i,\text{obs}}, X_i; \psi)
$$

So:

- $\theta$ is about the Alzheimer outcome process,
- $\psi$ is about the visit attendance/dropout process.

And when people say “$\theta$ and $\psi$ are orthogonal” (sometimes
written $\phi$ instead of $\psi$), they mean: **you can estimate
$\theta$ from the outcome model using the observed data likelihood
without having to correctly model $\psi$**, *provided* MAR + other
regularity conditions hold.

In human terms: *if missingness only depends on stuff you already
observed, you can focus on modeling the disease trajectory and still get
valid inference for treatment/time effects.*

On that slide, the **likelihood** is doing a very specific trick:

It **fits the joint model for $(Y_1, Y_2)$** (here: bivariate Normal),
but it only uses what you actually observed:

- For **completers** ($i=1,\dots,R$): you observed $(y_{i1}, y_{i2})$,
  so they contribute the **joint density**:

  $$
  f(y_{i1}, y_{i2} \mid \theta).
  $$

- For **incompleters** ($i=R+1,\dots,N$): you only observed $y_{i1}$
  (because $y_{i2}$ is missing), so they contribute the **marginal
  density** of $Y_1$:

  $$
  f(y_{i1} \mid \theta) = \int f(y_{i1}, y_{i2} \mid \theta) \, dy_{i2}.
  $$

So the **observed-data likelihood** is:

$$
L(\theta) = \prod_{i=1}^R f(y_{i1}, y_{i2} \mid \theta) \cdot \prod_{i=R+1}^N f(y_{i1} \mid \theta).
$$

### What that *means* for estimating $\mu_2$

If you do the naive “frequentist” thing in the table, you estimate:

$$
\hat{\mu}_2^{\text{CC}} = \frac{1}{R} \sum_{i=1}^R y_{i2},
$$

i.e., **use only completers**. That’s only safe if missingness is
basically MCAR (missing completely at random).

The likelihood/MLE instead uses the fact that in a bivariate Normal:

$$
E(Y_2 \mid Y_1 = y_1) = \mu_2 + \frac{\sigma_{21}}{\sigma_{11}}(y_1 - \mu_1),
$$

(which the slide writes as $\beta_0 + \beta_1 y_1$).

So for people missing $y_{i2}$, the likelihood “fills in” their
contribution **not with a guess**, but with the **model-implied
conditional mean** based on their $y_{i1}$:

$$
\hat{\mu}_2^{\text{ML}}
= \frac{1}{N} \left[
\sum_{i=1}^R y_{i2}
+
\sum_{i=R+1}^N \Big[\bar{y}_2 + \hat{\beta}_1 (y_{i1} - \bar{y}_1)\Big]
\right].
$$

That bracket term is basically:

- Start at the average $\bar{y}_2$ (from those with $y_2$ observed),
- Then **adjust up/down** depending on whether the person’s $y_{i1}$ is
  above/below $\bar{y}_1$,
- Scaled by $\hat{\beta}_1$ (the estimated $Y_2$–$Y_1$ relationship).

### Intuition in one sentence

**The likelihood uses everyone’s $Y_1$ values plus the observed
$Y_1$–$Y_2$ relationship to “borrow strength” and correct the bias you
get if incompleters systematically differ in $Y_1$.**

This is essentially what the EM algorithm would do behind the scenes:
**E-step** computes $E(Y_2 \mid Y_1)$ for missing cases; **M-step**
updates $\mu_2$ as the mean of observed $Y_2$ plus those conditional
expectations.

### Why the covariance structure doesn’t affect the means (with complete data)

In that slide, you’re estimating a **“saturated” mean model** on
**complete (balanced) data** — i.e., you’re literally estimating the
mean at each time point (Age 8 and Age 10) with no extra constraints. In
this special setup, the **MLEs of the means are just the sample means**,
and they **do not depend on the covariance model**.

### The key math (why the covariance drops out)

Write each boy’s two measurements as a vector:

$$
y_i =
\begin{pmatrix}
y_{i1} \\
y_{i2}
\end{pmatrix}
\sim N\!\left(
\mu =
\begin{pmatrix}
\mu_1 \\
\mu_2
\end{pmatrix},
\ \Sigma
\right).
$$

With **complete data**, the (negative) log-likelihood has the quadratic
form:

$$
\sum_{i=1}^N (y_i - \mu)^\top \Sigma^{-1} (y_i - \mu) + N\log|\Sigma|.
$$

Different covariance structures (Unstructured, CS, Independence) just
mean different choices for $\Sigma$. Now differentiate with respect to
$\mu$ and set to zero:

$$
\frac{\partial \ell}{\partial \mu}
\propto \sum_{i=1}^N \Sigma^{-1}(y_i - \mu) = 0
\quad\Rightarrow\quad
\Sigma^{-1}\sum_{i=1}^N (y_i - \mu) = 0.
$$

Multiply by $\Sigma$:

$$
\sum_{i=1}^N (y_i - \mu) = 0
\quad\Rightarrow\quad
\hat{\mu} = \frac{1}{N}\sum_{i=1}^N y_i.
$$

Boom: **$\hat{\mu}$ is the sample mean vector**. No $\Sigma$ anywhere.  
So the estimated means at Age 8 and Age 10 stay the same whether you
assume Unstructured, CS, or Independence.

### What *does* change with covariance structure (even when balanced)?

Mainly the **standard errors** and **tests/CI widths**, because those
depend on $\Sigma$ (through the information matrix). So the point
estimates stay fixed, but your uncertainty about them can change.

### Why it changes once the data become unbalanced (missing)

When some boys are missing Age 10, each subject contributes only the
part they observed. The score equations become weighted by the **inverse
of the observed sub-covariance**:

$$
\sum_i \Sigma_{i,\text{obs}}^{-1} (y_{i,\text{obs}} - \mu_{i,\text{obs}}) = 0.
$$

Now those weights **do depend on the covariance assumptions**, and the
model starts “borrowing strength” across time points differently. So
estimates can shift under Unstructured vs CS vs Independence.

### One important caveat

This invariance is **exact** here because the mean model is essentially
“two free means.” If instead you fit a constrained mean structure (e.g.,
a linear time trend, treatment×time interaction, covariates), covariance
assumptions can sometimes influence point estimates even with complete
data (though often not dramatically).

### Summary

That slide is basically shouting: **balanced + saturated mean ⇒
covariance structure doesn’t move the means; missingness breaks that
nice decoupling.**

### Motivation for using a linear mixed model for BPRS (rather than a pure repeated-measures covariance model)

The repeated-measures model written above is a **marginal mean model**
with a flexible within-subject covariance matrix ($\Sigma_i$). This
framework is appropriate when the primary goal is estimation of
population-average mean profiles and when the covariance structure can
be estimated reliably. However, in the Alzheimer cohort, two features
make a linear mixed model (LMM) a more robust primary analysis for BPRS.

#### 1) Unbalanced follow-up and dropout: estimation stability

BPRS is measured at baseline and planned annual visits up to year 6, but
follow-up is typically **unbalanced** because patients miss visits,
withdraw, become institutionalised, or die. Let $n_i$ denote the number
of observed BPRS measurements for subject $i$. In a repeated-measures
covariance model, the likelihood contribution for subject $i$ depends on
the inverse of the observed submatrix ($\Sigma_{i,\text{obs}}$):

$$
\ell_i(\beta, \Sigma) \propto
-\frac{1}{2}\log|\Sigma_{i,\text{obs}}|
-\frac{1}{2}(y_{i,\text{obs}}-\mu_{i,\text{obs}})^\top
\Sigma_{i,\text{obs}}^{-1}
(y_{i,\text{obs}}-\mu_{i,\text{obs}}).
$$

When many subjects have few late follow-up measurements, the most
flexible choices (e.g., unstructured correlation with visit-specific
variances) can become weakly identified, leading to unstable estimation
of $\Sigma$ and, by propagation, unstable inference for $\beta$.

In contrast, an LMM represents within-subject dependence through a
low-dimensional random-effects structure, which remains estimable even
when follow-up is sparse at later years.

#### 2) Parsimonious correlation through random effects matches clinical heterogeneity

Alzheimer progression is heterogeneous: individuals differ in baseline
symptom burden and in their rate of change. The LMM directly represents
this heterogeneity using subject-specific random effects. For example,
with random intercept and slope:

$$
Y_{ij} = x_{ij}^\top\beta + z_{ij}^\top b_i + \varepsilon_{ij},
$$

$$
b_i \sim N(0, D),
\qquad
\varepsilon_{ij} \sim N(0, \sigma^2),
\qquad
b_i \perp \varepsilon_{ij}.
$$

This induces a structured marginal covariance:

$$
\text{Var}(Y_i \mid X_i) = Z_i D Z_i^\top + \sigma^2 I,
$$

which captures the main sources of correlation (shared baseline level
and shared trajectory) without requiring estimation of an unconstrained
$\Sigma_i$ across all visit pairs. This is a natural representation for
BPRS in an Alzheimer cohort where persistent between-subject differences
are expected.

#### 3) Primary inference under MAR: likelihood-based LMM and observed-data likelihood

Given the observed dropout, inference is conducted under a working
Missing At Random (MAR) assumption, where the probability of observing
BPRS at a visit may depend on observed history but not on the unobserved
current value, conditional on the observed information.

Let $r_{ij}$ indicate whether $Y_{ij}$ is observed. MAR can be expressed
as:

$$
f(r_i \mid X_i, y_i^{o}, y_i^{m}) = f(r_i \mid X_i, y_i^{o}).
$$

Under MAR and distinct parameters, likelihood-based inference for
$\theta$ (the outcome model parameters) can be based on the
observed-data likelihood without specifying the full missingness model.
For the LMM this corresponds to maximizing:

$$
L(\theta; y^o) =
\prod_{i=1}^N
\int
f(y_i^{o} \mid b_i, X_i, \theta)
f(b_i \mid \theta)
db_i.
$$

This provides a principled way to use all observed BPRS measurements
while accommodating unbalanced designs. The repeated-measures covariance
model can also be fit under MAR via likelihood, but the practical
performance depends more strongly on the estimability and stability of
$\Sigma$ when late follow-up is sparse.

#### 4) Interpretation: subject-specific trajectories versus purely marginal means

The repeated-measures model targets population-average mean differences.
The LMM additionally yields subject-specific trajectories through the
conditional mean:

$$
E(Y_{ij} \mid b_i, X_i) = x_{ij}^\top\beta + z_{ij}^\top b_i,
$$

which aligns with the scientific setting in which individuals vary in
baseline symptom level and progression rate. This is particularly
relevant for BPRS where baseline severity and progression are expected
to differ between patients.

### Practical role of the repeated-measures model in the analysis plan

The flexible repeated-measures model remains useful as a diagnostic
starting point for exploring mean structure and for checking whether
observed within-subject variability patterns are compatible with the
simplified model chosen for reporting. However, for the primary
longitudinal analysis of BPRS in the presence of dropout and unbalanced
follow-up, the LMM provides a parsimonious dependence structure, stable
estimation, and likelihood-based inference under MAR, while retaining a
clinically interpretable representation of individual trajectories.

### A punchy defense answer (10 seconds)

“We started with a flexible repeated-measures covariance model as a
screening tool, but because follow-up is unbalanced and late visits are
sparse, an unstructured $\Sigma$ can be unstable. The LMM captures the
dependence through random intercepts/slopes, which is parsimonious,
clinically interpretable, and gives valid likelihood inference under MAR
using all observed BPRS data.”

### Inverse-probability weights for monotone dropout (what the slide is doing)

#### Setup and notation

Suppose each patient is scheduled for four follow-up visits (after
baseline). Define:

- $R_{ij} = 1$ if patient $i$ is **observed** at visit $j$,
- $R_{ij} = 0$ if patient $i$ is **missing** at visit $j$.

For monotone dropout, once a 0 appears, all later visits are 0 as well.

A convenient way to model dropout is through **discrete-time dropout
hazards**. Let:

- $p_2$: probability of becoming missing at visit 2, **given observed at
  visit 1**,
- $p_3$: probability of becoming missing at visit 3, **given observed at
  visits 1–2**,
- $p_4$: probability of becoming missing at visit 4, **given observed at
  visits 1–3**.

So $1-p_j$ is the conditional probability of **still being observed** at
visit $j$, given the person has remained in follow-up up to the previous
visit.

This is exactly what the green terms on your slide are multiplying.

#### Pattern probabilities and inverse weights (the table on the slide)

The slide lists monotone patterns, computes their probabilities
($\pi_k$), then takes weights ($w_k = 1/\pi_k$):

| Pattern (last observed visit) | $R_{i1}$ | $R_{i2}$ | $R_{i3}$ | $R_{i4}$ | Pattern probability ($\pi_k$) | Inverse weight ($w_k$)             |
|------------------------------:|:--------:|:--------:|:--------:|:--------:|-------------------------------|------------------------------------|
|                  5 (complete) |    1     |    1     |    1     |    1     | $(1-p_2)(1-p_3)(1-p_4)$       | $\dfrac{1}{(1-p_2)(1-p_3)(1-p_4)}$ |
|                 4 (drop at 4) |    1     |    1     |    1     |    0     | $(1-p_2)(1-p_3)p_4$           | $\dfrac{1}{(1-p_2)(1-p_3)p_4}$     |
|                 3 (drop at 3) |    1     |    1     |    0     |    0     | $(1-p_2)p_3$                  | $\dfrac{1}{(1-p_2)p_3}$            |
|                 2 (drop at 2) |    1     |    0     |    0     |    0     | $p_2$                         | $\dfrac{1}{p_2}$                   |

#### How to read one row (example: “drop at 4”)

To have pattern $1110$, the person must:

1.  stay observed at visit 2: $(1-p_2)$,
2.  stay observed at visit 3: $(1-p_3)$,
3.  then become missing at visit 4: $p_4$.

So:

$$
\pi_4 = (1-p_2)(1-p_3)p_4,
$$

and therefore:

$$
w_4 = \frac{1}{\pi_4}.
$$

#### What these weights are trying to achieve (conceptually)

The “unweighted” observed data at later visits is typically a **selected
subset** of patients (often healthier / less severe / more adherent).

Inverse-probability weighting “rebuilds” the target population by giving
more influence to patients whose observed pattern is **rare**:

- If a patient had a **low probability** of being observed up to a
  certain visit (small $\pi$), then $1/\pi$ is **large**, so they count
  more.
- If a patient had a **high probability** of being observed up to that
  visit (large $\pi$), then $1/\pi$ is smaller.

So weighted analyses aim to correct bias from dropout **under MAR**,
where dropout depends on observed history.

#### Visit-specific (more standard) IPW formulation

In practice, especially for **weighted GEE**, we often weight each
**observed outcome** at visit $j$ by the inverse probability of being
observed at that visit, given the observed past.

Let:

$$
\pi_{ij} = \Pr(R_{ij}=1 \mid \text{observed history up to } j-1),
$$

then the basic weight is:

$$
w_{ij} = \frac{1}{\pi_{ij}}.
$$

For monotone dropout, $\pi_{ij}$ can be written as a product of
“survival” terms. For example, the probability of being observed at
visit 4 is:

$$
\Pr(R_{i4}=1) = (1-p_2)(1-p_3)(1-p_4).
$$

That matches the “complete case” pattern probability in the table.

#### How the $p_j$ are estimated in the Alzheimer/BPRS setting

In your Alzheimer cohort, dropout is rarely random: it can be related to
observed severity/progression (e.g., higher prior BPRS, worsening
function, institutionalisation signals).

A typical discrete-time dropout model is a logistic regression for being
observed at the next visit, conditional on the observed past:

$$
\Pr(R_{i,j}=1 \mid \mathcal{H}_{i,j-1})
=
\operatorname{logit}^{-1}\Big(
\alpha_0
+
\alpha_1 \text{BPRS}_{i,j-1}
+
\alpha_2^\top X_i
+
\alpha_3 j
\Big),
$$

where $\mathcal{H}_{i,j-1}$ includes prior observed outcomes and
baseline covariates.

From this fitted model you get $\widehat{\pi}_{ij}$, and then compute
weights $w_{ij} = 1/\widehat{\pi}_{ij}$ (often with stabilization and
truncation, see below).

#### Two practical details you’ll usually mention (because statisticians will ask)

1.  **Stabilized weights (reduce variance):**  
    Instead of $1/\widehat{\pi}_{ij}$, use:

    $$
    w_{ij}^{\text{stab}}
    =
    \frac{\Pr(R_{ij}=1 \mid \text{baseline only})}{\Pr(R_{ij}=1 \mid \text{history})}.
    $$

    This keeps the same “bias correction” logic but avoids extremely
    large weights as often.

2.  **Truncation (avoid a few subjects dominating):**  
    If some weights are huge, it can inflate variance and make results
    unstable. A common approach is truncation at, say, the 99th
    percentile:

    $$
    w_{ij}^{\text{trunc}} = \min\big(w_{ij}, c\big),
    $$

    for some chosen cap $c$.

#### One-line summary you can use in a methodology section

Inverse-probability weights are constructed from a model for the
probability of being observed at each visit; each observed measurement
is then weighted by the inverse of that probability so that, under MAR,
the weighted estimating equations target the population that would have
been observed in the absence of dropout.

### First-stage MAR dropout model for inverse-probability weights (IPW)

That slide is **not** MNAR sensitivity analysis. It’s the **first-stage
MAR dropout model** you fit to build **inverse-probability weights**
(IPW) for a **weighted GEE** (or related “MAR methods” on your overview
slide).

In the taxonomy from your lecture:

- **Mechanism assumed:** MAR (dropout depends on observed history, like
  the last observed outcome).
- **Method class:** **weighted GEE** (a marginal model approach under
  MAR), sometimes described as a **selection-model factorisation** where
  you model $R$ given observed history, but you’re using it for
  weighting—not for MNAR.
- **Purpose:** correct for bias due to informative dropout **assuming
  MAR**.

This is typically done **early** as part of the main analysis *if* you
choose weighted GEE as your MAR-valid marginal model for the binary
outcome.

### What exactly is being modeled on the slide?

The slide defines a dropout time $D_i$ (the visit where subject $i$
drops). The term:

$$
\Pr(D_i = j \mid D_i \ge j, \cdot)
$$

is a **discrete-time hazard**: the probability of dropping **at time
$j$** given the subject has not dropped before $j$. The model:

$$
\operatorname{logit}\Big[\Pr(D_i = j \mid D_i \ge j, \cdot)\Big]
=
\psi_0
+
\psi_{11} I(\text{GSA}_{i,j-1} = 1)
+
\psi_{12} I(\text{GSA}_{i,j-1} = 2)
+
\psi_{13} I(\text{GSA}_{i,j-1} = 3)
+
\psi_{14} I(\text{GSA}_{i,j-1} = 4)
+
\psi_2 \text{PCA0}_i
+
\psi_3 \text{PF}_i
+
\psi_4 \text{GD}_i.
$$

This says: dropout at visit $j$ depends on **the previous observed
outcome** plus baseline covariates. This is a textbook
**MAR-consistent** idea because it conditions on observed history rather
than unobserved current $Y_{ij}$.

In your Alzheimer/BPRS context, the analogous model would use things
like:

- previous BPRS (or previous binary outcome status),
- baseline markers (age, sex, baseline severity, ADL, living situation),
- time $j$ (visit number),
- possibly interactions if clinically motivated.

### How this connects to weights

From the fitted dropout hazard, you get predicted hazards
$\widehat{p}_{ij}$. Convert hazards to “survival” probabilities (being
still observed up to time $j$). For monotone dropout:

$$
\widehat{\pi}_{ij}
=
\Pr(R_{ij} = 1 \mid \text{history})
=
\prod_{k=2}^{j} \big(1 - \widehat{p}_{ik}\big).
$$

Then the basic inverse probability weight is:

$$
w_{ij} = \frac{1}{\widehat{\pi}_{ij}}.
$$

These weights are then used in **weighted GEE** for your marginal model.

### Where does *sensitivity analysis* start?

Sensitivity analysis begins when you move beyond MAR and ask:

- “What if dropout still depends on the unobserved outcomes?” (MNAR)

Typical MNAR sensitivity tools include:

- **pattern-mixture $\delta$-adjustments** (shift post-dropout
  outcomes),
- **selection models with an MNAR parameter** (include current
  unobserved $Y_{ij}$ in dropout model via a sensitivity parameter),
- **shared-parameter (joint) models** linking outcome and dropout via
  latent random effects.

The slide you posted is **not** doing that: it’s building a dropout
model that depends only on **observed** information, which is exactly
the MAR/IPW setup.

### How to defend this in 10 minutes (fast, clean)

You can say something like:

- We first modeled dropout using a discrete-time logistic hazard
  depending on the **previous observed outcome and baseline
  covariates**.
- This provides estimated observation probabilities and
  inverse-probability weights.
- We used these weights in a weighted GEE to obtain marginal estimates
  valid under **MAR**.
- Since MNAR cannot be tested from observed data, robustness was
  assessed separately via an MNAR sensitivity analysis (e.g.,
  $\delta$-shift / pattern-mixture).

That’s coherent, standard, and matches your lecture’s “MAR methods vs
MNAR sensitivity” overview.

### If you need exact RMarkdown methodology

If your marginal model for the binary outcome is weighted GEE, here’s
how you’d write it:

#### Step 1: Dropout model (discrete-time hazard)

We modeled the probability of dropout at visit $j$ given the subject has
not dropped before $j$ as:

$$
\operatorname{logit}\Big[\Pr(D_i = j \mid D_i \ge j, \mathcal{H}_{i,j-1})\Big]
=
\psi_0
+
\psi_1 \text{BPRS}_{i,j-1}
+
\psi_2^\top X_i
+
\psi_3 j,
$$

where $\mathcal{H}_{i,j-1}$ includes prior observed outcomes and
baseline covariates.

#### Step 2: Observation probabilities and weights

From the fitted dropout model, we computed the probability of being
observed at visit $j$:

$$
\widehat{\pi}_{ij}
=
\prod_{k=2}^{j} \big(1 - \widehat{p}_{ik}\big),
$$

where $\widehat{p}_{ij}$ is the predicted dropout hazard at visit $j$.
The inverse-probability weight for each observed outcome was:

$$
w_{ij} = \frac{1}{\widehat{\pi}_{ij}}.
$$

#### Step 3: Weighted GEE

We fit a weighted GEE for the binary outcome using the weights $w_{ij}$
to account for informative dropout under MAR. The GEE targeted
population-average effects, and the working correlation structure was
\[insert your choice, e.g., exchangeable\].

### Weighted GEE (WGEE) for Alzheimer Data: RMarkdown Implementation

This section explains and implements **Weighted GEE (WGEE)** for your
Alzheimer data, using the same **GEE setup from Homework 2**
(`geepack::geeglm`, `corstr = "ar1"`, robust SEs), but adding
**inverse-probability-of-observation weights** derived from a **logistic
dropout model**.

### Why We Need Weights (Link to Homework 2 GEE)

Your Homework 2 marginal GEE fit was:

``` r
gee_model <- geeglm(
  CDRSBbin ~ time + sex + age + bmi + adl + abpet + taupet + 
    time:abpet + time:taupet,
  id = patid,
  data = alzheimer_long,
  family = binomial(link = "logit"),
  corstr = "ar1",
  std.err = "san.se"
)
```

This **unweighted** GEE uses only the observed responses and solves
estimating equations of the form:

$$
\sum_{i=1}^N D_i^\top V_i^{-1} \big(Y_i - \mu_i(\beta)\big) = 0,
$$

but if the probability of being observed depends on **observed past
outcomes/covariates** (a typical MAR situation in Alzheimer follow-up),
the observed sample at later visits is no longer representative of the
full cohort. WGEE corrects this by re-weighting contributions so that
individuals with high dropout risk “count more.”

### Step 1: Define the Observation/Dropout Indicators

Let $r_{ij}=1$ if the outcome at visit $j$ is observed and $0$
otherwise:

$$
r_{ij} =
\begin{cases}
1, & \text{if } Y_{ij} \text{ is observed}, \\
0, & \text{if } Y_{ij} \text{ is missing}.
\end{cases}
$$

For monotone dropout, define the discrete-time dropout event between $j$
and $j+1$ among those still observed at $j$:

$$
D_{i,j+1} = 1 \quad \Longleftrightarrow \quad r_{ij}=1 \text{ and } r_{i,j+1}=0.
$$

This is exactly what the lecturer’s “logistic regression for dropout
indicator” is doing: a **hazard model** for dropout at the next visit,
conditional on being present now.

### Step 2: Fit a Logistic Dropout (Hazard) Model

A typical MAR-consistent dropout model conditions on **observed
history**, e.g., current outcome and covariates:

$$
\Pr(D_{i,j+1}=1 \mid \mathcal{H}_{ij}) =
\operatorname{logit}^{-1}\Big(
\psi_0
+
\psi_1 Y_{ij}
+
\psi_2^\top X_i
+
\psi_3 j
\Big),
$$

where $\mathcal{H}_{ij}$ denotes the observed history up to visit $j$
(what you already know at visit $j$).

In your Alzheimer setting, this can include things like current/previous
cognition status and biomarkers (because those are observed up to time
$j$ for those still in follow-up).

### Step 3: Convert Predicted Dropout Risks into Inverse-Probability Weights

Let:

$$
\widehat{p}_{i,j+1} = \widehat{\Pr}(D_{i,j+1}=1 \mid \mathcal{H}_{ij}).
$$

Then the estimated probability of **staying observed** from $j$ to $j+1$
is:

$$
\widehat{s}_{i,j+1} = 1 - \widehat{p}_{i,j+1}.
$$

The estimated probability of being observed at visit $j$ is the product
of survival probabilities up to $j$:

$$
\widehat{\pi}_{ij} = \prod_{k=0}^{j-1} \widehat{s}_{i,k+1}.
$$

The (unstabilized) inverse-probability weight is:

$$
w_{ij} = \frac{1}{\widehat{\pi}_{ij}}.
$$

Often, **stabilized weights** are used to reduce variance:

$$
w_{ij}^{\text{stab}} =
\frac{\widehat{\pi}_{ij}^{\text{(num)}}}{\widehat{\pi}_{ij}^{\text{(den)}}},
$$

where the numerator comes from a simpler dropout model (e.g., baseline
covariates only), and the denominator is the full dropout model.

### R Code: Derive Weights and Fit WGEE

Below is the R code to derive weights and fit WGEE, aligned with your
Homework 2 setup:

``` r
library(data.table)
library(geepack)

setDT(alzheimer_long)
setorder(alzheimer_long, patid, time)

# 1) Create a complete subject x time grid to mark missing visits explicitly
times <- sort(unique(alzheimer_long$time))
ids   <- unique(alzheimer_long$patid)

grid <- CJ(patid = ids, time = times)
dat_full <- merge(grid, alzheimer_long, by = c("patid", "time"), all.x = TRUE)

# Observation indicator for the binary outcome used in GEE
dat_full[, R := as.integer(!is.na(CDRSBbin))]

# Next-visit observation indicator (within subject)
setorder(dat_full, patid, time)
dat_full[, R_next := shift(R, type = "lead"), by = patid]

# Dropout event between j and j+1
max_time <- max(dat_full$time, na.rm = TRUE)

dat_full[, dropout_next := fifelse(
  R == 1 & time < max_time & is.na(R_next), NA_integer_,
  fifelse(R == 1 & time < max_time & R_next == 0, 1L,
          fifelse(R == 1 & time < max_time & R_next == 1, 0L, NA_integer_))
)]

# Risk set for the hazard model
risk_dat <- dat_full[R == 1 & time < max_time & !is.na(dropout_next)]

# 2) Denominator model (full): dropout depends on observed history
drop_den <- glm(
  dropout_next ~ time + sex + age + bmi + adl + abpet + taupet + CDRSBbin,
  family = binomial(),
  data = risk_dat
)

risk_dat[, p_drop_den := predict(drop_den, type = "response")]
risk_dat[, p_stay_den := 1 - p_drop_den]

# 3) Numerator model (simpler): stabilize weights
drop_num <- glm(
  dropout_next ~ time + sex + age + bmi,
  family = binomial(),
  data = risk_dat
)

risk_dat[, p_drop_num := predict(drop_num, type = "response")]
risk_dat[, p_stay_num := 1 - p_drop_num]

# 4) Build cumulative probabilities pi_ij for being observed at each visit j
keep_cols <- risk_dat[, .(patid, time, p_stay_den, p_stay_num)]
dat_full <- merge(dat_full, keep_cols, by = c("patid", "time"), all.x = TRUE)

dat_full[is.na(p_stay_den), p_stay_den := 1]
dat_full[is.na(p_stay_num), p_stay_num := 1]

setorder(dat_full, patid, time)
dat_full[, pi_den := cumprod(shift(p_stay_den, fill = 1)), by = patid]
dat_full[, pi_num := cumprod(shift(p_stay_num, fill = 1)), by = patid]

dat_full[, w_stab := pi_num / pi_den]

# Keep weights only where the outcome is observed
w_dat <- dat_full[R == 1, .(patid, time, w_stab)]

# Merge weights back to the analysis dataset
alzheimer_wgee <- merge(alzheimer_long, w_dat, by = c("patid", "time"), all.x = TRUE)

# 5) Fit Weighted GEE
wgee_model <- geeglm(
  CDRSBbin ~ time + sex + age + bmi + adl + abpet + taupet +
    time:abpet + time:taupet,
  id = patid,
  data = alzheimer_wgee,
  family = binomial(link = "logit"),
  corstr = "ar1",
  weights = w_stab,
  std.err = "san.se"
)

# summary(wgee_model)
# tidy_with_ci(wgee_model, exponentiate = TRUE)
```

### What to Report in Your Methodology

- The dropout model is a **MAR model for missingness**, because it
  conditions on **observed history** (e.g., current cognition status,
  biomarkers, covariates).
- The resulting weights estimate $\widehat{\pi}_{ij}$, the probability
  that an observation at visit $j$ is observed.
- WGEE then solves a weighted version of the GEE estimating equations:

$$
\sum_{i=1}^N D_i^\top V_i^{-1} W_i \big(Y_i - \mu_i(\beta)\big) = 0,
$$

where $W_i$ is diagonal with entries $w_{ij}$.

### Where This Fits in Your Workflow

What the lecturer shows (logistic dropout model → weights → WGEE) is
**not** MNAR sensitivity analysis. It is a **primary MAR method** for
marginal inference.

- **GLMM** (likelihood-based) is typically valid under MAR when the mean
  model and random effects structure are correctly specified.
- **Unweighted GEE** is safest under MCAR; under MAR you generally move
  to **WGEE** (or MI).

MNAR sensitivity is a separate step (e.g., $\delta$-shift PMM for BPRS),
where you deliberately vary untestable MNAR assumptions.

### Observation-level vs Subject-level Weights in WGEE (and Why Observation-level is Preferred)

In a longitudinal Alzheimer cohort, the probability that a measurement
is observed typically **changes over follow-up** (later visits are more
likely to be missing) and can depend on **observed history** (e.g.,
prior CDR-SB/BPRS, ADL, biomarkers). Weighted GEE (WGEE) handles this
under **MAR** by reweighting the observed data to represent the full
target cohort.

### 1) Notation for the Missingness / Dropout Process

Let $Y_{ij}$ be the binary outcome (e.g.,
$Y_{ij} = \text{CDRSBbin}_{ij}$) for subject $i$ at visit $j$
($j=0,\dots,J$). Define the observation indicator:

$$
R_{ij} =
\begin{cases}
1, & \text{if } Y_{ij}\ \text{is observed}, \\
0, & \text{if } Y_{ij}\ \text{is missing}.
\end{cases}
$$

For monotone dropout, define “at-risk” (still under follow-up) as:

$$
A_{ij} = \prod_{k=0}^{j-1} R_{ik},
\qquad j \geq 1,
$$

so $A_{ij} = 1$ means the subject has not dropped out before visit $j$.

A standard MAR dropout/observation model uses the observed history up to
$j-1$, denoted $\mathcal{H}_{i,j-1}$ (baseline covariates and previously
observed outcomes/biomarkers), and specifies:

$$
p_{ij} = \Pr(R_{ij}=1 \mid A_{ij}=1,\ \mathcal{H}_{i,j-1}).
$$

### 2) Observation-level Weights

#### Definition (Monotone Dropout)

The probability that subject $i$ is observed **up to and including**
visit $j$ is:

$$
\pi_{ij} = \Pr(R_{i1}=1,\dots,R_{ij}=1 \mid \mathcal{H}_{i0})
= \prod_{k=1}^{j} p_{ik}.
$$

The **inverse probability weight for the measurement at visit $j$** is:

$$
w_{ij} = \frac{1}{\pi_{ij}}
= \frac{1}{\prod_{k=1}^{j} p_{ik}}.
$$

This produces a *different* weight at each time point within a subject
(usually increasing with $j$), which matches the reality that “being
observed at year 6” is rarer than “being observed at year 1.”

### Stabilized Weights (Often Better Behaved)

To reduce extreme weights, use stabilized weights:

$$
sw_{ij} =
\frac{\prod_{k=1}^{j} p^{(N)}_{ik}}
{\prod_{k=1}^{j} p^{(D)}_{ik}},
$$

where:

- $p^{(D)}_{ik}$ is from a richer denominator model (history +
  covariates),
- $p^{(N)}_{ik}$ is from a simpler numerator model (often baseline
  covariates + time only).

### 3) Subject-level Weights (and Why They’re Less Precise)

A **subject-level** weight assigns one constant weight to the entire
subject, commonly based on the probability of being observed through the
last planned visit:

$$
w_i = \frac{1}{\pi_{iJ}}
= \frac{1}{\prod_{k=1}^{J} p_{ik}}.
$$

Then *every* observed row for subject $i$ is weighted by $w_i$.

### Why This is Less Precise in Your Setting:

- It treats a subject who drops out at year 2 and one who drops out at
  year 5 too similarly (both get a single “overall” weight).
- It cannot adapt to visit-specific selection, which is exactly what
  happens in Alzheimer follow-up.
- It typically inflates variance (less efficient) because it uses a
  coarser correction than necessary.

That’s why your lecturer is pushing **observation-level weights**: they
correct selection **at the level where missingness happens (each
visit)**.

### 4) How Weights Enter WGEE (What Changes vs Your Current GEE)

Your current (unweighted) GEE solves an estimating equation of the form:

$$
U(\beta) = \sum_{i=1}^N D_i^\top V_i^{-1}(Y_i - \mu_i) = 0,
$$

where $D_i = \partial \mu_i / \partial \beta^\top$ and $V_i$ is the
working covariance.

WGEE replaces this with a weighted version (one common form):

$$
U_w(\beta) = \sum_{i=1}^N D_i^\top V_i^{-1} W_i (Y_i - \mu_i) = 0,
$$

where $W_i$ is diagonal with entries $w_{ij}$ (or $sw_{ij}$). This
restores validity under **MAR**, provided the weight model for $p_{ij}$
is correctly specified (or at least captures the key predictors of
missingness).

### 5) Correlation Structures in GEE (Quickly, for Your Write-up)

- **Exchangeable Correlation:** Constant within-subject correlation
  (another name: **equicorrelated**; in Gaussian settings it corresponds
  to **compound symmetry (CS)**):

  $$
  \text{Corr}(Y_{ij},Y_{ik}) = \rho \quad \text{for all } j \neq k.
  $$

- **AR(1) Correlation:** Correlation decays with time separation:

  $$
  \text{Corr}(Y_{ij},Y_{ik}) = \rho^{|j-k|}.
  $$

- **Independence:** Working independence (often OK if you rely on robust
  SEs, but can be inefficient).

### R Code: Deriving Observation-level Weights and Fitting WGEE

Below is a template that matches your variables (`patid`, `time`,
`CDRSBbin`, `sex`, `age`, `bmi`, `adl`, `abpet`, `taupet`) and your
existing `geeglm()` call.

``` r
library(data.table)
library(geepack)

setDT(alzheimer_long)
setorder(alzheimer_long, patid, time)

# 1) Create a complete person-period grid to represent missing visits explicitly
times <- sort(unique(alzheimer_long$time))

full_dt <- CJ(patid = unique(alzheimer_long$patid), time = times)
full_dt <- merge(full_dt, alzheimer_long,
                 by = c("patid", "time"),
                 all.x = TRUE)

setorder(full_dt, patid, time)

# 2) Define observation indicator for the outcome
full_dt[, R := as.integer(!is.na(CDRSBbin))]

# Optional: at-risk indicator for monotone dropout
full_dt[, A := cumprod(shift(R, 1, fill = 1L)), by = patid]

# 3) Build lagged predictors (observed history)
full_dt[, lag_y      := shift(CDRSBbin, 1), by = patid]
full_dt[, lag_abpet  := shift(abpet,    1), by = patid]
full_dt[, lag_taupet := shift(taupet,   1), by = patid]

# Keep rows that are "at risk" of being observed (monotone dropout case)
wdat <- full_dt[A == 1 & time > min(times)]

# 4) Denominator model for p_ij = P(R_ij = 1 | history)
den_fit <- glm(
  R ~ time + sex + age + bmi + adl + lag_y + lag_abpet + lag_taupet,
  family = binomial(link = "logit"),
  data = wdat
)

wdat[, p_den := predict(den_fit, type = "response")]

# 5) Numerator model for stabilized weights
num_fit <- glm(
  R ~ time + sex + age + bmi + adl,
  family = binomial(link = "logit"),
  data = wdat
)

wdat[, p_num := predict(num_fit, type = "response")]

# Put predictions back into full_dt
full_dt[wdat, `:=`(p_den = i.p_den, p_num = i.p_num), on = .(patid, time)]

# For baseline or non-at-risk rows, set probabilities to 1
full_dt[is.na(p_den), p_den := 1]
full_dt[is.na(p_num), p_num := 1]

# 6) Compute observation-level stabilized weights:
full_dt[, pi_den := cumprod(p_den), by = patid]
full_dt[, pi_num := cumprod(p_num), by = patid]

full_dt[, sw := pi_num / pi_den]

# Optional: truncate extreme weights
cap <- quantile(full_dt$sw[is.finite(full_dt$sw)], probs = 0.99, na.rm = TRUE)
full_dt[, sw_trunc := pmin(sw, cap)]

# 7) Fit WGEE using only observed outcomes
wgee_model <- geeglm(
  CDRSBbin ~ time + sex + age + bmi +
    adl + abpet + taupet +
    time:abpet + time:taupet,
  id = patid,
  data = full_dt,
  family = binomial(link = "logit"),
  corstr = "ar1",
  weights = sw_trunc,
  std.err = "san.se"
)

summary(wgee_model)
```

### What is happening in this **PROC GEE WGEE** slide?

This slide shows a **weighted GEE (WGEE)** analysis for a **binary
longitudinal outcome** under **MAR dropout**, where the weights come
from a **logistic model for being observed** at each visit.

### 1) The data step: creating predictors for the missingness (weight) model

``` sas
data help;
  set armdwgee;
  by subject;

  prevbindif = lag(bindif);
  if first.id then prevbindif = 1;

  time2 = 0; if time = 2 then time2 = 1;
  time3 = 0; if time = 3 then time3 = 1;
run;
```

### What these lines are doing:

- `prevbindif = lag(bindif);`  
  Creates the **previous-visit outcome** ($Y_{i,j-1}$). This is used
  because MAR dropout is often driven by **observed history**,
  especially the previous response.

- `if first.id then prevbindif = 1;`  
  Sets the baseline “previous outcome” to a fixed value for the first
  record per subject (since there is no prior visit). The specific value
  is just a coding convenience; the goal is to define something at
  baseline so the model runs.

- `time2`, `time3`  
  Create **time indicators** so the dropout/observation model can differ
  by visit. This allows the model to capture different missingness risks
  at visit 2 and visit 3 compared to visit 1.

In your Alzheimer setting, this is the same idea as using **observed
history** (e.g., previous symptom status or previous BPRS level) plus
baseline covariates to explain who remains observed.

### 2) The PROC GEE block: outcome model + correlation + missingness model (weights)

``` sas
proc gee data=help;
  class time treat subject lesion;

  model bindif = time treat*time / noint dist=binomial;

  repeated subject=subject / withinsubject=time type=exch corrw modelse;

  missmodel prevbindif treat lesion time2 time3 / type=obslevel;
run;
```

### A) The **outcome (mean) model**

`model bindif = time treat*time / noint dist=binomial;`

- `bindif` is the **binary outcome** ($Y_{ij} \in \{0,1\}$).
- `dist=binomial` specifies a **logistic mean model**.

Mathematically:

$$
Y_{ij} \mid X_{ij} \sim \text{Bernoulli}(\mu_{ij}),
$$

$$
\text{logit}(\mu_{ij}) = X_{ij}^\top \beta.
$$

Here:

- `time` is treated as a factor (class variable).
- `treat*time` allows **treatment effects to vary by visit**.
- `/ noint` excludes a separate intercept, so the `time` indicators act
  like visit-specific intercepts (a “cell means” parameterization).

### B) The **working correlation** (GEE “repeated” statement)

`repeated subject=subject / withinsubject=time type=exch ...;`

This specifies that repeated measures within the same `subject` are
correlated, using a **working correlation structure**.

- `type=exch`: **exchangeable correlation**, also called **compound
  symmetry (CS)** or **constant correlation**.

Mathematically:

$$
\text{Corr}(Y_{ij}, Y_{ik}) = \rho \quad \text{for all } j \neq k.
$$

Interpretation: any two observations from the same subject have the
**same correlation**, regardless of time separation.

- `corrw`: Requests printing/estimating the working correlation.
- `modelse`: Requests **model-based SEs** (naïve SEs). The default GEE
  SE is the robust “sandwich” SE, which is typically preferred for
  reporting.

### 3) The key WGEE piece: the `MISSModel` statement (weights)

`missmodel prevbindif treat lesion time2 time3 / type=obslevel;`

This is the heart of **Weighted GEE** under MAR.

### A) What is being modeled?

Define a missingness/observation indicator:

$$
R_{ij} =
\begin{cases}
1, & \text{if } Y_{ij} \text{ is observed}, \\
0, & \text{if } Y_{ij} \text{ is missing}.
\end{cases}
$$

WGEE fits a model for the probability of being observed, conditional on
**observed history** and baseline covariates:

$$
\pi_{ij} = \Pr(R_{ij}=1 \mid \mathcal{H}_{i,j-1}),
$$

where $\mathcal{H}_{i,j-1}$ includes observed history up to visit $j-1$.
A typical logistic version (matching the code) is:

$$
\text{logit}(\pi_{ij}) =
\gamma_0
+ \gamma_1 Y_{i,j-1}
+ \gamma_2^\top \text{treat}_i
+ \gamma_3^\top \text{lesion}_i
+ \gamma_4 I(j=2)
+ \gamma_5 I(j=3).
$$

Here:

- $Y_{i,j-1}$ corresponds to `prevbindif`.
- $I(j=2), I(j=3)$ correspond to `time2`, `time3`.

This is **not** MNAR modeling. It’s a **MAR-consistent** approach:
explaining observation using what is already observed.

### B) How are weights formed?

Once $\widehat{\pi}_{ij}$ is estimated, the inverse-probability weight
is:

$$
w_{ij} = \frac{1}{\widehat{\pi}_{ij}}.
$$

Intuition: If an observation is “rare to be observed” (small
$\widehat{\pi}_{ij}$), it gets a **large weight**, because it represents
many similar observations that were lost to missingness.

### 4) Why `type=obslevel` (observation-level) weights?

Your lecturer is correct to emphasize observation-level weights for
longitudinal dropout.

### Observation-level weights

A **different weight for each visit** within a subject:

$$
w_{ij} = \frac{1}{\widehat{\pi}_{ij}} \quad \text{(can change over } j\text{)}.
$$

This is appropriate when the probability of being observed depends on
**time-varying history** (e.g., previous outcome), which is exactly what
`prevbindif` is modeling.

### Subject-level weights

A **single weight per subject**, typically based on the probability of
being fully observed to the end. This is simpler but less precise
because it ignores visit-specific selection.

For Alzheimer data, where missingness increases with follow-up time and
depends on evolving clinical status, observation-level weights are more
appropriate.

### 5) Where this fits in your workflow (so you don’t mislabel it)

- Fitting a dropout/observation model like this is **part of a MAR-valid
  approach** for **marginal models**.
- It is **not sensitivity analysis** by itself.
- Sensitivity analysis starts when you deliberately move beyond MAR
  (e.g., $\delta$-adjustment PMM, selection/shared-parameter with MNAR
  parameters).

### Summary

**Logistic model for observation + WGEE = MAR-based primary analysis for
the marginal model.**

### Weighted GEE (WGEE) for Alzheimer Data: Observation-Level Weights

This section explains and implements **Weighted GEE (WGEE)** for your
Alzheimer data, using **observation-level inverse-probability weights**
derived from a **logistic dropout model**. The approach aligns with your
Homework 2 GEE setup (`geepack::geeglm`, `corstr = "ar1"`, robust SEs).

### 1) Outcome and Covariates (Your Homework 2 Model)

For subject $i = 1, \dots, N$ at visit $t = 0, 1, \dots, J$, define the
binary outcome:

$$
Y_{it} =
\begin{cases}
1, & \text{if } \texttt{CDRSBbin}_{it} = 1 \ (\text{e.g., CDRSB} > 10), \\
0, & \text{otherwise.}
\end{cases}
$$

Let the covariate vector match your Homework 2 GEE:

$$
X_{it} =
\Big(
t,\ \text{sex}_i,\ \text{age}_i,\ \text{bmi}_i,\ \text{adl}_i,\ \text{abpet}_{it},\ \text{taupet}_{it},\ t \cdot \text{abpet}_{it},\ t \cdot \text{taupet}_{it}
\Big).
$$

The marginal mean model is:

$$
\text{logit}\big(\Pr(Y_{it}=1 \mid X_{it})\big) = X_{it}^\top \beta.
$$

### 2) Observation Indicator and Probability of Being Observed

Define the observation indicator:

$$
R_{it} =
\begin{cases}
1, & \text{if } Y_{it} \text{ is observed}, \\
0, & \text{if } Y_{it} \text{ is missing}.
\end{cases}
$$

For **WGEE under MAR**, we model the probability of being observed at
each visit, conditional on **observed history**. This ensures that
subjects/visits with a higher likelihood of missingness are
appropriately upweighted.

#### Conditional Probability of Being Observed

Let $\mathcal{H}_{i,t-1}$ denote the observed history up to visit $t-1$
(e.g., $Y_{i,t-1}$, biomarkers at $t-1$, baseline covariates, and time).
The conditional probability of being observed at visit $t$ is modeled
as:

$$
\Pr(R_{it}=1 \mid R_{i,t-1}=1,\ \mathcal{H}_{i,t-1}) =
\text{logit}^{-1}\Big(
\psi_0
+ \psi_1 Y_{i,t-1}
+ \psi_2^\top X_i
+ \psi_3 t
+ \psi_4 \text{abpet}_{i,t-1}
+ \psi_5 \text{taupet}_{i,t-1}
\Big).
$$

Let the fitted conditional probability be:

$$
\hat{p}_{it} = \widehat{\Pr}(R_{it}=1 \mid R_{i,t-1}=1,\ \mathcal{H}_{i,t-1}).
$$

#### Cumulative Probability of Being Observed

The cumulative probability of being observed up to visit $t$ is:

$$
\hat{\pi}_{it} = \prod_{s=1}^t \hat{p}_{is}, \quad \hat{\pi}_{i0} = 1.
$$

#### Observation-Level Weights

The **inverse-probability weight** for visit $t$ is:

$$
w_{it} = \frac{1}{\hat{\pi}_{it}}.
$$

To reduce extreme weights, we use **stabilized weights**:

$$
w_{it}^{\text{stab}} =
\frac{\hat{\pi}_{it}^{\text{num}}}{\hat{\pi}_{it}^{\text{den}}},
$$

where:

- $\hat{\pi}_{it}^{\text{num}}$: Cumulative probability from a simpler
  model (e.g., baseline covariates + time).
- $\hat{\pi}_{it}^{\text{den}}$: Cumulative probability from the full
  history model.

### 3) R Code: Estimate Weights and Fit WGEE

Below is the R implementation for deriving observation-level weights and
fitting WGEE.

``` r
library(data.table)
library(geepack)

# Assume your long data is called alzheimer_long and has:
# patid, time, CDRSBbin, sex, age, bmi, adl, abpet, taupet

dt <- as.data.table(alzheimer_long)
setorder(dt, patid, time)

# 1) Observation indicator
dt[, R := as.integer(!is.na(CDRSBbin))]

# 2) Build lagged (observed-history) predictors
dt[, R_lag := shift(R, 1L), by = patid]
dt[, Y_lag := shift(CDRSBbin, 1L), by = patid]
dt[, abpet_lag := shift(abpet, 1L), by = patid]
dt[, taupet_lag := shift(taupet, 1L), by = patid]

# At baseline (time == 0), set "at risk" = 0 (no prior visit)
dt[time == 0, `:=`(R_lag = NA_integer_, Y_lag = NA_integer_,
                   abpet_lag = NA_real_, taupet_lag = NA_real_)]

# 3) Define the risk set for the conditional model:
#    only model R_it among those observed at previous visit (R_{i,t-1} = 1)
dt[, at_risk := as.integer(!is.na(R_lag) & (R_lag == 1L) & (time > 0))]

# Keep rows where the conditional probability is defined
wt_dat <- dt[at_risk == 1]

# 4) Denominator model (uses observed history) for P(R_it = 1 | history)
fit_den <- glm(
  R ~ factor(time) + sex + age + bmi + adl + Y_lag + abpet_lag + taupet_lag,
  data = wt_dat,
  family = binomial()
)

wt_dat[, p_den := predict(fit_den, type = "response")]

# 5) Numerator model (baseline-only + time), for stabilized weights
fit_num <- glm(
  R ~ factor(time) + sex + age + bmi + adl,
  data = wt_dat,
  family = binomial()
)

wt_dat[, p_num := predict(fit_num, type = "response")]

# 6) Merge predicted probabilities back
dt[at_risk == 1, `:=`(p_den = wt_dat$p_den, p_num = wt_dat$p_num)]

# For non-risk rows, set p_* = 1 so cumulative products stay unchanged
dt[is.na(p_den), p_den := 1]
dt[is.na(p_num), p_num := 1]

# 7) Cumulative probabilities pi_it and stabilized weights
dt[, pi_den := cumprod(p_den), by = patid]
dt[, pi_num := cumprod(p_num), by = patid]

dt[, w_stab := pi_num / pi_den]

# We only use weights for observed outcomes
analysis_dat <- dt[R == 1]

# (Optional but recommended) cap extreme weights
cap <- quantile(analysis_dat$w_stab, probs = 0.99, na.rm = TRUE)
analysis_dat[w_stab > cap, w_stab := cap]

# 8) Fit WGEE using the same mean model you used before
wgee_model <- geeglm(
  CDRSBbin ~ time + sex + age + bmi +
    adl + abpet + taupet +
    time:abpet + time:taupet,
  id = patid,
  data = analysis_dat,
  family = binomial(link = "logit"),
  corstr = "ar1",
  weights = w_stab,
  std.err = "san.se"
)

summary(wgee_model)
```

### Multiple Imputation (MI) and the Impact of Heavy Missingness

### Rubin’s Rules for Combining Estimates

With $M$ imputations, the combined estimate of the parameter of interest
is:

$$
\hat\beta^{*} = \frac{1}{M} \sum_{m=1}^{M} \hat\beta^{(m)}.
$$

The total variance is:

$$
V = W + \left(\frac{M+1}{M}\right)B,
$$

where:

- $W$: Within-imputation variance (average of the variances from each
  imputed dataset),

  $$
  W = \frac{1}{M} \sum_{m=1}^{M} U^{(m)},
  $$

- $B$: Between-imputation variance (how much the estimates
  $\hat\beta^{(m)}$ differ across imputations),

  $$
  B = \frac{1}{M-1} \sum_{m=1}^{M} \left(\hat\beta^{(m)} - \hat\beta^{*}\right)
  \left(\hat\beta^{(m)} - \hat\beta^{*}\right)^{\top}.
  $$

### What Happens When Missing Values Are Many?

#### 1) Between-Imputation Variance ($B$) Increases

When a large proportion of values are missing, the imputer has less real
information and must rely more on model-based guesses. This leads to
greater variability in the completed datasets, causing the estimates
$\hat\beta^{(m)}$ to differ more across imputations. As a result:

- The uncertainty about the missing data increases.
- The between-imputation variance ($B$) rises.

#### 2) Within-Imputation Variance ($W$) Often Increases

Each imputed dataset is analyzed as if it were complete, but the imputed
values behave like noisy measurements because they are drawn from a
predictive distribution. In many cases:

- Effective information is reduced.
- Standard errors within each dataset increase.
- The average variance $U^{(m)}$ (and thus $W$) increases.

However, $W$ does not always increase dramatically, as each dataset is
treated as “complete” after imputation. The primary penalty for
missingness is often reflected in $B$.

#### 3) Total Variance ($V$) Increases

The total variance is:

$$
V = W + \left(\frac{M+1}{M}\right)B.
$$

As $B$ increases with heavy missingness, the total variance $V$ grows,
leading to wider confidence intervals and larger $p$-values. This
weakens inference.

### Intuition for the Alzheimer Context (BPRS / Longitudinal Data)

In the Alzheimer cohort, later follow-up years often have fewer observed
outcomes due to dropout or illness progression. As the proportion of
missing BPRS values increases in later years, imputations rely more on
model-based predictions and less on direct observations. This increases
the variability of imputed trajectories across imputations, inflating
the between-imputation variance ($B$) and, consequently, the total
variance ($V$) used for inference.

### Diagnostic: Fraction of Missing Information (FMI)

A useful diagnostic for assessing the impact of missingness is the
**fraction of missing information (FMI)**. One component of FMI is the
relative increase in variance:

$$
r = \left(1 + \frac{1}{M}\right) \frac{\text{tr}(B)}{\text{tr}(W)}.
$$

When missingness is heavy:

- The ratio $B/W$ increases.
- The relative increase in variance ($r$) grows.
- FMI rises, leading to wider confidence intervals.

### Summary

As missingness increases:

- The between-imputation variance ($B$) grows because imputations
  disagree more.
- The total variance ($V$) increases, making inference less precise.
- Even though each completed dataset appears “complete,” the uncertainty
  from missing data propagates through the imputation process.

1.  **Simple ad-hoc fixes** (often taught first, often problematic),
2.  **Model-based methods under MAR** (including **multiple
    imputation**),
3.  **MNAR tools** (mostly for **sensitivity analysis**, not for “the
    one true fit”).

### Summary of Imputation Approaches (from the Class Notes)

#### 1) “MCAR/Simple” Approaches (Risky Defaults)

These are the **simple** options (e.g., **CC** and **LOCF**) but come
with clear warnings (bias/inefficiency).

**(a) Complete Case (CC)**  
Use only patients with no missing values for the outcome history being
analyzed.

- Works cleanly only in very special situations (e.g., missingness
  unrelated to outcomes after conditioning appropriately).
- In longitudinal Alzheimer data, CC can quietly shift the sample toward
  “healthier completers,” as later follow-ups are often missing due to
  worsening, institutionalization, or death.

**(b) Last Observation Carried Forward (LOCF)**  
Replace a missing value by the last observed value.

- In progressive diseases, LOCF builds in a *“no further change”*
  assumption after dropout, which is usually clinically implausible.
- LOCF is **not automatically conservative**; it can bias effects and
  distort uncertainty.

**BPRS Translation:**  
Because BPRS evolves over time and missingness increases at later years,
methods that freeze trajectories after dropout (LOCF) or discard partial
histories (CC) are not aligned with the data-generation story in an
Alzheimer follow-up cohort.

#### 2) Single Imputation (Fills Once, Then Pretends It Was Observed)

Single imputation is conceptually straightforward but statistically
fragile because it typically **understates uncertainty**: you plug in
one guess and proceed as if it were truth.

**Common Single-Imputation Methods:**

- **Mean/Median Imputation:** Crude; distorts variance and
  relationships.
- **Regression Imputation:** Predict missing values from other
  variables; over-smooths because predictions sit too neatly on the
  regression surface.
- **Hot Deck:** Replace with an observed value from a “similar” donor;
  better for preserving realistic values but requires careful donor
  definition.
- **Predictive Mean Matching (PMM):** Realistic-value imputation, but
  becomes principled only in the **multiple imputation** framework.

**BPRS Translation:**  
A single “best guess” for missing BPRS at year 5 is rarely believable as
*the* value. What matters is acknowledging uncertainty about the missing
value and propagating it into standard errors and conclusions.

#### 3) Multiple Imputation (MI): The Main Framework (MAR)

The MI workflow:

1.  **Impute** missing values **M** times to create **M** completed
    datasets.
2.  **Analyze** each completed dataset using the chosen model.
3.  **Pool** estimates and uncertainty across the **M** analyses.

##### (a) MI Methods Depend on the Missingness Pattern

**Monotone Missingness:** Common with dropout (e.g., patients drop out
and remain missing).  
For monotone patterns, the notes suggest:

- **Regression** for continuous variables
- **Logistic Regression** for binary variables
- **Discriminant Methods** for categorical variables
- **Regression + Predictive Mean Matching (REGPMM)**
- **Propensity Score Methods**

**Non-Monotone Missingness:** Missingness can appear and disappear
across time/variables.  
The notes suggest:

- **MCMC Imputation:** Multivariate model + iterative simulation.
- **FCS/Chained Equations:** Specify a conditional model for each
  variable with missingness and iterate through them.

**BPRS Translation:**  
- If BPRS missingness is mainly dropout, monotone MI methods are a
natural match.  
- If there’s a mix of intermittent missing items + dropout, **FCS** is
the flexible default because the dataset contains a mix of variable
types (continuous scores, binary indicators, categorical factors).

### Pooling Rules (Within vs Between)

With **M** imputations, let $\hat{\beta}^{(m)}$ be the estimate from
imputed dataset $m$, and let $U^{(m)}$ be its estimated variance.

**Pooled Estimate:**

$$
\hat{\beta}^{\ast} = \frac{1}{M} \sum_{m=1}^{M} \hat{\beta}^{(m)}.
$$

**Within-Imputation Variance:**

$$
W = \frac{1}{M} \sum_{m=1}^{M} U^{(m)}.
$$

**Between-Imputation Variance:**

$$
B = \frac{1}{M-1} \sum_{m=1}^{M} \left(\hat{\beta}^{(m)} - \hat{\beta}^{\ast}\right) \left(\hat{\beta}^{(m)} - \hat{\beta}^{\ast}\right)^{\top}.
$$

**Total Variance:**

$$
V = W + \left(\frac{M+1}{M}\right)B.
$$

### What Happens to Variance When Missingness is High?

- **Many Missing Values:** Less observed signal pushes **$W$** upward
  (less precision within each dataset).
- **More Imputation Influence:** Imputed datasets differ more,
  increasing **$B$**.

**Practical Implication:**

$$
V = W + \left(\frac{M+1}{M}\right)B \quad \Rightarrow \quad V \text{ increases as missingness increases.}
$$

**BPRS Translation:**  
If later-year BPRS is frequently missing, uncertainty includes both
**sampling variability** ($W$) and **imputation uncertainty** ($B$).
With heavier dropout, $B$ becomes a non-trivial part of the final
standard errors.

### One Clean Sentence to Link This to Your BPRS Methodology

In a cohort with progressive disease and dropout, multiple imputation
under MAR provides a principled way to replace missing BPRS values while
carrying forward both sources of uncertainty—model-fit uncertainty
within each imputed dataset ($W$) and imputation-to-imputation
variability ($B$)—into the final inference through Rubin’s pooling
rules.

### Choosing an appropriate method for longitudinal data with missingness

This section summarizes when different modeling strategies are
appropriate in longitudinal studies with incomplete data, and how they
relate to the missing-data mechanism. The discussion is framed around
repeated outcomes such as yearly BPRS measurements in an Alzheimer
cohort.

### Generalized Estimating Equations (GEE)

GEE targets **population-average effects** and treats the within-subject
correlation as a nuisance, specified through a *working correlation
structure*.

#### When standard GEE is appropriate

Standard (unweighted) GEE is appropriate when:

- The outcome is **fully observed** (complete data), or
- Missingness is **MCAR** (Missing Completely At Random).

Under MCAR, the estimating equations remain unbiased because the
observed data form a random subsample of the full data.

#### Limitations under MAR

When missingness is **MAR**, standard GEE is **not valid** in general:

- GEE conditions only on the mean model  
  $$E(Y_{ij} \mid X_{ij}) = \mu_{ij},$$  
  and does **not model the outcome distribution**.
- If dropout depends on observed history (e.g., earlier BPRS values),
  the estimating equations are biased unless this dependence is
  explicitly corrected.

Hence, under MAR, **plain GEE is insufficient**.

### Weighted GEE (WGEE) under MAR

To extend GEE to MAR settings, **inverse probability weights** are
introduced.

#### Idea

Each observed measurement is weighted by the inverse probability of
being observed:

$$
w_{ij} = \frac{1}{\hat{\pi}_{ij}},
$$

where

$$
\pi_{ij} = P(R_{ij} = 1 \mid \text{observed history}).
$$

The probabilities $\pi_{ij}$ are typically estimated using a **logistic
regression for dropout**, fitted to observed data only.

#### Why WGEE works under MAR

- If missingness depends only on observed history, reweighting
  reconstructs the population that would have been observed without
  dropout.
- The weighted estimating equations are unbiased for population-average
  effects.

#### Observation-level vs subject-level weights

The lecturer explicitly recommends **observation-level weights**:

- **Subject-level weights**: One weight per individual; simpler but less
  precise.
- **Observation-level weights**: A separate weight for each $(i,j)$;
  better reflects time-varying dropout risk.

In longitudinal Alzheimer data, dropout risk typically **changes over
time** and depends on evolving symptoms, so observation-level weights
are more appropriate.

### Mixed Models (LMM / GLMM)

Mixed models specify a **full likelihood** for the outcome process,
including subject-specific random effects.

#### Key property under MAR

Under MAR and **parameter distinctness**, the likelihood factorizes as:

$$
f(Y^{o}, R \mid X, \theta, \psi) = f(Y^{o} \mid X, \theta) f(R \mid Y^{o}, X, \psi),
$$

so inference for $\theta$ can be based on the observed-data likelihood
alone.

#### Implication

- **LMM (Gaussian outcomes)** and **GLMM (non-Gaussian outcomes)** are
  **valid under MAR without modeling the missingness process
  explicitly**.
- This is why your lecturer states that **under MAR, GLMM is OK**.

#### When mixed models are particularly attractive

- Continuous outcomes (e.g., BPRS total score).
- Interest in **subject-specific trajectories**.
- Unbalanced data with dropout.
- Correct specification of the mean structure is more important than the
  random-effects covariance.

### Multiple Imputation (MI)

Multiple imputation replaces each missing value with multiple plausible
draws, creating $M$ completed datasets.

#### When MI is appropriate

MI is appropriate when:

- Missingness is **MAR**.
- There are **incomplete covariates** (which GEE and GLMM cannot
  automatically handle).
- One wants a **general framework** applicable across different models.

#### Workflow

1.  Impute missing values $M$ times using models compatible with the
    data.
2.  Fit the analysis model (e.g., GEE, GLMM) in each dataset.
3.  Pool estimates using Rubin’s rules.

#### MI + GEE

Combining MI with GEE (**MI-GEE**) is a valid MAR approach:

- MI handles missing covariates and outcomes.
- GEE provides population-average inference.
- Variability across imputations is reflected in standard errors.

### MNAR and Sensitivity Analysis

When missingness may depend on **unobserved outcomes** (MNAR), no method
can identify the truth from the data alone.

In this case:

- MNAR models are used for **sensitivity analysis**, not as the primary
  analysis.
- Examples:
  - $\delta$-adjusted pattern-mixture models,
  - Selection models with sensitivity parameters,
  - Shared-parameter models.

The goal is to assess **robustness**, not to “fix” MNAR.

### Summary: Choosing a Method

| Situation                            | Recommended Approach  |
|--------------------------------------|-----------------------|
| Complete data or MCAR                | GEE                   |
| MAR, population-average effects      | WGEE or MI-GEE        |
| MAR, subject-specific interpretation | LMM / GLMM            |
| Incomplete covariates                | MI (with GLMM or GEE) |
| Suspected MNAR                       | Sensitivity analysis  |

### Practical Interpretation for the Alzheimer BPRS Study

- Dropout is common and plausibly depends on observed disease history →
  **MAR is a reasonable working assumption**.
- For primary inference:
  - **LMM / GLMM** are valid likelihood-based tools.
  - **WGEE with observation-level weights** provides a robust
    population-average alternative.
- **MI** allows consistent handling of missing covariates and
  complements both approaches.
- MNAR methods are reserved for **sensitivity analysis**, aligning
  statistical modeling with the clinical reality of disease progression.

### Local Influence (Local Sensitivity) for Alzheimer BPRS Dropout

#### Why This Exists (The Short Story)

In the Alzheimer cohort, **later BPRS measurements are missing more
often** due to dropout caused by progression, institutionalization,
death, or loss to follow-up. If dropout depends only on **observed
history**, MAR (Missing At Random) can be reasonable. However, if
dropout depends on the **unobserved current BPRS** (e.g.,
worse-than-expected symptoms), that is MNAR (Missing Not At Random).

**Local influence** is a way to ask:

> “If we make a *tiny* MNAR deviation from MAR, which patients (or which
> parts of the model) cause the biggest change in our conclusions?”

It is a **sensitivity analysis tool**, not your primary analysis.

### 1) Setup in Alzheimer Notation

Let $Y_{ij}$ be the BPRS for subject $i$ at visit $j$, and let $D_i$
denote the dropout time (e.g., the first visit after which they stop
contributing data).

A common selection-model framing is:

$$
f(Y_i, D_i \mid X_i) =
f(Y_i \mid X_i, \theta) \cdot f(D_i \mid Y_i, X_i, \psi),
$$

where:

- $f(Y_i \mid X_i, \theta)$: the **longitudinal outcome model** (e.g.,
  an LMM for continuous BPRS),
- $f(D_i \mid Y_i, X_i, \psi)$: the **dropout model** (e.g., logistic
  regression or hazard-style).

### 2) The “Perturbed MAR” Idea

#### Baseline MAR Dropout Model (Uses Observed History Only)

A standard MAR dropout model uses **past observed BPRS** (and
covariates), e.g.,

$$
\operatorname{logit}\Big(P(D_i = j \mid D_i \geq j, Y_{i,0:(j-1)}^{o}, X_i)\Big) =
\psi_0 + \psi_1 Y_{i,j-1} + \psi_2^\top X_i.
$$

This is MAR-style because dropout depends on observed history (e.g.,
$Y_{i,j-1}$), not on the missing $Y_{ij}$.

#### Local MNAR Perturbation (One Extra “Forbidden” Term)

Local influence adds a small MNAR deviation:

$$
\operatorname{logit}\Big(P(D_i = j \mid D_i \geq j, Y_{i,0:(j-1)}^{o}, X_i, Y_{ij})\Big) =
\psi_0 + \psi_1 Y_{i,j-1} + \psi_2^\top X_i + \omega Y_{ij}.
$$

- $\omega = 0$: **MAR**.
- $\omega \neq 0$: **MNAR**, because dropout now depends on the
  **current (possibly missing) BPRS**.

**Local sensitivity** asks what happens when $\omega$ is perturbed
*slightly away from 0*.

### 3) Likelihood Displacement: What the Slide is Measuring

Let $L_{\omega}(\theta, \psi)$ be the (observed-data) likelihood under
perturbation $\omega$.

1.  Fit the **MAR model** at $\omega = 0$, giving
    $(\widehat{\theta}, \widehat{\psi})$.
2.  Fit the perturbed model at $\omega \neq 0$, giving
    $(\widehat{\theta}_\omega, \widehat{\psi}_\omega)$.

The slide’s **likelihood displacement** is:

$$
LD(\omega) =
2 \Big[
\ell_{\omega=0}(\widehat{\theta}, \widehat{\psi}) -
\ell_{\omega=0}(\widehat{\theta}_\omega, \widehat{\psi}_\omega)
\Big] \geq 0,
$$

where $\ell$ is the log-likelihood.

#### Interpretation for BPRS:

- If $LD(\omega)$ grows quickly even for tiny $\omega$, your inference
  about BPRS trajectories is **fragile** to small MNAR deviations.
- If $LD(\omega)$ stays flat near $\omega = 0$, your inference is
  **robust** (at least locally).

### 4) Normal Curvature ($C_h$): What “Local Influence Direction” Means

The slide’s geometry is:

- $\omega$: a direction in which we “nudge” the model away from MAR.
- $LD(\omega)$: the “height” of a surface above $\omega = 0$.
- The **curvature** tells you how sharply the surface rises.

A common expression is:

$$
C_h = 2 \left| h^\top \Delta^\top \mathcal{L}^{-1} \Delta h \right|,
$$

where:

- $h$: a perturbation direction (e.g., “perturb subject $i$” or “perturb
  MNAR term”),
- $\mathcal{L}$: an information-like matrix (how sharply the likelihood
  identifies parameters),
- $\Delta$: links perturbations to the score/information.

#### Practical Explanation:

> It ranks which perturbations (or which subjects) have the biggest
> impact on the fitted model and the key conclusions.

### 5) “Local Influence” Translated to Your Alzheimer BPRS Analysis

#### What is the Practical Output?

You typically end up with:

1.  **Which subjects are influential for missingness sensitivity?**
    - Example: Patients with **rapid BPRS worsening before dropout** may
      dominate the sensitivity.
    - Patients with unusually “good” observed BPRS who still drop out
      early could also drive sensitivity (suggesting dropout is not only
      severity-driven).
2.  **Which conclusions are sensitive?**
    - The **overall time trend** in BPRS.
    - Whether baseline severity markers (e.g., CDR-SB, ABPET) relate to
      level and/or slope.
    - Key interactions (if any remain).

Local influence says:

> “Do those estimated effects change a lot under tiny MNAR deviations?”

### 6) The Questions You Should Ask Yourself

#### A. Clinical Plausibility Questions (Before Touching Math)

- **Why do people drop out?** Progression, death, institutionalization,
  side effects, logistics?
- **Is dropout plausibly related to *current unobserved BPRS* even after
  conditioning on past BPRS?**
  - If yes → MNAR is plausible → sensitivity analysis is not optional.

#### B. Statistical Identification Questions (The Honest Limitations)

- Can $\omega$ be learned from observed data?
  - Usually **no** (MNAR parameter is not identified without
    assumptions).
  - Local influence is about **robustness**, not “finding the true
    MNAR.”

#### C. Model-Structure Questions (So Your Sensitivity is Meaningful)

- Does your dropout model include the **key observed predictors** of
  dropout (prior BPRS, time, baseline markers)?
- Are you mixing different missingness causes (death vs missed
  appointment) into one dropout indicator?
  - If death is present, it can behave like a different mechanism (and
    should be discussed explicitly).

#### D. Reporting Questions (What the Marker/Examiner Will Actually Look For)

- What is your *primary* assumption (MAR) and why is it reasonable?
- What is your *sensitivity scenario* (small MNAR deviation) and what
  does it represent clinically?
- Which conclusions are stable vs fragile?

### 7) How to “Handle This” in a Homework/Presentation Setting

Given your **10 minutes total**, local influence is best used as:

1.  **One slide framing:**
    - “Primary analysis assumes MAR (LMM / likelihood).”
    - “We assess robustness to small MNAR deviations using a perturbed
      dropout model.”
2.  **One slide result:**
    - “Key fixed effects (time trend, baseline severity association)
      were / were not robust.”
    - “Sensitivity was driven mainly by a small subset of early dropouts
      (top influential subjects).”

Keep the heavy machinery in an appendix.

### 8) Minimal Implementation Roadmap (Conceptual + R Skeleton)

#### Steps:

1.  Fit outcome model for BPRS (LMM): $\widehat{\theta}$.
2.  Fit dropout model under MAR: $\widehat{\psi}$.
3.  Add the MNAR perturbation term ($\omega Y_{ij}$) and examine small
    $\omega$ (local).
4.  Measure impact on:
    - $\widehat{\beta}(\omega)$ for key mean parameters,
    - Likelihood displacement $LD(\omega)$,
    - Subject-wise influence (case perturbations).

#### R Skeleton:

``` r
# 1) Outcome model (MAR primary): LMM for BPRS
# lmm <- lme4::lmer(BPRS ~ time + age + sex + cdrsb0 + abpet0 + ... + (1 + time | id), data = dat)

# 2) Dropout indicator per visit (example):
# R_ij = 1 if observed at j, 0 if missing at j
# You can define a discrete-time dropout hazard among those still present.
# dat$at_risk <- with(dat, ave(R, id, FUN = function(x) cumprod(c(1, head(x, -1)))))

# 3) MAR dropout model: depends on observed history (lagged BPRS etc.)
# drop_mar <- glm(dropout_event ~ lag_BPRS + time + age + sex + cdrsb0 + ..., family = binomial, data = dat_at_risk)

# 4) Local MNAR perturbation:
# You can't directly use Y_ij when it's missing. The selection-model approach uses the joint model.
# In practice for homework, local sensitivity is often approximated by:
#   - delta-shifting imputed values, or
#   - perturbing weights / dropout probabilities, or
#   - case-weight perturbations and checking changes in beta_hat.
#
# Then quantify change in key betas over small perturbations.
```

### Sensitivity Analysis for Missing BPRS in the Alzheimer Cohort

#### Why Sensitivity Analysis is Needed

In the Alzheimer cohort, follow-up BPRS measurements are often
incomplete due to dropout caused by worsening health,
institutionalization, or death. This creates a realistic risk that later
BPRS values are observed primarily among “healthier survivors.” While
the primary analysis assumes **MAR (Missing At Random)**, **MNAR
(Missing Not At Random)** cannot be ruled out. Sensitivity analysis
evaluates whether the main conclusions remain robust under plausible
MNAR departures.

### 1) Notation for Outcome and Dropout

Let $y_i = (y_{i0}, y_{i1}, \dots, y_{iJ})$ represent the BPRS
measurements for subject $i$ over visits $j = 0, \dots, J$. Define the
observation indicator:

$$
r_{ij} =
\begin{cases}
1, & \text{if } y_{ij} \text{ is observed}, \\
0, & \text{if } y_{ij} \text{ is missing}.
\end{cases}
$$

Split $y_i$ into observed and missing components:
$y_i = (y_i^o, y_i^m)$. The **dropout time** (last observed visit) is:

$$
D_i = \max \{ j : r_{ij} = 1 \}.
$$

In this setting, $D_i$ is informative if subjects with worse symptoms
tend to drop out earlier.

### 2) Primary Analysis: LMM Under MAR

The primary analysis uses a **linear mixed model (LMM)** under the MAR
assumption:

$$
y_{ij} = x_{ij}^\top \beta + z_{ij}^\top b_i + \varepsilon_{ij},
\qquad
b_i \sim N(0, G),
\qquad
\varepsilon_{ij} \sim N(0, \sigma^2),
$$

where $b_i$ captures subject-specific baseline levels and progression.
This model:

- Uses all available repeated measures without discarding incomplete
  subjects.
- Provides valid likelihood-based inference under MAR.

#### MAR Assumption in the Alzheimer Context

The MAR assumption implies that missingness at visit $j$ depends only on
observed history (e.g., previous BPRS, baseline severity):

$$
f(r_i \mid X_i, y_i^o, y_i^m) = f(r_i \mid X_i, y_i^o).
$$

This assumption is untestable but serves as a **working assumption** for
the primary analysis.

### 3) Sensitivity Analysis: MNAR Concerns

#### Key Questions for Sensitivity Analysis

1.  **Direction:** Are dropouts likely to have *higher* BPRS (worse
    symptoms) than predicted under MAR?
2.  **Magnitude:** How much worse would unobserved BPRS need to be to
    change conclusions?
3.  **Timing:** Does MNAR primarily affect late follow-up (e.g., years
    3–6)?
4.  **Robustness Target:** Which inference (e.g., time trend,
    group-by-time interaction) is most sensitive to MNAR?

### 4) MNAR Sensitivity Option A: Pattern-Mixture ($\delta$-Adjustment)

#### Step A1: Fit the MAR LMM

Fit the primary LMM under MAR and obtain MAR-based predictions
$\widehat{y}_{ij}^{MAR}$.

#### Step A2: Adjust Missing Outcomes

For $j > D_i$, define MNAR-adjusted outcomes:

$$
y_{ij}^m(\delta) = \widehat{y}_{ij}^{MAR} + \delta.
$$

- $\delta = 0$: MAR (no departure).
- $\delta > 0$: Dropouts have worse symptoms than predicted under MAR.
- $\delta < 0$: Dropouts have better symptoms (less plausible
  clinically).

#### Step A3: Refit and Evaluate

Refit the model for a grid of $\delta$ values (e.g.,
$\delta \in \{0, 0.5, 1, 2, 3\}$) and track how key estimates (e.g.,
time trend, group-by-time interaction) change.

#### Reporting

- Present a table/plot showing how estimates vary with $\delta$.
- Example conclusion: “Findings remain robust up to $\delta = 2$, but
  attenuate beyond $\delta = 3$.”

### 5) MNAR Sensitivity Option B: Selection Model

#### Joint Model for Outcome and Dropout

The joint model factorizes as:

$$
f(y_i, r_i \mid X_i, \theta, \psi) =
f(y_i \mid X_i, \theta) f(r_i \mid X_i, y_i, \psi).
$$

The dropout model includes an MNAR term:

$$
\operatorname{logit}\Big(
\Pr(D_i = j \mid D_i \geq j, \mathcal{H}_{i,j-1}, y_{ij})
\Big) =
\psi_0 + \psi_1 y_{i,j-1} + \psi_2 y_{ij} + \psi_3^\top X_i.
$$

- $\psi_2$: MNAR parameter linking dropout to unobserved $y_{ij}$.

#### Sensitivity Analysis

- Fix $\psi_2$ to plausible values (or vary over a range).
- Refit the joint model and observe how inference on $\beta$ changes.

### 6) Summary: Recommended Approach

- **Primary Analysis:** LMM under MAR.
- **Sensitivity Analysis:** Pattern-mixture ($\delta$-adjustment) for
  transparency and clinical interpretability.
- **Optional:** Selection model for formal MNAR modeling.

### 7) Q&A Defense

“We assume MAR for the primary LMM because dropout plausibly depends on
observed history like prior BPRS and baseline severity. However, MNAR
cannot be excluded, so we conducted a sensitivity analysis to evaluate
how much worse unobserved BPRS would need to be to change conclusions.
The $\delta$-adjustment approach is transparent, clinically
interpretable, and shows that our findings are robust to plausible MNAR
departures.”

### Sensitivity Analysis for Missing BPRS in the Alzheimer Cohort

#### Why Sensitivity Analysis is Needed

In the Alzheimer cohort, follow-up BPRS measurements are often
incomplete due to dropout caused by worsening health,
institutionalization, or death. This creates a realistic risk that later
BPRS values are observed primarily among “healthier survivors.” While
the primary analysis assumes **MAR (Missing At Random)**, **MNAR
(Missing Not At Random)** cannot be ruled out. Sensitivity analysis
evaluates whether the main conclusions remain robust under plausible
MNAR departures.

### 1) Notation for Outcome and Dropout

Let $y_i = (y_{i0}, y_{i1}, \dots, y_{iJ})$ represent the BPRS
measurements for subject $i$ over visits $j = 0, \dots, J$. Define the
observation indicator:

$$
r_{ij} =
\begin{cases}
1, & \text{if } y_{ij} \text{ is observed}, \\
0, & \text{if } y_{ij} \text{ is missing}.
\end{cases}
$$

Split $y_i$ into observed and missing components:
$y_i = (y_i^o, y_i^m)$. The **dropout time** (last observed visit) is:

$$
D_i = \max \{ j : r_{ij} = 1 \}.
$$

In this setting, $D_i$ is informative if subjects with worse symptoms
tend to drop out earlier.

### 2) Primary Analysis: Likelihood-Based LMM Under MAR

The primary analysis uses a **linear mixed model (LMM)** under the MAR
assumption:

$$
y_{ij} = x_{ij}^\top \beta + z_{ij}^\top b_i + \varepsilon_{ij},
\qquad
b_i \sim N(0, G),
\qquad
\varepsilon_{ij} \sim N(0, \sigma^2),
$$

where $b_i$ captures subject-specific baseline levels and progression.
This model:

- Uses all available repeated measures without discarding incomplete
  subjects.
- Directly models subject heterogeneity (random intercept/slope).
- Provides valid likelihood-based inference under MAR.

#### MAR Assumption in the Alzheimer Context

The MAR assumption implies that missingness at visit $j$ depends only on
observed history (e.g., previous BPRS, baseline severity):

$$
f(r_i \mid X_i, y_i^o, y_i^m) = f(r_i \mid X_i, y_i^o).
$$

This assumption is untestable but serves as a **working assumption** for
the primary analysis.

### 3) What Question Are We Asking in Sensitivity Analysis?

A sensitivity analysis starts by making the MNAR concern concrete in the
Alzheimer context. These are the kinds of questions to explicitly
address:

1.  **Direction:** Are participants who drop out likely to have *higher*
    BPRS than predicted from their observed history (worse symptoms), or
    *lower* BPRS (less plausible clinically)?
2.  **Magnitude:** If dropouts are worse, how much worse would their
    unobserved BPRS need to be to materially change conclusions (e.g.,
    time trend, group differences)?
3.  **Timing:** Does MNAR primarily affect late follow-up (e.g., years
    3–6), where missingness is heavier?
4.  **Whose Dropout:** Is dropout related to observed markers like
    baseline severity, ADL, cognition, prior BPRS? (If yes, MAR becomes
    more plausible after adjustment, but MNAR may still remain.)
5.  **Robustness Target:** Which inference do we want to
    stress-test—mean slope over time, treatment/group-by-time
    interaction, or baseline covariate associations?

### 4) MNAR Sensitivity Option A: Pattern-Mixture ($\delta$-Adjustment)

This is the easiest to explain and defend for a short presentation, and
it links naturally to BPRS.

#### Step A1: Fit the MAR LMM and Get MAR-Based Predictions

Let $\widehat{y}_{ij}^{MAR}$ be the fitted/predicted mean for a missing
value under the MAR LMM.

#### Step A2: Shift the Unobserved Outcomes After Dropout

For $j > D_i$, define an MNAR-adjusted value:

$$
y_{ij}^m(\delta) = \widehat{y}_{ij}^{MAR} + \delta.
$$

**Interpretation for BPRS:**

- $\delta = 0$: MAR (no departure).
- $\delta > 0$: Unobserved BPRS is systematically higher than MAR
  predicts (worse symptoms among dropouts).
- $\delta < 0$: Unobserved BPRS is lower than MAR predicts (less
  plausible clinically, but can be included as a check).

#### Step A3: Refit and Track the Conclusions Across $\delta$

Repeat the analysis for a small grid of clinically plausible shifts,
e.g.,

$$
\delta \in \{0, 0.5, 1, 2, 3\},
$$

in BPRS units (choose the scale based on what “clinically meaningful
change” looks like for your data).

**What to Report:**

- A table/figure showing how the key estimates (e.g., time effect,
  covariate-by-time interaction) change as $\delta$ increases.
- A short conclusion like: “Main findings remain unchanged up to
  $\delta = 2$, but beyond $\delta = 3$ the time trend attenuates.”

### 5) MNAR Sensitivity Option B: Selection Model (Diggle–Kenward Style)

This approach models the longitudinal outcome and the dropout mechanism
jointly, allowing dropout to depend on the (possibly missing) current
outcome.

#### Joint Factorization

$$
f(y_i, r_i \mid X_i, \theta, \psi) =
f(y_i \mid X_i, \theta) f(r_i \mid X_i, y_i, \psi).
$$

Likelihood contribution uses integration over missing outcomes:

$$
L(\theta, \psi) =
\prod_{i=1}^N
\int
f(y_i^o, y_i^m \mid X_i, \theta)
f(r_i \mid X_i, y_i^o, y_i^m, \psi)
dy_i^m.
$$

#### Discrete-Time Dropout Hazard with MNAR Term

A common choice is:

$$
\operatorname{logit}\Big(
\Pr(D_i = j \mid D_i \geq j, \mathcal{H}_{i,j-1}, y_{ij})
\Big) =
\psi_0 + \psi_1 y_{i,j-1} + \psi_2 y_{ij} + \psi_3^\top X_i.
$$

- $y_{i,j-1}$: Observed history (MAR-consistent part).
- $y_{ij}$: Current BPRS, which may be unobserved if dropout occurs at
  $j$.
- $\psi_2$: The key MNAR parameter encoding how dropout depends on
  unobserved current severity.

#### Why This Becomes Sensitivity Analysis in Practice

$\psi_2$ is not identified from observed data alone. Sensitivity
analysis involves:

- Fixing $\psi_2$ to plausible values (or varying it over a range),
- Refitting the joint model,
- Observing how inference on $\beta$ changes.

### 6) Summary: Recommended Approach

- **Primary Analysis:** LMM under MAR.
- **Sensitivity Analysis:** Pattern-mixture ($\delta$-adjustment) for
  transparency and clinical interpretability.
- **Optional:** Selection model for formal MNAR modeling.

### 7) Q&A Defense

“We assume MAR for the primary LMM because dropout plausibly depends on
observed history like prior BPRS and baseline severity. However, MNAR
cannot be excluded, so we conducted a sensitivity analysis to evaluate
how much worse unobserved BPRS would need to be to change conclusions.
The $\delta$-adjustment approach is transparent, clinically
interpretable, and shows that our findings are robust to plausible MNAR
departures.”

### MAR vs MNAR: What You Can and Cannot Do with BPRS Data

#### Why You Can’t Test MAR vs MNAR

It’s important to acknowledge that **MNAR (Missing Not At Random) is not
identifiable** from the observed data alone. This means you cannot
formally test:

$$
H_0: \text{MAR} \quad \text{vs} \quad H_1: \text{MNAR}.
$$

Any “test” would rely on unverifiable modeling assumptions. Instead, you
can:

1.  **Support MAR** by showing that dropout is explainable by observed
    history.
2.  **Run MNAR sensitivity analyses** to assess how robust your
    conclusions are under plausible MNAR departures.

### Supporting MAR: Diagnostics and Evidence

#### 1) Model the Dropout Indicator Using Observed History

Define a dropout indicator:

$$
H_{i,j} = \mathbb{1}(D_i = j \mid D_i \geq j),
$$

where $D_i$ is the dropout time for subject $i$. Fit a logistic model
for the dropout hazard using **only observed information up to visit
$j-1$**:

$$
\operatorname{logit}\Big(\Pr(H_{i,j} = 1 \mid \mathcal{H}_{i,j-1})\Big) =
\alpha_0 + \alpha_1 y_{i,j-1} + \alpha_2^\top X_i + \alpha_3 j.
$$

- **Interpretation:** If dropout is strongly predicted by **previous
  BPRS** and baseline severity markers, this supports the idea that
  missingness depends on observed history, making MAR plausible *after
  conditioning on that history*.

#### 2) Extend the Observed History

Expand $\mathcal{H}_{i,j-1}$ to include:

- The last two BPRS values (or the slope between the last two visits),
- ADL/cognition history (if available),
- Baseline severity markers,
- Time.

If the model’s predictive ability improves and the patterns align with
clinical expectations, MAR becomes more credible.

#### 3) Compare Pre-Dropout Trajectories vs Completers

Group subjects by their dropout time ($D_i$) and plot the mean BPRS
trajectories for each group. If early dropouts already have **worse
observed BPRS before dropout**, this suggests informative dropout (not
MCAR). While this doesn’t distinguish MAR from MNAR, it shows that
missingness is *not random*.

### What You Cannot Do

You cannot formally test MAR vs MNAR using only $(y^o, r, X)$. Instead,
you must rely on sensitivity analyses to explore the impact of MNAR.

### MNAR Sensitivity Analyses: What You Should Do

#### Option A: $\delta$-Shift Pattern-Mixture Model

1.  **Fit the MAR LMM** and predict missing values under MAR:

    $$
    \widehat{y}_{ij}^{MAR}.
    $$

2.  **Shift the post-dropout values** by a constant $\delta$:

    $$
    y_{ij}^m(\delta) = \widehat{y}_{ij}^{MAR} + \delta, \quad j > D_i.
    $$

3.  **Sweep $\delta$ over plausible values** (e.g.,
    $\delta \in \{0, 0.5, 1, 2, 3\}$) and refit the model or recompute
    the key estimand.

4.  **Report results:** For example, “Results remain unchanged up to
    $\delta = 2$, but conclusions change beyond $\delta = 3$.”

#### Option B: Selection Model with an MNAR Parameter

1.  **Specify a dropout hazard model** that depends on the current
    (possibly missing) BPRS:

    $$
    \operatorname{logit}\Big(\Pr(H_{i,j} = 1 \mid \mathcal{H}_{i,j-1}, y_{ij})\Big) =
    \psi_0 + \psi_1 y_{i,j-1} + \psi_2 y_{ij} + \psi_3^\top X_i.
    $$

2.  **Treat $\psi_2$ as a sensitivity parameter**:

    - Fix $\psi_2$ to plausible values (e.g., $\psi_2 > 0$ implies
      dropout is more likely for worse unobserved BPRS).
    - Refit the model and track how the key inference changes.

### The “5-Second Defense” for Your Presentation

“We can’t test MAR vs MNAR from the observed BPRS data. We support MAR
by showing that dropout is explained by observed history (previous BPRS,
baseline severity, time). We then perform MNAR sensitivity analyses
using a $\delta$-shift pattern-mixture model (and optionally a selection
model) to quantify how strong MNAR would need to be to change
conclusions.”

### Next Steps: Tailoring the Dropout Model

If you provide details about your time coding and the covariates
measured longitudinally (e.g., ADL each year vs only at baseline), I can
help you write the exact dropout model formula and suggest specific
plots to include.

### Notation (BPRS Scenario)

Let patient $i = 1, \dots, N$ have planned visits $j = 0, 1, \dots, J$
(baseline to year $J$).

- **Outcome vector:** $Y_i = (Y_{i0}, \dots, Y_{iJ})^\top$, where
  $Y_{ij}$ is **BPRS**.
- **Missingness indicator:** $R_{ij} = 1$ if $Y_{ij}$ is observed, $0$
  if missing.
- **Split outcome into observed/missing parts:**
  $Y_i = (Y_i^{o}, Y_i^{m})$.
- **Dropout time:** $D_i = \max \{j : R_{ij} = 1\}$, the last observed
  visit (monotone dropout).

------------------------------------------------------------------------

### The “MAR in 3 Frameworks” Slide: What It Means

The joint distribution of **(BPRS outcomes, missingness pattern)** can
always be written as:

$$
f(Y_i, R_i) = f(Y_i) f(R_i \mid Y_i) = f(R_i) f(Y_i \mid R_i).
$$

This means the same joint model can be factorized in different ways,
leading to different modeling approaches.

### 1) Selection Models (Diggle–Kenward Style)

**Factorization:**

$$
f(Y_i, R_i) = f(Y_i \mid \theta) f(R_i \mid Y_i, \psi),
$$

where:

- $f(Y_i \mid \theta)$: the **BPRS longitudinal model** (e.g., LMM).
- $f(R_i \mid Y_i, \psi)$: the **dropout model**.

**MAR Statement:**

$$
f(r_i \mid y_i, \psi) = f(r_i \mid y_i^{o}, \psi).
$$

This means that once you condition on observed BPRS history ($y_i^{o}$),
dropout does **not** depend on the unobserved BPRS values.

#### Logistic Dropout Model for BPRS (Slide Example)

For monotone dropout, a discrete-time hazard model is:

$$
\text{logit}\Big(P(D_i = j \mid D_i \geq j, Y_{ij}, Y_{i,j-1})\Big) =
\psi_0 + \psi_1 Y_{i,j-1} + \psi_2 Y_{ij}.
$$

**Interpretation:**

- $\psi_1$: Dropout depends on **past observed symptom severity** (e.g.,
  higher BPRS at the last visit increases dropout risk).

- $\psi_2$: Dropout depends on the **current BPRS value** ($Y_{ij}$),
  which is unobserved at dropout.

  - If $\psi_2 \neq 0$: Dropout depends on unobserved outcomes →
    **MNAR**.
  - If $\psi_2 = 0$: Dropout depends only on observed history → **MAR**.
  - If $\psi_1 = \psi_2 = 0$: Dropout is unrelated to BPRS → **MCAR**.

**Implementation in Your Project:**

- **MAR Analysis (Primary):** Fit an **LMM** for BPRS and ignore the
  dropout model (valid under MAR + distinct parameters).
- **Sensitivity (MNAR):** Allow a nonzero $\psi_2$ or introduce an
  explicit sensitivity parameter.

------------------------------------------------------------------------

### 2) Pattern-Mixture Models

**Factorization:**

$$
f(Y_i, R_i) = f(R_i \mid \psi) f(Y_i \mid R_i, \theta).
$$

Here, you condition on the **missingness pattern** ($R_i$) or dropout
time ($D_i$) and model outcomes within patterns.

**MAR Statement:**

$$
f(y_i^{m} \mid y_i^{o}, r_i, \theta) = f(y_i^{m} \mid y_i^{o}, \theta).
$$

This means that after conditioning on observed BPRS history, the
distribution of missing BPRS does **not** depend on the dropout pattern.

#### How This Helps for BPRS

Pattern-mixture models are useful for **sensitivity analysis** because
they allow you to specify how missing outcomes differ from MAR
predictions. A common approach is **delta-adjustment**:

1.  Fit an MAR model to predict missing $Y_{ij}$.

2.  Adjust the predictions:

    $$
    Y_{ij}^{m} = \widehat{Y}_{ij}^{MAR} + \delta, \quad j > D_i,
    $$

    where $\delta$ is a sensitivity parameter (e.g., $\delta > 0$ means
    “dropouts would have had higher BPRS than MAR predicts”).

This approach is defensible in Alzheimer data, where dropout is often
linked to worsening health or institutionalization.

### 3) Shared-Parameter Models

**Factorization:**

$$
f(Y_i, R_i) = \int f(Y_i \mid b_i, \theta) f(R_i \mid b_i, \psi) f(b_i) db_i.
$$

Here, a latent subject-specific process ($b_i$) drives both:

- The BPRS trajectory ($f(Y_i \mid b_i, \theta)$).
- The dropout propensity ($f(R_i \mid b_i, \psi)$).

#### BPRS Interpretation

Patients have an unobserved “disease progression tendency” ($b_i$) that
affects both:

- Their BPRS trajectory over time, and
- Their likelihood of dropping out.

This framework is natural if dropout is connected to deterioration or
death. It is often implemented as a **joint model** (longitudinal
outcome + time-to-dropout) with shared random effects.

### Recommended Approach for Your BPRS Report

#### Primary Analysis

**LMM under MAR** is the default choice:

$$
Y_{ij} = x_{ij}^\top \beta + z_{ij}^\top b_i + \varepsilon_{ij},
$$

where:

- $b_i \sim N(0, G)$,
- $\varepsilon_{ij} \sim N(0, \sigma^2)$.

**Why This Is Appropriate:**

- BPRS is continuous and repeatedly measured → LMM fits the data
  structure.
- Under **MAR + ignorability**, you don’t need to model dropout
  explicitly to get valid inference for $\beta$.
- In Alzheimer follow-up, dropout is rarely MCAR; MAR is a reasonable
  working assumption if you include strong predictors of missingness
  (e.g., previous BPRS, age, baseline severity, ADL).

#### Sensitivity Analysis

Choose one sensitivity approach:

1.  **Pattern-Mixture ($\delta$-Adjustment):**

    - Assume unobserved BPRS is $\delta$ points higher than MAR
      predicts.
    - Vary $\delta$ over a plausible range.

2.  **Selection Model (MNAR):**

    - Allow dropout to depend on the current (possibly unobserved) BPRS
      via $\psi_2$.
    - Examine how conclusions change as $\psi_2$ moves away from 0.

3.  **Shared-Parameter Model:**

    - Fit a joint model where dropout depends on the same random effects
      as BPRS.
    - Compare results to the MAR-LMM.

For a concise presentation, **(1)** is the easiest to explain.

### Key Questions for BPRS Sensitivity Analysis

1.  **What is the plausible dropout story in Alzheimer follow-up?**

    - Worsening symptoms? Caregiver burden? Institutionalization? Death?
      (If yes, MNAR is plausible.)

2.  **Does dropout depend on observed history?**

    - If higher last-visit BPRS predicts dropout, MAR is plausible
      *conditional on that history*.

3.  **If MNAR is plausible, in what direction would missing BPRS differ
    from MAR predictions?**

    - Typically “worse than predicted” → motivates $\delta > 0$ or
      $\psi_2 > 0$.

4.  **Do your main conclusions (time trend, covariate effects) survive
    reasonable perturbations?**

    - Sensitivity analysis tests robustness, not “proves MAR.”

### Missing Data Methods: Overview and Application to BPRS in Alzheimer Cohort

#### Missing Data Toolkit from Course Notes

The course notes organize missing data methods by **missingness
mechanism** (MCAR → MAR → MNAR) and **modeling framework** (selection,
pattern-mixture, shared-parameter). Below is a summary of the methods
and their relevance to your BPRS analysis.

### 1) Complete-Case (CC) Analysis

- **What it is:** Analyze only subjects/visits with no missing data
  (discard incomplete cases).
- **When it’s appropriate:** Only under **MCAR** (Missing Completely At
  Random).
- **Why it’s discouraged:** The notes explicitly warn that **CC is
  biased outside MCAR** and inefficient because it discards partial
  data.

### 2) Last Observation Carried Forward (LOCF)

- **What it is:** Replace missing future outcomes with the subject’s
  last observed value.

- **When it’s appropriate:** Rarely appropriate; the notes discourage
  LOCF as a default.

- **Why it’s discouraged:**

  - **Bias:** LOCF assumes no change after dropout, which is implausible
    in progressive diseases like Alzheimer’s.
  - **Inefficiency:** LOCF underestimates variability and produces
    confidence intervals that are too narrow.

### 3) Direct Likelihood Methods (Primary MAR Approach)

- **What it is:** Fit a likelihood-based model directly to the observed
  data without imputing missing values.
- **When it’s appropriate:** Under **MAR** (Missing At Random).
- **Why it’s recommended:**
  - The notes describe direct likelihood as **“easy to conduct”** under
    MAR.
  - Your BPRS **LMM/GLMM** is a direct-likelihood approach, making it
    the default for MAR-valid inference.

### 4) Direct Bayesian Methods

- **What it is:** Bayesian analog of likelihood-based analysis.
- **When it’s appropriate:** Also valid under **MAR**.
- **Why it’s useful:** Incorporates prior information and naturally
  handles uncertainty about missing data.

### 5) Inverse Probability Weighting (IPW) / Weighted GEE (W-GEE)

- **What it is:** Reweight observations by the inverse probability of
  being observed, so subjects likely to drop out contribute more.
- **When it’s appropriate:** Under **MAR**, the notes recommend
  **weighted GEE** as an alternative to standard GEE, which is valid
  only under MCAR.
- **Why it’s useful:** Corrects for selection bias due to dropout by
  modeling the dropout process.

### 6) Multiple Imputation (MI)

- **What it is:** Replace missing values with multiple plausible draws,
  analyze each completed dataset, and pool results.
- **When it’s appropriate:** Primarily under **MAR**.
- **Why it’s useful:** Handles incomplete covariates and outcomes, and
  provides valid inference by propagating uncertainty about missing
  data.

### 7) MNAR Methods (Sensitivity Analysis)

- **What it is:** Methods for **MNAR (Missing Not At Random)**, where
  missingness depends on unobserved outcomes.

- **When it’s appropriate:** The notes emphasize MNAR methods for
  **sensitivity analysis**, not as the primary analysis.

- **Examples:**

  - **Selection models** (e.g., Diggle–Kenward).
  - **Pattern-mixture models** (e.g., $\delta$-adjustment).
  - **Shared-parameter models**.
  - **Local influence**.

### Application to BPRS in Alzheimer Cohort

#### Primary Analysis: MAR Assumption

Your primary analysis assumes **MAR**, which is reasonable if dropout
depends on observed history (e.g., prior BPRS, baseline severity). The
recommended methods under MAR are:

1.  **Direct Likelihood (LMM/GLMM):**

    - Fit a linear mixed model (LMM) for BPRS under MAR.
    - This is the default approach for continuous outcomes in
      longitudinal studies.

2.  **Weighted GEE (W-GEE):**

    - Use observation-level inverse-probability weights to correct for
      dropout.
    - Provides population-average effects under MAR.

3.  **Multiple Imputation (MI):**

    - Impute missing BPRS values multiple times and pool results.
    - Useful for handling incomplete covariates or outcomes.

#### Sensitivity Analysis: MNAR Assumption

Since MNAR cannot be ruled out, sensitivity analysis is essential. The
notes recommend the following MNAR methods:

1.  **Pattern-Mixture Models ($\delta$-Adjustment):**

    - Shift post-dropout BPRS values by a constant $\delta$ to simulate
      worse outcomes than MAR predicts.
    - Example: $y_{ij}^m(\delta) = \widehat{y}_{ij}^{MAR} + \delta$,
      where $\delta > 0$ represents MNAR departures.

2.  **Selection Models:**

    - Model dropout as a function of both observed and unobserved
      outcomes.
    - Example: Include the current (possibly missing) BPRS in the
      dropout model.

3.  **Shared-Parameter Models:**

    - Link the outcome and dropout processes through shared random
      effects (e.g., frailty).

4.  **Local Influence:**

    - Assess the impact of small MNAR deviations on key conclusions.

### Summary of Recommendations for BPRS Analysis

1.  **Primary Analysis (MAR):**

    - Use a likelihood-based LMM for BPRS.
    - Optionally, complement with W-GEE or MI for robustness under MAR.

2.  **Sensitivity Analysis (MNAR):**

    - Perform a $\delta$-adjustment pattern-mixture analysis.
    - Optionally, add a selection model or shared-parameter model for
      comparison.

3.  **Avoid LOCF/CC:**

    - LOCF and CC are discouraged due to bias and inefficiency.

4.  **Report Robustness:**

    - Highlight how key conclusions (e.g., time trend, group
      differences) change under MNAR assumptions.

This approach aligns with the course notes’ guidance and provides a
defensible framework for analyzing BPRS in the Alzheimer cohort.

### Sensitivity Analysis for Missingness in BPRS (Alzheimer Cohort)

#### Key Truth: You Cannot Prove or Disprove MAR

It is important to acknowledge that **MAR (Missing At Random)** cannot
be proven or disproven from the observed data alone. Instead, you can:

1.  **Show that MCAR is implausible**: Demonstrate that dropout is
    related to observed history.
2.  **Perform sensitivity analysis**: Assess how conclusions change
    under plausible MNAR (Missing Not At Random) scenarios.

This is why frameworks like **selection models**, **pattern-mixture
models**, and **shared-parameter models** are used—not to “prove” MAR,
but as structured tools for **what-if sensitivity analysis**.

------------------------------------------------------------------------

### 1) Primary MAR Analysis: Linear Mixed Model (LMM)

The primary analysis assumes MAR and uses a likelihood-based linear
mixed model (LMM) for BPRS. This approach is valid under MAR if
missingness depends only on observed history.

``` r
library(lme4)

# Fit the primary MAR model
fm_adj <- bprs ~ time + sex + age + bmi + job + adl + wzc + cdrsb0 +
  time:cdrsb0 + (time | patid) + (1 | trial)

fit_mar <- lmer(fm_adj, data = alzheimer_long, na.action = na.exclude)
summary(fit_mar)
```

This model includes:

- **Fixed effects**: Time trend, baseline predictors, and time×cdrsb0
  interaction.
- **Random effects**: Patient-specific intercepts/slopes and trial-level
  intercepts.

### 2) Questions to Ask About Missingness

#### A. What Kind of Missingness Exists?

- **Monotone dropout**: Once BPRS is missing, it stays missing (common
  in Alzheimer follow-up).
- **Intermittent missingness**: Patients skip a visit but return later.

This distinction matters because:

- **Dropout models** (e.g., hazard/logistic for “drop next”) are suited
  for monotone dropout.
- **Observation-level models** are needed for intermittent patterns.

#### B. Is MCAR Plausible?

Exploratory analysis shows that higher baseline BPRS is associated with
more missed visits, and missingness increases over time. This makes MCAR
implausible.

#### C. Is MAR Plausible?

MAR assumes that, after conditioning on observed history (e.g., prior
BPRS, time, baseline covariates), the probability of missingness does
**not** depend on unobserved BPRS values. While this cannot be fully
verified, you can check whether missingness is strongly explained by
observed history, which supports MAR as a working assumption.

### 3) Selection Model: Modeling Dropout Using Observed History

A **selection model** examines the dropout hazard as a function of past
and (possibly) current outcomes. For example:

$$
\text{logit}\Big(P(D_{ij} = 1 \mid D_i \geq j, y_{i,j-1}, y_{ij})\Big) =
\psi_0 + \psi_1 y_{i,j-1} + \psi_2 y_{ij}.
$$

- **MAR**: $\psi_2 = 0$ (dropout depends only on observed history).
- **MNAR**: $\psi_2 \neq 0$ (dropout depends on unobserved outcomes).

Since $\psi_2$ is not identifiable from observed data, this model is
used for **sensitivity analysis**.

#### Practical Implementation for BPRS

``` r
library(dplyr)

# Prepare data for dropout modeling
dat_miss <- alzheimer_long %>%
  arrange(patid, time) %>%
  group_by(patid) %>%
  mutate(
    r_t = as.integer(!is.na(bprs)),            # Observed at time t
    r_next = lead(r_t),                        # Observed at time t+1
    drop_next = if_else(r_t == 1 & r_next == 0, 1L, 0L),
    prev_bprs = lag(bprs)                      # Previous BPRS
  ) %>%
  ungroup()

# Filter rows where dropout is defined
drop_dat <- dat_miss %>%
  filter(r_t == 1, !is.na(r_next), !is.na(prev_bprs))

# Fit dropout model
drop_fit_mar <- glm(
  drop_next ~ time + prev_bprs + sex + age + bmi + job + adl + wzc + cdrsb0 + time:cdrsb0,
  family = binomial(),
  data = drop_dat
)

summary(drop_fit_mar)
```

**How to Report This:**

- This is a **selection-model-style dropout regression**.
- If dropout is explainable by **past observed BPRS and covariates**,
  MAR is more plausible.
- This does **not** prove MAR but shows whether MCAR is implausible and
  whether an observed-history mechanism is reasonable.

### 4) MAR Robustness Check: Inverse-Probability Weighting (IPW)

Inverse-probability weighting (IPW) corrects for informative dropout
under MAR by reweighting observations based on their probability of
being observed.

#### Compute Stabilized Weights

``` r
# Predict dropout probabilities
drop_dat <- drop_dat %>%
  mutate(p_drop = predict(drop_fit_mar, type = "response"),
         p_stay = 1 - p_drop)

# Cumulative probability of being observed up to time t
w_dat <- drop_dat %>%
  arrange(patid, time) %>%
  group_by(patid) %>%
  mutate(p_obs_cum = cumprod(p_stay),
         w_ipw = 1 / p_obs_cum) %>%
  ungroup()

# Merge weights back to the original data
alzheimer_w <- alzheimer_long %>%
  left_join(w_dat %>% select(patid, time, w_ipw), by = c("patid", "time")) %>%
  mutate(w_ipw = if_else(is.na(w_ipw), 1, w_ipw))  # Weight = 1 where undefined
```

#### Fit Weighted LMM

``` r
# Fit weighted LMM
fit_mar_ipw <- lmer(fm_adj, data = alzheimer_w, weights = w_ipw, na.action = na.exclude)

summary(fit_mar)
summary(fit_mar_ipw)
```

**Interpretation:**

- If key effects (e.g., time, time×cdrsb0) are stable between weighted
  and unweighted models, conclusions are robust under MAR.
- Large shifts suggest fragility even within MAR, requiring MNAR
  sensitivity analysis.

### 5) MNAR Sensitivity: Pattern-Mixture Models (PMM)

Pattern-mixture models re-factorize the joint distribution by dropout
pattern:

$$
f(y_i, r_i) = f(y_i \mid r_i) f(r_i).
$$

#### Delta Adjustment for MNAR

Assume unobserved BPRS after dropout is systematically worse than MAR
predicts:

$$
y_{ij}^{(MNAR)} = y_{ij}^{(MAR)} + \delta, \quad j > D_i.
$$

- $\delta = 0$: MAR.
- $\delta > 0$: Dropouts have worse symptoms than predicted under MAR.

#### Implementation Skeleton

1.  Define dropout time ($D_i$).
2.  Impute missing BPRS under MAR.
3.  Add $\delta$ to imputed values after $D_i$.
4.  Refit the model for a grid of $\delta$ values (e.g.,
    $\delta \in \{0, 2, 5, 10\}$).

### 6) Shared-Parameter Models: Latent MNAR Mechanism

Shared-parameter models assume that dropout depends on the same latent
process driving BPRS:

- **Longitudinal model**:
  $y_{ij} = x_{ij}^\top \beta + z_{ij}^\top b_i + \varepsilon_{ij}$,
- **Dropout hazard**:
  $h_i(t) = h_0(t) \exp(\gamma^\top w_i + \alpha m_i(t))$, where
  $m_i(t)$ is the latent trajectory.

This is a principled MNAR approach but requires specialized software
(e.g., `joineR`, `JMbayes`).

### Reporting Recommendations

1.  **Primary Analysis**: Fit an LMM under MAR.
2.  **MAR Robustness**:
    - Model dropout using observed history.
    - Perform IPW-weighted sensitivity analysis.
3.  **MNAR Sensitivity**:
    - Use pattern-mixture models ($\delta$-adjustment).
    - Optionally, fit shared-parameter models.

**Conclusion**: Report whether key effects (e.g., time, time×cdrsb0) are
robust across MAR and MNAR scenarios.

### MAR vs MNAR: What You Can and Cannot Do with BPRS Data

#### Why You Can’t Test MAR vs MNAR

It’s important to acknowledge that **MNAR (Missing Not At Random) is not
identifiable** from the observed data alone. This means you cannot
formally test:

$$
H_0: \text{MAR} \quad \text{vs} \quad H_1: \text{MNAR}.
$$

Any “test” would rely on unverifiable modeling assumptions. Instead, you
can:

1.  **Support MAR** by showing that dropout is explainable by observed
    history.
2.  **Run MNAR sensitivity analyses** to assess how robust your
    conclusions are under plausible MNAR departures.

### Supporting MAR: Diagnostics and Evidence

#### 1) Model the Dropout Indicator Using Observed History

Define a dropout indicator:

$$
H_{i,j} = \mathbb{1}(D_i = j \mid D_i \geq j),
$$

where $D_i$ is the dropout time for subject $i$. Fit a logistic model
for the dropout hazard using **only observed information up to visit
$j-1$**:

$$
\operatorname{logit}\Big(\Pr(H_{i,j} = 1 \mid \mathcal{H}_{i,j-1})\Big) =
\alpha_0 + \alpha_1 y_{i,j-1} + \alpha_2^\top X_i + \alpha_3 j.
$$

- **Interpretation:** If dropout is strongly predicted by **previous
  BPRS** and baseline severity markers, this supports the idea that
  missingness depends on observed history, making MAR plausible *after
  conditioning on that history*.

#### 2) Extend the Observed History

Expand $\mathcal{H}_{i,j-1}$ to include:

- The last two BPRS values (or the slope between the last two visits),
- ADL/cognition history (if available),
- Baseline severity markers,
- Time.

If the model’s predictive ability improves and the patterns align with
clinical expectations, MAR becomes more credible.

#### 3) Compare Pre-Dropout Trajectories vs Completers

Group subjects by their dropout time ($D_i$) and plot the mean BPRS
trajectories for each group. If early dropouts already have **worse
observed BPRS before dropout**, this suggests informative dropout (not
MCAR). While this doesn’t distinguish MAR from MNAR, it shows that
missingness is *not random*.

### What You Cannot Do

You cannot formally test MAR vs MNAR using only $(y^o, r, X)$. Instead,
you must rely on sensitivity analyses to explore the impact of MNAR.

### MNAR Sensitivity Analyses: What You Should Do

#### Option A: $\delta$-Shift Pattern-Mixture Model

1.  **Fit the MAR LMM** and predict missing values under MAR:

    $$
    \widehat{y}_{ij}^{MAR}.
    $$

2.  **Shift the post-dropout values** by a constant $\delta$:

    $$
    y_{ij}^m(\delta) = \widehat{y}_{ij}^{MAR} + \delta, \quad j > D_i.
    $$

3.  **Sweep $\delta$ over plausible values** (e.g.,
    $\delta \in \{0, 0.5, 1, 2, 3\}$) and refit the model or recompute
    the key estimand.

4.  **Report results:** For example, “Results remain unchanged up to
    $\delta = 2$, but conclusions change beyond $\delta = 3$.”

#### Option B: Selection Model with an MNAR Parameter

1.  **Specify a dropout hazard model** that depends on the current
    (possibly missing) BPRS:

    $$
    \operatorname{logit}\Big(\Pr(H_{i,j} = 1 \mid \mathcal{H}_{i,j-1}, y_{ij})\Big) =
    \psi_0 + \psi_1 y_{i,j-1} + \psi_2 y_{ij} + \psi_3^\top X_i.
    $$

2.  **Treat $\psi_2$ as a sensitivity parameter**:

    - Fix $\psi_2$ to plausible values (e.g., $\psi_2 > 0$ implies
      dropout is more likely for worse unobserved BPRS).
    - Refit the model and track how the key inference changes.

### The “5-Second Defense” for Your Presentation

“We can’t test MAR vs MNAR from the observed BPRS data. We support MAR
by showing that dropout is explained by observed history (previous BPRS,
baseline severity, time). We then perform MNAR sensitivity analyses
using a $\delta$-shift pattern-mixture model (and optionally a selection
model) to quantify how strong MNAR would need to be to change
conclusions.”

### Next Steps: Tailoring the Dropout Model

If you provide details about your time coding and the covariates
measured longitudinally (e.g., ADL each year vs only baseline), I can
help you write the exact dropout model formula and suggest specific
plots to include.

### Shared-Parameter Models (SPM) for MNAR Sensitivity in BPRS Analysis

#### Why Use Shared-Parameter Models?

In Alzheimer follow-up studies, missingness in BPRS data is often
related to latent clinical severity. Shared-parameter models (SPM) allow
missingness to depend on the same latent subject-specific effects that
drive the outcome trajectory. This provides a natural framework for
**MNAR (Missing Not At Random)** sensitivity analysis.

### 1) Translating the SPM Framework to BPRS

#### (A) Outcome Model (BPRS Trajectory)

The outcome model captures the longitudinal BPRS trajectory:

$$
Y_{ij} = x_{ij}^\top \beta + b_{0i} + b_{1i} t_{ij} + u_{\text{trial}(i)} + \varepsilon_{ij},
$$

where:

- $b_{0i}$: random intercept (subject-specific baseline severity),
- $b_{1i}$: random slope (subject-specific progression rate),
- $u_{\text{trial}(i)}$: trial-level random intercept,
- $\varepsilon_{ij}$: residual error.

This is equivalent to your MAR model:

``` r
fm_adj <- bprs ~ time + sex + age + bmi + job + adl + wzc + cdrsb0 +
  time:cdrsb0 + (time | patid) + (1 | trial)

fit_mar <- lmer(fm_adj, data = alzheimer_long, na.action = na.exclude)
```

``` rmd
### Shared-Parameter Models (SPM) for MNAR Sensitivity in BPRS Analysis

#### Why Use Shared-Parameter Models?

In Alzheimer follow-up studies, missingness in BPRS data is often related to latent clinical severity. Shared-parameter models (SPM) allow missingness to depend on the same latent subject-specific effects that drive the outcome trajectory. This provides a natural framework for **MNAR (Missing Not At Random)** sensitivity analysis.

---

### 1) Translating the SPM Framework to BPRS

#### (A) Outcome Model (BPRS Trajectory)

The outcome model captures the longitudinal BPRS trajectory:

$$
Y_{ij} = x_{ij}^\top \beta + b_{0i} + b_{1i} t_{ij} + u_{\text{trial}(i)} + \varepsilon_{ij},
$$

where:

- \(b_{0i}\): random intercept (subject-specific baseline severity),
- \(b_{1i}\): random slope (subject-specific progression rate),
- \(u_{\text{trial}(i)}\): trial-level random intercept,
- \(\varepsilon_{ij}\): residual error.

This is equivalent to your MAR model:

```r
fm_adj <- bprs ~ time + sex + age + bmi + job + adl + wzc + cdrsb0 +
  time:cdrsb0 + (time | patid) + (1 | trial)

fit_mar <- lmer(fm_adj, data = alzheimer_long, na.action = na.exclude)
```

#### (B) Missingness Model (Dropout Probability)

The missingness model links the probability of observing $Y_{ij}$ to the
latent random effects:

$$
\text{logit}\Big(P(R_{ij} = 1 \mid R_{i,j-1} = 1, x_{ij}, b_{0i}, b_{1i})\Big) =
\gamma_0 + \gamma^\top w_{ij} + \alpha_0 b_{0i} + \alpha_1 b_{1i}.
$$

- $R_{ij} = 1$: BPRS observed at visit $j$,
- $R_{i,j-1} = 1$: subject was still under follow-up at visit $j-1$
  (monotone dropout setup),
- $\alpha_0, \alpha_1$: MNAR parameters linking missingness to latent
  severity/progression.

------------------------------------------------------------------------

### 2) Practical Implementation in R

#### (A) Two-Stage Diagnostic/Sensitivity Approach

This approach is simpler than a full joint model and is often sufficient
for course-level sensitivity analysis.

##### Step 1: Fit the Outcome Model (Done)

Fit the MAR model and extract subject-level random effects (BLUPs):

``` r
re <- ranef(fit_mar)$patid
# Extract random intercepts (b0) and slopes (b1)
re_df <- tibble(
  patid = rownames(re),
  b0 = re[["(Intercept)"]],
  b1 = re[["time"]]
)
```

##### Step 2: Build the Response Indicator ($R_{ij}$)

Create a dataset with scheduled visits and observed/missing indicators:

``` r
library(dplyr)

datR <- alzheimer_long %>%
  distinct(patid, trial, sex, age, bmi, job, adl, wzc, cdrsb0, time) %>%
  left_join(alzheimer_long %>% select(patid, time, bprs), by = c("patid", "time")) %>%
  mutate(R = as.integer(!is.na(bprs))) %>%
  group_by(patid) %>%
  arrange(time) %>%
  mutate(R_prev = lag(R)) %>%
  ungroup()
```

##### Step 3: Add Random Effects to the Missingness Model

Merge the random effects into the dataset:

``` r
datR2 <- datR %>% left_join(re_df, by = "patid")
```

##### Step 4: Fit the Missingness Model

Fit a logistic regression for dropout probability, including the random
effects:

``` r
library(lme4)

miss_dat <- datR2 %>% filter(!is.na(R_prev), R_prev == 1)

fit_R <- glmer(
  R ~ time + sex + age + bmi + job + adl + wzc + cdrsb0 + time:cdrsb0 +
    b0 + b1 + (1 | trial),
  family = binomial(),
  data = miss_dat
)

summary(fit_R)
```

**Interpretation:**

- Significant coefficients for `b0` or `b1` suggest that missingness
  depends on latent severity/progression, indicating MNAR-type
  dependence.
- This motivates further MNAR sensitivity analysis.

#### (B) Best Practice: Joint Shared-Parameter Model

A full SPM is fit jointly, combining the outcome and missingness models.
This requires specialized software (e.g., `brms`, `JMbayes`).

### 3) Why This Framework Fits Alzheimer BPRS Data

- **Clinical Relevance:** Missingness in Alzheimer cohorts is often tied
  to worsening symptoms, caregiver burden, or
  institutionalization—factors naturally linked to latent severity.
- **Interpretability:** Shared-parameter models explicitly capture the
  relationship between disease progression and dropout.

### 4) Next Steps

- Confirm whether missingness is monotone or intermittent.
- If monotone, use the dropout model conditioning on $R_{i,j-1} = 1$.
- Report sensitivity results, emphasizing how MNAR assumptions affect
  key conclusions.

<!-- -->




    ### Pattern-Mixture Sensitivity Analysis for BPRS in Alzheimer Cohort

    #### Why Sensitivity Analysis is Needed

    In the Alzheimer cohort, follow-up BPRS measurements are often incomplete due to dropout caused by worsening health, institutionalization, or death. While the primary analysis assumes **MAR (Missing At Random)**, **MNAR (Missing Not At Random)** cannot be ruled out. Sensitivity analysis evaluates whether the main conclusions remain robust under plausible MNAR departures.

    ---

    ### 1) Translating the Slide’s Logic to BPRS

    The slide demonstrates a **pattern-mixture sensitivity analysis** for a binary outcome at two time points, split by **treatment arm** and **missingness pattern** (Completers vs Dropouts). For BPRS, we extend this to:

    - **Continuous outcome**: \(Y_{ij} = \text{BPRS}_{ij}\),
    - **Multiple follow-up times**: \(j = 0, 1, \dots, J\),
    - **Dropout patterns**: Defined by the last observed visit (\(D_i\)).

    ---

    ### 2) Defining Dropout Patterns

    Let \(R_{ij} = 1\) if BPRS is observed at visit \(j\), and \(R_{ij} = 0\) if missing. Define the **dropout time**:

    $$
    D_i = \max \{j : R_{ij} = 1\},
    $$

    where:

    - **Completers**: \(D_i = J\) (last planned visit),
    - **Dropouts**: \(D_i < J\) (missing after \(D_i\)).

    In R:

    ```r
    library(dplyr)

    # Define dropout patterns
    alzheimer_long <- alzheimer_long %>%
      group_by(patid) %>%
      mutate(
        R = as.integer(!is.na(bprs)),  # Observation indicator
        D_i = max(time[R == 1], na.rm = TRUE)  # Last observed visit
      ) %>%
      ungroup()

``` rmd
### Pattern-Mixture Sensitivity Analysis for BPRS in Alzheimer Cohort

#### Why Sensitivity Analysis is Needed

In the Alzheimer cohort, follow-up BPRS measurements are often incomplete due to dropout caused by worsening health, institutionalization, or death. While the primary analysis assumes **MAR (Missing At Random)**, **MNAR (Missing Not At Random)** cannot be ruled out. Sensitivity analysis evaluates whether the main conclusions remain robust under plausible MNAR departures.

---

### 1) Translating the Slide’s Logic to BPRS

The slide demonstrates a **pattern-mixture sensitivity analysis** for a binary outcome at two time points, split by **treatment arm** and **missingness pattern** (Completers vs Dropouts). For BPRS, we extend this to:

- **Continuous outcome**: \(Y_{ij} = \text{BPRS}_{ij}\),
- **Multiple follow-up times**: \(j = 0, 1, \dots, J\),
- **Dropout patterns**: Defined by the last observed visit (\(D_i\)).

---

### 2) Defining Dropout Patterns

Let \(R_{ij} = 1\) if BPRS is observed at visit \(j\), and \(R_{ij} = 0\) if missing. Define the **dropout time**:

$$
D_i = \max \{j : R_{ij} = 1\},
$$

where:

- **Completers**: \(D_i = J\) (last planned visit),
- **Dropouts**: \(D_i < J\) (missing after \(D_i\)).

In R:

```r
library(dplyr)

# Define dropout patterns
alzheimer_long <- alzheimer_long %>%
  group_by(patid) %>%
  mutate(
    R = as.integer(!is.na(bprs)),  # Observation indicator
    D_i = max(time[R == 1], na.rm = TRUE)  # Last observed visit
  ) %>%
  ungroup()
```

------------------------------------------------------------------------

### 3) Pattern-Mixture Model (PMM)

A **pattern-mixture model** allows the mean trajectory to differ by
dropout pattern. The mean structure is:

$$
E(Y_{ij} \mid D_i) = x_{ij}^\top \beta + \eta_{D_i} + \kappa_{D_i} t_{ij},
$$

where:

- $\eta_{D_i}$: Pattern-specific shift,
- $\kappa_{D_i}$: Pattern-specific slope adjustment.

In R:

``` r
# Fit a pattern-mixture model
fit_pmm <- lmer(
  bprs ~ time + sex + age + bmi + job + adl + wzc + cdrsb0 +
    time:cdrsb0 + factor(D_i) + time:factor(D_i) +
    (time | patid) + (1 | trial),
  data = alzheimer_long,
  na.action = na.exclude
)
summary(fit_pmm)
```

------------------------------------------------------------------------

### 4) Identifying Restrictions: MAR vs MNAR

Since outcomes after dropout are unobserved, we impose **identifying
restrictions** to complete the missing data. These restrictions define
the sensitivity scenarios:

#### (A) MAR (Baseline Assumption)

Under MAR, missing BPRS values follow the predicted trajectory from the
observed data:

$$
Y_{ij}^{\text{MAR}} = \widehat{Y}_{ij}^{\text{MAR}}, \quad j > D_i.
$$

#### (B) MNAR Scenarios

1.  **$\delta$-Shift Adjustment**: Add a constant shift to the
    MAR-predicted values after dropout:

    $$
    Y_{ij}^{\text{MNAR}} = Y_{ij}^{\text{MAR}} + \delta, \quad j > D_i.
    $$

2.  **Slope Change**: Accelerate deterioration after dropout:

    $$
    Y_{ij}^{\text{MNAR}} = Y_{ij}^{\text{MAR}} + \delta_{\text{slope}} (t_{ij} - t_{iD_i}), \quad j > D_i.
    $$

------------------------------------------------------------------------

### 5) Implementation: $\delta$-Shift Sensitivity Analysis

#### Step 1: Fit the MAR Model

Fit the primary MAR model:

``` r
fit_mar <- lmer(
  bprs ~ time + sex + age + bmi + job + adl + wzc + cdrsb0 +
    time:cdrsb0 + (time | patid) + (1 | trial),
  data = alzheimer_long,
  na.action = na.exclude
)
```

#### Step 2: Predict Missing Values Under MAR

Predict missing BPRS values:

``` r
# Predict MAR-based values
alzheimer_long <- alzheimer_long %>%
  mutate(
    bprs_mar = predict(fit_mar, newdata = alzheimer_long, allow.new.levels = TRUE)
  )
```

#### Step 3: Apply $\delta$-Shift Adjustments

Create MNAR-adjusted datasets:

``` r
# Define delta values
delta_values <- c(0, 0.5, 1, 2, 3)

# Apply delta shifts
mnar_datasets <- lapply(delta_values, function(delta) {
  alzheimer_long %>%
    mutate(
      bprs_mnar = if_else(time > D_i, bprs_mar + delta, bprs_mar)
    )
})
```

#### Step 4: Refit the Model for Each $\delta$

Refit the model for each MNAR scenario:

``` r
# Fit models for each delta
mnar_models <- lapply(mnar_datasets, function(data) {
  lmer(
    bprs_mnar ~ time + sex + age + bmi + job + adl + wzc + cdrsb0 +
      time:cdrsb0 + (time | patid) + (1 | trial),
    data = data,
    na.action = na.exclude
  )
})
```

------------------------------------------------------------------------

### 6) Reporting Results

#### Key Estimates to Compare

- **Time effect** ($\beta_{\text{time}}$),
- **Interaction** ($\beta_{\text{time:cdrsb0}}$).

#### Results Table Template

| $\delta$ | Time Effect ($\beta_{\text{time}}$) | Interaction ($\beta_{\text{time:cdrsb0}}$) |
|----------|-------------------------------------|--------------------------------------------|
| 0 (MAR)  | …                                   | …                                          |
| 0.5      | …                                   | …                                          |
| 1        | …                                   | …                                          |
| 2        | …                                   | …                                          |
| 3        | …                                   | …                                          |

------------------------------------------------------------------------

### 7) Interpretation

- **Robustness**: If estimates remain stable across $\delta$,
  conclusions are robust to MNAR.
- **Sensitivity**: If estimates change significantly, conclusions depend
  on the MNAR assumption.

------------------------------------------------------------------------

### 8) Conclusion

This pattern-mixture sensitivity analysis evaluates the robustness of
BPRS conclusions under plausible MNAR scenarios. The $\delta$-shift
approach is transparent, clinically interpretable, and highlights the
impact of untestable assumptions about missing data.

### Which other missingness models are in our **revision.md**, besides the MAR LMM?

In the revision notes, we essentially laid out **two MNAR sensitivity
families** (plus the MAR primary model):

1.  **Pattern-mixture model (PMM)** sensitivity (your main plan)

    - Define dropout pattern by the last observed visit (D_i).

    - Let the **mean BPRS trajectory depend on the dropout pattern**,
      e.g. a pattern-specific shift and/or slope:

      \[ E(Y\_{ij}D_i)=x\_{ij}^+*{D_i}+*{D_i}t\_{ij},\]

      and then impose *identifying restrictions* for the unobserved part
      (your ()-shift and slope-change ideas).

2.  **Selection model (Diggle–Kenward style)** MNAR sensitivity
    (optional in the notes)

    - Jointly model outcomes and dropout via

      \[ f(y_i,r_iX_i)=f(y_iX_i,),f(r_iX_i,y_i,),\]

      with a dropout hazard that can depend on the **current (possibly
      unobserved) outcome** through an MNAR parameter (your (\_2)).

    - This aligns with the course notes idea that in selection models,
      MNAR corresponds to allowing dropout to depend on the current
      (y\_{ij}) (not just observed history).

A third family is also covered in the course notes (even if we didn’t
build it out in the revision file): **shared-parameter models** (a.k.a.
joint models with shared random effects), where dropout depends on
latent subject-specific effects that also drive the longitudinal
outcome.

------------------------------------------------------------------------

### Which is better for *our* BPRS case, and why?

For your setting (Alzheimer follow-up, **very high dropout ~60% by later
years**, and the scientific goal is “how does BPRS evolve over time?”),
the best *practical* choice is:

**Primary:** MAR **direct likelihood** via the LMM (your Homework 1
model). **Sensitivity:** PMM with **()-adjustment / slope-change**
scenarios.

Here’s why I’d pick **PMM sensitivity** as the main MNAR check:

- **Clinician-facing interpretability:** a ()-shift is literally “after
  dropout, symptoms are on average () points worse (or better) than what
  MAR would predict,” and slope-change is “progression accelerates after
  dropout.” That maps cleanly to clinical stories (dropout because
  worsening, institutionalisation, etc.).

- **With ~60% missingness, extrapolation is unavoidable no matter
  what.** PMM makes that extrapolation explicit and tunable (you show
  *how large* an MNAR departure is needed to change the conclusion).
  That’s exactly the spirit of sensitivity analysis emphasized in the
  course material: MNAR methods require **strong, untestable
  assumptions** and are often most useful as sensitivity tools rather
  than “the one true model.”

- **Selection models are powerful but fragile in practice.** The
  Diggle–Kenward MNAR term (dependence on the unobserved current
  outcome) is not identified from observed data alone, so you end up
  fixing/varying (\_2) and checking robustness—i.e., sensitivity
  analysis again. Also, the Molenberghs et al. paper explicitly warns
  that it’s hard to interpret “evidence for MNAR” using only the
  observed data and notes that complicated MNAR fits can become
  computationally prohibitive.

- **Shared-parameter models** are elegant, but they typically add
  another layer of latent structure that’s harder to defend clinically
  (“dropout depends on an unobserved frailty/random effect”), and with
  heavy dropout they can become very model-dependent.

So: keep the **LMM under MAR** as the clean estimand-driven primary
analysis, and use **PMM ()/slope** as the sensitivity analysis you can
explain to both a statistician *and* a clinician without hand-waving.
That pairing is also consistent with the message from Molenberghs et
al. that likelihood-based ignorable methods are generally more
justifiable than crude rules like CC/LOCF, while MNAR work should be
handled carefully.
