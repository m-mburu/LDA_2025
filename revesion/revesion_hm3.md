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
