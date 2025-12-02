Differences between GEE and random-effects models for binary
longitudinal data
================

## 1. Big picture: what each model is answering

**GEE (Question 1)**

- Targets a **marginal / population-average** effect:

  > “On average in this Alzheimer population, how does time / tau / ADL
  > affect the probability that CDRSB (\> 10)?”

- Parameters (^{}) describe how the **average probability** changes with
  covariates.

**Random-effects GLMM (Question 2)**

- Targets a **subject-specific / conditional** effect:

  > “For a given patient, after accounting for their own random
  > intercept/slope, how does time / tau / ADL change their individual
  > risk of CDRSB (\> 10)?”

- Parameters (^{}) describe how the **individual log-odds** change,
  conditional on that patient’s random effect.

**Key message for the jury**

> For Gaussian outcomes, marginal and random-effects models can often be
> re-expressed into each other. For **binary** longitudinal data they
> are genuinely different models, targeting different parameters. That’s
> exactly the contrast our homework is illustrating.

------------------------------------------------------------------------

## 2. Interpretation: population-average vs subject-specific

For a **binary** outcome with a **logit** link:

- **GEE**: \[ = X\_{it}^\_M\] (\_M) is a **marginal log-odds ratio**.

- **Random-intercepts GLMM**: \[ = X\_{it}^\_S + b_i\] (\_S) is a
  **subject-specific log-odds ratio** (conditional on that patient’s
  (b_i)).

Important points you can say:

- Subject-specific effects are **larger in magnitude** than marginal
  ones (for non-Gaussian links). Intuitively, averaging over random
  effects **dampens** the effect at the population level.
- There is **no simple exact transformation** between (\_S) and (\_M)
  for logistic models; they answer different questions.
- In your results the two sets of estimates are actually quite similar
  because your random-effect variance is essentially zero (singular
  fit), so the GLMM is almost a marginal model in disguise.

------------------------------------------------------------------------

## 3. Modelling the association / correlation

**GEE**

- Uses a **working correlation** (e.g. AR(1)) to describe ((Y\_{it},
  Y\_{is} X)).

- This working correlation is treated as a **nuisance**. The strength of
  GEE1 is:

  - If the **mean model is correct**, (^{}) is **consistent** even when
    the working correlation is wrong, with robust (sandwich) SEs.

- It only uses **first and second moments**; it does *not* specify a
  full joint distribution for ((Y\_{i1},,Y\_{iT})).

**Random-effects GLMM**

- Specifies a **full hierarchical model**:

  - Conditional on (b_i), responses at different times are independent
    Bernoulli.
  - Correlation arises **through the random effect**.

- Once you choose the **distribution for (b_i)** (normal) and the link,
  you’ve fully specified the joint distribution for the repeated binary
  vector.

- Correlation structure is **implicit**:

  - Larger (\_b^2) ⇒ stronger within-patient correlation.

For your dataset:

- The GLMM estimated **almost zero random-effect variance**, so it
  implies **very weak extra correlation** beyond what time and
  covariates explain.
- Your GEE AR(1) correlation was modest (()), but you treated it as
  working correlation only; the main scientific focus stayed on the
  mean.

------------------------------------------------------------------------

## 4. Estimation methods: quasi vs full likelihood

### GEE (Question 1)

- Uses **quasi-likelihood / estimating equations**, not a full
  likelihood.

- Assumes:

  - Correct mean model,
  - Correct variance function (binomial: ((1-))),
  - Some “reasonable” working correlation.

- () solves a **set of estimating equations**; inference uses the
  **sandwich variance**.

### GLMM (Question 2)

- Fits a **full likelihood**, integrating over random effects:

  - Here: **Laplace approximation** ((nAGQ = 1)) via `glmer`.

- Alternatives mentioned in the notes:

  - **MQL/PQL** (marginal / penalised quasi-likelihood) – simpler but
    biased for binary outcomes with strong random effects and small
    cluster sizes.
  - **Adaptive Gauss–Hermite quadrature** – more accurate but
    computationally heavier.

- You correctly **chose Laplace** rather than MQL/PQL, because:

  - The notes stress that **PQL/MQL can underestimate variance and
    shrink effects** for binary data,
  - Laplace is a standard compromise for GLMMs, giving better estimates
    of both fixed and random effects.

For the defence, you can say:

> “For Question 2 we wanted a proper subject-specific model. In the
> notes, MQL/PQL are mainly described as approximations that can perform
> poorly, especially for binary data. We therefore used a GLMM fitted by
> Laplace approximation, which is closer to true maximum likelihood for
> the random-effects logistic model.”

------------------------------------------------------------------------

## 5. Missingness and dropout

This is a big theme in the notes.

**GEE**

- Standard GEE1 is valid under:

  - **MCAR**, or
  - Some forms of **MAR** where dropout depends only on observed
    covariates and past outcomes, and where the working correlation is
    chosen carefully.

- For more complex MAR dropout (especially if it depends on the
  unobserved outcome), the notes introduce:

  - **Weighted GEE (WGEE)**,
  - **Joint models** (shared parameter),
  - Or hierarchical (random-effects) models.

**GLMM**

- Under a correctly specified hierarchical model, GLMM can handle more
  general forms of **MAR** where:

  - Dropout depends on the random effects (e.g. sicker subjects,
    captured by (b_i), drop out earlier).

- If the random-effects distribution is misspecified, this protection
  can be lost.

For your defence:

- You already showed that **severe patients drop out earlier** (Homework
  1).
- For the GEE part, you assumed **MAR given observed history**, but
  didn’t implement WGEE.
- For the GLMM, you implicitly assumed that dropout is MAR given
  covariates and random effects; however, your random-effect variance is
  almost zero, so in practice both models are working under a fairly
  **weak MAR structure**.
- You can acknowledge:

> “In theory, a random-effects model can be more natural when dropout
> depends on latent subject frailty. In our data the estimated subject
> heterogeneity was almost zero, so the GLMM didn’t buy us much extra
> protection against informative dropout. A more advanced analysis could
> consider WGEE or joint modelling if dropout is strongly driven by
> unobserved worsening.”

------------------------------------------------------------------------

## 6. Efficiency and robustness (this is very ‘Chapter 26’)

**GEE**

- **Robust** to misspecification of the correlation (for ()), but:

  - If the working correlation is very wrong, you lose **efficiency**
    (SEs bigger than necessary).

- Uses only **first two moments** ⇒ robust, but not the most efficient
  if the full distribution is correctly specified.

**Random-effects GLMM**

- If the model is correctly specified (link, covariates, random-effects
  distribution):

  - **More efficient** than GEE, because it uses the **full
    likelihood**.

- But it’s **less robust**:

  - If the random-effects distribution is wrong, both estimates and SEs
    can be biased.
  - Also more sensitive to link mis-specification.

In your homework:

- GEE was a **safe choice** for question 1 (population-average target).

- GLMM for question 2 gave **very small random-effects variances**,
  hinting either:

  - True subject/trial heterogeneity is tiny, *or*
  - The model structure (random slope, interactions, scale of the
    predictors) is not capturing the real heterogeneity.

- Because the random-effects variance is almost zero, the **efficiency
  gain of GLMM over GEE is minimal** in your fitted models.

------------------------------------------------------------------------

## 7. When would you *prefer* one over the other (according to the notes)?

You can memorise and say something like this:

**Prefer GEE when…**

- Main interest is **population-average effects** (public health, policy
  evaluation).
- You have **many subjects, relatively few time points**.
- You want robustness to the exact correlation structure, and are happy
  with marginal odds ratios.
- Prediction at individual level is not the main goal.

**Prefer random-effects GLMM when…**

- You’re interested in **subject-specific trajectories** and
  **heterogeneity** (which patients are systematically worse/better).
- You need **individual predictions** or **empirical Bayes estimates**
  as part of the question (which you did in Q2).
- You want to incorporate **complex hierarchical structures** (patients
  within trials, random slopes, etc.)
- You are willing to specify and defend a full hierarchical model.

------------------------------------------------------------------------

## 8. Critique of your own random-effects analysis (what you can say)

This ties all the theory back to your results:

1.  **Random-effects structure**

    - You fitted ( ( ) ) but your reported variance table suggests
      near-zero variance and a singular fit.
    - In hindsight, a **simpler random intercept model** might have been
      more stable and closer to what the data support.
    - Alternatively, you could argue that the process is mainly driven
      by the fixed effect of time and tau, with little extra
      subject-level variation.

2.  **Trial random intercept**

    - The empirical Bayes plot shows **almost no between-trial
      heterogeneity** after adjustment.

    - That’s fine to report, but you might say:

      > “Given the almost zero variation, the random trial effect is not
      > essential; a fixed trial effect or even no trial effect might
      > have been sufficient.”

3.  **Approximate equivalence to GEE**

    - Because the random-effect variances are tiny, your GLMM
      **collapses to something very close to a marginal model**.
    - That’s why your fixed-effect estimates from GLMM and GEE look
      similar; this is consistent with the theory in the notes that as
      the random-effect variance shrinks, subject-specific and marginal
      effects move closer.

4.  **Missingness**

    - Neither model explicitly models dropout.

    - You can honestly say:

      > “We treated missing data as MAR given observed covariates and
      > outcome history. A more advanced analysis, especially if dropout
      > is strongly outcome-dependent, would require WGEE or a
      > shared-parameter/joint model, as discussed in the notes.”

------------------------------------------------------------------------

## 9. A short 3–4 sentence “defence summary” you can say out loud

If they ask, “How do GEE and random-effects models compare for this
binary longitudinal data?”, you can answer:

> “For binary longitudinal data, GEE and random-effects models are
> fundamentally different. GEE targets **population-average** effects
> and is robust to the exact correlation structure, which is why we used
> it in Question 1 to answer marginal questions about the average
> probability of CDRSB (\>10). The random-effects logistic model is a
> **subject-specific** model that specifies a full hierarchical
> structure; it’s more efficient when correctly specified and gives
> individual-level predictions, which we used in Question 2 to look at
> empirical Bayes random effects. In our data, the estimated
> random-effects variances were almost zero, so the GLMM behaved very
> similarly to the GEE, suggesting limited remaining patient or trial
> heterogeneity after adjustment. A more sophisticated treatment of
> dropout, like weighted GEE or joint modelling, would be needed if we
> wanted to fully address non-random dropout as discussed in Chapter
> 26.”
