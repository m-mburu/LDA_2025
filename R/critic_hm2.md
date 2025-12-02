Critical Appraisal of GLMMs for Longitudinal Data
================

Nice, this is a really solid piece of work already. I’ll go through it
systematically and point out **what is good**, **what is slightly off
compared to the notes**, and **what you could tighten for the defense**.

I’ll follow this order:

1.  Data prep & EDA
2.  Question 1 – GEE model
3.  Question 2 – GLMM / random-effects model
4.  Cross-cutting issues (missing data, correlation, estimation, etc.)
5.  Short “if I had 2–3 days to polish this” checklist

------------------------------------------------------------------------

## 1. Data preparation and EDA

### 1.1 Long-format creation

- **Good**

  - You correctly reshape to long using `melt()` with `patid`, `trial`,
    and baseline covariates as id variables.
  - You construct `CDRSBbin` exactly as required in the assignment.
  - You keep `time` as an integer index 0,…,6 which matches yearly
    visits.

- **Minor issues / possible clarifications**

  - In the notes they often distinguish **baseline covariates** vs
    **time-varying covariates** explicitly.

    - In your long data, `age`, `bmi`, `adl` are baseline (constant over
      `time`), while `abpet`, `taupet` are time-varying.
    - You *say* that for the GEE/GLMM age, BMI, ADL are baseline and
      ABPET/TAUPET time-varying, but this could be made a bit more
      explicit when you write the model (e.g. (X_i) vs (W\_{it})).

  - Check whether in the original SAS file ADL is repeated (adl0–adl6)
    or only baseline. If ADL is actually longitudinal in the notes, you
    should at least **acknowledge** that you only used baseline ADL.

### 1.2 Descriptive plots

- **What you did well**

  - You start with **outcome over time** (stacked bar + sex plot) before
    looking at covariates → that matches the logic in the book.
  - You give **percentages** (75–80% vs 20–25%; 65–70% etc.), not just
    eyeball statements.
  - You connect sex differences to time: small difference at baseline,
    convergence over time.
  - For age/BMI, you correctly emphasise that the boxplots are
    **cross-sectional among survivors**, not full-cohort summaries.
  - You bring in knowledge from Homework 1: older / more severe patients
    drop out → good continuity.

- **Gaps vs what the notes usually emphasise**

  - The assignment is about **non-Gaussian longitudinal data with
    dropout**. The notes usually push you to:

    - Show **overall retention by time** (not only for high baseline
      CDRSB). You did a nice retention plot for those with high CDRSB at
      baseline, but *not* for the full cohort. I would:

      - either add a second curve for “all” patients, or
      - write 1–2 sentences in the text summarising: “overall, only X%
        remain at visit 7”.

    - Make the missingness pattern more explicit: e.g. “dropout is
      mostly monotone; most patients leave permanently after their last
      observed visit”.

  - For the **semivariogram**, you compute it and plot it (good!), but
    you never *verbally* connect it back to the choice of working
    correlation in the main text, only in the “Motivation” section. It
    would help to say right under the plot:

    - “The semivariogram increases with time lag and then levels off,
      consistent with positive within-patient correlation that decays
      with lag. This motivates an AR(1)-type working correlation for the
      GEE model.”

------------------------------------------------------------------------

## 2. Question 1 – GEE (marginal) model

### 2.1 Model formulation

- **Very good**

  - You write down the **mean model clearly**: \[ (\_{it}) = \_0 + *1
    *{it} + \]
  - You state the **binomial variance** and the **AR(1) working
    correlation** explicitly.
  - You clearly state the **target**: population-average effects
    (marginal odds ratios).
  - You mention **sandwich standard errors** and the key robustness
    property – this matches the notes.

- **What could be sharpened**

  - In the notes, with **incomplete data**, there is a subtle point:

    - For some missingness mechanisms under MAR, **independence working
      correlation** is safer than AR(1) because the sandwich correction
      is derived under the “independence estimating equations”.
    - They often warn that with non-independence working correlations
      and missing data, **GEE1 isn’t automatically unbiased** unless
      missingness is quite benign.

  - In your homework you go straight to AR(1) and justify it via the
    semivariogram. That’s good for **efficiency**, but for the defense
    you should add a one-line caveat like:

    > “We chose AR(1) as a reasonable working model to gain efficiency,
    > but we know from the notes that with dropout this can be
    > sensitive; an independence working correlation would be more
    > conservative but less efficient.”

  - You include **time × ABPET** and **time × TAUPET** in the model, but
    you don’t comment on whether these are actually needed /
    significant. The notes usually expect a brief remark like:

    - “We included time-by-biomarker interactions to allow biomarker
      effects on the log-odds to vary over follow-up, but these were not
      strongly supported by the data (p-values …), so for simplicity one
      might also consider a model without interactions.”

### 2.2 Choice of method & estimation

- **Good**

  - You correctly describe GEE1 as a **semi-parametric** method that:

    - targets marginal mean,
    - treats correlation as nuisance,
    - is robust via sandwich SEs.

  - You mention **large-sample normality** and that you have many
    patients (~1250), which supports asymptotics.

- **Missing vs notes**

  - In the theory, they often emphasise **how to choose the working
    correlation**:

    - e.g. by comparing QIC, or checking residual plots.

  - You justify AR(1) by the semivariogram, which is fine, but you
    don’t:

    - compare AR(1) vs independence vs exchangeable, or
    - at least mention that you checked another structure and found
      similar β. For homework this is OK, but in an oral defense they
      can ask “How did you know AR(1) was not harming you?”

### 2.3 Assumptions and pros/cons

- **Very good**

  - You list correctly:

    - Correct mean model,
    - Correct variance function,
    - Independent clusters,
    - MAR/MCAR for missingness,
    - Possible impact of mis-specified working correlation on
      **efficiency** (not on consistency).

  - Your advantages/disadvantages list is basically textbook-level.

- **Two things the notes often stress that you could mention
  explicitly**

  1.  **GEE does *not* produce a likelihood**, so:

      - you cannot use AIC/BIC,
      - you must rely on **QIC** or similar for model comparison. You
        used AIC/BIC later for the GLMM but not for the GEE → worth
        explicitly stating “we did not use AIC/BIC for the GEE, as there
        is no full likelihood.”

  2.  **Interpretation of coefficients**:

      - For logistic GEE, the coefficients are **population-average
        log-odds ratios**, which are generally **smaller in magnitude**
        than subject-specific ones from GLMM.
      - You say “easy to interpret for clinicians”, which is true, but
        you don’t mention this “shrinkage towards 0” relative to GLMM.

### 2.4 Results & interpretation

- **Good**

  - You present:

    - baseline odds / probability from the intercept,
    - clear statement that time has a strong positive effect,
    - tau PET is associated with higher odds, others not so much.

  - You explicitly say **“population-average longitudinal model”** –
    nice.

- **Potential improvements**

  - The notes like **odds-ratio phrasing** rather than raw coefficients.
    You computed ORs in the table but in the text you mostly talk about
    coefficients:

    - E.g. say: “Each additional year is associated with about ((0.25)
      )-fold higher odds of CDRSB \> 10, on average in the population.”

  - You might want to **comment on effect sizes** more: for TAU, OR ~
    1.5 per unit – is that clinically large given the scale of TAUPET?

  - You never discuss the **AR(1) parameter estimate** () in words; you
    could add:

    > “The working correlation estimate α ≈ 0.31 indicates moderate
    > within-patient correlation at lag 1, which decays as lag
    > increases.”

------------------------------------------------------------------------

## 3. Question 2 – GLMM (random-effects model)

Here there are some **important inconsistencies** between what you
claimed and what the code actually fits. This is exactly what an
examiner will pick up.

### 3.1 Model specification vs code

- In the text you say:

  - “We fitted a GLMM with **random intercepts for both patients and
    trials**.”
  - The displayed model equation has: ((*{ijt}) = X*{ijt}^+ b\_{0i} +
    u\_{0j}).

- But the `glmer` call you show is:

  ``` r
  glmm_model <- glmer(
    CDRSBbin ~ time + sex + age + bmi + abpet + taupet+
     adl + time:abpet + time:taupet + 1 + (time | patid),
    ...
  )
  ```

  - This includes a **random intercept and random slope for time per
    patient** `(time|patid)`.
  - It does **not** include `(1|trial)` at all.

- Later you compute:

  ``` r
  ranef_trial <- as.data.frame(r_effs$trial)
  ```

  which would actually **fail** for the shown model, because there is no
  `trial` random effect.

👉 So something changed between the model you actually ran and the code
you pasted. For the exam you must make them **consistent**:

- Either:

  - use `(1 | patid) + (1 | trial)` and describe it as such (random
    intercepts only), **or**
  - keep `(time | patid)` and explicitly say “random intercept and
    random slope for time per patient; trial entered only as a fixed
    effect or not at all.”

Right now the text, code, and random-effects section contradict each
other.

### 3.2 Estimation method

- **Good**

  - You clearly say you used **maximum likelihood with Laplace
    approximation** (`nAGQ = 1`), and that this is preferable to PQL/MQL
    for binary data – this aligns perfectly with the notes.

- **What you could add**

  - The notes sometimes say that for **logistic GLMMs with many
    clusters**, Laplace is usually OK, but if you want more accurate
    variance component estimates you could use `nAGQ > 1` at extra
    computational cost. You can just mention: “We used Laplace (AGQ=1)
    for speed; in principle, AGQ\>1 would give slightly more accurate
    random-effects estimates.”

### 3.3 Random-effects variance / singular fit

- In your “Results” you say:

  > “The estimated random-intercept variances were very small …
  > indicating a boundary (singular) fit.”

  and later:

  > “Empirical Bayes predictions show all random intercepts essentially
  > zero.”

- **This is exactly what the notes would warn about:**

  - When random-effect variances are virtually zero, the data **do not
    support** that random effect.
  - With a singular fit, the GLMM is effectively collapsing back to a
    **fixed-effects logistic regression** (almost no extra clustering
    structure).

- In that situation you should **not talk** about “patient-specific
  trajectories” or “heterogeneity between sites” as if they were
  substantial:

  - because your own results say there is almost no unexplained
    variation left.

- You do say “most of the variation is explained at the overall
  population level” – that’s good – but you don’t draw the obvious
  practical consequence:

  - that **simpler models** (e.g. only fixed effects, or random
    intercept only) may be sufficient.
  - and that the EB predictions are not very informative because they
    are almost all 0.

### 3.4 Residual variance table

- Your `tidy_re_varcorr()` function adds a “Residual” row via
  `sigma(model)` and reports a “residual variance” for the GLMM.
- For binomial GLMMs in `lme4`, `sigma(model)` is not really a
  **residual SD** in the usual sense (the Bernoulli variance is ((1-)),
  and on the logit scale the residual variance is fixed at (^2/3) for a
  logistic distribution).
- The **VarCorr** for a GLMM typically only lists the random-effects
  variances; there is no free residual variance parameter.
- So the “Residual variance” row in Table 6 is a bit misleading; it
  suggests an estimated residual variance like in a Gaussian LMM, which
  is not what you have here. It’s not fatal, but in a defense they could
  ask about it.

### 3.5 Interpretation

- **Good**

  - You correctly say that GLMM coefficients are **subject-specific**
    and not directly comparable to the GEE coefficients.
  - You highlight that **time** again has a strong positive effect,
    consistent with GEE.

- **Could be richer**

  - The notes emphasise the difference between **population-average vs
    subject-specific odds ratios**:

    - subject-specific ORs are typically **larger in magnitude**
      (further from 1) than marginal ORs, because conditioning reduces
      variability.

  - Your GLMM and GEE time coefficients are actually quite similar
    (~0.25), which is interesting:

    - This is partly because your random effects are almost zero.
    - You could mention this explicitly: “Because the random-effects
      variances are tiny, the subject-specific and marginal time effects
      are very similar; in a dataset with stronger heterogeneity we
      would expect larger subject-specific effects.”

### 3.6 Empirical Bayes predictions

- Conceptually, what you do is exactly what the assignment asked:
  extract `ranef()`, plot histograms, comment.

- But because of the **singular fit**, the EB predictions are almost
  trivial: everything near 0.

- That’s fine, but for the defense I would add one sentence connecting
  to the theory:

  > “The near-zero empirical Bayes random effects are consistent with
  > the almost zero variance components, and they show that the GLMM is
  > effectively a fixed-effects model for this outcome. In practice it
  > means that, once we adjust for the covariates, there is little
  > evidence of extra clustering by patient or trial.”

------------------------------------------------------------------------

## 4. Cross-cutting issues

### 4.1 Missing data / dropout

- You **acknowledge** MAR/MCAR assumptions in both GEE and GLMM
  sections.
- You show a retention plot for high-CDRSB patients at baseline – that’s
  nice and specific.

What the notes would still want:

- Some **explicit statement** about the likely missingness type:

  - “Given that dropout is higher in older and more severe patients
    (from Homework 1 results), missingness is plausibly MAR given
    baseline age, severity and previous outcomes, but could still be
    MNAR.”

- And a short remark about **methods if MAR fails**:

  - e.g. “If dropout were strongly MNAR, both GEE and GLMM could be
    biased, and we would need joint models or pattern-mixture models,
    which are beyond this homework.”

### 4.2 Correlation structure vs semi-variogram

- You computed a **semi-variogram** for the binary outcome in the GEE
  part and saw increasing semivariance with lag.

- In the GLMM, the correlation is induced by the random intercept and
  slope; you do not explore it.

- The notes often like to see that **your chosen model is consistent
  with the exploratory dependence**:

  - You can just say: “A random intercept + slope structure produces
    decaying within-patient correlation, which is in line with the
    semivariogram we saw for CDRSBbin; so the GLMM’s implied correlation
    is qualitatively compatible with the data.”

------------------------------------------------------------------------

## 5. Polishing checklist (what I’d fix if I had time)

If you want to make this bullet-proof for the defense, I’d do:

1.  **Align model text and code for GLMM**

    - Decide: `(1|patid) + (1|trial)` *or* `(time|patid)` (with or
      without `(1|trial)`).
    - Rewrite both the formula and the prose to match exactly.
    - Re-run `VarCorr` and EB plots accordingly.

2.  **Acknowledge singular fit more strongly**

    - Explicitly state: “Because the random-effect variances are at the
      boundary, simpler models without random effects are likely
      adequate. We still show the GLMM as an illustration of
      subject-specific modelling, but interpret random effects with
      caution.”

3.  **Tweak GEE working-correlation discussion**

    - Add one line: “Independence working correlation is the safest
      option under dropout; we chose AR(1) for efficiency, but checked
      that the main conclusions are similar under independence.”

4.  **Minor theoretical points to add**

    - GEE has no likelihood → can’t use AIC/BIC; QIC would be the formal
      tool.
    - Subject-specific vs marginal OR difference, and why they are
      similar here (tiny RE variances).

5.  **Clean up the residual variance table for GLMM**

    - Either remove the “Residual” row or explain that for binomial GLMM
      there is no free residual variance like in Gaussian models.

------------------------------------------------------------------------

Overall: conceptually you are in very good shape. The main weaknesses
are **consistency and small theoretical caveats**, not huge conceptual
mistakes. If you fix the GLMM description + explicitly acknowledge the
singular fit and GEE–dropout nuance, you’ll be very hard to catch out in
the oral defense.
