# Longitudinal Data Analysis

## Assignment 1: Continuous Longitudinal Data

## 2025–

## 1 Introduction on Alzheimer data

A group of elderly patients, spread across 20 clinical centers, are being followed over time for Alzheimers
disease. Some patients still live at home, while others livein nursing homes (WZC). Patients are included
in the study upon first diagnosis of Alzheimers disease. Naturally, this occurs at different ages for different
individuals.
The longitudinal follow-up is based on four measurements. First, there are two cognitive scales (CDR-SB and
BPRS); second, there are two possible biomarkers (ABPET andTAUPET). Each of these four variables is
measured at baseline (upon inclusion, moment 0), and then annually over six years.
In addition to these longitudinal measurements, several background and clinical parameters are recorded at
baseline only. These include age (AGE), sex (SEX), education level (EDU), body mass index (BMI), whether
the person has a paid job (JOB), activities of daily living (ADL), and whether they live in a nursing home
(WZC).
Note that not all participants remain in the study until the end. Baseline values are recorded for everyone, but
from the first follow-up moment (year 1) onwards, some participants drop out and do not return to the study.

The objective of the study is to investigate how BPRS changesover time, and how it relates to the patient
characteristics at baseline as well as the baseline measurements of CDR-SB, ABPET, and TAUPET.

## 2 Data file

- SAS file Alzheimer25.sas7bdat
- Variables:
    1. TRIAL: number of the participating clinical center
    2. PATID: unique patient ID within the study (consisting of acenter code and a patient sequence
       number)
    3. SEX: 0 = male, 1 = female
    4. AGE: age at inclusion (in years)
    5. EDU: highest education level (1: Primary education; 2: Lower secondary education; 3: Upper
       secondary education; 4: Higher education)
    6. BMI: body mass index at baseline (kg/m^2 )
    7. INKOMEN: monthly net disposable income (Euro)
    8. JOB: 1 if the person has a paid job; 0 otherwise


9. ADL: “Activities of daily living” at baseline
10. WZC: 1 for nursing home residents; 0 otherwise
11. CDRSB0 – CDRSB6: Clinical Dementia Rating Sum of Boxes atbaseline (0) and years 1–
12. BPRS0 – BPRS6: Brief Psychiatric Rating Scale at baseline (0) and years 1–
13. ABPET0 – ABPET6: Amyloid- Positron Emission Tomographyat baseline (0) and years 1–
14. TAUPET0 – TAUPET6: Tau Positron Emission Tomography at baseline (0) and years 1–

## 3 Assignments

1. Describe the data, and use graphical techniques to explore the mean structure, the variance structure
    and the correlation structure. Summarize your conclusions. What are the implications with respect to
    statistical modeling?
2. What summary statistics are appropriate for the analysisof these data? Why? Do they yield the same
    results? Summarize your conclusions.
3. Fit a multivariate model, and find the most parsimonious mean structure which can be used to describe
    the average evolutions in the data. What covariance structures are applicable in this case? What is the
    most parsimonious structure you can find?
4. Use an explicit two-stage analysis to get an initial impression about trends and the effect of covariates
    on those trends.
5. Formulate a plausible random-effects model. Fit your model, and compare the results with those from
    the multivariate model. Check the appropriateness of your random-effects model. Calculate the subject-
    specific intercepts/slopes and compare them to the ones you obtained from a two-stage analysis. What
    do you conclude?

## 4 General remarks

- For each question, motivate your choice of techniques, estimation methods, assumptions you make, and
    describe possible advantages/disadvantages, problems.
- For each of the above questions, summarize your conclusionsand report them to a clinician.
- Carefully reflect on the parameterization of your models.
- Do you have any recommendations with respect to future similar experiments?
- The deadline for submitting your report isNovember 21, 2025, 9am. Submit your report by mailing
    it toboth instructors(geert.verbeke@kuleuven.be, geert.molenberghs@uhasselt.be).
- As soon as your team composition is final, submit the composition toboth instructors, and indicate for
    each member of your team whether he/she takes the course on campus or online, and whether he/she
    will present on November 21 or on December 4. Do thisno later than November 14, 2025.


