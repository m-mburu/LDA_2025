
Introduction to Longitudinal Data Analysis
Geert Molenberghs
Geert Verbeke
geert.molenberghs@uhasselt.be
geert.verbeke@kuleuven.be
geert.molenberghs@kuleuven.be
Interuniversity Institute for Biostatistics and statistical Bioinformatics (I-BioStat)
Universiteit Hasselt & KU Leuven, Belgium
www.ibiostat.be
Interuniversity Institute for Biostatistics
and statistical Bioinformatics
Master in Statistics, UHasselt & KU Leuven
Contents
1
Related References . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
1
I Continuous Longitudinal Data
10
2
3
4
5
6
Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 11
Cross-sectional versus Longitudinal Data . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 42
Simple Methods . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 55
The Multivariate Regression Model . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 65
A Model for Longitudinal Data . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 93
Introduction to Longitudinal Data Analysis
i
7
8
9
Exploratory Data Analysis . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 124
Estimation of the Marginal Model . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 156
Inference for the Marginal Model
·	. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 182
10
Inference for the Random Eﬀects . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 237
11 General Guidelines for Model Building . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 261
12 Power Analyses under Linear Mixed Models . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 294
II Marginal Models for Non-Gaussian Longitudinal Data
313
13 The Toenail Data . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 314
14 The Analgesic Trial . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 317
15 The National Toxicology Program (NTP) Data . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 320
16 Generalized Linear Models . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 323
17 Parametric Modeling Families . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 341
Introduction to Longitudinal Data Analysis
ii
18 Conditional Models . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 345
19 Full Marginal Models . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 361
20 Generalized Estimating Equations . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 373
21 A Family of GEE Methods . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 385
III Generalized Linear Mixed Models for Non-Gaussian Longitudinal Data
405
22 The Beta-binomial Model . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 406
23 Generalized Linear Mixed Models (GLMM) . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 412
24 Fitting GLMM’s in SAS . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 430
IV Marginal Versus Random-eﬀects Models and Case Studies
440
25 Marginal Versus Random-eﬀects Models . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 441
26 Case Study: The NTP Data . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 454
27 Case Study: Binary Analysis of Analgesic Trial . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 476
Introduction to Longitudinal Data Analysis
iii
28 Case Study: Ordinal Analysis of Analgesic Trial . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 490
29 Count Data: The Epilepsy Study . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 505
V Incomplete Data
531
30 Setting The Scene . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 532
31 Proper Analysis of Incomplete Data . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 558
32 Analysis of the ARMD Trial . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 580
33 Weighted Generalized Estimating Equations . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 589
34 Multiple Imputation . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 618
35 Creating Monotone Missingness . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 641
VI Topics in Methods and Sensitivity Analysis for Incomplete Data
654
36 An MNAR Selection Model and Local Inﬂuence . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 655
37
Local Inﬂuence for the ARMD Trial . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 669
Introduction to Longitudinal Data Analysis
iv
38 Mechanism for Growth Data . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 680
39
Interval of Ignorance / Bodyguard . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 683
40 Pattern-mixture Models . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 726
41 PMM Analysis of the ARMD Trial
·	. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 741
42 Sensitivity Analysis Based on Multiple Imputation . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 742
43 Overview . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 748
Introduction to Longitudinal Data Analysis
v
Chapter 1
Related References
·	Aerts, M., Geys, H., Molenberghs, G., and Ryan, L.M. (2002). Topics in
Modelling of Clustered Data. London: Chapman & Hall.
·	Brown, H. and Prescott, R. (1999). Applied Mixed Models in Medicine. New
York: John Wiley & Sons.
·	Carpenter, J.R. and Kenward, M.G. (2013). Multiple Imputation and its
Application. New York: John Wiley & Sons.
·	Crowder, M.J. and Hand, D.J. (1990). Analysis of Repeated Measures. London:
Chapman & Hall.
Introduction to Longitudinal Data Analysis
1
·	Davidian, M. and Giltinan, D.M. (1995). Nonlinear Models For Repeated
Measurement Data. London: Chapman & Hall.
·	Davis, C.S. (2002). Statistical Methods for the Analysis of Repeated
Measurements. New York: Springer.
·	Demidenko, E. (2004). Mixed Models: Theory and Applications. New York: John
Wiley & Sons.
·	Diggle, P.J., Heagerty, P.J., Liang, K.Y. and Zeger, S.L. (2002). Analysis of
Longitudinal Data. (2nd edition). Oxford: Oxford University Press.
·	Fahrmeir, L. and Tutz, G. (2002). Multivariate Statistical Modelling Based on
Generalized Linear Models (2nd ed). New York: Springer.
·	Fitzmaurice, G.M., Davidian, M., Verbeke, G., and Molenberghs, G.(2009).
Longitudinal Data Analysis. Handbook. Hoboken, NJ: John Wiley & Sons.
Introduction to Longitudinal Data Analysis
2
·	Fitzmaurice, G.M., Laird, N.M., and Ware, J.H. (2004). Applied Longitudinal
Analysis. New York: John Wiley & Sons.
·	Ga lecki, A. and Burzykowski, T. (2013). Linear Mixed-Eﬀects Models Using R.
New York: Springer.
·	Goldstein, H. (1979). The Design and Analysis of Longitudinal Studies. London:
Academic Press.
o	Goldstein, H. (1995). Multilevel Statistical Models. London: Edward Arnold.
Hand, D.J. and Crowder, M.J. (1995). Practical Longitudinal Data Analysis.
London: Chapman & Hall.
·	Hedeker, D. and Gibbons, R.D. (2006). Longitudinal Data Analysis. New York:
John Wiley & Sons.
Hogan, J. and Daniels, M. (2008) Missing Data in Longitudinal Studies:
·	Introduction to Longitudinal Data Analysis
3
Strategies for Bayesian modelling and sensitivity analysis. Boca Raton:
CRC/Chapman & Hall.
·	Jones, B. and Kenward, M.G. (1989). Design and Analysis of Crossover Trials.
London: Chapman & Hall.
·	Kshirsagar, A.M. and Smith, W.B. (1995). Growth Curves. New York: Marcel
Dekker.
·	Leyland, A.H. and Goldstein, H. (2001) Multilevel Modelling of Health Statistics.
Chichester: John Wiley & Sons.
·	Lindsey, J.K. (1993). Models for Repeated Measurements. Oxford: Oxford
University Press.
·	Littell, R.C., Milliken, G.A., Stroup, W.W., Wolﬁnger, R.D., and Schabenberger,
O. (2005). SAS for Mixed Models (2nd ed.). Cary: SAS Press.
Introduction to Longitudinal Data Analysis
4
·	Little, R.J.A., D’Agostino, R., Dickerson, K., Emerson, S.S., Farrar, J.T.,
Frangakis, C., Hogan, J.W., Molenberghs, G., Murphy, S.A., Neaton, J.D.,
Rotnitzky, A., Scharfstein, D., Shih, W., Siegel, J.P., and Stern, H. National
Research Council (2010). The Prevention and Treatment of Missing Data in
Clinical Trials. Panel on Handling Missing Data in Clinical Trials. Committee on
National Statistics, Division of Behavioral and Social Sciences and Education.
Washington, D.C.: The National Academies Press.
·	Little, R.J.A. and Rubin, D.B. (1987, 2002, 2014). Statistical Analysis with
Missing Data (2nd ed.). New York: John Wiley & Sons.
·	Longford, N.T. (1993). Random Coeﬃcient Models. Oxford: Oxford University
Press.
·	Mallinckrodt, C.H. (2013). Preventing and Treating Missing Data in Longitudinal
Clinical Trials: A Practical Guide. New York: Cambridge University Press.
McCullagh, P. and Nelder, J.A. (1989). Generalized Linear Models (second
·	Introduction to Longitudinal Data Analysis
5
edition). London: Chapman & Hall.
·	McLachlan, G.J. and Krishnan, T. (1997). The EM Algorithm and Extensions.
New York: John Wiley & Sons.
·	Mallinckrodt, C.H. (2013). Preventing and Treating Missing Data in Longitudinal
Clinical Trials: A Practical Guide. New York: Cambridge University Press.
·	Molenberghs, G., Fitzmaurice, G., Kenward, M.G., Tsiatis, A.A., and Verbeke, G.
(2015). Handbook of Missing Data. Boca Raton: Chapman & Hall/CRC.
·	Molenberghs, G. and Kenward, M.G. (2007). Missing Data in Clinical Studies.
Chichester: John Wiley & Sons.
·	Molenberghs, G. and Verbeke, G. (2005). Models for Discrete Longitudinal Data.
New York: Springer.
O’Kelly, M. and Ratitch, B. (2014). A Guide to Planning for Missing Data. New
·	Introduction to Longitudinal Data Analysis
6
York: John Wiley & Sons.
·	Pinheiro, J.C. and Bates D.M. (2000). Mixed eﬀects models in S and S-Plus.
New York: Springer.
·	Raghunathan, T. (2017). Missing Data Analysis in Practice. London:
CRC/Chapman & Hall.
·	Rizopoulos, D. (2012). Joint Models for Longitudinal and Time-to-Event Data.
With Applications in R. Boca Raton: Chapman & Hall/CRC.
·	Rubin, D.B. (1987). Multiple Imputation for Nonresponse in Surveys. New York:
John Wiley & Sons.
·	Schafer J.L. (1997). Analysis of Incomplete Multivariate Data. London:
CRC/Chapman & Hall.
Searle, S.R., Casella, G., and McCulloch, C.E. (1992). Variance Components.
·	Introduction to Longitudinal Data Analysis
7
New-York: Wiley.
o	Senn, S.J. (1993). Cross-over Trials in Clinical Research. Chichester: Wiley.
Tan, M.T., Tian, G.-L., and Ng, K.W. (2010). Bayesian Missing Data Problems.
Boca Raton: Chapman & Hall/CRC.
·	van Buuren, S. (2012). Flexible Imputation of Missing Data. Boca Raton:
Chapman & Hall/CRC.
·	Verbeke, G. and Molenberghs, G. (1997). Linear Mixed Models In Practice: A
SAS Oriented Approach, Lecture Notes in Statistics 126. New-York: Springer.
·	Verbeke, G. and Molenberghs, G. (2000). Linear Mixed Models for Longitudinal
Data. Springer Series in Statistics. New-York: Springer.
·	Vonesh, E.F. and Chinchilli, V.M. (1997). Linear and Non-linear Models for the
Analysis of Repeated Measurements. Basel: Marcel Dekker.
Introduction to Longitudinal Data Analysis
8
o	Weiss, R.E. (2005). Modeling Longitudinal Data. New York: Springer.
West, B.T., Welch, K.B., and Ga lecki, A.T. (2007). Linear Mixed Models: A
Practical Guide Using Statistical Software. Boca Raton: Chapman & Hall/CRC.
·	Wu, H. and Zhang, J.-T. (2006). Nonparametric Regression Methods for
Longitudinal Data Analysis. New York: John Wiley & Sons.
·	Wu, L. (2010). Mixed Eﬀects Models for Complex Data. Boca Raton: Chapman
& Hall/CRC.
Introduction to Longitudinal Data Analysis
9
Part I
Continuous Longitudinal Data
Introduction to Longitudinal Data Analysis
10
Chapter 2
Introduction
·	Repeated Measures / Longitudinal data
·	Examples
Introduction to Longitudinal Data Analysis
11
2.1 Repeated Measures / Longitudinal Data
Repeated measures are obtained when a response
is measured repeatedly on a set of units
Units:
·	Subjects, patients, participants, . . .
·	Animals, plants, . . .
·	Clusters: families, towns, branches of a company,. . .
·	. . .
Special case: Longitudinal data
o	Introduction to Longitudinal Data Analysis
12
2.2 Captopril Data
·	Taken from Hand, Daly, Lunn,
McConway, & Ostrowski (1994)
o	15 patients with hypertension
The response of interest is the supine
blood pressure, before and after
treatment with CAPTOPRIL
Pati¨ent
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
Before
SBP DBP
130
210
122
169
124
187
104
160
112
167
101
176
121
185
124
206
115
173
102
146
98
174
119
201
106
198
107
148
100
154
After
SBP DBP
125
201
121
165
121
166
106
157
101
147
85
145
98
168
105
180
103
147
98
136
90
151
98
168
110
179
103
129
82
131
Introduction to Longitudinal Data Analysis
13
o	Research question:
How does treatment aﬀect BP ?
Remarks:
·	Paired observations:
Most simple example of longitudinal
data
·	Much variability between subjects
Introduction to Longitudinal Data Analysis
14
2.3 Growth Curves
§	Taken from Goldstein 1979
The height of 20 schoolgirls, with small, medium, or tall mothers, was measured
over a 4-year period:
Mothers height Children numbers
Small mothers
< 155 cm
Medium mothers [155cm; 164cm]
Tall mothers
164 cm
1
7
14
→
→
→
6
13
20
Is growth related to height of mother ?
Research question:
Introduction to Longitudinal Data Analysis
15
Individual proﬁles:
·	Introduction to Longitudinal Data Analysis
16
·	Remarks:
·	Almost perfect linear relation between Age and Height
·	Much variability between girls
·	Little variability within girls
·	Fixed number of measurements per subject
·	Measurements taken at ﬁxed time points
Introduction to Longitudinal Data Analysis
17
2.4 Growth Data
Taken from Potthoﬀ and Roy, Biometrika (1964)
The distance from the center of the pituitary to the maxillary ﬁssure was recorded
at ages 8, 10, 12, and 14, for 11 girls and 16 boys
Research question:
§	Is dental growth related to gender ?
Introduction to Longitudinal Data Analysis
18
Individual proﬁles:
·	Introduction to Longitudinal Data Analysis
19
·	Remarks:
·	Much variability between children
·	Considerable variability within children
·	Fixed number of measurements per subject
·	Measurements taken at ﬁxed time points
Introduction to Longitudinal Data Analysis
20
2.5 Rat Data
o	Research question (Dentistry, K.U.Leuven):
How does craniofacial growth
depend on testosteron production ?
Randomized experiment in which 50 male Wistar rats are randomized to:
·	Control (15 rats)
·	Low dose of Decapeptyl (18 rats)
·	High dose of Decapeptyl (17 rats)
Introduction to Longitudinal Data Analysis
21
·	Treatment starts at the age of 45 days; measurements taken every 10 days, from
day 50 on.
·	The responses are distances (pixels) between well deﬁned points on x-ray pictures
of the skull of each rat:
Introduction to Longitudinal Data Analysis
22
·	Measurements with respect to the roof, base and height of the skull. Here, we
consider only one response, reﬂecting the height of the skull.
Individual proﬁles:
·	Introduction to Longitudinal Data Analysis
23
o	Complication: Dropout due to anaesthesia (56%):
Age (days)
50
60
70
80
90
100
110
Observations
Control Low High
15
13
13
10
7
4
4
18
17
15
15
12
10
8
17
16
15
13
10
10
10
Total
50
46
43
38
29
24
22
Remarks:
·	Much variability between rats, much less variability within rats
·	Fixed number of measurements scheduled per subject, but not all
measurements available due to dropout, for known reason.
·	Measurements taken at ﬁxed time points
Introduction to Longitudinal Data Analysis
24
2.6 Toenail Data
·	Reference: De Backer, De Keyser, De Vroey, Lesaﬀre, British Journal of
Dermatology (1996).
·	Toenail Dermatophyte Onychomycosis: Common toenail infection, diﬃcult to
treat, aﬀecting more than 2% of population.
·	Classical treatments with antifungal compounds need to be administered until the
whole nail has grown out healthy.
o	New compounds have been developed which reduce treatment to 3 months
Randomized, double-blind, parallel group, multicenter study for the comparison of
two such new compounds (A and B) for oral treatment.
Introduction to Longitudinal Data Analysis
25
Research question:
Are both treatments equally eﬀective
for the treatment of TDO ?
189 patients randomized, 36 centers
2
×
48 weeks of total follow up (12 months)
12 weeks of treatment (3 months)
Measurements at months 0, 1, 2, 3, 6, 9, 12.
o	Introduction to Longitudinal Data Analysis
26
Response considered here: Unaﬀected nail length (mm):
·	Introduction to Longitudinal Data Analysis
27
As response is related to toe size, we restrict to patients with big toenail as target
nail =
150 and 148 subjects.
⇒
30 randomly selected proﬁles, in each group:
o	Introduction to Longitudinal Data Analysis
28
o	Complication: Dropout (24%):
Observations
Time (months)
Treatment A Treatment B
Total
0
1
2
3
6
9
12
150
149
146
140
131
120
118
148
142
138
131
124
109
108
298
291
284
271
255
229
226
Remarks:
·	Much variability between subjects
·	Much variability within subjects
·	Fixed number of measurements scheduled per subject, but not all
measurements available due to dropout, for unknown reason.
·	Measurements taken at ﬁxed time points
Introduction to Longitudinal Data Analysis
29
2.7 Mastitis in Dairy Cattle
Taken from Diggle & Kenward, Applied statistics (1994)
Mastitis : Infectious disease, typically reducing milk yields
Research question:
§	Are high yielding cows more susceptible ?
·	Hence, is the probability of occurrence of mastitis related to the yield that would
have been observed had mastitis not occured ?
Hypothesis cannot be tested directly since ‘covariate is missing for all events’
·	Introduction to Longitudinal Data Analysis
30
o	Individual proﬁles:
Remarks:
·	Paired observations: Most simple
example of longitudinal data
·	Much variability between cows
·	Missingness process itself is
of interest
Introduction to Longitudinal Data Analysis
31
2.8 The Baltimore Longitudinal Study of Aging (BLSA)
·	Reference: Shock, Greullich, Andres, Arenberg, Costa, Lakatta, & Tobin, National
Institutes of Health Publication, Washington, DC: National Institutes of Health
(1984).
·	BLSA: Ongoing, multidisciplinary observational study, started in 1958, with the
study of normal human aging as primary objective
·	Participants:
·	volunteers, predominantly white, well educated, and ﬁnancially comfortable
·	return approximately every 2 years for 3 days of biomedical and psychological
examinations
·	at ﬁrst only males (over 1500 by now), later also females
·	an average of almost 7 visits and 16 years of follow-up
Introduction to Longitudinal Data Analysis
32
o	The BLSA is a unique resource for rapidly evaluating longitudinal hypotheses:
·	data from repeated clinical examinations
·	a bank of frozen blood and urine samples
Drawbacks of such observational studies:
·	More complicated analyses needed (see later)
·	Observed evolutions may be highly inﬂuenced by many covariates which may or
may not be recorded in the study
Introduction to Longitudinal Data Analysis
33
2.8.1 Prostate Data
·	References:
·	Carter et al (1992, Cancer Research).
·	Carter et al (1992, Journal of the American Medical Association).
·	Morrell et al (1995, Journal of the American Statistical Association).
·	Pearson et al (1994, Statistics in Medicine).
·	Prostate disease is one of the most common and most costly medical problems in
the United States
o	Important to look for markers which can detect the disease at an early stage
Prostate-Speciﬁc Antigen is an enzyme produced by both normal and cancerous
prostate cells
Introduction to Longitudinal Data Analysis
34
o	PSA level is related to the volume of prostate tissue.
Problem: Patients with Benign Prostatic Hyperplasia also have an increased PSA
level
·	Overlap in PSA distribution for cancer and BPH cases seriously complicates the
detection of prostate cancer.
Research question (hypothesis based on clinical practice):
·	Can longitudinal PSA proﬁles be used to
detect prostate cancer in an early stage ?
Introduction to Longitudinal Data Analysis
35
·	A retrospective case-control study based on frozen serum samples:
·	16 control patients
·	20 BPH cases
·	14 local cancer cases
·	4 metastatic cancer cases
·	Complication: No perfect match for age at diagnosis and years of follow-up
possible
·	Hence, analyses will have to correct for these age diﬀerences between the
diagnostic groups.
Introduction to Longitudinal Data Analysis
36
Individual proﬁles:
·	Introduction to Longitudinal Data Analysis
37
·	Remarks:
·	Much variability between subjects
·	Little variability within subjects
·	Highly unbalanced data
Introduction to Longitudinal Data Analysis
38
2.8.2 Hearing Data
·	References:
·	Brant & Fozard, Journal of the Acoustic Society of America (1990).
·	Morrell & Brant, Statistics in Medicine (1991).
Hearing thresholds, by means of sound proof chamber and Bekesy audiometer
11 frequencies : 125
8000 Hz, both ears
→
Research question:
How does hearing depend on aging ?
Introduction to Longitudinal Data Analysis
39
Data considered here:
·	500 Hz
·	6170 observations (3089 left ear, 3081 right ear) from 681 males without any
otologic disease
·	followed for up to 22 years, with a maximum of 15 measurements/subject
30 randomly selected proﬁles, for each ear:
o	Introduction to Longitudinal Data Analysis
40
·	Remarks:
·	Much variability between subjects
·	Much variability within subjects
·	Highly unbalanced data
Introduction to Longitudinal Data Analysis
41
Chapter 3
Cross-sectional versus Longitudinal Data
·	Introduction
·	Paired verus unpaired t-test
·	Cross-sectional versus longitudinal data
Introduction to Longitudinal Data Analysis
42
3.1 Introduction
The examples have illustrated several aspects of longitudinal data structures:
·	Experimental and observational
·	Balanced and unbalanced
·	With or without missing data (dropout)
Often, there is far more variability between subjects than within subjects.
This is also reﬂected in correlation within units
§	Introduction to Longitudinal Data Analysis
43
·	For example, for the growth curves, the correlation matrix of the 5 repeated
measurements equals
1.00 0.95 0.96 0.93 0.87
0.95 1.00 0.97 0.96 0.89
0.96 0.97 1.00 0.98 0.94
0.93 0.96 0.98 1.00 0.98
0.87 0.89 0.94 0.98 1.00




























o	This correlation structure cannot be ignored in the analyses (Section 3.2)
The advantage however is that longitudinal data allow to study changes within
subjects (Section 3.3).
Introduction to Longitudinal Data Analysis
44
3.2 Paired versus Unpaired t-test
3.2.1 Paired t-test
The simplest case of longitudinal data are paired data
We re-consider the diastolic blood pressures from the Captopril data
The data can be summarized as:
§	Introduction to Longitudinal Data Analysis
45
There is an average decrease of more than 9 mmHG
The classical analysis of paired data is based on comparisons within subjects:
o	∆i = Yi1 −
Yi2,
i = 1, . . . , 15
·	A positive ∆i corresponds to a decrease of the BP, while a negative ∆i is
equivalent to an increase.
·	Testing for treatment eﬀect is now equivalent to testing whether the average
diﬀerence µ∆ equals zero.
Introduction to Longitudinal Data Analysis
46
Statistica output:
o	Hence, the average change in BP is statistically, signiﬁcantly diﬀerent from zero
(p = 0.001).
Introduction to Longitudinal Data Analysis
47
3.2.2 Unpaired, Two-sample, t-test
·	What if we had ignored the paired nature of the data ?
We then could have used a two-sample (unpaired) t-test to compare the average
BP of untreated patients (controls) with treated patiens.
We would still have found a signiﬁcant diﬀerence (p = 0.0366), but the p-value
larger compared to the one obtained using the
would have been more than 30
paired t-test (p = 0.001).
×
Conclusion:
15 × 2 6= 30 × 1
Introduction to Longitudinal Data Analysis
48
·	The two-sample t-test does not take into account the fact that the 30
measurements are not independent observations.
·	This illustrates that classical statistical models which assume independent
observations will not be valid for the analysis of longitudinal data
Introduction to Longitudinal Data Analysis
49
3.3 Cross-sectional versus Longitudinal Data
Suppose it is of interest to study the relation between some response Y and age
A cross-sectional study yields the following data:
The graph suggests a negative relation between Y and age.
§	Introduction to Longitudinal Data Analysis
50
·	Exactly the same observations could also have been obtained in a longitudinal
study, with 2 measurements per subject.
First case:
·	Are we now still inclined to conclude that there is a negative
relation between Y and Age ?
Introduction to Longitudinal Data Analysis
51
·	The graph suggests a negative cross-sectional relation but a positive longitudinal
trend.
Second case:
o	The graph now suggests the cross-sectional as well as longitudinal trend to be
negative.
Introduction to Longitudinal Data Analysis
52
Conclusion:
·	Longitudinal data allow to distinguish diﬀerences between
subjects from changes within subjects
Application: Growth curves for babies (next page)
·	Introduction to Longitudinal Data Analysis
53
Introduction to Longitudinal Data Analysis
54
Chapter 4
Simple Methods
·	Introduction
·	Overview of frequently used methods
·	Summary statistics
Introduction to Longitudinal Data Analysis
55
4.1 Introduction
·	The reason why classical statistical techniques fail in the context of longitudinal
data is that observations within subjects are correlated.
·	In many cases the correlation between two repeated measurements decreases as
the time span between those measurements increases.
o	A correct analysis should account for this
The paired t-test accounts for this by considering subject-speciﬁc diﬀerences
∆i = Yi1 −
Yi2.
·	This reduces the number of measurements to just one per subject, which implies
that classical techniques can be applied again.
Introduction to Longitudinal Data Analysis
56
o	In the case of more than 2 measurements per subject, similar simple techniques
are often applied to reduce the number of measurements for the ith subject, from
ni to 1.
Some examples:
·	Analysis at each time point separately
·	Analysis of Area Under the Curve (AUC)
·	Analysis of endpoints
·	Analysis of increments
·	Analysis of covariance
Introduction to Longitudinal Data Analysis
57
4.2 Overview of Frequently Used Methods
4.2.1 Analysis at Each Time Point
§	The data are analysed at each occasion separately.
Advantages:
·	Simple to interpret
·	Uses all available data
Disadvantages:
·	Does not consider ‘overall’ diﬀerences
·	Does not allow to study evolution diﬀerences
·	Problem of multiple testing
Introduction to Longitudinal Data Analysis
58
4.2.2 Analysis of Area Under the Curve
For each subject, the area under its curve is calculated :
AU Ci = (ti2 −
ti1)
(yi1 + yi2)/2 + (ti3 −
×
ti2)
×
(yi2 + yi3)/2 + . . .
Afterwards, these AU Ci are analyzed.
Advantages:
·	No problems of multiple testing
·	Does not explicitly assume balanced data
·	Compares ‘overall’ diﬀerences
Disadvantage: Uses only partial information : AU Ci
·	Introduction to Longitudinal Data Analysis
59
4.2.3 Analysis of Endpoints
·	In randomized studies, there are no systematic diﬀerences at baseline.
Hence, ‘treatment’ eﬀects can be assessed by only comparing the measurements
at the last occasion.
Advantages:
·	No problems of multiple testing
·	Does not explicitly assume balanced data
Disadvantages:
·	Uses only partial information : yini
·	Only valid for large data sets
Introduction to Longitudinal Data Analysis
60
4.2.4 Analysis of Increments
o	A simple method to compare evolutions between subjects, correcting for
diﬀerences at baseline, is to analyze the subject-speciﬁc changes yini −
yi1.
Advantages:
·	No problems of multiple testing
·	Does not explicitly assume balanced data
Disadvantage: Uses only partial information : yini −
·	yi1
Introduction to Longitudinal Data Analysis
61
4.2.5 Analysis of Covariance
§	Another way to analyse endpoints, correcting for diﬀerences at baseline, is to use
analysis of covariance techniques, where the ﬁrst measurement is included as
covariate in the model.
Advantages:
·	No problems of multiple testing
·	Does not explicitly assume balanced data
Disadvantages:
·	Uses only partial information : yi1 and yini
·	Does not take into account the variability of yi1
Introduction to Longitudinal Data Analysis
62
4.3 Summary Statistics
§	The AUC, endpoints and increments are examples of summary statistics
Such summary statistics summarize the vector of repeated measurements for each
subject separately.
This leads to the following general procedure :
·	Step 1 : Summarize data of each subject into one statistic, a summary
statistic
·	Step 2 : Analyze the summary statistics, e.g. analysis of covariance to
compare groups after correction for important covariates
·	This way, the analysis of longitudinal data is reduced to the analysis of
independent observations, for which classical statistical procedures are available.
Introduction to Longitudinal Data Analysis
63
o	However, all these methods have the disadvantage that (lots of) information is lost
Further, they often do not allow to draw conclusions about the way the endpoint
has been reached:
Introduction to Longitudinal Data Analysis
64
Chapter 5
The Multivariate Regression Model
·	The general multivariate model
·	Model ﬁtting with SAS
·	Model reduction
·	Remarks
Introduction to Longitudinal Data Analysis
65
5.1 The General Multivariate Model
We re-consider the growth data:
·	Introduction to Longitudinal Data Analysis
66
This is a completely balanced data set:
·	4 measurements for all subjects
·	measurements taken at exactly the same time points
Let Yi be the vector of n repeated measurements for the ith subject :
Yi = 
0
Yi1 Yi2 . . . Yin 


The general multivariate model assumes that Yi satisﬁes a regression model
§	Yi = Xiβ + εi with
Xi : matrix of covariates
β : vector of regression parameters
εi : vector of error components, εi



N (0, Σ)
∼
Introduction to Longitudinal Data Analysis
67
o	We then have the following distribution for Yi :
Yi
N (Xiβ, Σ)
∼
The mean structure Xiβ is modelled as in classical linear regression and ANOVA
models
·	Usually, Σ is just a general (n
×
However, special structures for Σ can be assumed (see later).
n) covariance matrix.
Assuming independence across individuals, β and the parameters in Σ can be
estimated by maximizing
·	LM L =
N
Yi=1



(2π)−
n/2
−
Σ
|
|
1
2 exp 
−


1
2
(yi
−
Xiβ)0 Σ−
1 (yi
−
Xiβ)





Introduction to Longitudinal Data Analysis
68
Inference is based on classical maximum likelihood theory:
·	LR tests
·	Asymptotic WALD tests
More details on inference will be discussed later
o	Introduction to Longitudinal Data Analysis
69
5.2 Model Fitting With SAS
5.2.1 Model Parameterization
·	As an example, we ﬁt a model with unstructured mean and unstructured
covariance matrix to the growth data (Model 1).
Let xi be equal to 0 for a boy, and equal to 1 for a girl
One possible parameterization of the model is
o	Yi1 = β0,8(1
Yi2 = β0,10(1
Yi3 = β0,12(1
Yi4 = β0,14(1
−
−
−
−
xi) + β1,8xi + εi1
xi) + β1,10xi + εi2
xi) + β1,12xi + εi3
xi) + β1,14xi + εi4
Introduction to Longitudinal Data Analysis
70
In matrix notation:
·	Y i = Xiβ + εi,
with
and with
Xi =
xi)
(1
(1
−
0
0
0








0
−
0
0
xi)
(1
0
0
−
0
xi)
0
0
0
(1
−
xi 0
0
0 xi 0
0
xi) 0
0
0
0 xi 0
0 xi
0








β = (β0,8, β0,10, β0,12, β0,14, β1,8, β1,10, β1,12, β1,14)0
Introduction to Longitudinal Data Analysis
71
5.2.2 SAS Program
SAS syntax:
·	proc mixed data = growth method = ml;
class idnr sex age;
model measure = age*sex / noint s;
repeated age / type = un subject = idnr;
run;
·	Data structure:
one record per observation:
idnr
age
sex
measure
1.0000 1.0000
1.0000 1.0000
2.0000 2.0000
8.0000 10.000
12.000 14.000
8.0000 10.000
1.0000 1.0000
1.0000 1.0000
1.0000 1.0000
21.000 20.000
21.500 23.000
21.000 21.500
......
......
......
......
26.000 26.000
27.000 27.000
27.000 27.000
12.000 14.000
8.0000 10.000
12.000 14.000
0.0000 0.0000
0.0000 0.0000
0.0000 0.0000
26.000 30.000
22.000 21.500
23.500 25.000
Introduction to Longitudinal Data Analysis
72
·	The mean is modeled in the MODEL statement, as in other SAS procedures for
linear models
·	The covariance matrix is modeled in the REPEATED statement:
·	option ‘type=’ speciﬁes covariance structure
·	option ‘subject=idnr’ speciﬁes the clusters in the data set
·	the variable ‘age’ is used to order measurements within clusters
Introduction to Longitudinal Data Analysis
73
5.2.3 Results
Maximized log-likelihood value: ` =
208.25 −
Estimates for parameters in mean structure, and implied ﬁtted averages:
o	Parameter
β0,8
β0,10
β0,12
β0,14
β1,8
β1,10
β1,12
β1,14
MLE
(s.e.)
22.8750 (0.5598)
23.8125 (0.4921)
25.7188 (0.6112)
27.4688 (0.5371)
21.1818 (0.6752)
22.2273 (0.5935)
23.0909 (0.7372)
24.0909 (0.6478)
Introduction to Longitudinal Data Analysis
74
Fitted covariance and correlation matrices:
·	Σ =
c





















5.0143 2.5156 3.6206 2.5095
2.5156 3.8748 2.7103 3.0714
3.6206 2.7103 5.9775 3.8248
2.5095 3.0714 3.8248 4.6164





















=
⇒
1.0000 0.5707 0.6613 0.5216
0.5707 1.0000 0.5632 0.7262
0.6613 0.5632 1.0000 0.7281
0.5216 0.7262 0.7281 1.0000










































Introduction to Longitudinal Data Analysis
75
5.3 Model Reduction
o	In many circumstances, one will be interested in reducing the model.
For the growth data for example, one may be interested in ﬁnding out whether the
ﬁtted average proﬁles can be well described by straight lines.
·	Also, the covariance matrix contained 10 parameters, not even of interest. If this
can be reduced, one may gain eﬃciency for the mean structure.
·	In practice, one therefore usually tries to reduce the mean and covariance
structures, yielding more parsimonious models
This is now illustrated using the growth data
·	Introduction to Longitudinal Data Analysis
76
5.3.1 Reduction of the Mean Structure
Model 2: Linear Average Trends
Linear average trend within each group, unstructured 4
4 covariance matrix Σ
×
Model 2 is given by (xi = 1 for girls):
Yij = β0 + β01xi + β10tj(1
xi) + β11tjxi + εij,
−
In matrix notation, this equals Y i = Xiβ + εi, with design matrix
§	Xi =
1 xi 8(1
xi)
8xi
−
1 xi 10(1
1 xi 12(1
1 xi 14(1
xi) 10xi
xi) 12xi
xi) 14xi
−
−
−


















·	

















Introduction to Longitudinal Data Analysis
77
Parameterization β = (β0, β01, β10, β11)0 :
·	β0 : intercept for boys
·	β0 + β01 : intercept for girls
·	β10 : slope for boys
·	β11 : slope for girls
SAS program :
proc mixed data = growth method = ml;
class idnr sex ageclss;
model measure = sex age*sex /
repeated ageclss / type = un subject = idnr;
run;
s;
The variable ageclss is a copy of the original variable age
§	Introduction to Longitudinal Data Analysis
78
LR test Model 2 versus Model 1:
Mean
Covar par
2`
Ref G2
df
p
1 unstr.
unstr.
−
18 416.509
2
= slopes unstr.
14 419.477
1 2.968 4
0.5632 Predicted trends:
girls : ˆYj = 17.43 + 0.4764tj
boys : ˆYj = 15.84 + 0.8268tj
o	Introduction to Longitudinal Data Analysis
79
6
Model 3: Parallel Average Proﬁles
Linear average trend within each sex group, the same slope for both groups
Unstructured 4
×
4 covariance matrix Σ
Model 3 is given by:
Yij = β0 + β01xi + β1tj + εij.
In matrix notation, this equals Y i = Xiβ + εi, with design matrix
·	Xi =
1 xi 8
1 xi 10
1 xi 12
1 xi 14




































Introduction to Longitudinal Data Analysis
80
·	Parameterization β = (β0, β01, β1)0 :
·	β0 : intercept for boys
·	β0 + β01 : intercept for girls
·	β1 : common slope for boys and girls
·	SAS program :
proc mixed data = growth method = ml;
class idnr sex ageclss;
model measure = sex age / s;
repeated ageclss / type = un subject = idnr;
run;
LR test:
·	Mean
Covar par
2`
Ref G2
df
p
1 unstr.
unstr.
−
18 416.509
2
= slopes unstr.
14 419.477
1 2.968 4
0.5632 3 = slopes unstr.
13 426.153
2 6.676 1
0.0098 Introduction to Longitudinal Data Analysis
81
6
Predicted trends: girls : ˆYj = 15.37 + 0.6747tj
boys : ˆYj = 17.42 + 0.6747tj
·	Introduction to Longitudinal Data Analysis
82
5.3.2 Reduction of the Covariance Structure
·	In order to reduce the number of parameters in the covariance structure, we can
now ﬁt models with more parsimonious structures
This often leads to more eﬃcient inferences for the mean parameters.
This is particularly useful when many repeated measurements are taken per
subject.
SAS includes a large variety of covariance structures (see SAS help function)
§	Introduction to Longitudinal Data Analysis
83
Some examples:
·	Structure
Example
Structure
Example
Unstructured
type=UN
Simple
type=SIMPLE
Compound
symmetry
type=CS
Banded
type=UN(2)
First-order
autoregressive
type=AR(1)




σ2
1 σ12 σ13
σ12 σ2
2 σ23
σ13 σ23 σ2
3
σ2
0
0
0 σ2 0
0 σ2
0




1 + σ2
σ2
σ2
1
σ2
1
σ2
1
1 + σ2
σ2
σ2
1
σ2
1
σ2
1
1 + σ2
σ2




σ2
1 σ12
0
σ12 σ2
2 σ23
0 σ23 σ2
3




ρσ2 ρ2σ2
σ2
ρσ2
σ2
ρσ2
σ2
ρ2σ2 ρσ2
























Toeplitz
type=TOEP
Toeplitz (1)
type=Toep(1)
Heterogeneous
compound
symmetry
type=CSH
Heterogeneous
ﬁrst-order
autoregressive
type=ARH(1)
Heterogeneous
Toeplitz
type=TOEPH




















Introduction to Longitudinal Data Analysis




σ2 σ12 σ13
σ12 σ2 σ12
σ13 σ12 σ2
σ2 0
0
0 σ2 0
0 σ2
0




σ2
1
ρσ1σ2
ρσ1σ3 ρσ2σ3
ρσ1σ2 ρσ1σ3
ρσ2σ3
σ2
3
σ2
2




σ2
1
ρσ1σ2
ρ2σ1σ3 ρσ2σ3
ρσ1σ2 ρ2σ1σ3
ρσ2σ3
σ2
3
σ2
2
σ2
1
ρ1σ1σ2
ρ2σ1σ3 ρ1σ2σ3
ρ1σ1σ2 ρ2σ1σ3
ρ1σ2σ3
σ2
3
σ2
2








84
Model 4: Toeplitz Covariance Structure
o	Linear average trend within each sex group
The estimated covariance matrix (s.e.) of the unstructured covariance matrix
under Model 2 equals:
5.12(1.42) 2.44(0.98) 3.61(1.28) 2.52(1.06)
2.44(0.98) 3.93(1.08) 2.72(1.07) 3.06(1.01)
3.61(1.28) 2.72(1.07) 5.98(1.63) 3.82(1.25)
2.52(1.06) 3.06(1.01) 3.82(1.25) 4.62(1.26)




































·	This suggests that a possible model reduction could consist of assuming equal
variances, and banded covariances.
Introduction to Longitudinal Data Analysis
85
·	This is the so-called Toeplitz covariance matrix Σ, with elements of the form
Σij = α
:
i
j
−
|
|
Σ =
α0 α1 α2 α3
α1 α0 α1 α2
α2 α1 α0 α1
α3 α2 α1 α0




































·	Note that this is only really meaningful when the time points at which
measurements are taken are equally spaced, as in the current example.
·	SAS program :
proc mixed data = growth method = ml;
class sex idnr ageclss;
model measure = sex age*sex / s;
repeated ageclss / type = toep subject = idnr;
run;
Introduction to Longitudinal Data Analysis
86
o	LR test Model 4 versus Model 2:
Mean
Covar
1 unstr.
unstr.
par
2`
Ref G2
df
p
−
18 416.509
2
4
= slopes unstr.
14 419.477
1 2.968 4
0.5632 = slopes banded
8 424.643
2 5.166 6
0.5227 Fitted covariance and correlation matrices:
Σ =
c





















4.9439 3.0507 3.4054 2.3421
3.0507 4.9439 3.0507 3.4054
3.4054 3.0507 4.9439 3.0507
2.3421 3.4054 3.0507 4.9439





















=
⇒
1.0000 0.6171 0.6888 0.4737
0.6171 1.0000 0.6171 0.6888
0.6888 0.6171 1.0000 0.6171
0.4737 0.6888 0.6171 1.0000





















Introduction to Longitudinal Data Analysis





















87
6
6
Model 5: AR(1) Covariance Structure
o	Linear average trend within each sex group
The AR(1) covariance structure assumes exponentially decaying correltions, i.e.,
j
elements of Σ of the form Σij = σ2ρ|
| :
−
i
Σ = σ2
1 ρ ρ2 ρ3
ρ 1 ρ ρ2
ρ2 ρ 1 ρ
ρ3 ρ2 ρ 1




































·	Note that this is also only really meaningful when the time points at which
measurements are taken are equally spaced.
Introduction to Longitudinal Data Analysis
88
SAS program:
proc mixed data = growth method = ml;
class sex idnr ageclss;
model measure = sex age*sex / s;
repeated ageclss / type = AR(1) subject = idnr;
run;
LR test Model 5 versus Models 2 and 4 :
o	Mean
Covar
1 unstr.
unstr.
par
2`
Ref G2
df
p
−
18 416.509
2
4
5
= slopes unstr.
14 419.477
= slopes banded
8 424.643
1
2
2.968 4
0.5632 5.166 6
0.5227 = slopes AR(1)
6 440.681
2 21.204 8
0.0066 4 16.038 2
0.0003 Introduction to Longitudinal Data Analysis
89
6
6
6
Fitted covariance and correlation matrices:
·	Σ =
c





















4.8903 2.9687 1.8021 1.0940
2.9687 4.8903 2.9687 1.8021
1.8021 2.9687 4.8903 2.9687
1.0940 1.8021 2.9687 4.8903





















=
⇒
1.0000 0.6070 0.3685 0.2237
0.6070 1.0000 0.6070 0.3685
0.3685 0.6070 1.0000 0.6070
0.2237 0.3685 0.6070 1.0000










































Introduction to Longitudinal Data Analysis
90
5.4 Remarks
·	The multivariate regression model is primarily suitable when measurements are
taken at a relatively small number of ﬁxed time points
o	Even if some measurements are missing, the multivariate regression model can be
applied, as long as the software allows for unequal numbers of measurements per
subject.
In the SAS procedure MIXED, this is taken care of in the REPEATED statement
repeated ageclss /
;
from which it can be derived which outcomes have been observed, and which ones
are missing.
Introduction to Longitudinal Data Analysis
91
o	In case of large numbers of repeated measurements:
·	Multivariate regression models can only be applied under very speciﬁc mean
and covariance structures, even in case of complete balance.
·	For example, unstructured means and/or unstructured covariances require
estimation of very many parameters
In case of highly unbalanced data:
·	Multivariate regression models can only be applied under very speciﬁc mean
and covariance structures.
·	For example, Toeplitz and AR(1) covariances are not meaningful since time
points are not equally spaced.
·	For example, compound symmetric covariances are meaningful, but based on
very strong assumptions.
Introduction to Longitudinal Data Analysis
92
Chapter 6
A Model for Longitudinal Data
·	Introduction
·	The 2-stage model formulation
·	Examples: Rat and prostate data
·	The general linear mixed-eﬀects model
·	Hierarchical versus marginal model
·	Examples: Rat and prostate data
·	A model for the residual covariance structure
Introduction to Longitudinal Data Analysis
93
6.1 Introduction
·	In practice: often unbalanced data:
·	unequal number of measurements per subject
·	measurements not taken at ﬁxed time points
Therefore, multivariate regression techniques are often not applicable
Often, subject-speciﬁc longitudinal proﬁles can be well approximated by linear
regression functions
This leads to a 2-stage model formulation:
·	Stage 1: Linear regression model for each subject separately
·	Stage 2: Explain variability in the subject-speciﬁc regression coeﬃcients using
known covariates
Introduction to Longitudinal Data Analysis
94
6.2 A 2-stage Model Formulation
6.2.1 Stage 1
§	Response Yij for ith subject, measured at time tij, i = 1, . . . , N , j = 1, . . . , ni
Response vector Yi for ith subject:
Yi = (Yi1, Yi2, . . . , Yini)0
Stage 1 model:
Yi = Ziβi + εi
Introduction to Longitudinal Data Analysis
95
Zi is a (ni ×
·	q) matrix of known covariates
βi is a q-dimensional vector of subject-speciﬁc regression coeﬃcients
εi
∼
N (0, Σi), often Σi = σ2Ini
Note that the above model describes the observed variability within subjects
§	Introduction to Longitudinal Data Analysis
96
6.2.2 Stage 2
·	Between-subject variability can now be studied from relating the βi to known
covariates
Stage 2 model:
βi = Kiβ + bi
Ki is a (q
×
p) matrix of known covariates
β is a p-dimensional vector of unknown regression parameters
N (0, D)
bi
∼
·	Introduction to Longitudinal Data Analysis
97
6.3 Example: The Rat Data
Individual proﬁles:
·	Introduction to Longitudinal Data Analysis
98
Transformation of the time scale to linearize the proﬁles:
·	Ageij −→
tij = ln[1 + (Ageij −
45)/10)]
·	Note that t = 0 corresponds to the start of the treatment (moment of
randomization)
Stage 1 model:
Yij = β1i + β2itij + εij,
j = 1, . . . , ni
Matrix notation:
o	Yi = Ziβi + εi
with
Zi =
1 ti1
1 ti2
...
...
1 tini










































Introduction to Longitudinal Data Analysis
99
·	In the second stage, the subject-speciﬁc intercepts and time eﬀects are related to
the treatment of the rats
Stage 2 model:
β1i = β0 + b1i,
β2i = β1Li + β2Hi + β3Ci + b2i,



Li, Hi, and Ci are indicator variables:
o	Li =



1 if low dose
0 otherwise
Hi =
1 if high dose
0 otherwise
Ci =



1 if control
0 otherwise



Introduction to Longitudinal Data Analysis
100
·	Parameter interpretation:
·	β0: average response at the start of the treatment (independent of treatment)
·	β1, β2, and β3: average time eﬀect for each treatment group
Introduction to Longitudinal Data Analysis
101
6.4 Example: The Prostate Data
Individual proﬁles:
·	Introduction to Longitudinal Data Analysis
102
Transformation of the response:
PSAij −→
Yij = ln(PSAij + 1)
Stage 1 model:
Yij = β1i + β2itij + β3it2
ij + εij,
j = 1, . . . , ni
Matrix notation:
§	Yi = Ziβi + εi
with
Zi =





















1 ti1
t2
i1
t2
i2
1 ti2
...
...
1 tini t2
ini





















·	In the second stage, the subject-speciﬁc intercepts and time eﬀects are related to
the age (at diagnosis) and disease status
Introduction to Longitudinal Data Analysis
103
Stage 2 model:
β1i = β1Agei + β2Ci + β3Bi + β4Li + β5Mi + b1i,
β2i = β6Agei + β7Ci + β8Bi + β9Li + β10Mi + b2i,
β3i = β11Agei + β12Ci + β13Bi + β14Li + β15Mi + b3i



Ci, Bi, Li and Mi are indicator variables:
o	1 if Control
0 otherwise
1 if L/R cancer case
0 otherwise
Ci =
Li =






1 if BPH case
0 otherwise
1 if Metastatic cancer case
0 otherwise
Bi =
Mi =






Introduction to Longitudinal Data Analysis
104
·	Parameter interpretation:
·	β2, β3, β4, and β5: average intercepts after correction for age
·	β7, β8, β9, and β10: average linear time eﬀects after correction for age.
·	β12, β13, β14, and β15: average quadratic time eﬀects after correction for age.
Introduction to Longitudinal Data Analysis
105
6.5 The General Linear Mixed-eﬀects Model
o	A 2-stage approach can be performed explicitly in the analysis
However, this is just another example of the use of summary statistics:
·	Yi is summarized by
βi
d
·	summary statistics
βi analysed in second stage
d
·	The associated drawbacks can be avoided by combining the two stages into one
model:
Yi = Ziβi + εi
βi = Kiβ + bi
=
⇒



Yi = ZiKi
β + Zibi + εi = Xiβ + Zibi + εi
Xi
{z
|
}
Introduction to Longitudinal Data Analysis
106
o	General linear mixed-eﬀects model:
Yi = Xiβ + Zibi + εi
N (0, D),
bi
∼
N (0, Σi),
εi
∼
b1, . . . , bN, ε1, . . . , εN independent



Terminology:
·	Fixed eﬀects: β
·	Random eﬀects: bi
·	Variance components: elements in D and Σi
Introduction to Longitudinal Data Analysis
107
6.6 Hierarchical versus Marginal Model
The general linear mixed model is given by:
Yi = Xiβ + Zibi + εi
N (0, D),
bi
∼
N (0, Σi),
εi
∼
b1, . . . , bN, ε1, . . . , εN independent



It can be rewritten as:
o	Yi
bi
|
∼
N (Xiβ + Zibi, Σi),
bi
∼
N (0, D)
Introduction to Longitudinal Data Analysis
108
It is therefore also called a hierarchical model:
·	A model for Yi given bi
·	A model for bi
Marginally, we have that Yi is distributed as:
o	Yi
∼
N (Xiβ, ZiDZ 0i + Σi)
·	Hence, very speciﬁc assumptions are made about the dependence of mean and
covariance onn the covariates Xi and Zi:
·	Implied mean : Xiβ
·	Implied covariance : Vi = ZiDZ 0i + Σi
Note that the hierarchical model implies the marginal one, NOT vice versa
·	Introduction to Longitudinal Data Analysis
109
6.7 Example: The Rat Data
Stage 1 model:
Yij = β1i + β2itij + εij,
j = 1, . . . , ni
Stage 2 model:
β1i = β0 + b1i,
β2i = β1Li + β2Hi + β3Ci + b2i,



Combined:
Yij = (β0 + b1i) + (β1Li + β2Hi + β3Ci + b2i)tij + εij
§	β0 + b1i + (β1 + b2i)tij + εij,
if low dose
β0 + b1i + (β2 + b2i)tij + εij,
if high dose
β0 + b1i + (β3 + b2i)tij + εij,
if control.
=



Introduction to Longitudinal Data Analysis
110
Implied marginal mean structure:
·	Linear average evolution in each group
·	Equal average intercepts
·	Diﬀerent average slopes
Implied marginal covariance structure (Σi = σ2Ini):
o	Cov(Yi(t1), Yi(t2)) = 
1 t1 
D


1
t2
















·	σ2δ
{
t1,t2}
= d22t1 t2 + d12(t1 + t2) + d11 + σ2δ
{
·	t1,t2}
·	Note that the model implicitly assumes that the variance function is quadratic
over time, with positive curvature d22.
Introduction to Longitudinal Data Analysis
111
·	A model which assumes that all variability in subject-speciﬁc slopes can be
ascribed to treatment diﬀerences can be obtained by omitting the random slopes
b2i from the above model:
Yij = (β0 + b1i) + (β1Li + β2Hi + β3Ci)tij + εij
β0 + b1i + β1tij + εij,
if low dose
β0 + b1i + β2tij + εij,
if high dose
β0 + b1i + β3tij + εij,
if control.
=



This is the so-called random-intercepts model
o	The same marginal mean structure is obtained as under the model with random
slopes
Introduction to Longitudinal Data Analysis
112
Implied marginal covariance structure (Σi = σ2Ini):
·	Cov(Yi(t1), Yi(t2)) = 
1 
D 
1 




·	σ2δ
{
t1,t2}
= d11 + σ2δ
{
·	t1,t2}
·	Hence, the implied covariance matrix is compound symmetry:
·	constant variance d11 + σ2
·	constant correlation ρI = d11/(d11 + σ2) between any two repeated
measurements within the same rat
Introduction to Longitudinal Data Analysis
113
6.8 Example: The Prostate Data
§	Stage 1 model:
Yij = β1i + β2itij + β3it2
ij + εij,
j = 1, . . . , ni
Stage 2 model:
β1i = β1Agei + β2Ci + β3Bi + β4Li + β5Mi + b1i,
β2i = β6Agei + β7Ci + β8Bi + β9Li + β10Mi + b2i,
β3i = β11Agei + β12Ci + β13Bi + β14Li + β15Mi + b3i,



Combined:
Yij = β1Agei + β2Ci + β3Bi + β4Li + β5Mi
·	(β6Agei + β7Ci + β8Bi + β9Li + β10Mi) tij
·	(β11Agei + β12Ci + β13Bi + β14Li + β15Mi) t2
ij
·	b1i + b2itij + b3it2
ij + εij.
Introduction to Longitudinal Data Analysis
114
Implied marginal mean structure:
·	Quadratic average evolution in each group
·	Average intercept and linear as well as quadratic slopes corrected for age
diﬀerences
Implied marginal covariance structure (Σi = σ2Ini):
Cov(Yi(t1), Yi(t2)) = 

1 t1 t2
1 

D
1
t2
t2
2




























·	σ2δ
{
t1,t2}
= d33t2
1 t2
2 + d23(t2
1 t2 + t1 t2
o	d22t1 t2
+d13(t2
1 + t2
o	d12(t1 + t2) + d11 + σ2δ
{
·	t1,t2}
The implied variance function is now a four-degree polynomial over time.
§	Introduction to Longitudinal Data Analysis
115
6.9 Example: Bivariate Observations
Balanced data, two measurements per subject (ni = 2), two models:
·	Model 1:
Random intercepts
heterogeneous errors
Model 2:
Uncorrelated intercepts and slopes
+
measurement error
(d) (1 1) + 
1

1






d + σ2
1
V = 






= 






d
d + σ2
2
d













σ2
1 0
0 σ2
2







V = 






= 






1 0
1 1














d1 + σ2
d1 0
0 d2
1 1
0 1














d1
d1
d1 + d2 + σ2














σ2 0
0 σ2







·	






Introduction to Longitudinal Data Analysis
116
o	Diﬀerent hierarchical models can produce the same marginal model
Hence, a good ﬁt of the marginal model cannot be interpreted as evidence for any
of the hierarchical models.
·	A satisfactory treatment of the hierarchical model is only possible within a
Bayesian context.
Introduction to Longitudinal Data Analysis
117
6.10 A Model for the Residual Covariance Structure
o	Often, Σi is taken equal to σ2Ini
We then obtain conditional independence:
Conditional on bi, the elements in Yi are independent
·	In the presence of no, or little, random eﬀects, conditional independence is often
unrealistic
·	For example, the random intercepts model not only implies constant variance, it
also implicitly assumes constant correlation between any two measurements within
subjects.
Introduction to Longitudinal Data Analysis
118
§	Hence, when there is no evidence for (additional) random eﬀects, or if they would
have no substantive meaning, the correlation structure in the data can be
accounted for in an appropriate model for Σi
Frequently used model:
Yi = Xiβ + Zibi + ε(1)i + ε(2)i
|
{z
↓εi
}
3 stochastic components:
·	bi: between-subject variability
·	ε(1)i: measurement error
·	ε(2)i: serial correlation component
Introduction to Longitudinal Data Analysis
119
ε(2)i represents the belief that part of an individual’s observed proﬁle is a response
to time-varying stochastic processes operating within that individual.
o	This results in a correlation between serial measurements, which is usually a
decreasing function of the time separation between these measurements.
o	The correlation matrix Hi of ε(2)i is assumed to have (j, k) element of the form
hijk = g(
|
) for some decreasing function g(
·
) with g(0) = 1
tij −
tik|
Frequently used functions g(
·
):
·	Exponential serial correlation: g(u) = exp(
φu)
·	Gaussian serial correlation: g(u) = exp(
−
−
φu2)
Introduction to Longitudinal Data Analysis
120
Graphically, for φ = 1:
·	Extreme cases:
·	∞
·	φ = +
: components in ε(2)i independent
·	φ = 0: components in ε(2)i perfectly correlated
Introduction to Longitudinal Data Analysis
121
In general, the smaller φ, the stronger is the serial correlation.
Resulting ﬁnal linear mixed model:
o	Yi = Xiβ + Zibi + ε(1)i + ε(2)i
N (0, D)
N (0, σ2Ini)
N (0, τ 2Hi)
bi
∼
ε(1)i
∼
ε(2)i
∼
independent



Introduction to Longitudinal Data Analysis
122
Graphical representation of all 4 components in the model:
·	Introduction to Longitudinal Data Analysis
123
Chapter 7
Exploratory Data Analysis
·	Introduction
·	Mean structure
·	Variance function
·	Correlation structure
·	Individual proﬁles
Introduction to Longitudinal Data Analysis
124
7.1 Introduction
·	A linear mixed model makes assumptions about:
·	mean structure: (non-)linear, covariates,. . .
·	variance function: constant, quadratic, . . .
·	correlation structure: constant, serial, . . .
·	subject-speciﬁc proﬁles: linear, quadratic, . . .
·	In practice, linear mixed models are often obtained from a two-stage model
formulation
However, this may or may not imply a valid marginal model
·	Introduction to Longitudinal Data Analysis
125
As an example, reconsider the growth curves:
·	Introduction to Longitudinal Data Analysis
126
o	The individual proﬁles support a random-intercepts model
However, the estimated covariance matrix suggests non-constant variance
function:
6.11 6.88
8.26 7.44
7.18 6.88 8.53
9.78 9.01
8.70 8.26 9.78 12.04 10.99 10.96
7.44 9.01 10.99 10.42 10.56
7.18 8.70 10.96 10.56 11.24



























·	


























·	Data exploration is therefore extremely helpful as additional tool in the selection
of appropriate models
Introduction to Longitudinal Data Analysis
127
7.2 Exploring the Mean Structure
·	For balanced data, averages can be calculated for each occasion separately, and
standard errors for the means can be added
·	Example: rat data:
·	SAS program:
filename fig1 ’d:\path\file.eps’;
goptions reset=all ftext=swiss device=psepsf gsfname=fig1 gsfmode=replace
rotate=landscape;
proc gplot data=test;
plot y*age / haxis=axis1 vaxis=axis2;
symbol c=red i=std1mjt w=2 mode=include;
axis1 label=(h=2 ’Age (days)’) value=(h=1.5) order=(40 to 120 by 10) minor=none;
axis2 label=(h=2 A=90 ’Response (pixels)’) value=(h=1.5) order=(70 to 85 by 5)
minor=none;
title h=3 ’Average evolution, with standard errors of means’;
run;quit;
Introduction to Longitudinal Data Analysis
128
·	SAS output:
·	Conclusion: non-linear average trend, increasing standard errors due to dropout
Introduction to Longitudinal Data Analysis
129
o	For unbalanced data:
·	Discretize the time scale and use simple averaging within intervals
·	Smoothing techniques to estimate the average evolution nonparametrically
Example: prostate data:
·	SAS program for loess smoothing:
proc loess data=test;
ods output scoreresults=out;
model lnpsa=time;
score data=test;
run;
proc sort data=out;
by time;
run;
filename fig1 ’d:\path\file.eps’;
goptions reset=all ftext=swiss device=psepsf
gsfname=fig1 gsfmode=replace rotate=landscape;
proc gplot data=out;
plot lnpsatime=1 p_lnpsatime=2
/ overlay haxis=axis1 vaxis=axis2;
symbol1 c=red v=dot h=0.2 mode=include;
symbol2 c=black i=join w=2 mode=include;
axis1 label=(h=2 ’Years before diagnosis’)
value=(h=1.5) order=(0 to 30 by 5) minor=none;
axis2 label=(h=2 A=90 ’ln(PSA+1)’) value=(h=1.5)
order=(0 to 4 by 1) minor=none;
title h=3 ’Loess smoothing’;
run;quit;
Introduction to Longitudinal Data Analysis
130
·	SAS output:
Introduction to Longitudinal Data Analysis
131
·	If (important) covariates or factors are known, similar plots can be constructed for
subgroups with diﬀerent values for these covariates or factors.
Example for the rat data:
·	Introduction to Longitudinal Data Analysis
132
Example for the prostate data:
·	Introduction to Longitudinal Data Analysis
133
7.3 Exploring the Variance Function
The variance function equals
·	σ2(t) = E[Y (t)
µ(t)]2
−
·	Hence, an estimate for σ2(t) can be obtained from applying any of the techniques
described for exploring the mean structure to squared residuals r2
ij
Introduction to Longitudinal Data Analysis
134
Example for the rat data (averages with standard deviations):
·	Introduction to Longitudinal Data Analysis
135
Example for the prostate data (based on group-speciﬁc smoothing of averages):
·	Introduction to Longitudinal Data Analysis
136
7.4 Exploring the Correlation Structure
7.4.1 Scatterplot and Correlation Matrix
·	For balanced longitudinal data, the correlation structure can be studied through
the correlation matrix, or a scatterplot matrix
Correlation matrix for the growth data:
·	1.00 0.63 0.71 0.60
0.63 1.00 0.63 0.76
0.71 0.63 1.00 0.80
0.60 0.76 0.80 1.00


















·	

















·	Graphically, pairwise scatterplots can be used for exploring the correlation between
any two repeated measurements
Introduction to Longitudinal Data Analysis
137
Scatterplot matrix for the growth data:
·	Introduction to Longitudinal Data Analysis
138
7.4.2 Semi-variogram
·	For unbalanced data, the same approach can be used, after discretizing the time
scale.
·	An alternative method, in case the variance function suggests constant variance is
the semi-variogram
Re-consider the general linear mixed model:
·	Yi = Xiβ + Zibi + ε(1)i + ε(2)i
N (0, D)
N (0, σ2Ini)
N (0, τ 2Hi)
bi
∼
ε(1)i
∼
ε(2)i
∼
independent



Introduction to Longitudinal Data Analysis
139
Based on a mean function exploration, residuals rij = yij −
·	µ(tij) can be obtained
These residuals are assumed to follow the model:
ri = Zibi + ε(1)i + ε(2)i
The semi-variogram assumes constant variance, which implies that the only
random eﬀects in the model will at most be intercepts, i.e., Zi = 
1 1
1 

· · ·
We will denote the variance of the random intercepts by ν2

The covariance matrix is then of the form
Vi = Var(Yi) = Var(ri) = ν2ZiZ 0i + σ2Ini + τ 2Hi
The residuals rij have constant variance ν2 + σ2 + τ 2
o	Introduction to Longitudinal Data Analysis
140
·	The correlation between any two residuals rij and rik from the same subject i is
given by
ρ(
|
tij −
tik|
) =
ν2 + τ 2 g(
tij −
|
ν2 + σ2 + τ 2
tik|
)
§	One can show that, for j
= k,
1
2
E (rij −
rik)2 = σ2 + τ 2 (1
= v(uijk)
g(
|
tij −
tik|
))
−
The function v(u) is called the semi-variogram, and it only depends on the time
points tij through the time lags uijk =
·	tij −
|
tik|
·	Decreasing serial correlation functions g(
·
with v(0) = σ2, which converge to σ2 + τ 2 as u grows to inﬁnity.
) yield increasing semi-variograms v(u),
Introduction to Longitudinal Data Analysis
141
6
·	Semi-variograms for exponential and Gaussian serial correlation functions g(
·
σ2 = 0.7, τ 2 = 1.3, and ν2 = 1, φ = 1:
),
Introduction to Longitudinal Data Analysis
142
·	Obviously, an estimate of v(u) can be used to explore the relative importance of
the stochastic components bi, ε(1)i, and ε(2)i, as well as the nature of the serial
correlation function g(
·
).
An estimate of v(u) is obtained from smoothing the scatter plot of the
N
ni(ni −
1)/2 half-squared diﬀerences vijk = (rij −
Xi=1
residuals within subjects versus the corresponding time lags uijk =
rik)2/2 between pairs of
·	tik|
tij −
|
One can also show that, for i
= k:
1
2E[rij −
rkl]2 = σ2 + τ 2 + ν2
Hence, the total variability in the data (assumed constant) can be estimated by
σ2 +
τ 2 +
ν2 =
c
c
c
1
ni
nl
2N ∗ Xi
=k
Xj=1
Xl=1
(rij −
rkl)2,
where N ∗ is the number of terms in the sum.
Introduction to Longitudinal Data Analysis
143
6
6
·	Example: prostate data
·	We now consider the control group only:
·	Assuming constant variability, the variogram can be constructed to explore the
3 stochastic components.
Introduction to Longitudinal Data Analysis
144
·	SAS program for loess smoothing:
/* Calculation of residuals, linear average trend */
proc glm data=prostate;
model lnpsa=time;
output out=out r=residual;
run;
/* Calculation of the variogram */
proc variogram data=out outpair=out;
coordinates xc=time yc=id;
compute robust novariogram;
var residual;
run;
data variogram;set out;
if y1=y2;vario=(v1-v2)**2/2; run;
data variance;set out;
if y1<y2; vario=(v1-v2)**2/2; run;
/* Calculation of the total variance (=0.148) */
proc means data=variance mean;
var vario;
run;
Introduction to Longitudinal Data Analysis
145
/* Loess smoothing of the variogram */
proc loess data=variogram;
ods output scoreresults=out;
model vario=distance;
score data=variogram;
run;
proc sort data=out;by distance;run;
filename fig1 ’d:\path\file.eps’;
goptions reset=all ftext=swiss device=psepsf gsfname=fig1
gsfmode=replace rotate=landscape;
proc gplot data=out;
plot variodistance=1 p_variodistance=2
/ overlay haxis=axis1 vaxis=axis2 vref=0.148 lvref=3;
symbol1 c=red v=dot h=0.2 mode=include;
symbol2 c=black i=join w=2 mode=include;
axis1 label=(h=2 ’Time lag’) value=(h=1.5)
order=(0 to 20 by 5) minor=none;
axis2 label=(h=2 A=90 ’v(u)’) value=(h=1.5)
order=(0 to 0.4 by 0.1) minor=none;
title h=3 ’Semi-variogram’;
run;quit;
Introduction to Longitudinal Data Analysis
146
·	SAS output:
·	The total variability is estimated to be 0.148
·	Random intercepts represent most of the variability, while there is very little
evidence for the presence of serial correlation.
Introduction to Longitudinal Data Analysis
147
7.5 Exploring the Individual Proﬁles
7.5.1 Introduction
·	As discussed before, linear mixed models are often obtained from a two-stage
model formulation
·	This is based on a good approximation of the subject-speciﬁc proﬁles by linear
regression models
This requires methods for the exploration of longitudinal proﬁles
·	Introduction to Longitudinal Data Analysis
148
7.5.2 Graphical Exploration
o	An natural way to explore longitudinal proﬁles is by plotting them
Example: Prostate data:
·	SAS program:
proc sort data=prostate;
by id time;
run;
filename fig1 ’d:\path\file.eps’;
goptions reset=all ftext=swiss device=psepsf gsfname=fig1
gsfmode=replace rotate=landscape i=join;
proc gplot data=test;
plot lnpsa*time=id / haxis=axis1 vaxis=axis2 nolegend;
axis1 label=(h=2 ’Years before diagnosis’) value=(h=1.5)
order=(0 to 30 by 5) minor=none;
axis2 label=(h=2 A=90 ’ln(PSA+1)’) value=(h=1.5)
order=(0 to 4 by 1) minor=none;
title h=3 ’Individual profiles’;
run;quit;
Introduction to Longitudinal Data Analysis
149
·	SAS output:
·	In case of large data sets:
·	Randomly select some proﬁles
·	Order subjects according to a speciﬁc proﬁle characteristic (mean,
variability,. . . ) and plot proﬁles for some proﬁles
Introduction to Longitudinal Data Analysis
150
7.5.3 Exploring Subject-speciﬁc Regression Model
Some ad hoc statistical procedures for checking the linear regression models
·	used in the ﬁrst stage of the model formulation.
Yi = Ziβi + εi
·	Extensions of classical linear regression techniques:
·	Coeﬃcient R2 of multiple determination
·	Formal test for the need of a model extension
Introduction to Longitudinal Data Analysis
151
Coeﬃcients of Multiple Determination
In linear regression:
R2 = SSTO
SSE
−
SSTO
Subject-speciﬁc coeﬃcients:
i = SSTOi −
R2
SSTOi
SSEi
Histogram of R2
i or scatterplot of R2
i versus ni
Overall R2:
N
R2
meta =
SAS macro available
Xi=1
(SSTOi −
SSTOi
N
Xi=1
SSEi)
,
o	Introduction to Longitudinal Data Analysis
152
Test for Model Extension
·	Test for the need to extend the linear regression model Y = Xβ + ε with
additional covariates in X ∗:
F =
(SSE(R)
SSE(F )/(N
−
SSE(F ))/p∗
p∗)
p
−
−
Overall test for the need to extend the stage 1 model:
·	Fmeta =
X
i:ni≥
p+p∗}
(SSEi(R)
SSEi(F ))

−
,


SSEi(F )

,




{

X
i:ni≥
p+p∗}
p∗

{

X
i:ni≥
p+p∗}
p
(ni −
−

p∗)


X
i:ni≥
p+p∗}
{


{




·	Null-distribution is F with
freedom
i:ni≥
p+p∗}
P{
p∗ and
i:ni≥
p+p∗}
P{
(ni −
p
−
p∗) degrees of
SAS macro available
·	Introduction to Longitudinal Data Analysis
153
Example: Prostate Data
Scatterplots of R2
i under linear and quadratic model:
·	Introduction to Longitudinal Data Analysis
154
·	Linear model:
·	R2
meta = 0.8188
·	F -test linear vs. quadratic: F54,301 = 6.2181 (p < 0.0001)
·	Quadratic model:
·	R2
meta = 0.9143
·	F -test quadratic vs. cubic: F54,247 = 1.2310 (p = 0.1484)
Introduction to Longitudinal Data Analysis
155
Chapter 8
Estimation of the Marginal Model
·	Introduction
·	Maximum likelihood estimation
·	Restricted maximum likelihood estimation
·	Fitting linear mixed models in SAS
·	Negative variance components
Introduction to Longitudinal Data Analysis
156
8.1 Introduction
§	Recall that the general linear mixed model equals
Yi = Xiβ + Zibi + εi
N (0, D)
N (0, Σi)
bi
εi
∼
∼



independent
The implied marginal model equals Yi
N (Xiβ, ZiDZ 0i + Σi)
∼
Note that inferences based on the marginal model do not explicitly assume the
presence of random eﬀects representing the natural heterogeneity between subjects
Introduction to Longitudinal Data Analysis
157
§	Notation:
·	β: vector of ﬁxed eﬀects (as before)
·	α: vector of all variance components in D and Σi
·	θ = (β0, α0)0: vector of all parameters in marginal model
Marginal likelihood function:
LML(θ) =
(2π)−
ni/2
Vi(α)
|
|
1
2 exp 
−


1
2
−
(Yi
−
N
Yi=1



1
Xiβ)0 V −
i
(α) (Yi
−
Xiβ)





If α were known, MLE of β equals
β(α) = 
N
−
X 0iWiXi
X 0iWiyi,
Xi=1


1 N
Xi=1


c
1
where Wi equals V −
i
·	Introduction to Longitudinal Data Analysis
158
o	In most cases, α is not known, and needs to be replaced by an estimate
α
d
Two frequently used estimation methods for α:
·	Maximum likelihood
·	Restricted maximum likelihood
Introduction to Longitudinal Data Analysis
159
8.2 Maximum Likelihood Estimation (ML)
αML obtained from maximizing
·	d
with respect to α
LML(α,
β(α))
c
The resulting estimate
β(
αML) for β will be denoted by
c
d
βML
c
·	αML and
with respect to α and β simultaneously.
c
·	d
βML can also be obtained from maximizing LML(θ) with respect to θ, i.e.,
Introduction to Longitudinal Data Analysis
160
8.3 Restricted Maximum Likelihood Estimation (REML)
8.3.1 Variance Estimation in Normal Populations
Consider a sample of N observations Y1, . . . , YN from N (µ, σ2)
For known µ, MLE of σ2 equals:
σ2 =
c
(Yi −
Xi
µ)2/N
o	σ2 is unbiased for σ2
·	c
When µ is not known, MLE of σ2 equals:
σ2 =
c
(Yi −
Xi
Y )2/N
Note that
σ2 is biased for σ2:
c
E
σ2
(cid:18)
(cid:19)
c
=
N
−
N
1
σ2
o	Introduction to Longitudinal Data Analysis
161
The bias expression tells us how to derive an unbiased estimate:
S2 =
(Yi −
Xi
Y )2/(N
−
Apparently, having to estimate µ introduces bias in MLE of σ2
How to estimate σ2, without estimating µ ﬁrst ?
The model for all data simultaneously:
·	Y =
Y1
...
YN




























N
∼
µ
...
µ










































, σ2IN














Introduction to Longitudinal Data Analysis
162
We transform Y such that µ vanishes from the likelihood:
U =
Y2
Y3
Y1 −
Y2 −
...
YN
−
2 −
YN
1
−
YN
−
1 −
YN






















































= A0Y
∼
N (0, σ2A0A)
MLE of σ2, based on U , equals:
S2 =
1
−
N
(Yi −
1 Xi
Y )2
A deﬁnes a set of N
−
1 linearly independent ‘error contrasts’
S2 is called the REML estimate of σ2, and S2 is independent of A
·	Introduction to Longitudinal Data Analysis
163
8.3.2 Estimation of Residual Variance in Linear Regression Model
§	Consider a sample of N observations Y1, . . . , YN from a linear regression model:
Y =
Y1
...
YN




























N (Xβ, σ2I)
∼
σ2 = (Y
c
X
β)0(Y
c
−
−
X
β)/N,
c
MLE of σ2:
Note that
σ2 is biased for σ2:
c
E
σ2
(cid:18)
(cid:19)
c
=
N
−
N
p
σ2
Introduction to Longitudinal Data Analysis
164
The bias expression tells us how to derive an unbiased estimate:
MSE = (Y
X
β)0(Y
c
−
−
X
β)/(N
c
p),
−
The MSE can also be obtained from transforming the data orthogonal to X:
U = A0Y
∼
N (0, σ2A0A)
The MLE of σ2, based on U , now equals the mean squared error, MSE
The MSE is again called the REML estimate of σ2
·	Introduction to Longitudinal Data Analysis
165
8.3.3 REML for the Linear Mixed Model
o	We ﬁrst combine all models
into one model
in which
Yi
∼
N (Xiβ, Vi)
Y
∼
N (Xβ, V )
Y =
Y1
...
YN




























, X =
X1
...
XN




























, V (α) =
V1 · · ·
...
·	. .
0
· · ·














0
...
VN














Again, the data are transformed orthogonal to X:
U = A0Y
∼
N (0, A0V (α)A)
Introduction to Longitudinal Data Analysis
166
The MLE of α, based on U is called the REML estimate, and is denoted by
αREML
d
The resulting estimate
β(
αREML) for β will be denoted by
c
d
βREML
c
o	αREML and
·	d
βREML can also be obtained from maximizing
c
LREML(θ) = (cid:12)
X 0iWi(α)Xi(cid:12)
N
Xi=1
(cid:12)
(cid:12)
(cid:12)
(cid:12)
(cid:12)
(cid:12)
1
2
−
LML(θ)
(cid:12)
(cid:12)
(cid:12)
(cid:12)
(cid:12)
(cid:12)
with respect to θ, i.e., with respect to α and β simultaneously.
α,
β(α)
LREML
REML likelihood function.
c
(cid:19)
(cid:18)
is the likelihood of the error contrasts U , and is often called the
Note that LREML(θ) is NOT the likelihood for our original data Y
o	Introduction to Longitudinal Data Analysis
167
8.4 Fitting Linear Mixed Models in SAS
o	Reconsider the model for the prostate data:
ln(PSAij + 1)
= β1Agei + β2Ci + β3Bi + β4Li + β5Mi
·	(β6Agei + β7Ci + β8Bi + β9Li + β10Mi) tij
·	(β11Agei + β12Ci + β13Bi + β14Li + β15Mi) t2
ij
·	b1i + b2itij + b3it2
ij + εij.
Factor group deﬁned by :
·	control : group = 1
·	BPH : group = 2
·	local cancer : group = 3
·	metastatic cancer : group = 4
Introduction to Longitudinal Data Analysis
168
We will assume Σi = σ2Ini
time and timeclss are time, expressed in decades before diagnosis
age is age at the time of diagnosis
lnpsa = ln(P SA + 1)
SAS program:
o	proc mixed data=prostate method=reml;
class id group timeclss;
model lnpsa = group age grouptime agetime grouptime2 agetime2 / noint solution;
random intercept time time2 / type=un subject=id g gcorr v vcorr;
repeated timeclss / type=simple subject=id r rcorr;
run;
Introduction to Longitudinal Data Analysis
169
§	PROC MIXED statement:
·	calls procedure MIXED
·	speciﬁes data-set (records correspond to occasions)
·	estimation method: ML, REML (default), . . .
CLASS statement: deﬁnition of the factors in the model
MODEL statement:
·	response variable
·	ﬁxed eﬀects
·	options similar to SAS regression procedures
Introduction to Longitudinal Data Analysis
170
o	RANDOM statement:
·	deﬁnition of random eﬀects (including intercepts !)
·	identiﬁcation of the ‘subjects’ : independence accross subjects
·	type of random-eﬀects covariance matrix D
·	options ‘g’ and ‘gcorr’ to print out D and corresponding correlation matrix
·	options ‘v’ and ‘vcorr’ to print out Vi and corresponding correlation matrix
REPEATED statement :
·	ordering of measurements within subjects
·	the eﬀect(s) speciﬁed must be of the factor-type
·	identiﬁcation of the ‘subjects’ : independence accross subjects
·	type of residual covariance matrix Σi
·	options ‘r’ and ‘rcorr’ to print out Σi and corresponding correlation matrix
Introduction to Longitudinal Data Analysis
171
·	Some frequently used covariance structures available in RANDOM and
REPEATED statements:
Structure
Example
Structure
Example
Unstructured
type=UN
Simple
type=SIMPLE
Compound
symmetry
type=CS
Banded
type=UN(2)
First-order
autoregressive
type=AR(1)
σ2
1 σ12 σ13
σ12 σ2
2 σ23
σ13 σ23 σ2
3



σ2
0
0
0 σ2 0
0 σ2
0



1 + σ2
σ2
σ2
1
σ2
1
σ2
1
1 + σ2
σ2
σ2
1
σ2
1
σ2
1
1 + σ2
σ2



σ2
1 σ12
0
σ12 σ2
2 σ23
0 σ23 σ2
3



ρσ2 ρ2σ2
σ2
ρσ2
ρσ2
σ2
σ2
ρ2σ2 ρσ2


















Toeplitz
type=TOEP
Toeplitz (1)
type=Toep(1)
Heterogeneous
compound
symmetry
type=CSH
Heterogeneous
ﬁrst-order
autoregressive
type=ARH(1)
Heterogeneous
Toeplitz
type=TOEPH















Introduction to Longitudinal Data Analysis



σ2 σ12 σ13
σ12 σ2 σ12
σ13 σ12 σ2
σ2 0
0
0 σ2 0
0 σ2
0



σ2
1
ρσ1σ2
ρσ1σ3 ρσ2σ3
ρσ1σ2 ρσ1σ3
ρσ2σ3
σ2
3
σ2
2



σ2
1
ρσ1σ2
ρ2σ1σ3 ρσ2σ3
ρσ1σ2 ρ2σ1σ3
ρσ2σ3
σ2
3
σ2
2
σ2
1
ρ1σ1σ2
ρ2σ1σ3 ρ1σ2σ3
ρ1σ1σ2 ρ2σ1σ3
ρ1σ2σ3
σ2
3
σ2
2






172
·	When serial correlation is to be ﬁtted, it should be speciﬁed in the REPEATED
statement, and the option ‘local’ can then be added to also include measurement
error, if required.
·	Some frequently used serial correlation structures available in RANDOM and
REPEATED statements:
Structure
Example
Power
type=SP(POW)(list)
Exponential
type=SP(EXP)(list)
Gaussian
type=SP(GAU)(list)
σ2
σ2
σ2









ρd12 ρd13
1
ρd23
ρd12
1
ρd13 ρd23
1



exp(
1
d12/ρ)
d13/ρ) exp(
d12/ρ) exp(
−
1
exp(
d23/ρ)
d13/ρ)
−
d23/ρ)
−
1

exp(
1
12/ρ2)
d2
13/ρ2) exp(
d2
12/ρ2) exp(
d2
−
exp(
1
23/ρ2)
d2
−
−


d2
13/ρ2)
−
23/ρ2)
d2
−
1



exp(
exp(
exp(
exp(
−
−
−
−
Introduction to Longitudinal Data Analysis
173
ML and REML estimates for ﬁxed eﬀects:
·	Eﬀect
Age eﬀect
Intercepts:
Control
BPH
L/R cancer
Met. cancer
time eﬀect
×
Age
Time eﬀects:
Control
BPH
L/R cancer
Met. cancer
time2 eﬀect
×
Age
Time2 eﬀects:
Control
BPH
L/R cancer
Met. cancer
Parameter
β1
MLE (s.e.)
0.026 (0.013)
REMLE (s.e.)
0.027 (0.014)
β2
β3
β4
β5
β6
β7
β8
β9
β10
β11
β12
β13
β14
β15
1.077 (0.919)
0.493 (1.026)
0.314 (0.997)
1.574 (1.022)
0.010 (0.020)
−
−
−
1.098 (0.976)
0.523 (1.090)
0.296 (1.059)
1.549 (1.086)
0.011 (0.021)
−
−
−
0.511 (1.359)
0.313 (1.511)
1.072 (1.469)
1.657 (1.499)
0.002 (0.008)
−
−
0.568 (1.473)
0.396 (1.638)
1.036 (1.593)
1.605 (1.626)
0.002 (0.009)
−
−
−
−
0.106 (0.549)
0.119 (0.604)
0.350 (0.590)
0.411 (0.598)
−
−
0.130 (0.610)
0.158 (0.672)
0.342 (0.656)
0.395 (0.666)
Introduction to Longitudinal Data Analysis
174
ML and REML estimates for variance components:
·	Eﬀect
Covariance of bi:
var(b1i)
var(b2i)
var(b3i)
cov(b1i, b2i)
cov(b2i, b3i)
cov(b3i, b1i)
Residual variance:
var(εij)
Log-likelihood
Parameter
MLE (s.e.)
REMLE (s.e.)
d11
d22
d33
d12 = d21 −
d23 = d32 −
d13 = d31
0.398 (0.083)
0.768 (0.187)
0.103 (0.032)
0.443 (0.113)
0.273 (0.076)
0.133 (0.043)
0.452 (0.098)
0.915 (0.230)
0.131 (0.041)
0.518 (0.136)
0.336 (0.095)
0.163 (0.053)
−
−
σ2
0.028 (0.002)
1.788 −
0.028 (0.002)
31.235 −
Introduction to Longitudinal Data Analysis
175
Fitted average proﬁles at median age at diagnosis:
·	Introduction to Longitudinal Data Analysis
176
8.5 Negative Variance Components
Reconsider the model for the rat data:
Yij = (β0 + b1i) + (β1Li + β2Hi + β3Ci + b2i)tij + εij
REML estimates obtained from SAS procedure MIXED:
o	Eﬀect
Intercept
Time eﬀects:
Low dose
High dose
Control
Covariance of bi:
var(b1i)
var(b2i)
cov(b1i, b2i)
Residual variance:
var(εij)
REML log-likelihood
Parameter
β0
REMLE (s.e.)
68.606 (0.325)
β1
β2
β3
d11
d22
d12 = d21
7.503 (0.228)
6.877 (0.231)
7.319 (0.285)
3.369 (1.123)
)
0.000 (
0.090 (0.381)
σ2
1.445 (0.145)
466.173 −
Introduction to Longitudinal Data Analysis
177
·	This suggests that the REML likelihood could be further increased by allowing
negative estimates for d22
·	In SAS, this can be done by adding the option ‘nobound’ to the PROC MIXED
statement.
Results:
·	Eﬀect
Intercept
Time eﬀects:
Low dose
High dose
Control
Covariance of bi:
var(b1i)
var(b2i)
cov(b1i, b2i)
Residual variance:
var(εij)
REML log-likelihood
Parameter restrictions for α
IR, σ2
IR
∈
REMLE (s.e.)
68.618 (0.313)
0, σ2
dii
0
REMLE (s.e.)
68.606 (0.325)
dii
≥
≥
∈
Parameter
β0
β1
β2
β3
7.503 (0.228)
6.877 (0.231)
7.319 (0.285)
7.475 (0.198)
6.890 (0.198)
7.284 (0.254)
d11
d22
d12 = d21
3.369 (1.123)
)
0.000 (
0.090 (0.381)
2.921 (1.019)
0.287 (0.169)
0.462 (0.357)
−
σ2
1.445 (0.145)
1.522 (0.165)
466.173 −
465.193 −
Introduction to Longitudinal Data Analysis
178
·	Note that the REML log-likelihood value has been further increased and a
negative estimate for d22 is obtained.
Brown & Prescott (1999, p. 237) :
·	Introduction to Longitudinal Data Analysis
179
·	Meaning of negative variance component ?
·	Fitted variance function:
Var(Yi(t)) = 
1 t 


D
d
1
t
















σ2
c
d11 +
·	The suggested negative curvature in the variance function is supported by the
σ2 = −0.287t2 + 0.924t + 4.443
d22t2 + 2
d12t +
=
d
d
d
c
sample variance function:
Introduction to Longitudinal Data Analysis
180
·	This again shows that the hierarchical and marginal models are not equivalent:
·	The marginal model allows negative variance components, as long as the
marginal covariances Vi = ZiDZ 0i + σ2Ini are positive deﬁnite
·	The hierarchical interpretation of the model does not allow negative variance
components
Introduction to Longitudinal Data Analysis
181
Chapter 9
Inference for the Marginal Model
·	Inference for ﬁxed eﬀects:
Wald test
t-test and F -test
Robust inference
LR test
∗
∗
∗
∗
·	Inference for variance components:
Wald test
LR test
∗
∗
·	Information criteria
Introduction to Longitudinal Data Analysis
182
9.1 Inference for the Fixed Eﬀects
Estimate for β:
β(α) = 
c


N
Xi=1
−
X 0iWiXi
1 N
Xi=1


X 0iWiyi
with α replaced by its ML or REML estimate
Conditional on α,
β(α) is multivariate normal with mean β and covariance
c
N
Var(
β) = 
c
N



Xi=1
X 0iWivar(Yi)WiXi


1
−
X 0iWiXi


1
−
X 0iWiXi


N



Xi=1
1
−
X 0iWiXi


Xi=1
N
Xi=1




= 
In practice one again replaces α by its ML or REML estimate
§	Introduction to Longitudinal Data Analysis
183
9.1.1 Approximate Wald Test
For any known matrix L, consider testing
H0 : Lβ = 0,
versus HA : Lβ
= 0
Wald test statistic:
G =
β0L0 
L 
c





N
Xi=1
1
X 0iV −
i
α)Xi
(
d


1
−
1
−
L0



L
β
c
Asymptotic null distribution of G is χ2 with rank(L) degrees of freedom
§	Introduction to Longitudinal Data Analysis
184
6
9.1.2 Approximate t-test and F -test
Wald test based on
·	Var(
β) = 
c
N
Xi=1


1
−
X 0iWi(α)Xi


·	Variability introduced from replacing α by some estimate is not taken into
account in Wald tests
§	Therefore, Wald tests will only provide valid inferences in suﬃciently large samples
In practice, this is often resolved by replacing the χ2 distribution by an appropriate
F -distribution (are the normal by a t).
For any known matrix L, consider testing
H0 : Lβ = 0,
versus HA : Lβ
= 0
Introduction to Longitudinal Data Analysis
185
6
F test statistic:
·	β0L0 
L 
F = c





N
Xi=1
1
X 0iV −
i
α)Xi
(
d
rank(L)


1
−
1
−
L0



L
β
c
o	Approximate null-distribution of F is F with numerator degrees of freedom equal
to rank(L)
·	Denominator degrees of freedom to be estimated from the data:
·	Containment method
·	Satterthwaite approximation
·	Kenward and Roger approximation
·	. . .
Introduction to Longitudinal Data Analysis
186
·	In the context of longitudinal data, all methods typically lead to large numbers of
degrees of freedom, and therefore also to very similar p-values.
For univariate hypotheses (rank(L) = 1) the F -test reduces to a t-test
·	Introduction to Longitudinal Data Analysis
187
9.1.3 Example: The Prostate Data
o	Linear hypotheses of the form
H0 : Lβ = 0,
versus HA : Lβ
= 0
can be tested in SAS using a CONTRAST statement
As an example, reconsider the model for the prostate data:
ln(PSAij + 1) = β1Agei + β2Ci + β3Bi + β4Li + β5Mi
·	(β6Agei + β7Ci + β8Bi + β9Li + β10Mi) tij
·	(β11Agei + β12Ci + β13Bi + β14Li + β15Mi) t2
ij
·	b1i + b2itij + b3it2
ij + εij.
·	We now test whether the local cancer cases evolve diﬀerent from the metastatic
cancer cases.
Introduction to Longitudinal Data Analysis
188
6
o	The null-hypothesis is speciﬁed by
H0 :
β4 = β5
β9 = β10
β14 = β15,



This is equivalent with testing
H0 :














0 0 0 1
−
1 0 0 0 0
0 0 0 0 0
0 0 0 0
0 0 0 0 1
1 0 0 0 0
−
0 0 0 0
0 0 0 0 0
0 0 0 0 1
0
0
1














−
β = 0,
which is of the form Lβ = 0
Introduction to Longitudinal Data Analysis
189
o	Related statements in SAS:
model lnpsa = group age grouptime agetime
grouptime2 agetime2
/ noint ddfm=satterth;
contrast ’L/R can = Met can’ group 0 0 1 -1,
grouptime 0 0 1 -1,
grouptime2 0
0 1 -1 / chisq;
Remarks:
·	The Satterthwaite approximation is used for the denominator degrees of
freedom
·	The option ‘chisq’ in CONTRAST statement is needed in order also to obtain
a Wald test
Introduction to Longitudinal Data Analysis
190
·	Additional table in the output :
CONTRAST Statement Results
Source
NDF DDF ChiSq
F Pr > ChiSq Pr > F
L/R can = Met can
3 24.4 17.57 5.86
0.0005 0.0037
·	Several CONTRAST statements can now be used to reduce the model, in a
stepwise procedure
·	This leads to the following simpliﬁcations :
·	no interaction age
·	no interaction age
time2
time
×
×
·	quadratic time eﬀect the same for both cancer groups
·	the quadratic time eﬀect is not signiﬁcant for the non-cancer groups
·	the linear time eﬀect is not signiﬁcant for the controls
Introduction to Longitudinal Data Analysis
191
o	Simultaneous testing of all these hypotheses is testing the null hypothesis
(no age
(no age by time interaction)
(no linear time eﬀect for controls)
β6 = 0
β7 = 0
β11 = 0
β12 = 0
β13 = 0
β14 = β15 (equal quadratic time eﬀect for both cancer groups).
time2 interaction)
(no quadratic time eﬀect for controls)
(no quadratic time eﬀect for BPH)
×
H0 :



This hypothesis is of the form
H0 :
0 0 0 0 0 1 0 0 0 0 0 0 0 0
0 0 0 0 0 0 1 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 1 0 0 0
0 0 0 0 0 0 0 0 0 0 0 1 0 0
0 0 0 0 0 0 0 0 0 0 0 0 1 0
0 0 0 0 0 0 0 0 0 0 0 0 0 1





















0
0
0
0
0
1





















−
β = 0
Introduction to Longitudinal Data Analysis
192
o	The hypothesis can be tested with the following statements:
model lnpsa = group age grouptime agetime grouptime2 agetime2
/ noint ddfm=satterth;
contrast ’Final model’ age*time 1,
grouptime 1 0 0 0,
agetime2 1,
grouptime2 1 0 0 0,
grouptime2 0 1 0 0,
group*time2 0 0 1 -1 / chisq;
This results in the following table in the output (Satterthwaite approximation):
CONTRAST Statement Results
Source
NDF
DDF ChiSq
F Pr > ChiSq Pr > F
Final model
6
46.7 3.39 0.56
0.7587 0.7561
Introduction to Longitudinal Data Analysis
193
The simpliﬁed model is now given by:
·	ln(PSAij + 1) = β1Agei + β2Ci + β3Bi + β4Li + β5Mi
·	(β8Bi + β9Li + β10Mi) tij
·	β14 (Li + Mi) t2
ij
·	b1i + b2itij + b3it2
ij + εij,
·	SAS procedure MIXED also allows using an ESTIMATE statement to estimate
and test linear combinations of the elements of β
·	Using similar arguments as for the approximate Wald-test, t-test, and F -test,
approximate conﬁdence intervals can be obtained for such linear combinations,
also implemented in the ESTIMATE statement.
·	Speciﬁcation of L remains the same as for the CONTRAST statement, but L can
now only contain one row.
Introduction to Longitudinal Data Analysis
194
9.1.4 Robust Inference
§	Estimate for β:
β(α) = 
N
−
X 0iWiXi
c
Xi=1


1 N
Xi=1


with α replaced by its ML or REML estimate
X 0iWiYi
Conditional on α,
β has mean
c
β(α)
E
(cid:20)
(cid:21)
c
N
Xi=1
N
Xi=1
= 


= 


−
X 0iWiXi
−
X 0iWiXi
1 N
Xi=1
1 N
Xi=1




X 0iWiE(Yi)
X 0iWiXiβ = β
provided that E(Yi) = Xiβ
Hence, in order for
is correctly speciﬁed.
c
β to be unbiased, it is suﬃcient that the mean of the response
Introduction to Longitudinal Data Analysis
195
·	Conditional on α,
β has covariance
c
N
Var(
β) = 
c
1
−
N



Xi=1
X 0iWiVar(Yi)WiXi





N
Xi=1
1
−
X 0iWiXi


1
−
X 0iWiXi


X 0iWiXi


Xi=1


N
Xi=1


= 
·	Note that this assumes that the covariance matrix Var(Yi) is correctly modelled as
Vi = ZiDZ 0i + Σi
o	This covariance estimate is therefore often called the ‘naive’ estimate.
The so-called ‘robust’ estimate for Var(
matrix to be correctly speciﬁed is obtained from replacing Var(Yi) by
Yi
0 rather than Vi
Xi
Xi
Yi
β
β
c
β), which does not assume the covariance
(cid:20)
−
(cid:21) (cid:20)
c
−
(cid:21)
c
Introduction to Longitudinal Data Analysis
196
·	Yi
The only condition for
the mean is again correctly speciﬁed.
Xi
Yi
−
β
(cid:21) (cid:20)
c
(cid:20)
Xi
β
c
(cid:21)
−
0 to be unbiased for Var(Yi) is that
·	The so-obtained estimate is called the ‘robust’ variance estimate, also called the
sandwich estimate:
Var(
β) = 
c
N
Xi=1


X 0iWiXi


1
−
N



Xi=1
X 0iWiVar(Yi)WiXi


N



Xi=1
1
−
X 0iWiXi


|
{z
}
|
↓
BREAD
{z
↓
MEAT
}
|
{z
}
↓
BREAD
·	Based on this sandwich estimate, robust versions of the Wald test as well as of
the approximate t-test and F -test can be obtained.
Introduction to Longitudinal Data Analysis
197
§	Note that this suggests that as long as interest is only in inferences for the mean
structure, little eﬀort should be spent in modeling the covariance structure,
provided that the data set is suﬃciently large
Extreme point of view: OLS with robust standard errors
Appropriate covariance modeling may still be of interest:
·	for the interpretation of random variation in data
·	for gaining eﬃciency
·	in presence of missing data, robust inference only valid under very severe
assumptions about the underlying missingness process (see later)
Introduction to Longitudinal Data Analysis
198
9.1.5 Example: Prostate Data
We reconsider the reduced model for the prostate data:
·	ln(PSAij + 1)
= β1Agei + β2Ci + β3Bi + β4Li + β5Mi
·	(β8Bi + β9Li + β10Mi) tij
·	β14 (Li + Mi) t2
ij
·	b1i + b2itij + b3it2
ij + εij,
·	Robust inferences for the ﬁxed eﬀects can be obtained from adding the option
‘empirical’ to the PROC MIXED statement:
proc mixed data=prostate method=reml empirical;
Introduction to Longitudinal Data Analysis
199
Comparison of naive and robust standard errors (only ﬁxed eﬀects !):
·	Eﬀect
Age eﬀect
Intercepts:
Control
BPH
L/R cancer
Met. cancer
Time eﬀects:
BPH
L/R cancer
Met. cancer
Time2 eﬀects:
Cancer
Parameter Estimate (s.e.(1),s.e.(2))
0.016 (0.006;0.006)
β1
β2
β3
β4
β5
β8
β9
β10
−
0.564 (0.428;0.404)
0.275 (0.488;0.486)
1.099 (0.486;0.499)
2.284 (0.531;0.507)
0.410 (0.068;0.067)
1.870 (0.233;0.360)
2.303 (0.262;0.391)
−
−
−
β14 = β15
0.510 (0.088;0.128)
s.e.(1): Naive, s.e.(2): Robust
·	For some parameters, the robust standard error is smaller than the naive,
model-based one. For other parameters, the opposite is true.
Introduction to Longitudinal Data Analysis
200
9.1.6 Example: Growth Data
·	Comparison of naive and robust standard errors under Model 1 (unstructured
mean as well as covariance), for the orthodontic growth data:
Parameter
MLE
(naive s.e.)
(robust s.e.)
β0,8
β0,10
β0,12
β0,14
β1,8
β1,10
β1,12
β1,14
22.8750 23.8125
(0.5598)
(0.4921)
25.7188 (0.6112)
27.4688 21.1818
(0.5371)
(0.6752)
22.2273 (0.5935)
23.0909 24.0909
(0.7372)
(0.6478)
(0.5938)
(0.5170)
(0.6419)
(0.5048)
(0.6108)
(0.5468)
(0.6797)
(0.7007)
How could the covariance structure be improved ?
·	Introduction to Longitudinal Data Analysis
201
o	We ﬁt a model with a separate covariance structure for each group (Model 0)
SAS program:
proc mixed data=test method=ml ;
class idnr sex age;
model measure = age*sex / noint s;
repeated age / type=un subject=idnr r rcorr group=sex;
run;
LR test for Model 1 versus Model 0 : p = 0.0082
The ﬁxed-eﬀects estimates remain unchanged.
The naive standard errors under Model 0 are exactly the same as the sandwich
estimated standard errors under Model 1.
Introduction to Longitudinal Data Analysis
202
9.1.7 Likelihood Ratio Test
·	Comparison of nested models with diﬀerent mean structures, but equal covariance
structure
Null hypothesis of interest equals H0 : β
∈
parameter space Θβ of the ﬁxed eﬀects β.
·	Θβ,0, for some subspace Θβ,0 of the
·	Notation:
·	LML: ML likelihood function
·	θML,0: MLE under H0
c
·	θML: MLE under general model
c
Introduction to Longitudinal Data Analysis
203
Test statistic:
·	2 ln λN =
−
−
2 ln 



LML(
LML(
θML,0)
θML)
c
c




Asymptotic null distribution: χ2 with d.f. equal to the diﬀerence in dimension of
Θβ and Θβ,0.
·	Introduction to Longitudinal Data Analysis
204
9.1.8 Example: Prostate Data
§	We reconsider the reduced model:
ln(PSAij + 1)
= β1Agei + β2Ci + β3Bi + β4Li + β5Mi + (β8Bi + β9Li + β10Mi) tij
ij + b1i + b2itij + b3it2
·	β14 (Li + Mi) t2
ij + εij,
Testing for the need of age correction, i.e., H0 : β1 = 0
Results under ML estimation:
ML estimation
Under β1 ∈
Under H0 : β1 = 0 LML =
LML =
IR
−
3.575 6.876
2 ln λN
−
degrees of freedom
p-value
−
6.602 1
0.010 Introduction to Longitudinal Data Analysis
205
Results under REML estimation:
·	REML estimation
Under β1 ∈
IR
LREML =
Under H0 : β1 = 0 LREML =
−
20.165 19.003
−
−2.324
2 ln λN
−
degrees of freedom
p-value
Negative LR test statistic !
Introduction to Longitudinal Data Analysis
206
9.1.9 LR Test for Fixed Eﬀects Under REML
§	How can the negative LR test statistic be explained ?
Under REML, the response Y is transformed into error contrasts U = A0Y , for
some matrix A with A0X = 0.
Afterwards, ML estimation is performed based on the error contrasts
The reported likelihood value, LREML(
contrasts U
c
θ) is the likelihood at maximum for the error
Models with diﬀerent mean structures lead to diﬀerent sets of error contrasts
Hence, the corresponding REML likelihoods are based on diﬀerent observations,
which makes them no longer comparable
Introduction to Longitudinal Data Analysis
207
Conclusion:
·	LR tests for the mean structure are not valid under REML
Introduction to Longitudinal Data Analysis
208
9.2 Inference for the Variance Components
o	Inference for the mean structure is usually of primary interest.
However, inferences for the covariance structure is of interest as well:
·	interpretation of the random variation in the data
·	overparameterized covariance structures lead to ineﬃcient inferences for mean
·	too restrictive models invalidate inferences for the mean structure
Introduction to Longitudinal Data Analysis
209
9.2.1 Approximate Wald Test
·	Asymptotically, ML and REML estimates of α are normally distributed with
correct mean and inverse Fisher information matrix as covariance
Hence approximate s.e.’s and Wald tests can easily be obtained
·	Introduction to Longitudinal Data Analysis
210
9.2.2 Example: Prostate Data
We reconsider the reduced model:
·	ln(PSAij + 1)
= β1Agei + β2Ci + β3Bi + β4Li + β5Mi + (β8Bi + β9Li + β10Mi) tij
ij + b1i + b2itij + b3it2
·	β14 (Li + Mi) t2
ij + εij,
·	Standard errors and approximate Wald tests for variance components can be
obtained in PROC MIXED from adding the option ‘covtest’ to the PROC MIXED
statement:
proc mixed data=prostate method=reml covtest;
Introduction to Longitudinal Data Analysis
211
Related output:
·	Covariance Parameter Estimates
Cov Parm
Subject
Estimate
UN(1,1)
UN(2,1)
UN(2,2)
UN(3,1)
UN(3,2)
UN(3,3)
timeclss
XRAY
XRAY
XRAY
XRAY
XRAY
XRAY
XRAY
0.4432 -0.4903
0.8416 0.1480
-0.3000
0.1142 0.02837
Standard
Error
0.09349 0.1239
0.2033 0.04702
0.08195 0.03454
0.002276 Z
Value
4.74 -3.96
4.14 3.15
-3.66
3.31 12.47
Pr Z
<.0001
<.0001
<.0001
0.0017 0.0003
0.0005 <.0001
The reported p-values often do not test meaningful hypotheses
The reported p-values are often wrong
There are corrections available in SAS 9,4 (stating with Level 1M3), using the
‘nobound’ option
§	Introduction to Longitudinal Data Analysis
212
9.2.3 Caution with Wald Tests for Variance Components
Marginal versus Hierarchical Model
One of the Wald tests for the variance components in the reduced model for the
prostate data was
·	Cov Parm
Subject
Estimate
Standard
Error
Z
Value
Pr Z
UN(3,3)
XRAY
0.1142 0.03454
3.31 0.0005
o	This presents a Wald test for H0 : d33 = 0
However, under the hierarchical model interpretation, this null-hypothesis is not of
any interest, as d23 and d13 should also equal zero whenever d33 = 0.
·	Hence, the test is meaningful under the marginal model only, i.e., when no
underlying random eﬀects structure is believed to describe the data.
Introduction to Longitudinal Data Analysis
213
Boundary Problems
·	The quality of the normal approximation for the ML or REML estimates strongly
depends on the true value α
·	Poor normal approximation if α is relatively close to the boundary of the
parameter space
If α is a boundary value, the normal approximation completely fails
One of the Wald tests for the variance components in the reduced model for the
prostate data was
Cov Parm
Subject
Estimate
Standard
Error
Z
Value
Pr Z
UN(3,3)
XRAY
0.1142 0.03454
3.31 0.0005
This presents a Wald test for H0 : d33 = 0
§	Introduction to Longitudinal Data Analysis
214
o	Under the hierarchical model interpretation, d33 = 0 is a boundary value, implying
the the calculation of the above p-value is based on an incorrect null-distribution
for the Wald test statistic.
Indeed, how could ever, under H0,
is estimated under the restriction d33 ≥
d
0 ?
d33 be normally distributed with mean 0, if d33
·	Hence, the test is only correct, when the null-hypothesis is not a boundary value
(e.g., H0 : d33 = 0.1).
·	Note that, even under the hierarchical model interpretation, a classical Wald test
is valid for testing H0 : d23 = 0.
Introduction to Longitudinal Data Analysis
215
9.2.4 Likelihood Ratio Test
·	Comparison of nested models with equal mean structures, but diﬀerent covariance
structure
§	Null hypothesis of interest equals H0 : α
parameter space Θα of the variance components α.
∈
Θα,0, for some subspace Θα,0 of the
Notation:
·	LML: ML likelihood function
·	θML,0: MLE under H0
c
·	θML: MLE under general model
c
Test statistic:
2 ln λN =
−
−
2 ln 



LML(
LML(
θML,0)
θML)
c
c




Introduction to Longitudinal Data Analysis
216
Asymptotic null distribution: χ2 with d.f. equal to the diﬀerence in dimension of
Θα and Θα,0.
o	Note that, as long as models are compared with the same mean structure, a valid
LR test can be obtained under REML as well.
·	Indeed, both models can be ﬁtted using the same error contrasts, making the
likelihoods comparable.
Note that, if H0 is a boundary value, the classical χ2 approximation may not be
valid.
o	For some very speciﬁc null-hypotheses on the boundary, the correct asymptotic
null-distribution has been derived
Introduction to Longitudinal Data Analysis
217
9.2.5 Marginal Testing for the Need of Random Eﬀects
o	Under a hierarchical model interpretation, the asymptotic null-distribution for the
LR test statistic for testing signiﬁcance of all variance components related to one
or multiple random eﬀects, can be derived.
Example: for the prostate model, testing whether the variance components
associated to the quadratic random time eﬀect are equal to zero, is equivalent to
testing
H0 : d13 = d23 = d33 = 0
·	Note that, under the hierarchical interpretation of the model, H0 is on the
boundary of the parameter space
Introduction to Longitudinal Data Analysis
218
Case 1: No Random Eﬀects versus one Random Eﬀect
Hypothesis of interest:
·	H0 : D = 0
versus HA : D = d11
for some non-negative scalar d11
·	Asymptotic null-distribution equals
with equal weights 0.5:
2 ln λN −→
−
0:1, the mixture of χ2
χ2
0 and χ2
1
Introduction to Longitudinal Data Analysis
219
Under H0,
−
2 ln λN equals 0 in 50% of the cases
Intuitive explanation:
·	consider the extended parameter space IR for d11
d11 will be negative in 50% of the cases
·	under H0,
·	under the restriction d11 ≥
θML,0) = LML(
·	hence, LML(
d
c
c
Graphically (τ 2 = d11):
0, these cases lead to
d11 = 0
θML) in 50% of the cases
d
§	Introduction to Longitudinal Data Analysis
220
Case 2: One versus two Random Eﬀects
Hypothesis of interest:
·	H0 : D =
d11 0


0 0







,







for d11 > 0, versus HA that D is (2
×
2.	positive semideﬁnite
·	Asymptotic null-distribution:
equal weights 0.5:
2 ln λN −→
−
1:2, the mixture of χ2
χ2
1 and χ2
2 with
Introduction to Longitudinal Data Analysis
221
Case 3: q versus q + 1 Random Eﬀects
Hypothesis of interest:
·	H0 : D =
D11 0
00
0








,








for D11 (q
×
semideﬁnite.
q) positive deﬁnite, versus HA that D is ((q + 1)
(q + 1)) positive
×
·	Asymptotic null-distribution:
with equal weights 0.5.
2 ln λN −→
−
q:q+1, the mixture of χ2
χ2
q and χ2
q+1
Introduction to Longitudinal Data Analysis
222
Case 4: q versus q + k Random Eﬀects
Hypothesis of interest:
H0 : D =
D11 0
0
0








,








for D11 (q
×
semideﬁnite.
q) positive deﬁnite, versus HA that D is ((q + k)
(q + k)) positive
×
Simulations needed to derive asymptotic null distribution
o	Introduction to Longitudinal Data Analysis
223
Conclusions
o	Correcting for the boundary problem reduces p-values
Thus, ignoring the boundary problem too often leads to over-simpliﬁed covariance
structures
·	Hence, ignoring the boundary problem may invalidate inferences, even for the
mean structure
Introduction to Longitudinal Data Analysis
224
9.2.6 Example: Rat Data
We reconsider the model with random intercepts and slopes for the rat data:
·	Yij = (β0 + b1i) + (β1Li + β2Hi + β3Ci + b2i)tij + εij
in which tij equals ln[1 + (Ageij −
45)/10)]
·	The marginal model assumes linear average trends with common intercept for the
3 groups, and covariance structure:
Cov(Yi(t1), Yi(t2)) = 
1 t1 
D


1
t2
















·	σ2δ
{
t1,t2}
= d22t1 t2 + d12(t1 + t2) + d11 + σ2δ
{
·	t1,t2}
Introduction to Longitudinal Data Analysis
225
Exploring the variance function yields:
o	This suggested earlier that the above random-eﬀects model might not be valid, as
it does not allow negative curvature in the variance function
·	It is therefore of interest to test whether the random slopes b2i may be left out of
the model.
·	Interpretation:
·	On hierarchical level: all rats receiving the same treatment have the same slope
·	On marginal level: constant variance, constant correlation
Introduction to Longitudinal Data Analysis
226
o	Null-hypothesis to be tested: H0 : d12 = d22 = 0
REML estimates under hierarchical and marginal interpretation, as well as under
H0:
Eﬀect
Intercept
Time eﬀects:
Low dose
High dose
Control
Covariance of bi:
var(b1i)
var(b2i)
cov(b1i, b2i)
Residual variance:
var(εij)
REML log-likelihood
Parameter restrictions for α
0, σ2
dii
REMLE (s.e.)
68.606 (0.325)
IR, σ2
IR
∈
REMLE (s.e.)
68.618 (0.313)
0 dii
≥
≥
∈
Under H0
REMLE (s.e.)
68.607 (0.331)
Parameter
β0
β1
β2
β3
7.503 (0.228)
6.877 (0.231)
7.319 (0.285)
7.475 (0.198)
6.890 (0.198)
7.284 (0.254)
7.507 (0.225)
6.871 (0.228)
7.507 (0.225)
d11
d22
d12 = d21
3.369 (1.123)
0.000 (
)
0.090 (0.381)
2.921 (1.019)
0.287 (0.169)
0.462 (0.357)
−
3.565 (0.808)
)
(
)
(
σ2
1.445 (0.145)
1.522 (0.165)
1.445 (0.145)
466.173 −
465.193 −
466.202 −
Introduction to Longitudinal Data Analysis
227
Test Under Marginal Interpretation
o	Unrestricted parameter space for α, no boundary problem
Wald test:
·	Test statistic:
d12

 d
d22 
d

Var(
d12)
d
d
Cov(
d12,
d22)
d
d
d
Cov(
d12,
d22)
d
d
d
Var(
d22)
d
d








1
−








d12
d
d22
d
















= 
0.462 
0.287 

−
0.127 0.038
−
0.038 0.029








−
1
−








0.462 0.287








−








= 2.936,
·	p-value:
P (χ2
2 ≥
2.936 |
H0) = 0.2304
Introduction to Longitudinal Data Analysis
228
·	LR test:
·	Test statistic:
·	p-value:
2 ln λN =
2(
−
−
−
466.202 + 465.193) = 2.018
P (χ2
2 ≥
2.018 |
H0) = 0.3646
Introduction to Longitudinal Data Analysis
229
Test Under Hierarchical Interpretation
§	Restricted parameter space for α (positive semi-deﬁnite D), boundary problem !
LR test statistic:
2 ln λN =
2(
−
−
−
466.202 + 466.173) = 0.058
p-value:
P (χ2
1:2 ≥
0.058 H0)
|
= 0.5 P (χ2
1 ≥
0.058 |
H0) + 0.5 P (χ2
2 ≥
0.058 |
H0) = 0.8906
·	Note that the naive p-value, obtained from ignoring the boundary problem is
indeed larger:
P (χ2
2 ≥
0.058 |
H0) = 0.9714
Introduction to Longitudinal Data Analysis
230
Reduced Model
§	Under both model interpretations, H0 was accepted, leading to the reduced model:
Yij = (β0 + b1i) + (β1Li + β2Hi + β3Ci)tij + εij
Marginal interpretation:
·	linear average trends with common intercept for the 3 groups
·	constant variance estimated to be
d11 +
σ2 = 3.565 + 1.445 = 5.010
·	constant (intraclass) correlation
c
d
ρI =
d
d11
d11 +
d
d
σ2 = 0.712
c
The hierarchical interpretation, possible since
heterogeneity between rats is restricted to diﬀerences in starting values, not slopes.
d11 = 3.565 > 0, is that
d
Introduction to Longitudinal Data Analysis
231
9.3 Information Criteria
9.3.1 Deﬁnition of Information Criteria
LR tests can only be used to compare nested models
How to compare non-nested models ?
The general idea behind the LR test for comparing model A to a more extensive
model B is to select model A if the increase in likelihood under model B is small
compared to increase in complexity
A similar argument can be used to compare non-nested models A and B
·	Introduction to Longitudinal Data Analysis
232
·	One then selects the model with the largest (log-)likelihood provided it is not
(too) complex
The model is selected with the highest penalized log-likelihood `
some function
) of the number #θ of parameters in the model.
− F
(
·
F
(#θ) for
Diﬀerent functions
(
·
F
) lead to diﬀerent criteria:
o	Criterion
Akaike (AIC)
Schwarz (SBC)
Hannan and Quinn (HQIC)
Bozdogan (CAIC)
?: n∗ = n =
?: n∗ = n
N
i=1 ni under ML
P
p under REML
−
Deﬁnition of
)?
(
·
F
(#θ) = #θ
(#θ) = (#θ ln n∗)/2
(#θ) = #θ ln(ln n∗)
(#θ) = #θ (ln n∗ + 1)/2
F
F
F
F
Introduction to Longitudinal Data Analysis
233
o	Information criteria are no formal testing procedures !
For the comparison of models with diﬀerent mean structures, information criteria
should be based on ML rather than REML, as otherwise the likelihood values
would be based on diﬀerent sets of error contrasts, and therefore would no longer
be comparable.
Introduction to Longitudinal Data Analysis
234
9.3.2 Example: Rat Data
Consider the random-intercepts model for the rat data:
·	Yij = (β0 + b1i) + (β1Li + β2Hi + β3Ci)tij + εij
in which tij equals ln[1 + (Ageij −
45)/10)]
·	We now want to compare this model with a model which assumes common
average slope for the 3 treatments.
·	Information criteria can be obtained in SAS from adding the option ‘ic’ to the
PROC MIXED statement:
proc mixed data=rats method=ml ic;
Introduction to Longitudinal Data Analysis
235
o	Summary of results:
Mean structure
`ML
#θ
AIC
SBC
Separate average slopes
Common average slope
464.326 466.622
−
−
6
4
470.326 470.622
−
−
480.914 477.681
−
−
Selected models:
·	AIC: model with separate slopes
·	SBC: model with common slopes
·	Based on Wald test, the average slopes are found not to be signiﬁcantly diﬀerent
from each other (p = 0.0987)
Introduction to Longitudinal Data Analysis
236
Chapter 10
Inference for the Random Eﬀects
·	Empirical Bayes inference
·	Best linear unbiased prediction
·	Example: Prostate data
·	Shrinkage
·	Example: Random-intercepts model
·	Example: Prostate data
·	Normality assumption for random eﬀects
Introduction to Longitudinal Data Analysis
237
10.1 Empirical Bayes Inference
·	Random eﬀects bi reﬂect how the evolution for the ith subject deviates from the
expected evolution Xiβ.
Estimation of the bi helpful for detecting outlying proﬁles
This is only meaningful under the hierarchical model interpretation:
Yi
bi
|
∼
N (Xiβ + Zibi, Σi)
bi
∼
N (0, D)
Since the bi are random, it is most natural to use Bayesian methods
Terminology: prior distribution N (0, D) for bi
·	Introduction to Longitudinal Data Analysis
238
Posterior density:
·	f (bi
yi)
|
f (bi
|
Yi = yi) =
Z
f (yi
|
bi) f (bi)
·	. .
≡
∝
∝
f (yi
f (yi
bi) f (bi)
|
bi) f (bi) dbi
|
exp 

−
∝
1
2
(bi
−
DZ 0iWi(yi
−
1
Xiβ))0 Λ−
i
(bi
DZ 0iWi(yi
Xiβ))

−
−

for some positive deﬁnite matrix Λi.
Posterior distribution:
·	yi
bi
|
∼
N (DZ 0iWi(yi
Xiβ), Λi)
−
Introduction to Longitudinal Data Analysis

239
Posterior mean as estimate for bi:
·	bi(θ) = E [bi
d
|
Yi = yi] =
bi f (bi
|
Z
yi) dbi = DZ 0iWi(α)(yi
Xiβ)
−
bi(θ) is normally distributed with covariance matrix
·	d
var(
bi(θ)) = DZ 0i 

d
Wi −

WiXi 
N
Xi=1


1
−
X 0iWiXi


X 0iWi

ZiD

Note that inference for bi should account for the variability in bi
Therefore, inference for bi is usually based on
o	var(
bi(θ)
d
−
bi) = D
var(
bi(θ))
d
−
Introduction to Longitudinal Data Analysis
240
o	Wald tests can be derived
Parameters in θ are replaced by their ML or REML estimates, obtained from
ﬁtting the marginal model.
bi =
·	d
bi(
θ) is called the Empirical Bayes estimate of bi.
d
c
·	Approximate t- and F -tests to account for the variability introduced by replacing
θ by
θ, similar to tests for ﬁxed eﬀects.
c
Introduction to Longitudinal Data Analysis
241
10.2 Best Linear Unbiased Prediction (BLUP)
·	Often, parameters of interest are linear combinations of ﬁxed eﬀects in β and
random eﬀects in bi
§	For example, a subject-speciﬁc slope is the sum of the average slope for subjects
with the same covariate values, and the subject-speciﬁc random slope for that
subject.
In general, suppose u = λ0ββ + λ0bbi is of interest
Conditionally on α,
u = λ0β
β + λ0b
bi is BLUP:
·	linear in the observations Yi
c
c
d
·	unbiased for u
·	minimum variance among all unbiased linear estimators
Introduction to Longitudinal Data Analysis
242
10.3 Example: Prostate Data
We reconsider the reduced model:
·	ln(PSAij + 1)
= β1Agei + β2Ci + β3Bi + β4Li + β5Mi + (β8Bi + β9Li + β10Mi) tij
ij + b1i + b2itij + b3it2
·	β14 (Li + Mi) t2
ij + εij
·	In SAS the estimates can be obtained from adding the option ‘solution’ to the
random statement:
random intercept time time2
/ type=un subject=id solution;
ods listing exclude solutionr;
ods output solutionr=out;

```r
# R translation: EB (BLUP) random effects analogous to SAS 'solution' on RANDOM
library(here)
library(haven)
library(lme4)
library(dplyr)
library(ggplot2)
library(tibble)

# Load prostate data and construct outcome (aligns with Chapter 9 reduced model)
prostate <- read_sas(here("year2/Longitudinal_Data_Analysis", "data", "prostate.sas7bdat")) %>%
	mutate(lpsa = log(psa + 1))

# Fit reduced model with unstructured random-effects covariance for (Intercept, time, time^2)
# In lme4, (1 + time + I(time^2) | id) estimates a full (UN) D matrix
m_reduced <- lmer(
	lpsa ~ age + cap + bon + lun + met + (bon + lun + met):time + I(lun + met):I(time^2) +
		(1 + time + I(time^2) | id),
	data = prostate,
	REML = TRUE
)

# Empirical Bayes (BLUP) random effects: like SAS solutionR written to a dataset
re_list <- ranef(m_reduced, condVar = TRUE)  # attach conditional var-cov via postVar attribute
re_id <- re_list$id                          # random effects for grouping factor 'id'

# Tidy into a data frame
re_df <- re_id %>%
	as.data.frame() %>%
	rownames_to_column(var = "id")

# Extract conditional standard errors for each random-effect component
pv <- attr(re_id, "postVar")                # array [p x p x n_subjects]
se_intercept <- sqrt(pv[1, 1, ])
se_time      <- sqrt(pv[2, 2, ])
se_time2     <- sqrt(pv[3, 3, ])

re_df <- re_df %>%
	mutate(
		se_intercept = se_intercept,
		se_time = se_time,
		se_time2 = se_time2,
		z_intercept = `(Intercept)` / se_intercept,
		z_time      = time / se_time,
		z_time2     = `I(time^2)` / se_time2
	)

# Examples: histograms and scatterplots to detect exceptional evolutions over time
ggplot(re_df, aes(`(Intercept)`)) + geom_histogram(bins = 30) + labs(title = "EB random intercepts")
ggplot(re_df, aes(time)) + geom_histogram(bins = 30) + labs(title = "EB random slopes (time)")
ggplot(re_df, aes(`I(time^2)`)) + geom_histogram(bins = 30) + labs(title = "EB random slopes (time^2)")

ggplot(re_df, aes(time, `I(time^2)`)) +
	geom_point(alpha = 0.7) +
	labs(title = "EB slope vs quadratic term", x = "random slope for time", y = "random slope for time^2")

# Optional: subject-specific (conditional) coefficients (fixed + random)
coef_df <- coef(m_reduced)$id %>%
	as.data.frame() %>%
	rownames_to_column(var = "id")
```
Introduction to Longitudinal Data Analysis
243
·	The ODS statements are used to write the EB estimates into a SAS output data
set, and to prevent SAS from printing them in the output window.
·	In practice, histograms and scatterplots of certain components of
detect model deviations or subjects with ‘exceptional’ evolutions over time
d
bi are used to
Introduction to Longitudinal Data Analysis
244
I
n
t
r
o
d
u
c
t
i
o
n
t
o
L
o
n
g
i
t
u
d
n
a
i
l
D
a
t
a
A
n
a
l
y
s
i
s
2
4
5
·	Strong negative correlations in agreement with correlation matrix corresponding to
ﬁtted D:
Dcorr =
d
1.000 −
0.803 0.658
0.803 1.000
0.968 −
0.658 −
0.968 1.000














−














o	Histograms and scatterplots show outliers
Subjects #22, #28, #39, and #45, have highest four slopes for time2 and
smallest four slopes for time, i.e., with the strongest (quadratic) growth.
·	Subjects #22, #28 and #39 have been further examined and have been shown to
be metastatic cancer cases which were misclassiﬁed as local cancer cases.
Subject #45 is the metastatic cancer case with the strongest growth
·	Introduction to Longitudinal Data Analysis
246
10.4 Shrinkage Estimators
bi
d
o	Consider the prediction of the evolution of the ith subject:
Yi
d
≡
Xi
β + Zi
bi
c
d
= Xi
1
β + ZiDZ 0iV −
i
(yi
c
Xi
β)
c
−
=
Ini −
(cid:18)
1
ZiDZ 0iV −
i
(cid:19)
Xi
1
i yi
β + ZiDZ 0iV −
c
1
i Xi
= ΣiV −
β +
c
Ini −
(cid:18)
1
ΣiV −
i
yi,
(cid:19)
Hence,
observed data yi, with weights
Yi is a weighted mean of the population-averaged proﬁle Xi
1
respectively.
V −
i
d
Σi
and Ini − c
1
V −
i
d
Σi
d
c
c
β and the
Introduction to Longitudinal Data Analysis
247
o	Note that Xi
to the total variability.
c
β gets much weight if the residual variability is ‘large’ in comparison
This phenomenon is usually called shrinkage :
The observed data are shrunk towards the prior average
proﬁle Xiβ.
·	This is also reﬂected in the fact that for any linear combination λ0bi of random
eﬀects,
var(λ0
bi)
d
≤
var(λ0bi).
Introduction to Longitudinal Data Analysis
248
10.5 Example: Random-intercepts Model
o	Consider the random-intercepts model, without serial correlation:
·	Zi = 1ni, vector of ones
·	D = σ2
b , scalar
·	Σi = σ2Ini
The EB estimate for the random intercept bi then equals
bi = σ2
b
1ni0
c
σ2
b
1ni
(cid:18)
1ni0 + σ2Ini (cid:19)
1
−
(yi
Xiβ)
−
=
=
σ2
b
σ2
1ni0 



Ini −
1
ni
niσ2
b
σ2 + niσ2
b
σ2
b
σ2 + niσ2
b
1ni
1ni0



(yi
−
Xiβ)
ni
(yij −
Xj=1
X [j]
i β)
Introduction to Longitudinal Data Analysis
249
·	Remarks:
·	bi is weighted average of 0 (prior mean) and the average residual for subject i
c
·	less shrinkage the larger ni
·	less shrinkage the smaller σ2 relative to σ2
b
Introduction to Longitudinal Data Analysis
250
10.6 Example: Prostate Data
·	Comparison of predicted, average, and observed proﬁles for the subjects #15 and
#28, obtained under the reduced model:
Introduction to Longitudinal Data Analysis
251
Illustration of the shrinkage eﬀect :
·	Var(
bi) =
d
d
0.403 −
0.440 0.131
0.440 0.729
0.253 −
0.131 −
0.253 0.092














−
,














D =
d
0.443 −
0.490 0.148
0.490 0.842
0.300 −
0.148 −
0.300 0.114














−














Introduction to Longitudinal Data Analysis
252
10.7 The Normality Assumption for Random Eﬀects
·	In practice, histograms of EB estimates are often used to check the normality
assumption for the random eﬀects
However, since
·	bi = DZ 0iWi(yi
d
Xiβ)
−
var(
bi) = DZ 0i 

d
Wi −

WiXi 
N
Xi=1


1
−
X 0iWiXi


X 0iWi

ZiD

one should at least ﬁrst standardize the EB estimates
·	Further, due to the shrinkage property the EB estimates do not fully reﬂect the
heterogeneity in the data.
Introduction to Longitudinal Data Analysis
253
·	Small simulation example:
·	1000 proﬁles with 5 measurements, balanced
·	1000 random intercepts sampled from
1
2
N (
−
2, 1) +
1
2
N (2, 1)
·	Σi = σ2Ini, σ2 = 30
·	Data analysed assuming normality for the intercepts
Introduction to Longitudinal Data Analysis
254
·	Histogram of sampled intercepts and empirical Bayes estimates:
·	Clearly, severe shrinkage forces the estimates
bi to satisfy the normality
assumption
b
Introduction to Longitudinal Data Analysis
255
Conclusion:
·	EB estimates obtained under normality
cannot be used to check normality
·	This suggests that the only possibility to check the normality assumption is to ﬁt
a more general model, with the classical linear mixed model as special case, and to
compare both models using LR methods
Introduction to Longitudinal Data Analysis
256
10.8 The Heterogeneity Model
·	One possible extension of the linear mixed model is to assume a ﬁnite mixture as
random-eﬀects distribution:
bi
g
∼
Xj=1
pjN (µj, D),
with
g
Xj=1
pj = 1 and
g
Xj=1
pjµj = 0
Interpretation:
·	Population consists of g subpopulations
·	Each subpopulation contains fraction pj of total population
·	In each subpopulation, a linear mixed model holds
The classical model is a special case: g = 1
o	Introduction to Longitudinal Data Analysis
257
Very ﬂexible class of parametric models for random-eﬀects distribution:
·	Introduction to Longitudinal Data Analysis
258
·	Fitting of the model is based on the EM algorithm
SAS macro available
EB estimates can be calculated under the heterogeneity model
Small simulation example:
·	1000 proﬁles with 5 measurements, balanced
·	1000 random intercepts sampled from
1
2
N (
−
2, 1) +
1
2
N (2, 1)
·	Σi = σ2Ini, σ2 = 30
·	Data analysed under heterogeneity model
Introduction to Longitudinal Data Analysis
259
·	Histogram of sampled intercepts and empirical Bayes estimates:
·	The correct random-eﬀects distribution is (much) better reﬂected, than before
under the assumption of normality
Introduction to Longitudinal Data Analysis
260
Chapter 11
General Guidelines for Model Building
·	Introduction
·	General strategy
·	Example: The prostate data
Introduction to Longitudinal Data Analysis
261
11.1 Introduction
Marginal linear mixed model:
·	Yi
∼
N (Xiβ, ZiDZ 0i + σ2Ini + τ 2Hi)
·	Fitting a linear mixed model requires speciﬁcation of a mean structure, as well as
covariance structure
·	Mean structure:
·	Covariates
·	Time eﬀects
·	Interactions
·	Covariance structure:
·	Random eﬀects
·	Serial correlation
Introduction to Longitudinal Data Analysis
262
Both components aﬀect each other:
'
$
·	Mean structure Xiβ
Covariance structure Vi
&
(cid:27)
?
Estimation of θ
Covariance matrix for
t-tests and F-tests
Conﬁdence intervals
Eﬃciency
Prediction
θ
c
&
%
?
%
Introduction to Longitudinal Data Analysis
263
·	When most variability is due to between-subject variability, the two-stage
approach will often lead to acceptable marginal models
·	In the presence of a lot within-subject variability, the two-stage approach is less
straightforward
Also, a two-stage approach may imply unrealistic marginal models
·	Introduction to Longitudinal Data Analysis
264
·	For example, reconsider the growth curves:
·	Individual proﬁles:
·	A random-intercepts model seems reasonable
Introduction to Longitudinal Data Analysis
265
·	However, the covariance matrix equals
6.11 6.88
8.26 7.44
7.18 6.88 8.53
9.78 9.01
8.70 8.26 9.78 12.04 10.99 10.96
7.44 9.01 10.99 10.42 10.56
7.18 8.70 10.96 10.56 11.24



























·	


























The aim of this chapter is to discuss some general guidelines for model building.
·	Introduction to Longitudinal Data Analysis
266
11.2 General Strategy
Yi = Xiβ + Zibi + εi
1.	Preliminary mean structure Xiβ
2.	Preliminary random-eﬀects structure Zibi
3.	Residual covariance structure Σi
4.	Reduction of the random-eﬀects structure Zibi
5.	Reduction of the mean structure Xiβ
Introduction to Longitudinal Data Analysis
267
11.3 Preliminary Mean Structure
11.3.1 Strategy
o	Remove all systematic trends from the data, by calculating OLS residual proﬁles :
ri = yi
Xi
βOLS ≈
c
−
Zibi + εi
For balanced designs with few covariates :
Saturated mean structure
Introduction to Longitudinal Data Analysis
268
For balanced designs with many covariates, or for highly unbalanced data sets :
·	The most elaborate model one is prepared
to consider for the mean structure
·	Selection of preliminary mean structures will be based on exploratory tools for the
mean.
·	Note that the calculation of
obtained in any regression module
c
βOLS ignores the longitudinal structure, and can be
·	Provided the preliminary mean structure is ‘suﬃciently richt’, consistency of
follows from the theory on robust inference for the ﬁxed eﬀects.
βOLS
c
Introduction to Longitudinal Data Analysis
269
11.3.2 Example: Prostate Data
Smoothed average trend within each group:
·	Introduction to Longitudinal Data Analysis
270
§	Quadratic function over time, within each diagnostic group
Correction for age, via the inclusion of age, age
time and age
×
time2.
×
Note that this yields the same model as the model originally obtained from a
two-stage approach, containing 15 ﬁxed eﬀects
Introduction to Longitudinal Data Analysis
271
11.4 Preliminary Random-eﬀects Structure
11.4.1 Stragegy
ri
≈
Zibi + εi
Explore the residual proﬁles
Any structure left, may indicate the presence of subject-speciﬁc regression
coeﬃcients
Try to describe the each residual proﬁle with a (relatively) simple model.
§	Introduction to Longitudinal Data Analysis
272
·	Do not include covariates in Zi which are not included in Xi. Otherwise, it is not
justiﬁed to assume E(bi) = 0.
·	Use ‘well-formulated’ models: Do not include higher-order terms unless all
lower-order terms are included as well.
·	Compare implied variance and covariance functions with results from exploratory
tools for covariance structure
Introduction to Longitudinal Data Analysis
273
11.4.2 Example: Prostate Data
§	OLS residual proﬁles and smoothed average of squared OLS residuals:
We assume a quadratic function for each residual proﬁle
This results in a model with random intercepts, and random slopes for the linear
as well as quadratic time eﬀect.
Introduction to Longitudinal Data Analysis
274
Variance function:
·	1 t t2
D




·	σ2
1
t
t2




























·	Comparison of smoothed average of squared OLS residuals and ﬁtted variance
function:
Introduction to Longitudinal Data Analysis
275
·	Possible explanation for observed diﬀerences:
·	Small t: some subjects have extremely large responses close to diagnosis. This
may have inﬂated the ﬁtted variance
·	Large t: few observations available: only 24 out of 463 measurements taken
earlier than 20 years prior to diagnosis.
Introduction to Longitudinal Data Analysis
276
11.5 Residual Covariance Structure
11.5.1 Strategy
ri
≈
Zibi + εi
Which covariance matrix Σi for εi ?
In many applications, random eﬀects explain most of the variability
Therefore, in the presence of random eﬀects other than intercepts, often
Σi = σ2Ini is assumed
§	Introduction to Longitudinal Data Analysis
277
However, many other covariance structures can be speciﬁed as well
·	Introduction to Longitudinal Data Analysis
278
·	A special class of parametric models for Σi is obtained from splitting εi into a
measurement error component ε(1)i and a serial correlation component ε(2)i:
Yi = Xiβ + Zibi + ε(1)i + ε(2)i
N (0, D)
N (0, σ2Ini)
N (0, τ 2Hi)
bi
∼
ε(1)i
∼
ε(2)i
∼
independent



o	Only the correlation matrix Hi then still needs to be speciﬁed
Hi is assumed to have (j, k) element of the form hijk = g(
|
) with g(0) = 1
decreasing function g(
·
tij −
tik|
) for some
Introduction to Longitudinal Data Analysis
279
Frequently used functions g(
·
):
·	Exponential serial correlation: g(u) = exp(
φu)
·	Gaussian serial correlation: g(u) = exp(
−
Graphical representation (φ = 1):
−
φu2)
o	Introduction to Longitudinal Data Analysis
280
·	When only random intercepts are included, the semi-variogram can be used to
explore the presence and the nature of serial correlation
·	When other random eﬀects are present as well, an extension of the variogram is
needed.
Also, a variety of serial correlation functions can be ﬁtted and compared.
·	Introduction to Longitudinal Data Analysis
281
11.5.2 Example: Prostate Data
·	Based on the preliminary mean and random-eﬀects structures, several serial
correlation functions can be ﬁtted.
·	For example, a model with Gaussian serial correlation can be ﬁtted in SAS using
the following program:
proc mixed data=prostate method=reml;
class id group timeclss;
model lnpsa = group age grouptime agetime grouptime2 agetime2 / noint solution;
random intercept time time2 / type=un subject=id g gcorr v vcorr;
repeated timeclss / type=sp(gau)(time) local subject=id r rcorr;
run;
·	REPEATED statement:
·	the serial correlation model is speciﬁed in the ‘type’ option
·	‘local’ is added to include measurement error
Introduction to Longitudinal Data Analysis
282
§	Summary of model ﬁts:
Residual covariance structure
REML log-likelihood
Measurement error
Measurement error + Gaussian
Measurement error + exponential
31.235 24.787
24.266 −
−
−
The presence of serial correlation is clearly detected
However, there seems to be little information in the data to distinguish between
diﬀerent serial correlation structures
·	Practical experience suggests that including serial correlation, if present, is far
more important than correctly specifying the serial correlation function.
Introduction to Longitudinal Data Analysis
283
Variance function:
·	1 t t2
D




·	σ2 + τ 2
1
t
t2




























·	Comparison of smoothed average of squared OLS residuals and ﬁtted variance
function:
Introduction to Longitudinal Data Analysis
284
·	Inclusion of serial correlation leads to diﬀerent estimates for the variance
components in D
·	Therefore, the ﬁtted variance function diﬀers from the one obtained before
without serial correlation
·	The deviation for small values of t remains, but the functions coincide better for
large t.
Introduction to Longitudinal Data Analysis
285
11.6 Reduction of Preliminary Random-eﬀects Structure
·	Once an appropriate residual covariance model is obtained, one can try to reduce
the number of random eﬀects in the preliminary random-eﬀects structure
This is done based on inferential tools for variance components
·	Introduction to Longitudinal Data Analysis
286
11.7 Reduction of Preliminary Mean Structure
·	Once an appropriate covariance model is obtained, one can try to reduce the
number of covariates in the preliminary mean structure
o	This is done based on inferential tools for ﬁxed eﬀects
In case there is still some doubt about the validity of the marginal covariance
structure, robust inference can be used to still obtain correct inferences.
Introduction to Longitudinal Data Analysis
287
11.8 Example: Prostate Data
·	Fixed eﬀects estimates from the ﬁnal model, under Gaussian serial correlation, and
without serial correlation:
Eﬀect
Age eﬀect
Intercepts:
Control
BPH
L/R cancer
Met. cancer
Time eﬀects:
BPH
L/R cancer
Met. cancer
Time2 eﬀects:
Cancer
Serial corr.
No serial corr.
Parameter Estimate (s.e.) Estimate (s.e.)
β1
β2
β3
β4
β5
β8
β9
β10
0.015 (0.006)
0.016 (0.006)
−
0.496 (0.411)
0.320 (0.470)
−
0.564 (0.428)
0.275 (0.488)
1.216 (0.469)
1.099 (0.486)
2.353 (0.518)
2.284 (0.531)
0.376 (0.070)
0.410 (0.068)
1.877 (0.210)
2.274 (0.244)
1.870 (0.233)
2.303 (0.262)
−
−
−
−
−
−
β14 = β15
0.484 (0.073)
0.510 (0.088)
Introduction to Longitudinal Data Analysis
288
·	Variance components estimates from the ﬁnal model, under Gaussian serial
correlation, and without serial correlation:
Eﬀect
Covariance of bi:
var(b1i)
var(b2i)
var(b3i)
cov(b1i, b2i)
cov(b2i, b3i)
cov(b3i, b1i)
Measurement error variance:
var(ε(1)ij)
Gaussian serial correlation:
Serial corr.
No serial corr.
Parameter
Estimate (s.e.)
Estimate (s.e.)
d11
d22
d33
0.393 (0.093)
0.550 (0.187)
0.443 (0.093)
0.842 (0.203)
0.056 (0.028)
0.114 (0.035)
d12 = d21 −
d23 = d32 −
d13 = d31
0.382 (0.114)
0.170 (0.070)
0.490 (0.124)
0.300 (0.082)
−
−
0.098 (0.039)
0.148 (0.047)
σ2
0.023 (0.002)
0.028 (0.002)
var(ε(2)ij)
Rate of exponential decrease
τ 2
1/√φ
0.029 (0.018)
0.599 (0.192)
(
(
)
)
REML log-likelihood
13.704 −
20.165 −
Introduction to Longitudinal Data Analysis
289
·	Many standard errors are smaller under the model which includes the Gaussian
serial correlation component
·	Hence, adding the serial correlation leads to more eﬃcient inferences for most
parameters in the marginal model.
Introduction to Longitudinal Data Analysis
290
11.9 Random-eﬀects Structure versus Residual Covariance
Structure
The marginal covariance structue equals
·	Vi = ZiDZ 0i + Σi
·	Hence, the residual covariance Σi models all variation not yet been accounted for
by random eﬀects
·	In practice, one therefore often observes strong competition between these two
sources of stochastic variation
·	This is also reﬂected in substantial correlations between the variance components
estimates
Introduction to Longitudinal Data Analysis
291
·	As an example, consider the ﬁnal model for the prostate data, with Gaussian serial
correlation
Estimated correlation matrix for variance components estimates:
·	Corr 
d11,
d12,
d22,
d13,
d23,
d33,
τ 2, 1/s
φ,
σ2
d
d
d
d
d
d
1.00 0.87
0.62 0.70
0.49 0.39
0.18 0.10
0.00 −
−
−
−
−
0.87 1.00
0.85 0.94
0.75 0.63
0.21 0.08
0.03 −
−
−
−
−
0.62 0.85
1.00 0.88
0.97 0.91
−
−
0.70 0.94
0.88 1.00
0.82 0.72
−
−
−
−
0.46 0.29
0.02 −
−
0.22 0.06
0.05 0.49
0.75 0.97
0.82 1.00
0.97 0.51
0.33 0.02
−
−
−
−
−





























c
0.39 0.63
0.91 0.72
0.97 1.00
−
−
−
−
0.57 0.38
0.01 c
c
0.18 −
−
0.10 0.21
0.08 0.46
0.22 −
−
0.29 0.06
−
−
0.51 0.33
0.57 −
−
0.38 1.00
0.81 0.04
0.81 1.00
0.32 =


−
−
−
0.00 0.03
0.02 0.05
0.02 0.01
0.04 0.32
1.00 
·	


























Introduction to Longitudinal Data Analysis
292
·	Relatively large correlations between
parameters in D
τ 2 and the estimates of some of the
c
Small correlations between
σ2 and the other estimates, except for 1/
r
c
φ.
c
Indeed, the serial correlation component vanishes for φ becoming inﬁnitely large.
o	Introduction to Longitudinal Data Analysis
293
Chapter 12
Power Analyses under Linear Mixed Models
·	F test for ﬁxed eﬀects
·	Calculation in SAS
·	Examples
Introduction to Longitudinal Data Analysis
294
12.1 F Statistics for Fixed Eﬀects
o	Consider a general linear hypothesis
H0 : Lβ = 0,
versus HA : Lβ
= 0
F test statistic:
β0L0 
L 
F = c





N
Xi=1
1
X 0iV −
i
α)Xi
(
d
rank(L)


1
−
1
−
L0



L
β
c
o	Approximate null-distribution of F is F with numerator degrees of freedom equal
to rank(L)
Introduction to Longitudinal Data Analysis
295
6
·	Denominator degrees of freedom to be estimated from the data:
·	Containment method
·	Satterthwaite approximation
·	Kenward and Roger approximation
·	. . .
·	In general (not necessarily under H0), F is approximately F distributed with the
same numbers of degrees of freedom, but with non-centrality parameter
φ = β0L0 
L 
N
Xi=1


1
X 0iV −
i
1
−
α)Xi
(
d


L0



1
−
Lβ
which equals 0 under H0.



·	This can be used to calculate powers under a variety of models, and under a
variety of alternative hypotheses
Introduction to Longitudinal Data Analysis
296
Note that φ is equal to rank(L)
F , and with
β replaced by β
c
×
The SAS procedure MIXED can therefore be used for the calculation of φ and the
related numbers of degrees of freedom.
o	Introduction to Longitudinal Data Analysis
297
12.2 Calculation in SAS
·	Construct a data set of the same dimension and with the same covariates and
factor values as the design for which power is to be calculated
Use as responses yi the average values Xiβ under the alternative model
The ﬁxed eﬀects estimate will then be equal to
1 N
N
β(α) = 
c
−
X 0iWi(α)Xi
Xi=1


N
Xi=1


= 
−
X 0iWi(α)Xi
X 0iWi(α)yi
X 0iWi(α)Xiβ = β
Xi=1
1 N
Xi=1




Hence, the F -statistic reported by SAS will equal φ/rank(L)
§	Introduction to Longitudinal Data Analysis
298
·	This calculated F value, and the associated numbers of degrees of freedom can be
saved and used afterwards for calculation of the power.
·	Note that this requires keeping the variance components in α ﬁxed, equal to the
assumed population values.
Steps in calculations:
·	Use PROC MIXED to calculate φ, and degrees of freedom ν1 and ν2
·	Calculate critical value Fc:
P (Fν1,ν2,0 > Fc) = level of signiﬁcance
·	Calculate power:
power = P (Fν1,ν2,φ > Fc)
The SAS functions ‘ﬁnv’ and ‘probf’ are used to calculated Fc and the power
o	Introduction to Longitudinal Data Analysis
299
12.3 Example 1
o	Re-consider the random-intercepts model previously discussed for the rat data:
Yij = (β0 + b1i) + (β1Li + β2Hi + β3Ci)tij + εij
in which tij equals ln[1 + (Ageij −
45)/10)]
This model is ﬁtted in SAS as follows:
proc mixed data = test;
class treat rat;
model y = treatt / solution ddfm=kr;
random intercept / subject=rat;
contrast ’Equal slopes’ treatt 1 -1 0, treat*t 1 0 -1;
run;
Introduction to Longitudinal Data Analysis
300
o	The CONTRAST statement is added to test equality of the average slopes.
Suppose a new experiment is to be designed, to test the above hypothesis, when
the true parameter values are given by:
Eﬀect
Intercept
Time eﬀects:
Low dose
High dose
Control
Covariance of bi:
var(b1i)
Residual variance:
var(εij)
Parameter True value
β0
β1
β2
β3
d11
σ2
68
7
7.5 6.5
3.6 1.4
Introduction to Longitudinal Data Analysis
301
·	The power of a design with 10 rats per treatment group is calculated as follows:
·	Construction of data set with expected averages as response values:
data power;
do treat=1 to 3;
do rat=1 to 10;
do age=50 to 110 by 10;
t=log(1+(age-45)/10);
if treat=1 then y=68 + 7.5t;
if treat=2 then y=68 + 7.0t;
if treat=3 then y=68 + 6.5*t;
output;
end;
end;
end;
Introduction to Longitudinal Data Analysis
302
·	Fit model, keeping the variance components equal to their true values:
proc mixed data = power noprofile;
class treat rat;
model y = treatt ;
random intercept / subject=rat(treat);
parms (3.6) (1.4) / noiter;
contrast ’Equal slopes’ treatt 1 -1 0,
treat*t 1 0 -1;
ods output contrasts=c;
run;
·	PARMS statement to specify starting values for the variance components.
·	The ‘noiter’ and ‘noproﬁle’ options request that no iterations be performed and
that inferences are based on the speciﬁed values.
·	ODS statement needed to save F , ν1 and ν2.
Introduction to Longitudinal Data Analysis
303
·	Calculation of φ, Fc and power:
data power;
set c;
alpha=0.05;
ncparm=numdf*fvalue;
fc=finv(1-alpha,numdf,dendf,0);
power=1-probf(fc,numdf,dendf,ncparm);
run;
proc print;run;
Output:
·	Label
Num
DF
Den
DF
FValue
ProbF
alpha
ncparm
fc
power
Equal slopes
2
177
4.73 0.0100
0.05 9.46367 3.04701
0.78515 Introduction to Longitudinal Data Analysis
304
·	Hence, there is a power of 78.5% to detect the prespeciﬁed diﬀerences at the 5%
level of signiﬁcance.
Increasing the number of rats yields the following powers:
·	Group size
10
11
12
13
14
15
20
Power
78.5%
82.5%
85.9%
88.7%
91.0%
92.9%
97.9%
Introduction to Longitudinal Data Analysis
305
12.4 Example 2
·	We continue the previous random-intercepts model and study the eﬀect of varying
the variance components values
o	Results (10 rats per group):
d11
3.6 3.2
4.0 σ2
1.0 1.4
1.8 89.3%
88.5%
87.9%
79.8%
78.5%
77.4%
71.9%
70.3%
68.9%
Conclusions:
·	The power decreases as the total variance increases
·	Keeping the total variance constant, the power increases as the intraclass
correlation ρI = d11/(d11 + σ2) increases
Introduction to Longitudinal Data Analysis
306
12.5 Example 3
12.5.1 Introduction
Experiment for the comparison of two treatments A and B
A total of N general practitioners (GP’s) involved
Each GP treats n subjects
Yij is the response for subject j treated by GP i
The analysis should account for the variability between GP’s
o	Introduction to Longitudinal Data Analysis
307
·	We use the following random-intercepts model, where the random intercepts
reﬂect random GP eﬀects:
β1 + b1i + εij
if treatment A
β2 + b1i + εij
if treatment B
Yij =



Assumed true parameter values:
·	Eﬀect
Fixed eﬀects:
Parameter True value
Average treatment A
Average treatment B
β1
β2
Variance components:
var(b1i)
var(εij)
d11
σ2
d11 + σ2
1
2
?
?
4
Introduction to Longitudinal Data Analysis
308
·	Hence, the individual variance components are unknown. Only the total variability
is known to equal 4.
·	Power analyses will be performed for several values for the intraclass correlation
ρI = d11/(d11 + σ2)
Introduction to Longitudinal Data Analysis
309
12.5.2 Case 1: Treatments Assigned to GP’s
·	We now consider the situation in which the treatments will be randomly assigned
to GP’s, and all subjects with the same GP will be treated identically.
Powers for 2
·	×
25 = 50 GP’s, each treating 10 subjects (α = 0.05):
ρI
0.25 0.50
0.75 Power
86%
65%
50%
The power decreases as the intraclass correlation increases
·	Introduction to Longitudinal Data Analysis
310
12.5.3 Case 2: Treatments Assigned to Subjects
We now consider the situation in which the treatments will be randomly assigned
to subjects within GP’s, with the same number n/2 of subjects assigned to both
treatments
Powers for 2
×
5 = 10 subjects within 10 GP’s (α = 0.05):
ρI
0.25 0.50
0.75 Power
81%
94%
100%
The power increases as the intraclass correlation increases
Note also that Case 2 requires many less observations than Case 1
·	Introduction to Longitudinal Data Analysis
311
12.5.4 Conclusion
Within-‘subject’ correlation
increases power for inferences on within-‘subject’ eﬀects,
but decreases power for inferences on between-‘subject’ eﬀects
Introduction to Longitudinal Data Analysis
312
Part II
Marginal Models for Non-Gaussian Longitudinal Data
Introduction to Longitudinal Data Analysis
313
Chapter 13
The Toenail Data
·	Toenail Dermatophyte Onychomycosis: Common toenail infection, diﬃcult to
treat, aﬀecting more than 2% of population.
·	Classical treatments with antifungal compounds need to be administered until the
whole nail has grown out healthy.
o	New compounds have been developed which reduce treatment to 3 months
Randomized, double-blind, parallel group, multicenter study for the comparison of
two such new compounds (A and B) for oral treatment.
Introduction to Longitudinal Data Analysis
314
Research question:
Severity relative to treatment of TDO ?
189 patients randomized, 36 centers
2
×
48 weeks of total follow up (12 months)
12 weeks of treatment (3 months)
measurements at months 0, 1, 2, 3, 6, 9, 12.
o	Introduction to Longitudinal Data Analysis
315
Frequencies at each visit (both treatments):
·	Introduction to Longitudinal Data Analysis
316
Chapter 14
The Analgesic Trial
single-arm trial with 530 patients recruited (491 selected for analysis)
analgesic treatment for pain caused by chronic nonmalignant disease
treatment was to be administered for 12 months
we will focus on Global Satisfaction Assessment (GSA)
GSA scale goes from 1=very good to 5=very bad
GSA was rated by each subject 4 times during the trial, at months 3, 6, 9, and 12.
§	Introduction to Longitudinal Data Analysis
317
Research questions:
·	Evolution over time
·	Relation with baseline covariates: age, sex, duration of the pain, type of pain,
disease progression, Pain Control Assessment (PCA), . . .
·	Investigation of dropout
Frequencies:
o	GSA Month 3
Month 6
Month 9
Month 12
1
2
3
4
5
55 14.3% 38 12.6% 40 17.6% 30 13.5%
112 29.1% 84 27.8% 67 29.5% 66 29.6%
151 39.2% 115 38.1% 76 33.5% 97 43.5%
52 13.5% 51 16.9% 33 14.5% 27 12.1%
15
3.9% 14
4.6% 11
4.9%
3
1.4%
Tot
385
302
227
223
Introduction to Longitudinal Data Analysis
318
Missingness:
·	Measurement occasion
Month 3 Month 6 Month 9 Month 12 Number
%
Completers
O
Dropouts
O
M
M
O
M
M
M
Non-monotone missingness
M
O
O
M
O
O
M
M
O
O
M
O
O
M
O
M
O
O
O
M
O
M
M
M
O
O
O
O
163
41.2 51
51
63
30
7
2
18
2
1
1
3
12.91 12.91
15.95 7.59
1.77 0.51
4.56 0.51
0.25 0.25
0.76 O
O
O
O
O
O
O
O
M
M
M
M
Introduction to Longitudinal Data Analysis
319
Chapter 15
The National Toxicology Program (NTP) Data
Developmental Toxicity Studies
o	Research Triangle Institute
The eﬀect in mice of 3 chemicals:
·	DEHP: di(2-ethyhexyl)-phtalate
·	EG: ethylene glycol
·	DYME: diethylene glycol dimethyl ether
Introduction to Longitudinal Data Analysis
320
·	Implanted fetuses:
·	death/resorbed
·	viable:
∗
∗
weight
malformations: visceral,
skeletal, external
Data structure:
·	dam
(cid:0)
(cid:0)
@
@
(cid:0)
(cid:0)
@
@
(cid:0)
(cid:0)(cid:9)
?
@
@R
·	. .implant (mi). . .
(cid:0)
@
@
(cid:0)
(cid:0)
(cid:0)
(cid:0)
(cid:0)(cid:9)
@
@
@
@R
viable (ni) non-viable (ri)
(cid:1)
A
(cid:1)
(cid:1)
A
A
(cid:1)
(cid:1)
(cid:1)(cid:11)
A
A
AU
(cid:1)
A
(cid:1)
(cid:1)
A
A
(cid:1)
(cid:1)
(cid:1)(cid:11)
A
A
AU
malf. (zi)weight death resorption
(cid:1)
A
A
(cid:1)
A
A
(cid:1)
(cid:1)
(cid:1)
A
AU
(cid:1)(cid:11)
?
1 . . . K
Introduction to Longitudinal Data Analysis
321
# Dams,
1
≥
Exposure Dose Impl. Viab.
EG
DEHP
0
750
1500
3000
0
44
91
191
292
DYME
0
62.5 125
250
500
25
24
23
23
30
26
26
24
25
21
20
24
23
22
25
24
22
23
30
26
26
17
9
21
20
24
23
22
Litter
Size
Malformations
Live (mean) Ext. Visc. Skel.
297
276
229
226
330
288
277
137
50
282
225
290
261
141
11.9 11.5
10.4 9.8
13.2 11.1
10.7 8.1
5.6 13.4
11.3 12.1
11.3 0.0
1.1 1.7
7.1 0.0
1.0 5.4
0.0 0.0
0.9 4.0
1.5 0.4
7.2 0.3
8.7 36.7
55.8 1.2
0.4 4.3
17.5 15.3
18.3 54.0 50.0
48.0 0.0
0.0 1.0
2.7 0.0
0.0 0.0
0.1 0.0
0.0 1.0
20.0 6.1
66.0 19.9
79.4 Introduction to Longitudinal Data Analysis
322
Chapter 16
Generalized Linear Models
·	The model
·	Maximum likelihood estimation
·	Examples
·	McCullagh and Nelder (1989)
Introduction to Longitudinal Data Analysis
323
16.1 The Generalized Linear Model
Suppose a sample Y1, . . . , YN of independent observations is available
·	All Yi have densities f (yi|
·	θi, φ) which belong to the exponential family:
φ−
(cid:26)
1[yθi −
ψ(θi)] + c(y, φ)
(cid:27)
θi, φ) = exp
f (y
|
θi the natural parameter
Linear predictor: θi = xi0β
θ is the scale parameter (overdispersion parameter)
ψ(.) is a function to be discussed next
·	Introduction to Longitudinal Data Analysis
324
16.2 Mean and Variance
o	We start from the following general propterty:
f (y
|
Z
θ, φ)dy
=
Z
exp
φ−
1[yθ
(cid:26)
ψ(θ)] + c(y, φ)
(cid:27)
dy = 1
−
Taking ﬁrst and second-order derivatives with respect to θ yields
∂
∂θ Z
f (y
|
θ, φ) dy = 0
∂2
∂θ2
Z
f (y
|
θ, φ) dy = 0



Introduction to Longitudinal Data Analysis
325
[y
−
ψ0(θ)] f (y
|
θ, φ) dy = 0
[φ−
1(y
−
ψ0(θ))2
−
ψ00(θ)] f (y
|
θ, φ) dy = 0
Z
Z
E(Y ) = ψ0(θ)
Var(Y ) = φψ00(θ)
⇐⇒
⇐⇒






Note that, in general, the mean µ and the variance are related:
·	Var(Y ) = φψ00 "ψ0−
1(µ)# = φv(µ)
Introduction to Longitudinal Data Analysis
326
The function v(µ) is called the variance function.
The function ψ0−
1 which expresses θ as function of µ is called the link function.
ψ0 is the inverse link function
§	Introduction to Longitudinal Data Analysis
327
16.3 Examples
16.3.1 The Normal Model
Model:
Density function:
o	N (µ, σ2)
Y
∼
θ, φ) =
f (y
|
1
√2πσ2
exp 

−

1
σ2(y
−
µ)2



= exp 


1
σ2

yµ



−
µ2
2




·	



ln(2πσ2)
2
y2
2σ2
−







Introduction to Longitudinal Data Analysis
328
Exponential family:
·	θ = µ
·	φ = σ2
·	ψ(θ) = θ2/2
·	c(y, φ) = ln(2πφ)
2 −
y2
2φ
Mean and variance function:
·	µ = θ
·	v(µ) = 1
Note that, under this normal model, the mean and variance are not related:
φv(µ) = σ2
The link function is here the identity function: θ = µ
·	Introduction to Longitudinal Data Analysis
329
16.3.2 The Bernoulli Model
Model:
Density function:
o	Bernoulli(π)
Y
∼
θ, φ) = πy(1
f (y
|
y
π)1
−
−
= exp
y ln π + (1
{
−
y) ln(1
π)
}
−
= exp 


y ln 


π
−
1
π 


·	ln(1
−
π)


Introduction to Longitudinal Data Analysis
330
·	Exponential family:
·	θ = ln
·	φ = 1
π
π !
1
−
·	ψ(θ) = ln(1
·	c(y, φ) = 0
−
π) = ln(1 + exp(θ))
·	Mean and variance function:
·	µ = exp θ
1+exp θ = π
·	v(µ) = exp θ
(1+exp θ)2 = π(1
π)
−
Note that, under this model, the mean and variance are related:
φv(µ) = µ(1
µ)
−
The link function here is the logit link: θ = ln
µ
µ !
1
−
o	Introduction to Longitudinal Data Analysis
331
16.3.3 The Poisson Model
Model:
Density function:
o	Poisson(λ)
Y
∼
θ, φ) =
f (y
|
e−
λλy
y!
= exp
y ln λ
{
λ
ln y!
}
−
−
Introduction to Longitudinal Data Analysis
332
Exponential family:
·	θ = ln λ
·	φ = 1
·	ψ(θ) = λ = exp θ
·	c(y, φ) =
ln y!
−
Mean and variance function:
·	µ = exp θ = λ
·	v(µ) = exp θ = λ
Note that, under this model, the mean and variance are related:
φv(µ) = µ
The link function is here the log link: θ = ln µ
·	Introduction to Longitudinal Data Analysis
333
16.4 Generalized Linear Models (GLM)
Suppose a sample Y1, . . . , YN of independent observations is available
·	All Yi have densities f (yi|
·	θi, φ) which belong to the exponential family
·	In GLM’s, it is believed that the diﬀerences between the θi can be explained
through a linear function of known covariates:
θi = xi0β
xi is a vector of p known covariates
β is the corresponding vector of unknown regression parameters, to be estimated
from the data.
o	Introduction to Longitudinal Data Analysis
334
16.5 Maximum Likelihood Estimation
Log-likelihood:
`(β, φ) =
1
φ Xi
[yiθi −
ψ(θi)] +
c(yi, φ)
Xi
First order derivative with respect to β:
∂`(β, φ)
∂β
=
1
φ Xi
∂θi
∂β
[yi −
ψ0(θi)]
The score equations for β to be solved:
§	S(β) =
∂θi
∂β
Xi
[yi −
ψ0(θi)] = 0
Introduction to Longitudinal Data Analysis
335
Since µi = ψ0(θi) and vi = v(µi) = ψ00(θi), we have that
∂µi
β
= ψ00(θi)
∂θi
∂β
= vi
∂θi
∂β
The score equations now become
o	S(β) =
∂µi
∂β
Xi
1
v−
i
(yi −
µi) = 0
·	Note that the estimation of β depends on the density only through the means µi
and the variance functions vi = v(µi).
Introduction to Longitudinal Data Analysis
336
o	The score equations need to be solved numerically:
·	iterative (re-)weighted least squares
·	Newton-Raphson
·	Fisher scoring
Inference for β is based on classical maximum likelihood theory:
·	asymptotic Wald tests
·	likelihood ratio tests
·	score tests
Introduction to Longitudinal Data Analysis
337
·	In some cases, φ is a known constant, in other examples, estimation of φ may be
required to estimate the standard errors of the elements in β
Estimation can be based on Var(Yi) = φvi:
φ =
c
N
1
−
(yi −
p Xi
µi)2/vi(
µi)
d
d
For example, under the normal model, this would yield:
o	σ2 =
1
xi0
β)2,
(yi −
p Xi
the mean squared error used in linear regression models to estimate the residual
variance.
N
−
c
c
Introduction to Longitudinal Data Analysis
338
16.6 Illustration: The Analgesic Trial
o	Early dropout (did the subject drop out after the ﬁrst or the second visit) ?
Binary response
PROC GENMOD can ﬁt GLMs in general
PROC LOGISTIC can ﬁt models for binary (and ordered) responses
SAS code for logit link:
proc genmod data=earlydrp;
model earlydrp = pca0 weight psychiat physfct / dist=b;
run;
proc logistic data=earlydrp descending;
model earlydrp = pca0 weight psychiat physfct;
run;
Introduction to Longitudinal Data Analysis
339
o	SAS code for probit link:
proc genmod data=earlydrp;
model earlydrp = pca0 weight psychiat physfct / dist=b link=probit;
run;
proc logistic data=earlydrp descending;
model earlydrp = pca0 weight psychiat physfct / link=probit;
run;
Selected output:
Analysis Of Parameter Estimates
Parameter DF Estimate
Standard
Error
Wald 95%
Confidence Limits
Chi-
Square Pr > ChiSq
Intercept
PCA0
WEIGHT
PSYCHIAT
PHYSFCT
Scale
1
1
1
1
1
0
-1.0673
0.3981 -0.0211
0.7169 0.0121
1.0000 0.7328
0.1343 0.0072
0.2871 0.0050
0.0000 -2.5037
0.1349 -0.0353
0.1541 0.0024
1.0000 0.3690
0.6614 -0.0070
1.2796 0.0219
1.0000 2.12
8.79 8.55
6.23 5.97
0.1453 0.0030
0.0034 0.0125
0.0145 NOTE: The scale parameter was held fixed.
Introduction to Longitudinal Data Analysis
340
Chapter 17
Parametric Modeling Families
·	Continuous outcomes
·	Longitudinal generalized linear models
·	Notation
Introduction to Longitudinal Data Analysis
341
17.1 Continuous Outcomes
Marginal Models:
Random-Eﬀects Models:
E(Yij|
xij) = x0ijβ
E(Yij|
bi, xij) = x0ijβ + z0ijbi
Transition Models:
§	E(Yij|
Yi,j
1, . . . , Yi1, xij) = x0ijβ + αYi,j
−
1
−
Introduction to Longitudinal Data Analysis
342
17.2 Longitudinal Generalized Linear Models
§	Normal case: easy transfer between models
Also non-normal data can be measured repeatedly (over time)
Lack of key distribution such as the normal [=
]
⇒
·	A lot of modeling options
·	Introduction of non-linearity
·	No easy transfer between model families
cross-
sectional
longitudinal
normal outcome
linear model
LMM
non-normal outcome
GLM
?
Introduction to Longitudinal Data Analysis
343
17.3 Notation
Let the outcomes for subject i = 1, . . . , N be denoted as (Yi1, . . . , Yini).
Group into a vector Y i:
·	Binary data: each component is either 0 or 1.
·	(Binary data: each component is either 1 or 2.)
·	(Binary data: each component is either
·	(Categorical data: Yij ∈ {
The corresponding covariate vector is xij.
.)
1, . . . , c
}
1 or +1.)
−
It is convenient to use (binary 0/1 data):
·	E(Yij) = Pr(Yij = 1) = µij
and
µijk = E(YijYik) = Pr(Yij = 1, Yik = 1)
Introduction to Longitudinal Data Analysis
344
Chapter 18
Conditional Models
·	A log-linear model
·	Quadratic version of the model
·	Linear version of the model
·	Clustered-data versions of the model
·	Transition models
Introduction to Longitudinal Data Analysis
345
18.1 A Log-linear Model
·	Cox (1972)
Joint distribution of Y i in terms of a multivariate exponential family:
f (yi, θi) = exp 


n
Xj=1
θijyij +
θij1j2yij1yij2 + . . . + θi1...nyi1 . . . yin −
Xj1<j2
A(θi)


= c(θi) exp 


n
Xj=1
θijyij +
Xj1<j2
θij1j2yij1yij2 + . . . + θi1...nyi1 . . . yin


A(θi) [equivalently, c(θi)] is the normalizing constant
θi is the canonical parameter, consisting of ﬁrst, second, up to nth order
components.
Introduction to Longitudinal Data Analysis
346
Interpretation of Parameters:
o	The parameters have a conditional interpretation:
θij = ln 



Pr(Yij = 1
|
Pr(Yij = 0
|
Yik = 0; k
Yik = 0; k
= j)
= j)




·	the ﬁrst order parameters (main eﬀects) are interpreted as conditional
⇒
logits.
·	Similarly,
θijk = ln 


Pr(Yij = 1, Yik = 1
Pr(Yij = 1, Yik = 0
|
|
Yi = 0; k, j Yi = 0; k, j
= )Pr(Yij = 0, Yik = 0 = )Pr(Yij = 0, Yik = 1
Yi = 0; k, j Yi = 0; k, j
|
|
= ) = ) 


·	These are conditional log odds ratios.
Introduction to Longitudinal Data Analysis
347
6
6
6
6
6
6
o	Advantages:
·	The parameter vector is not constrained. All values of θ
nonnegative probabilities.
IR yield
∈
·	Calculation of the joint probabilities is fairly straightforward:
ignore the normalizing constant
evaluate the density for all possible sequences y
1
sum all terms to yield c(θ)−
∗
∗
∗
Drawbacks:
·	Due to above conditional interpretation, the models are less useful for
regression.
The dependence of E(Yij) on covariates involves all parameters, not only the main eﬀects.
·	The interpretation of the parameters depends on the length ni of a sequence.
These drawbacks make marginal models or models that combine marginal and
conditional features better suited.
Introduction to Longitudinal Data Analysis
348
18.2 Quadratic and Linear Versions
·	Cox (1972) and others suggest that often the higher order interactions can be
neglected. This claim is supported by empirical evidence.
The quadratic exponential model:
f (yi, θi) = exp 


n
Xj=1
θijyij +
θij1j2yij1yij2 −
Xj1<j2
A(θi)


= c(θi) exp 
n
θijyij +
Xj1<j2
θij1j2yij1yij2 
·	

Xj=1


The linear exponential model:
f (yi, θi) = exp 
n
Xj=1


θijyij −
A(θi)


then this model reﬂects the assumption of independence.
The linear model equals logistic regression.
§	Introduction to Longitudinal Data Analysis
349
18.3 A Version for Clustered Binary Data
o	NTP data: Yij is malformation indicator for fetus j in litter i
Code Yij as
1 or 1
−
di is dose level at which litter i is exposed
Simpliﬁcation:
θij = θi = β0 + βddi
and
θij1j2 = βa
Using
we obtain
f (zi|
Zi =
ni
Xj=1
Yij
exp
θizi + βazi(ni −
{
zi)
A(θi)
}
−





θi, βa) = 
ni
zi




Introduction to Longitudinal Data Analysis
350
18.4 Transition Models
Molenberghs and Verbeke (2005, Section 11.5)
Outcome Yij or error term εij is a function of history hij = (Yi1, . . . , Yi,j
−
Order of transition model: # of previous outcomes in regression
Stationary model: functional form of dependence independent of occurrence time
A stationary ﬁrst-order autoregressive model for continuous data is:
o	Yi1 = x0i1β + εi1
Yij = x0ijβ + αYi,j
1 + εij
−
Introduction to Longitudinal Data Analysis
351
Assume
·	then
=
⇒
N (0, σ2)
and
εi1 ∼
εij ∼
j
j0−
N (0, σ2(1
α2))
−
|σ2
cov(Yij, Yij0) = α|
a marginal multivariate normal model with AR(1) variance-covariance matrix.
·	For non-Gaussian outcomes, ﬁrst write
Yij = µc
ij + εc
ij
and then
µc
ij = E(Yij|
ij) = var(Yij|
hij)
hij)
φvc(µc
Introduction to Longitudinal Data Analysis
352
Example of a linear predictor:
ηij(µc
ij) = x0ijβ + κ(hij, β, α)
κ is a function of the history.
This model is easy to ﬁt since it leads to independent GLM contributions:
§	f (yi1, . . . , yini) = f (yi1)
= f (yi1)
·
·
f (yi2|
yi1)
ni
Yj=2
f (yij|
yi1, yi2)
f (yini|
f (yi3|
·
hij) = f (yi1, . . . , yiq)
·
yi1, . . . , yi,ni−
ni
hij)
f (yij|
Yj=q+1
·
This product yields ni −
·	q independent univariate GLM contributions.
A separate model may need to be considered for the ﬁrst q measurements.
·	Introduction to Longitudinal Data Analysis
353
§	A logistic-regression type example:
logit[P (Yij = 1
|
xij, Yi,j
1 = yi,j
−
1, β, α)] = x0ijβ + αyi,j
−
−
The marginal means and variances do not follow easily, except in the normal case.
Recursive formulas are:
µij = µc
ij(0)[1
vij = [µc
ij(1)
µi,j
1] + µc
−
−
µc
ij(0)]2vi,j
ij(1)µi,j
1 + vc
−
1
−
ij(0)[1
−
µi,j
1] + vc
−
ij(1)µi,j
1
−
−
Introduction to Longitudinal Data Analysis
354
18.4.1 Analysis of the Toenail Data
§	Formulate a transition model (Model I):
Yij ∼
Bernoulli(µij)
logit 



µij
µij
1
−




= β0 + β1Ti + β2tij + β3Titij + α1yi,j
1
−
To account for unequal spacing (Model II):
·	α1 describes the transition eﬀect for the later measurements
·	α1a is the ‘excess’ during the ﬁrst quarter
·	hence: autoregressive eﬀect at months 1, 2, and 3 is α1 + α1a
Alternatively: dependence on two prior occasions:
logit 



µij
µij
1
−




= β0 + β1Ti + β2tij + β3Titij + α1yi,j
1 + α2yi,j
−
2
−
Introduction to Longitudinal Data Analysis
355
Fitted models:
·	First order
Eﬀect
Intercept
Ti
tij
Ti ·
Dep. on Yi,j
tij
Dep. on Yi,j
Dep. on Yi,j
1
−
1
−
2
−
Par.
I
II
Second order
β0
β1
β2
β3
α1
α1a
α2
-3.14 (0.27)
-3.77 (0.34)
-3.28 (0.34)
0.00 (0.31)
-0.08 (0.32)
0.13 (0.39)
-0.09 (0.04)
0.03 (0.05)
-0.05 (0.04)
-0.08 (0.06)
-0.06 (0.06)
-0.09 (0.07)
4.48 (0.22)
3.59 (0.29)
4.01 (0.39)
1.56 (0.35)
0.25 (0.38)
Introduction to Longitudinal Data Analysis
356
Two separate models, depending on the level of the previous outcome:
logit 



µij
µij
1
−




= (β00 + β10Ti + β20tij + β30Titij)IYi,j
1=0
−
+(β01 + β11Ti + β21tij + β31Titij)IYi,j
1=1
−
Fitted model:
o	Yi,j
1 = 0
−
Yi,j
1 = 1
−
Eﬀect
Par. Estimate (s.e.) Par. Estimate (s.e.)
Intercept
Ti
tij
Ti ·
tij
β00
β10
β20
β30
-3.92 (0.56)
0.45 (0.70)
-0.06 (0.09)
0.07 (0.10)
β01
β11
β21
β31
1.56 (1.26)
-0.01 (0.37)
-0.20 (0.06)
0.04 (0.07)
Introduction to Longitudinal Data Analysis
357
18.4.2 Transition Model in SAS
Prepare models so that the previous outcome can be used as a covariate (using
the same code as used to ﬁt a model for dropout – see Part V)
·	%dropout(data=test,id=idnum,time=time,response=onyresp,out=test2);
data test2a;
set test2;
prev1=prev;
drop prev;
run;
%dropout(data=test2a,id=idnum,time=time,response=prev1,out=test3);
data test3a;
set test3;
prev2=prev;
drop prev;
run;
Introduction to Longitudinal Data Analysis
358
o	The result for the ﬁrst subject is
Obs
idnum
time
treatn onyresp prev1 prev2
1
2
3
4
5
6
7
1
1
1
1
1
1
1
0
1
2
3
6
9
12
1
1
1
1
1
1
1
1
1
1
0
0
0
0
·	1
1
1
0
0
0
o	1
1
1
0
0
Code to ﬁt a transition model:
proc genmod data=test3a descending;
model onyresp = treatn time treatn*time prev1 / dist=binomial;
run;
Introduction to Longitudinal Data Analysis
359
When both predecessors are used, one merely adds ‘prev2’ to MODEL statement:
model onyresp = prev1 treatnprev1 timeprev1 treatntimeprev1
/ noint dist=binomial;
To ﬁt Model II, an additional variable ‘prev1a’ needs to be created:
o	data test3b;
set test3a;
prev1a=prev1;
if time>3 then prev1a=0;
run;
which is then added to the logistic regression, next to ‘prev1.’
Introduction to Longitudinal Data Analysis
360
Chapter 19
Full Marginal Models
·	Introduction
·	Link functions
·	Associations
·	Bahadur model
·	Multivariate Probit model
·	Example: POPS data
Introduction to Longitudinal Data Analysis
361
19.1 Introduction
o	Choices to make:
·	Description of mean proﬁles (univariate parameters) and of association
(bivariate and higher order parameters)
·	Degree of modeling:
joint distribution fully speciﬁed
⇒
only a limited number of moments
∗
∗
likelihood procedures
e.g., generalized estimating equations
⇒
Minimally, one speciﬁes:
·	ηi(µi) =
ηi1(µi1), . . . , ηin(µin)
{
}
·	E(Y i) = µi
·	var(Y i) = φv(µi) where v(.) is a known variance function
·	corr(Y i) = R(α)
ηi(µi) = X iβ
and
Introduction to Longitudinal Data Analysis
362
19.2 Univariate Link Functions
The marginal logit link:
ηij = ln(µij)
ln(1
−
−
µij) = logit(µij).
The probit link:
1
1 (µij).
ηij = Φ−
The complementary log-log link
·	· · ·
Introduction to Longitudinal Data Analysis
363
19.3 Pairwise Association
·	Success probability approach. (Ekholm 1991)
Logit link for two-way probabilities
ηijk = ln(µijk)
ln(1
−
−
µijk) = logit(µijk),
Marginal correlation coeﬃcient. (Bahadur model)
·	ρijk =
r
µijµik
µij)µik(1
µijk −
−
−
µij(1
µik)
ηijk = ln(1 + ρijk)
ln(1
−
−
ρijk)
(Fisher’s z transform)
Introduction to Longitudinal Data Analysis
364
Marginal odds ratio. (Dale model)
ψijk =
µik + µijk)
µijk)
µij −
(µijk)(1
−
µijk)(µij −
(µik −
Pr(Yij = 1, Yik = 1)Pr(Yij = 0, Yik = 0)
Pr(Yij = 0, Yik = 1)Pr(Yij = 1, Yik = 0)




= 



ηijk = ln(ψijk)
(log odds ratio)
Higher order association deﬁned similarly
Calculations can become combersome
§	Introduction to Longitudinal Data Analysis
365
19.4 The Bahadur Model
Univariate: E(Yij) = P (Yij = 1)
πij.
≡
Bivariate: E(YijYik) = P (Yij = 1, Yik = 1)
πijk.
≡
Correlation structure:
Corr(Yij, Yik)
ρijk =
≡
[πij(1
πijπik
πijk −
πij)πik(1
−
−
πik)]1/2.
This yields expression for pairwise probabilities:
πijk = πijπik + ρijk[πij(1
πij)πik(1
−
πik)]1/2.
−
Similarly for the full joint distribution f (y).
o	Introduction to Longitudinal Data Analysis
366
Let
·	and
εij =
Yij −
πij(1
−
πij
πij)
r
and
eij =
yij −
πij(1
−
πij
πij)
r
,
ρijk = E(εijεik),
ρijk = E(εijεikεi),
...
ρi12...ni = E(εi1εi2 . . . εini).
Introduction to Longitudinal Data Analysis
367
A general expression:
·	f (yi) = f1(yi)c(yi),
with
and
f1(yi) =
ni
Yj=1
yij
ij (1
π
−
πij)1
−
yij
c(yi) = 1 +
ρijkeijeik +
Xj<k
Xj<k<`
ρijkeijeikei + . . . + ρi12...niei1ei2 . . . eini.
Introduction to Longitudinal Data Analysis
368
19.5 The Multivariate Probit Model
·	3 categorical out-
E.g., 4
come arises from underlying
bivariate normal
×
·	Covariate eﬀects
cut oﬀ points
≡
shift of
·	Correlation = polychoric cor-
relation: allowed to depend on
covariates
Introduction to Longitudinal Data Analysis
369
19.6 The POPS Data
Project On Preterm and Small for Gestational Age Infants
1530 Dutch children (1983)
Collected data:
§	Perinatal information:
Ability scores at the age of 2:
·	Bilirubin value
·	Are the child’s movements natural ?
·	Neonatal seizures
·	Can the child pile three bricks ?
·	Congenital malformations
·	Can the child put a ball
in a boxed
when asked to ?
Introduction to Longitudinal Data Analysis
370
19.7 Application to POPS Data
Bahad
Probit
Dale-Norm
Dale-Logist
First Ability Score
Intercept
3.67(0.49)
2.01(0.26)
2.03(0.27)
3.68(0.52)
Neonatal seiz.
-1.94(0.42)
-1.12(0.26)
-1.16(0.26)
-2.06(0.44)
Congenital malf.
-1.21(0.31)
-0.61(0.18)
-0.62(0.18)
-1.17(0.33)
100
×
Bilirubin
-0.69(0.25)
-0.32(0.14)
-0.32(0.14)
-0.64(0.27)
Second Ability Score
Intercept
4.03(0.51)
2.19(0.27)
2.21(0.27)
4.01(0.54)
Neonatal seiz.
-2.26(0.43)
-1.27(0.26)
-1.29(0.26)
-2.28(0.44)
Congenital malf.
-1.08(0.32)
-0.56(0.19)
-0.59(0.19)
-1.11(0.34)
100
×
Bilirubin
-0.85(0.26)
-0.42(0.14)
-0.41(0.14)
-0.80(0.27)
Third Ability Score
Intercept
3.32(0.50)
1.84(0.27)
1.91(0.27)
3.49(0.54)
Neonatal seiz.
-1.55(0.44)
-0.88(0.27)
-0.93(0.27)
-1.70(0.46)
Congenital malf.
-0.96(0.32)
-0.47(0.19)
-0.49(0.19)
-0.96(0.35)
100
×
Bilirubin
-0.44(0.26)
-0.21(0.14)
-0.24(0.14)
-0.49(0.28)
Introduction to Longitudinal Data Analysis
371
Bahad
Probit
Dale-Norm
Dale-Logist
Association parameters
ρ
ρ
ψ
ψ
(1,2): ρ or ψ
0.27(0.05)
0.73(0.05)
17.37(5.19)
17.35(5.19)
(1,2): z(ρ) or ln ψ
0.55(0.11)
1.85(0.23)
2.85(0.30)
2.85(0.30)
(1,3): ρ or ψ
0.39(0.05)
0.81(0.04)
30.64(9.78)
30.61(9.78)
(1,3): z(ρ) or ln ψ
0.83(0.12)
2.27(0.25)
3.42(0.32)
3.42(0.32)
(2,3): ρ or ψ
0.23(0.05)
0.72(0.05)
17.70(5.47)
17.65(5.47)
(2,3): z(ρ) or ln ψ
0.47(0.10)
1.83(0.23)
2.87(0.31)
2.87(0.31)
(1,2,3): ρ or ψ
(1,2,3): z(ρ) or ln ψ
—
—
—
—
0.91(0.69)
0.92(0.69)
-0.09(0.76)
-0.09(0.76)
Log-likelihood
-598.44
-570.69
-567.11
-567.09
Introduction to Longitudinal Data Analysis
372
Chapter 20
Generalized Estimating Equations
·	General idea
·	Asymptotic properties
·	Working correlation
·	Special case and application
·	SAS code and output
Introduction to Longitudinal Data Analysis
373
20.1 General Idea
o	Univariate GLM, score function of the form (scalar Yi):
S(β) =
N
Xi=1
∂µi
∂β
1
v−
i
(yi −
µi) = 0
with
vi = Var(Yi)
In longitudinal setting: Y = (Y 1, . . . , Y N):
∂µij
∂β
(yij −
S(β) =
µij) =
1
v−
ij
Xi Xj
N
Xi=1
D0i [Vi(α)]−
1 (yi −
µi) = 0
where
p matrix with (i, j)th elements
∂µij
·	Di is an ni ×
∂β
·	yi and µi are ni-vectors with elements yij and µij
ni diagonal or more complex?
·	Is Vi ni ×
Introduction to Longitudinal Data Analysis
374
·	Vi = Var(Y i) is more complex since it involves a set of nuisance parameters α,
determining the covariance structure of Y i:
Vi(β, α) = φA1/2
i
(β)Ri(α)A1/2
i
(β)
in which

r
A1/2
i
(β) =
vi1(µi1(β)) . . .
·	. .
...
0
...
vini(µini(β))
and Ri(α) is the correlation matrix of Yi, parameterized by α.
·	. .
0
r



























·	Same form as for full likelihood procedure, but we restrict speciﬁcation to the ﬁrst
moment only
Liang and Zeger (1986)
·	Introduction to Longitudinal Data Analysis
375
20.2 Large Sample Properties
As N
→ ∞
where
√N( ˆβ
β)
∼
N (0, I −
1
0 )
D0i[Vi(α)]−
1Di
−
N
Xi=1
I0 =
(Unrealistic) Conditions:
·	α is known
·	the parametric form for Vi(α) is known
This is the naive
≡
purely model based variance estimator
Solution: working correlation matrix
§	Introduction to Longitudinal Data Analysis
376
20.3 Unknown Covariance Structure
Keep the score equations
S(β) =
N
Xi=1
[Di]0 [Vi(α)]−
1 (yi −
µi) = 0
BUT
·	suppose Vi(.) is not the true variance of Y i but only a plausible guess, a so-called
working correlation matrix
·	specify correlations and not covariances, because the variances follow from the
mean structure
the score equations are solved as before
·	Introduction to Longitudinal Data Analysis
377
The asymptotic normality results change to
√N( ˆβ
β)
−
∼
N (0, I −
1
1
0 )
0 I1I −
I0 =
I1 =
N
Xi=1
N
Xi=1
D0i[Vi(α)]−
1Di
D0i[Vi(α)]−
1Var(Y i)[Vi(α)]−
1Di.
This is the robust
·	I0 is the bread
≡
empirically corrected
≡
sandwich variance estimator
·	I1 is the ﬁlling (ham or cheese)
Correct guess =
⇒
likelihood variance
§	Introduction to Longitudinal Data Analysis
378
The estimators ˆβ are consistent even if the working correlation matrix is incorrect
An estimate is found by replacing the unknown variance matrix Var(Y i) by
o	(Y i −
ˆµi)(Y i −
ˆµi)0.
Even if this estimator is bad for Var(Y i) it leads to a good estimate of I1,
provided that:
o	replication in the data is suﬃciently large
·	same model for µi is ﬁtted to groups of subjects
·	observation times do not vary too much between subjects
A bad choice of working correlation matrix can aﬀect the eﬃciency of ˆβ
Care needed with incomplete data (see Part V)
o	Introduction to Longitudinal Data Analysis
379
20.4 The Working Correlation Matrix
Vi = Vi(β, α, φ) = φA1/2
i
(β)Ri(α)A1/2
i
(β)
·	Variance function: Ai is (ni ×
GLM variance function.
ni) diagonal with elements v(µij), the known
Working correlation: Ri(α) possibly depends on a diﬀerent set of parameters
α.
§	Overdispersion parameter: φ, assumed 1 or estimated from the data.
The unknown quantities are expressed in terms of the Pearson residuals
Note that eij depends on β.
eij =
µij
yij −
v(µij)
r
·	Introduction to Longitudinal Data Analysis
380
20.5 Estimation of Working Correlation
Liang and Zeger (1986) proposed moment-based estimates for the working correlation.
Corr(Yij, Yik)
Independence
Exchangeable
AR(1)
0
α
α
Estimate
—
ˆα = 1
N
N
i=1
1
ni(ni
P
−
j
P
=k eijeik
ˆα = 1
N
N
i=1
1
ni
−
1
P
ni
j
≤
−
P
1 eijei,j+1
Dispersion parameter:
ˆφ =
1
N
N
Xi=1
1
ni
ni
Xj=1
e2
ij.
Unstructured
αjk
ˆαjk = 1
N
N
i=1 eijeik
P
Introduction to Longitudinal Data Analysis
381
6
20.6 Fitting GEE
The standard procedure, implemented in the SAS procedure GENMOD.
1.	Compute initial estimates for β, using a univariate GLM (i.e., assuming
independence).
2.	. Compute Pearson residuals eij.
·	Compute estimates for α and φ.
·	Compute Ri(α) and Vi(β, α) = φA1/2
i
(β)Ri(α)A1/2
i
(β).
3.	Update estimate for β:
β(t+1) = β(t)
N
−



Xi=1
1
i Di
D0iV −


1
−
N



Xi=1
1
D0iV −
i
(yi −
µi)
·	

4.	Iterate 2.–3. until convergence.
1
Estimates of precision by means of I −
0
1
and/or I −
1
0 .
0 I1I −
Introduction to Longitudinal Data Analysis
382
20.7 Special Case: Linear Mixed Models
§	Estimate for β:
β(α) = 
N
−
X 0iWiXi
c
Xi=1


1 N
Xi=1


with α replaced by its ML or REML estimate
X 0iWiYi
Conditional on α,
β has mean
c
β(α)
E
(cid:20)
(cid:21)
c
N
Xi=1
= 


−
X 0iWiXi
1 N
Xi=1


X 0iWiXiβ = β
provided that E(Yi) = Xiβ
Hence, in order for
is correctly speciﬁed.
c
β to be unbiased, it is suﬃcient that the mean of the response
Introduction to Longitudinal Data Analysis
383
·	Conditional on α,
β has covariance
c
−
X 0iWiXi

1
N


Xi=1
X 0iWiVar(Yi)WiXi



N
Xi=1
1
−
X 0iWiXi

N
= 
Xi=1

1
−
X 0iWiXi

N
Xi=1

Var(
β) = 
c
·	Note that this model-based version assumes that the covariance matrix
Var(Yi) is correctly modelled as Vi = ZiDZ 0i + Σi.
An empirically corrected version is:
·	Var(
β) = 
c
N
Xi=1


X 0iWiXi


1
−
N



Xi=1
X 0iWiVar(Yi)WiXi


N



Xi=1
1
−
X 0iWiXi


|
{z
}
|
↓
BREAD
{z
↓
MEAT
}
|
{z
}
↓
BREAD
Introduction to Longitudinal Data Analysis
384
Chapter 21
A Family of GEE Methods
·	Classical approach
·	Prentice’s two sets of GEE
·	Linearization-based version
·	GEE2
·	Alternating logistic regressions
Introduction to Longitudinal Data Analysis
385
21.1 Prentice’s GEE
N
Xi=1
1
D0iV −
i
(Y i −
µi) = 0,
N
Xi=1
1
E0iW −
i
(Z i −
δi) = 0
where
Zijk =
r
µij)(Yik −
(Yij −
µij)µik(1
µij(1
−
−
µik)
µik)
,
δijk = E(Zijk)
The joint asymptotic distribution of √N ( ˆβ
variance-covariance matrix consistently estimated by
−
β) and √N ( ˆα
α) normal with
−
A 0
B C
N
























Λ11 Λ12
Λ21 Λ22
A B0
0 C
























Introduction to Longitudinal Data Analysis
386
where
N
Xi=1
N
Xi=1
N
Xi=1
A = 

B = 

C = 

and
1
−
1
i Di
D0iV −

1
−
1
−
1
i Ei
E0iW −

1
i Ei
E0iW −

Λ11 =
Λ12 =
N
Xi=1
N
Xi=1
1
D0iV −
1
i Cov(Y i)V −
i Di,
1
D0iV −
1
i Cov(Y i, Zi)W −
i Ei,
Λ21 = Λ12,
Λ22 =
N
Xi=1
1
E0iW −
1
i Cov(Zi)W −
i Ei,
,
N


Xi=1
1
E0iW −
i
∂Z i
∂β 

N


Xi=1
1
i Di
D0iV −

1
−
,
,
Statistic
Var(Y i)
(Y i −
Cov(Y i, Zi) (Y i −
(Z i −
Var(Z i)
Estimator
µi)(Y i −
µi)(Z i −
δi)(Z i −
µi)0
δi)0
δi)0
Introduction to Longitudinal Data Analysis
387
21.2 GEE Based on Linearization
21.2.1 Model formulation
Previous version of GEE are formulated directly in terms of binary outcomes
This approach is based on a linearization:
yi = µi + εi
with
ηi = g(µi),
ηi = X iβ,
Var(yi) = Var(εi) = Σi.
ηi is a vector of linear predictors,
g(.) is the (vector) link function.
·	Introduction to Longitudinal Data Analysis
388
21.2.2 Estimation (Nelder and Wedderburn 1972)
o	Solve iteratively:
where
N
Xi=1
X 0iWiXiβ =
N
Xi=1
XiWiy∗i ,
1
i Fi,
Wi = F 0i Σ−
y∗i = ˆηi + (yi −
1
ˆµi)F −
i
,
Fi =
∂µi
∂ηi
,
Σi = Var(ε),
µi = E(yi).
Remarks:
·	y∗i is called ‘working variable’ or ‘pseudo data’.
·	Basis for SAS macro and procedure GLIMMIX
·	For linear models, Di = Ini and standard linear regression follows.
Introduction to Longitudinal Data Analysis
389
21.2.3 The Variance Structure
Σi = φA1/2
i
(β)Ri(α)A1/2
i
(β)
§	φ is a scale (overdispersion) parameter,
Ai = v(µi), expressing the mean-variance relation (this is a function of β),
Ri(α) describes the correlation structure:
·	If independence is assumed then Ri(α) = Ini.
·	Other structures, such as compound symmetry, AR(1),. . . can be assumed as
well.
Introduction to Longitudinal Data Analysis
390
21.3 GEE2
§	Model:
·	Marginal mean structure
·	Pairwise association:
Odds ratios
Correlations
∗
∗
Working assumptions: Third and fourth moments
Estimation:
·	Second-order estimating equations
·	Likelihood (assuming 3rd and 4th moments are correctly speciﬁed)
Introduction to Longitudinal Data Analysis
391
21.4 Alternating Logistic Regression
Diggle, Heagerty, Liang, and Zeger (2002) and Molenberghs and Verbeke (2005)
When marginal odds ratios are used to model association, α can be estimated
using ALR, which is
·	almost as eﬃcient as GEE2
·	almost as easy (computationally) than GEE1
µijk as before and αijk = ln(ψijk) the marginal log odds ratio:
§	logit Pr(Yij = 1
|
xij) = xijβ
logit Pr(Yij = 1
|
Yik = yik) = αijkyik + ln 



µij −
µij −
µijk
µik + µijk




1
−
Introduction to Longitudinal Data Analysis
392
αijk can be modelled in terms of predictors
the second term is treated as an oﬀset
the estimating equations for β and α are solved in turn, and the ‘alternating’
between both sets is repeated until convergence.
this is needed because the oﬀset clearly depends on β.
·	Introduction to Longitudinal Data Analysis
393
21.5 Application to the Toenail Data
21.5.1 The model
Consider the model:
Yij ∼
Bernoulli(µij),
log 



µij
µij
1
−




= β0 + β1Ti + β2tij + β3Titij
Yij: severe infection (yes/no) at occasion j for patient i
tij: measurement time for occasion j
Ti: treatment group
·	Introduction to Longitudinal Data Analysis
394
21.5.2 Standard GEE
o	SAS Code:
proc genmod data=test descending;
class idnum timeclss;
model onyresp = treatn time treatn*time
/ dist=binomial;
repeated subject=idnum / withinsubject=timeclss
type=exch covb corrw modelse;
run;
SAS statements:
·	The REPEATED statements deﬁnes the GEE character of the model.
·	‘type=’: working correlation speciﬁcation (UN, AR(1), EXCH, IND,. . . )
·	‘modelse’: model-based s.e.’s on top of default empirically corrected s.e.’s
·	‘corrw’: printout of working correlation matrix
·	‘withinsubject=’: speciﬁcation of the ordering within subjects
Introduction to Longitudinal Data Analysis
395
·	Selected output:
·	Regression parameters:
Analysis Of Initial Parameter Estimates
Parameter
DF Estimate
Standard
Error
Wald 95%
Confidence Limits
Chi-
Square
Intercept
treatn
time
treatn*time
Scale
1
1
1
1
0
-0.5571
0.0240 -0.1769
-0.0783
1.0000 0.1090
0.1565 0.0246
0.0394 0.0000
-0.7708
-0.2827
-0.2251
-0.1556
1.0000 -0.3433
0.3307 -0.1288
-0.0010
1.0000 26.10
0.02 51.91
3.95 . Estimates from ﬁtting the model, ignoring the correlation structure, i.e.,
from ﬁtting a classical GLM to the data, using proc GENMOD.
·	The reported log-likelihood also corresponds to this model, and therefore
should not be interpreted.
·	The reported estimates are used as starting values in the iterative estimation
procedure for ﬁtting the GEE’s.
Introduction to Longitudinal Data Analysis
396
Analysis Of GEE Parameter Estimates
Empirical Standard Error Estimates
Parameter
Estimate
Standard
Error
95% Confidence
Limits
Z Pr > |Z|
-0.5840
Intercept
0.0120 treatn
time
-0.1770
treatn*time -0.0886
0.1734 -0.9238 -0.2441
0.5241 0.2613 -0.5001
0.0311 -0.2380 -0.1161
0.0233 0.0571 -0.2006
-3.37
0.05 -5.69
-1.55
0.0008 0.9633
<.0001
0.1208 Analysis Of GEE Parameter Estimates
Model-Based Standard Error Estimates
Parameter
Estimate
Standard
Error
95% Confidence
Limits
Z Pr > |Z|
-0.5840
Intercept
0.0120 treatn
time
-0.1770
treatn*time -0.0886
0.1344 -0.8475 -0.3204
0.1866 -0.3537
0.3777 0.0209 -0.2180 -0.1361
0.0362 -0.1596 -0.0177
-4.34
0.06 -8.47
-2.45
<.0001
0.9486 <.0001
0.0143 . The working correlation:
Exchangeable Working Correlation
Correlation
0.420259237 Introduction to Longitudinal Data Analysis
397
21.5.3 Alternating Logistic Regression
‘type=exch’
−→
‘logor=exch’
Note that α now is a genuine parameter
o	Introduction to Longitudinal Data Analysis
398
Selected output:
·	Analysis Of GEE Parameter Estimates
Empirical Standard Error Estimates
Parameter
Estimate
Standard
Error
95% Confidence
Limits
Z Pr > |Z|
-0.5244
Intercept
0.0168 treatn
time
-0.1781
treatn*time -0.0837
3.2218 Alpha1
0.1686 -0.8548 -0.1940
0.2432 -0.4599
0.4935 0.0296 -0.2361 -0.1200
0.0182 0.0520 -0.1856
3.7917 2.6519
0.2908 -3.11
0.07 -6.01
-1.61
11.08 0.0019
0.9448 <.0001
0.1076 <.0001
Analysis Of GEE Parameter Estimates
Model-Based Standard Error Estimates
Parameter
Estimate
Standard
Error
95% Confidence
Limits
Z Pr > |Z|
-0.5244
Intercept
0.0168 treatn
-0.1781
time
treatn*time -0.0837
0.1567 -0.8315 -0.2173
0.2220 -0.4182
0.4519 0.0233 -0.2238 -0.1323
0.0392 -0.1606 -0.0068
-3.35
0.08 -7.63
-2.13
0.0008 0.9395
<.0001
0.0329 Introduction to Longitudinal Data Analysis
399
21.5.4 Linearization Based Method
GLIMMIX macro:
%glimmix(data=test, procopt=%str(method=ml empirical),
stmts=%str(
class idnum timeclss;
model onyresp = treatn time treatn*time / solution;
repeated timeclss / subject=idnum type=cs rcorr;
),
error=binomial,
link=logit);
GLIMMIX procedure:
o	proc glimmix data=test method=RSPL empirical;
class idnum;
model onyresp (event=’1’) = treatn time treatn*time
/ dist=binary solution;
random residual / subject=idnum type=cs;
run;
Introduction to Longitudinal Data Analysis
400
Both produce the same results
The GLIMMIX macro is a MIXED core, with GLM-type surrounding statements
The GLIMMIX procedure does not call MIXED, it has its own engine
PROC GLIMMIX combines elements of MIXED and of GENMOD
RANDOM residual
is the PROC GLIMMIX way to specify residual correlation
o	Introduction to Longitudinal Data Analysis
401
21.5.5 Results of Models Fitted to Toenail Data
Eﬀect Par.
IND
EXCH
GEE1
UN
Int.
Ti
tij
Ti
tij
·
Int.
Ti
tij
tij
Ti
·
Ass.
Int.
Ti
tij
Ti
tij
·
β0
β1
β2
β3
β0
β1
β2
β3
α
β0
β1
β2
β3
-0.557(0.109;0.171)
-0.584(0.134;0.173)
-0.720(0.166;0.173)
0.024(0.157;0.251)
0.012(0.187;0.261)
0.072(0.235;0.246)
-0.177(0.025;0.030)
-0.177(0.021;0.031)
-0.141(0.028;0.029)
-0.078(0.039;0.055)
-0.089(0.036;0.057)
-0.114(0.047;0.052)
ALR
-0.524(0.157;0.169)
0.017(0.222;0.243)
-0.178(0.023;0.030)
-0.084(0.039;0.052)
3.222(
;0.291)
Linearization based method
-0.557(0.112;0.171)
-0.585(0.142;0.174)
-0.630(0.171;0.172)
0.024(0.160;0.251)
0.011(0.196;0.262)
0.036(0.242;0.242)
-0.177(0.025;0.030)
-0.177(0.022;0.031)
-0.204(0.038;0.034)
-0.078(0.040:0.055)
-0.089(0.038;0.057)
-0.106(0.058;0.058)
estimate (model-based s.e.; empirical s.e.)
Introduction to Longitudinal Data Analysis
402
21.5.6 Discussion
·	GEE1: All empirical standard errors are correct, but the eﬃciency is higher for the
tij eﬀect:
more complex working correlation structure, as seen in p-values for Ti ·
Structure p-value
IND
0.1515 EXCH
0.1208 UN
0.0275 Thus, opting for reasonably adequate correlation assumptions still pays oﬀ, in
spite of the fact that all are consistent and asymptotically normal
Similar conclusions for linearization-based method
·	Introduction to Longitudinal Data Analysis
403
Model-based s.e. and empirically corrected s.e. in reasonable agreement for UN
Typically, the model-based standard errors are much too small as they are based
on the assumption that all observations in the data set are independent, hereby
overestimating the amount of available information, hence also overestimating the
precision of the estimates.
ALR: similar inferences but now also α part of the inferences
§	Introduction to Longitudinal Data Analysis
404
Part III
Generalized Linear Mixed Models for Non-Gaussian
Longitudinal Data
Introduction to Longitudinal Data Analysis
405
Chapter 22
The Beta-binomial Model
·	Genesis of the model
·	Implied marginal distribution
Introduction to Longitudinal Data Analysis
406
22.1 Genesis of the Beta-binomial Model
Skellam (1948), Kleinman (1973)
Let Yi be a ni-dimensional vector of Bernoulli-distributed outcomes, with success
probability bi.
Assume the elements in Yi to be independent, conditionally on bi
Then, he conditional density of Yi, given bi is proportional to the density of
Zi =
ni
Xj=1
Yij
The density of Zi, given bi is binomial with ni trials and success probability bi.
o	Introduction to Longitudinal Data Analysis
407
·	The beta-binomial model assumes the bi to come from a beta distribution with
parameters α and β:
α, β) =
f (bi|
bα
1
−
i
bi)β
(1
−
B(α, β)
1
−
B(., .): the beta function
·	α and β can depend on covariates, but this dependence is temporarily dropped
from notation
Introduction to Longitudinal Data Analysis
408
22.2 Implied Marginal Model
The marginal density of Zi is the so-called beta-binomial density:
·	α, β) =
fi(zi|
=
Z








ni
zi
















bzi
i (1
bi)ni−
zif (bi|
−
α, β)dbi
B(zi + α, ni −
B(α, β)
zi + β)
ni
zi








Introduction to Longitudinal Data Analysis
409
Useful moments and relationships (π = µi/ni):
·	1
α = π(ρ−
−
1
π)(ρ−
β = (1
−
Mean
µi = E(Zi) = ni
−
α
α + β
Correlation
ρ = Corr(Yij, Yik) =
1
α + β + 1
Variance
Var(Zi) = niπ(1
π)[1 + (ni −
−
1)ρ]
Introduction to Longitudinal Data Analysis
410
The density can now be written as:
·	π, ρ) =
fi(zi|
ni
zi
















1
B[zi + π(ρ−
−
1
B[π(ρ−
1), ni −
1), (1
−
zi + (1
−
1
π)(ρ−
−
1)]
−
1
π)(ρ−
1)]
−
·	When there are covariates (e.g., sub-populations, dose groups), rewrite π and/or
ρ as πi and/or ρi, respectively.
It is then easy to formulate a model through the marginal parameters πi and ρi:
·	πi can be modeled through, e.g., a logit link
·	ρi can be modeled through, e.g., Fisher’s z transformation
In Part IV, the NTP data will be analyzed using the beta-binomial model
o	Introduction to Longitudinal Data Analysis
411
Chapter 23
Generalized Linear Mixed Models (GLMM)
·	Introduction: LMM Revisited
·	Generalized Linear Mixed Models (GLMM)
·	Fitting Algorithms
·	Example
Introduction to Longitudinal Data Analysis
412
23.1 Introduction: LMM Revisited
§	We re-consider the linear mixed model:
Yi
bi
|
∼
N (Xiβ + Zibi, Σi),
bi
∼
N (0, D)
The implied marginal model equals Yi
N (Xiβ, ZiDZ 0i + Σi)
∼
Hence, even under conditional independence, i.e., all Σi equal to σ2Ini, a marginal
association structure is implied through the random eﬀects.
·	The same ideas can now be applied in the context of GLM’s to model association
between discrete repeated measures.
Introduction to Longitudinal Data Analysis
413
23.2 Generalized Linear Mixed Models (GLMM)
·	Given a vector bi of random eﬀects for cluster i, it is assumed that all responses
Yij are independent, with density
f (yij|
θij, φ) = exp
(cid:26)
φ−
1[yijθij −
ψ(θij)] + c(yij, φ)
(cid:27)
§	θij is now modelled as θij = xij0β + zij0bi
As before, it is assumed that bi
N (0, D)
∼
Let fij(yij|
density of Yi equals
bi, β, φ) denote the conditional density of Yij given bi, the conditional
fi(yi
|
bi, β, φ) =
ni
Yj=1
fij(yij|
bi, β, φ)
Introduction to Longitudinal Data Analysis
414
The marginal distribution of Yi is given by
fi(yi
|
β, D, φ) =
fi(yi
Z
|
bi, β, φ) f (bi
D) dbi
|
=
ni
Z
Yj=1
fij(yij|
bi, β, φ) f (bi
D) dbi
|
where f (bi
|
D) is the density of the N (0, D) distribution.
The likelihood function for β, D, and φ now equals
o	L(β, D, φ) =
=
N
Yi=1
N
fi(yi
|
β, D, φ)
ni
Yi=1 Z
Yj=1
fij(yij|
bi, β, φ) f (bi
D) dbi
|
Introduction to Longitudinal Data Analysis
415
§	Under the normal linear model, the integral can be worked out analytically.
In general, approximations are required:
·	Approximation of integrand
·	Approximation of data
·	Approximation of integral
Predictions of random eﬀects can be based on the posterior distribution
f (bi
|
Yi = yi)
·	‘Empirical Bayes (EB) estimate’:
Posterior mode, with unknown parameters replaced by their MLE
Introduction to Longitudinal Data Analysis
416
23.3 Laplace Approximation of Integrand
Integrals in L(β, D, φ) can be written in the form I =
eQ(b)db
Z
Second-order Taylor expansion of Q(b) around the mode yields
Q(b)
Q(
b) +
c
≈
1
2
(b
b)0Q00(
b)(b
− c
c
b),
− c
Quadratic term leads to re-scaled normal density. Hence,
b).
(2π)q/2
Q00(
eQ(
b)
1/2
I
−
≈
−
(cid:12)
(cid:12)
(cid:12)
(cid:12)
c
(cid:12)
(cid:12)
(cid:12)
(cid:12)
c
Exact approximation in case of normal kernels
Good approximation in case of many repeated measures per subject
o	Introduction to Longitudinal Data Analysis
417
23.4 Approximation of Data
23.4.1 General Idea
o	Re-write GLMM as:
Yij = µij + εij = h(x0ijβ + z0ijbi) + εij
with variance for errors equal to Var(Yij|
bi) = φv(µij)
Linear Taylor expansion for µij:
·	Penalized quasi-likelihood (PQL): Around current
·	Marginal quasi-likelihood (MQL): Around current
bi
β and
β and bi = 0
d
c
c
Introduction to Longitudinal Data Analysis
418
23.4.2 Penalized quasi-likelihood (PQL)
§	Linear Taylor expansion around current
β and
c
bi:
d
Yij ≈
h(x0ij
β + z0ij
bi) + h0(x0ij
β + z0ij
bi)x0ij(β
c
c
c
c
β) + h0(x0ij
β + z0ij
bi)z0ij(bi
c
c
c
bi) + εij
c
−
−
≈
µij + v(
µij)x0ij(β
c
c
−
β) + v(
c
c
µij)z0ij(bi
bi) + εij
c
−
In vector notation: Yi
µi +
≈ d
ViXi(β
d
β) +
− c
ViZi(bi
d
bi) + εi
− d
Re-ordering terms yields:
Yi∗
1
V −
i
≡ d
(Yi
µi) + Xi
− d
β + Zi
bi
c
d
≈
Xiβ + Zibi + ε∗
i ,
·	Model ﬁtting by iterating between updating the pseudo responses Yi∗ and ﬁtting
the above linear mixed model to them.
Introduction to Longitudinal Data Analysis
419
23.4.3 Marginal quasi-likelihood (MQL)
§	Linear Taylor expansion around current
β and bi = 0:
c
Yij ≈
h(x0ij
β) + h0(x0ij
β)x0ij(β
c
β) + h0(x0ij
β)z0ijbi + εij
c
c
−
≈
µij + v(
µij)x0ij(β
c
c
µij)z0ijbi + εij
β) + v(
c
c
c
−
In vector notation: Yi
µi +
≈ d
ViXi(β
d
β) +
− c
ViZibi + εi
d
Re-ordering terms yields:
Yi∗
1
V −
i
≡ d
(Yi
µi) + Xi
− d
β
c
≈
Xiβ + Zibi + ε∗
i
·	Model ﬁtting by iterating between updating the pseudo responses Yi∗ and ﬁtting
the above linear mixed model to them.
Introduction to Longitudinal Data Analysis
420
23.4.4 PQL versus MQL
MQL only performs reasonably well if random-eﬀects variance is (very) small
Both perform bad for binary outcomes with few repeated measurements per cluster
With increasing number of measurements per subject:
·	MQL remains biased
·	PQL consistent
Improvements possible with higher-order Taylor expansions
·	Introduction to Longitudinal Data Analysis
421
23.5 Approximation of Integral
The likelihood contribution of every subject is of the form
f (z)φ(z)dz
Z
where φ(z) is the density of the (multivariate) normal distribution
Gaussian quadrature methods replace the integral by a weighted sum:
o	f (z)φ(z)dz
Z
Q
≈
Xq=1
wqf (zq)
·	Q is the order of the approximation. The higher Q the more accurate the
approximation will be
Introduction to Longitudinal Data Analysis
422
·	The nodes (or quadrature points) zq are solutions to the Qth order Hermite
polynomial
o	The wq are well-chosen weights
The nodes zq and weights wq are reported in tables. Alternatively, an algorithm is
available for calculating all zq and wq for any value Q.
·	With Gaussian quadrature, the nodes and weights are ﬁxed, independent of
f (z)φ(z).
·	With adaptive Gaussian quadrature, the nodes and weights are adapted to
the ‘support’ of f (z)φ(z).
Introduction to Longitudinal Data Analysis
423
Graphically (Q = 10):
·	Introduction to Longitudinal Data Analysis
424
·	Typically, adaptive Gaussian quadrature needs (much) less quadrature points than
classical Gaussian quadrature.
On the other hand, adaptive Gaussian quadrature is much more time consuming.
Adaptive Gaussian quadrature of order one is equivalent to Laplace transformation.
Ample detail can be found in Molenberghs and Verbeke (2005, Sections 14.3–14.5)
§	Introduction to Longitudinal Data Analysis
425
23.6 Example: Toenail Data
Yij is binary severity indicator for subject i at visit j.
Model:
Yij|
bi ∼
Bernoulli(πij),
log 



πij
πij
1
−




Notation:
·	Ti: treatment indicator for subject i
= β0 + bi + β1Ti + β2tij + β3Titij
·	tij: time point at which jth measurement is taken for ith subject
Adaptive as well as non-adaptive Gaussian quadrature, for various Q.
·	Introduction to Longitudinal Data Analysis
426
Results:
·	Gaussian quadrature
Q = 3
Q = 5
Q = 10
Q = 20
Q = 50
-1.52 (0.31)
-2.49 (0.39)
-0.99 (0.32)
-1.54 (0.69)
-1.65 (0.43)
-0.39 (0.38)
0.19 (0.36)
0.47 (0.36)
-0.43 (0.80)
-0.09 (0.57)
-0.32 (0.03)
-0.38 (0.04)
-0.38 (0.05)
-0.40 (0.05)
-0.40 (0.05)
-0.09 (0.05)
-0.12 (0.07)
-0.15 (0.07)
-0.14 (0.07)
-0.16 (0.07)
2.26 (0.12)
3.09 (0.21)
4.53 (0.39)
3.86 (0.33)
4.04 (0.39)
1344.1 1259.6
1254.4 1249.6
1247.7 Adaptive Gaussian quadrature
Q = 3
Q = 5
Q = 10
Q = 20
Q = 50
-2.05 (0.59)
-1.47 (0.40)
-1.65 (0.45)
-1.63 (0.43)
-1.63 (0.44)
-0.16 (0.64)
-0.09 (0.54)
-0.12 (0.59)
-0.11 (0.59)
-0.11 (0.59)
-0.42 (0.05)
-0.40 (0.04)
-0.41 (0.05)
-0.40 (0.05)
-0.40 (0.05)
-0.17 (0.07)
-0.16 (0.07)
-0.16 (0.07)
-0.16 (0.07)
-0.16 (0.07)
4.51 (0.62)
3.70 (0.34)
4.07 (0.43)
4.01 (0.38)
4.02 (0.38)
1259.1 1257.1
1248.2 1247.8
1247.8 β0
β1
β2
β3
σ
2`
−
β0
β1
β2
β3
σ
2`
−
Introduction to Longitudinal Data Analysis
427
·	Conclusions:
·	(Log-)likelihoods are not comparable
·	Diﬀerent Q can lead to considerable diﬀerences in estimates and standard
errors
·	For example, using non-adaptive quadrature, with Q = 3, we found no
diﬀerence in time eﬀect between both treatment groups
(t =
0.09/0.05, p = 0.0833).
−
·	Using adaptive quadrature, with Q = 50, we ﬁnd a signiﬁcant interaction
0.16/0.07, p = 0.0255).
between the time eﬀect and the treatment (t =
−
·	Assuming that Q = 50 is suﬃcient, the ‘ﬁnal’ results are well approximated
with smaller Q under adaptive quadrature, but not under non-adaptive
quadrature.
Introduction to Longitudinal Data Analysis
428
Comparison of ﬁtting algorithms:
·	Adaptive Gaussian Quadrature, Q = 50
·	MQL and PQL
Summary of results:
Parameter
QUAD
PQL
MQL
Intercept group A
Intercept group B
Slope group A
Slope group B
Var. random intercepts (τ 2)
1.63 (0.44)
−
1.75 (0.45)
−
0.40 (0.05)
−
0.57 (0.06)
−
15.99 (3.02)
0.72 (0.24)
−
0.72 (0.24)
−
0.29 (0.03)
−
0.40 (0.04)
−
0.56 (0.17)
−
0.53 (0.17)
−
0.17 (0.02)
−
0.26 (0.03)
−
4.71 (0.60)
2.49 (0.29)
Severe diﬀerences between QUAD (gold standard ?) and MQL/PQL.
MQL/PQL may yield (very) biased results, especially for binary data.
·	Introduction to Longitudinal Data Analysis
429
Chapter 24
Fitting GLMM’s in SAS
·	Proc GLIMMIX for adaptive Gaussian quadrature, Laplace, PQL, and MQL
·	Proc NLMIXED for adaptive and non-adaptive Gaussian quadrature
Introduction to Longitudinal Data Analysis
430
24.1 Procedure GLIMMIX for PQL and MQL
o	Re-consider logistic model with random intercepts for toenail data
SAS code (adaptive Gaussian quadrature; 50 quadrature points):
proc glimmix data=test method=gauss(q=50);
class idnum;
model onyresp (event=’1’) = treatn time treatn*time
/ dist=binary solution;
random intercept / subject=idnum;
run;
PQL obtained with option ‘method=RSPL’
MQL obtained with option ‘method=RMPL’
Inclusion of random slopes:
random intercept time / subject=idnum type=un;
Introduction to Longitudinal Data Analysis
431
Selected SAS output (PQL):
·	Covariance Parameter Estimates
Cov Parm
Subject
Estimate
Standard
Error
Intercept
idnum
4.7095 0.6024
Solutions for Fixed Effects
Effect
Estimate
Intercept
treatn
time
treatn*time
-0.7204
-0.02594
-0.2782
-0.09583
Standard
Error
0.2370 0.3360
0.03222 0.05105
DF
t Value
Pr > |t|
292
1612
1612
1612
-3.04
-0.08
-8.64
-1.88
0.0026 0.9385
<.0001
0.0607 Introduction to Longitudinal Data Analysis
432
24.2 Procedure NLMIXED for Gaussian Quadrature
Re-consider logistic model with random intercepts for toenail data
SAS program (non-adaptive, Q = 3):
proc nlmixed data=test noad qpoints=3;
parms beta0=-1.6 beta1=0 beta2=-0.4 beta3=-0.5 sigma=3.9;
teta = beta0 + b + beta1treatn + beta2time + beta3*timetr;
expteta = exp(teta);
p = expteta/(1+expteta);
model onyresp  binary(p);
random b  normal(0,sigma**2) subject=idnum;
run;
Adaptive Gaussian quadrature obtained by omitting option ‘noad’
§	Introduction to Longitudinal Data Analysis
433
Automatic search for ‘optimal’ value of Q in case of no option ‘qpoints=’
Selected SAS output (non-adaptive, Q = 3):
Parameter Estimates
Parameter Estimate
Standard
Error
DF t Value Pr > |t| Alpha
Lower
Upper
Gradient
beta0
beta1
beta2
beta3
sigma
-1.5311
-0.4294
-0.3107
-0.07539
2.2681 0.2961 293
0.3728 293
0.03373 293
0.04998 293
0.1220 293
-5.17
-1.15
-9.21
-1.51
18.58 <.0001
0.2503 <.0001
0.1325 <.0001
0.05 -2.1139 -0.9483
0.05 -1.1631
0.3043 0.05 -0.3771 -0.2443
0.05 -0.1738 0.02298
2.5083 0.05
2.0279 2.879E-7
-2.11E-6
-0.00003
-0.00003
-3.6E-6
Good starting values needed !
§	Introduction to Longitudinal Data Analysis
434
The inclusion of random slopes can be speciﬁed as follows:
·	proc nlmixed data=test noad qpoints=3;
parms beta0=-1.6 beta1=0 beta2=-0.4 beta3=-0.5
d11=3.9 d12=0 d22=0.1;
teta = beta0 + b1 + beta1treatn + beta2time
·	b2time + beta3timetr;
expteta = exp(teta);
p = expteta/(1+expteta);
model onyresp  binary(p);
random b1 b2  normal([0, 0] , [d11, d12, d22])
subject=idnum;
run;
Introduction to Longitudinal Data Analysis
435
24.2.1 Some Comments on the NLMIXED Procedure
·	Diﬀerent optimization algorithms are available to carry out the maximization of
the likelihood.
Constraints on parameters are also allowed in the optimization process.
The conditional distribution (given the random eﬀects) can be speciﬁed as
Normal, Binomial, Poisson, or as any distribution for which you can specify the
likelihood by programming statements.
E-B estimates of the random eﬀects can be obtained.
Only one RANDOM statement can be speciﬁed.
Only normal random eﬀects are allowed.
o	Introduction to Longitudinal Data Analysis
436
§	Does not calculate automatic initial values.
Make sure your data set is sorted by cluster ID!
PROC NLMIXED can perform Gaussian quadrature by using the options NOAD
and NOADSCALE. The number of quadrature points can be speciﬁed with the
option QPOINTS=m.
·	PROC NLMIXED can maximize the marginal likelihood using the
Newton-Raphson algorithm by specifying the option TECHNIQUE=NEWRAP.
Introduction to Longitudinal Data Analysis
437
24.2.2 The Main Statements
o	NLMIXED statement:
·	option ‘noad’ to request no adaptive quadrature
·	by default, adaptive Gaussian quadrature is used
·	the option ‘qpoints’ speciﬁes the number of quadrature points
·	by default, the number of quadrature points is selected adaptively by
evaluating the log-likelihood function at the starting values of the parameters
until two successive evaluations show suﬃciently small relative change.
PARMS statement:
·	starting values for all parameters in the model
·	by default, parameters not listed in the PARMS statement are given an initial
value of 1
Introduction to Longitudinal Data Analysis
438
o	MODEL statement:
·	conditional distribution of the data, given the random eﬀects
·	valid distributions:
normal(m,v): Normal with mean m and variance v
binary(p): Bernoullie with probability p
binomial(n,p): Binomial with count n and probability p
poisson(m): Poisson with mean m
general(ll): General model with log-likelihood ll
∗
∗
∗
∗
∗
·	since no factors can be deﬁned, explicit creation of dummies is required
RANDOM statement:
·	speciﬁcation of the random eﬀects
·	the procedure requires the data to be ordered by subject !
·	empirical Bayes estimates can be obtained by adding out=eb
Introduction to Longitudinal Data Analysis
439
Part IV
Marginal Versus Random-eﬀects Models and Case Studies
Introduction to Longitudinal Data Analysis
440
Chapter 25
Marginal Versus Random-eﬀects Models
·	Interpretation of GLMM parameters
·	Marginalization of GLMM
·	Conclusion
Introduction to Longitudinal Data Analysis
441
25.1 Interpretation of GLMM Parameters: Toenail Data
·	We compare our GLMM results for the toenail data with those from ﬁtting GEE’s
(unstructured working correlation):
GLMM
GEE
Parameter
Estimate (s.e.)
Estimate (s.e.)
Intercept group A
Intercept group B
Slope group A
Slope group B
1.6308 (0.4356)
1.7454 (0.4478)
0.4043 (0.0460)
0.5657 (0.0601)
−
−
−
−
−
−
−
−
0.7219 (0.1656)
0.6493 (0.1671)
0.1409 (0.0277)
0.2548 (0.0380)
Introduction to Longitudinal Data Analysis
442
·	The strong diﬀerences can be explained as follows:
·	Consider the following GLMM:
Yij|
Bernoulli(πij),
bi ∼
·	The conditional means E(Yij|
log 
1



πij
πij




= β0 + bi + β1tij
−
bi), as functions of tij, are given by
E(Yij|
bi)
=
exp(β0 + bi + β1tij)
1 + exp(β0 + bi + β1tij)
Introduction to Longitudinal Data Analysis
443
·	The marginal average evolution is now obtained from averaging over the
random eﬀects:
E(Yij) = E[E(Yij|
bi)] = E 



exp(β0 + bi + β1tij)
1 + exp(β0 + bi + β1tij)
exp(β0 + β1tij)
1 + exp(β0 + β1tij)




=
Introduction to Longitudinal Data Analysis
444
6
Hence, the parameter vector β in the GEE model needs to be interpreted
completely diﬀerent from the parameter vector β in the GLMM:
o	GEE: marginal interpretation
·	GLMM: conditional interpretation, conditionally upon level of random eﬀects
·	In general, the model for the marginal average is not of the same parametric form
as the conditional average in the GLMM.
·	For logistic mixed models, with normally distributed random random intercepts, it
can be shown that the marginal model can be well approximated by again a
logistic model, but with parameters approximately satisfying
RE
β
M
c
β
c
= √c2σ2 + 1 > 1,
σ2 = variance random intercepts
c = 16√3/(15π)
Introduction to Longitudinal Data Analysis
445
·	For the toenail application, σ was estimated as 4.0164, such that the ratio equals
√c2σ2 + 1 = 2.5649.
o	The ratio’s between the GLMM and GEE estimates are:
GLMM
GEE
Parameter
Estimate (s.e.)
Estimate (s.e.) Ratio
Intercept group A
Intercept group B
Slope group A
Slope group B
1.6308 (0.4356)
−
1.7454 (0.4478)
−
0.4043 (0.0460)
−
0.5657 (0.0601)
−
0.7219 (0.1656) 2.2590
−
0.6493 (0.1671) 2.6881
−
0.1409 (0.0277) 2.8694
−
0.2548 (0.0380) 2.2202
−
Note that this problem does not occur in linear mixed models:
·	Conditional mean: E(Yi
bi) = Xiβ + Zibi
·	Speciﬁcally: E(Yi
bi = 0) = Xiβ
·	Marginal mean: E(Yi) = Xiβ
|
|
Introduction to Longitudinal Data Analysis
446
The problem arises from the fact that, in general,
·	E[g(Y )]
= g[E(Y )]
·	So, whenever the random eﬀects enter the conditional mean in a non-linear way,
the regression parameters in the marginal model need to be interpreted diﬀerently
from the regression parameters in the mixed model.
·	In practice, the marginal mean can be derived from the GLMM output by
integrating out the random eﬀects.
·	This can be done numerically via Gaussian quadrature, or based on sampling
methods.
Introduction to Longitudinal Data Analysis
447
6
25.2 Marginalization of GLMM: Toenail Data
·	As an example, we plot the average evolutions based on the GLMM output
obtained in the toenail example:
P (Yij = 1)
=
E 


E 





exp(
−
1 + exp(
1.6308 + bi −
1.6308 + bi −
−
0.4043tij)
,
0.4043tij) 


exp(
−
1 + exp(
1.7454 + bi −
1.7454 + bi −
−
0.5657tij)
,
0.5657tij) 


Introduction to Longitudinal Data Analysis
448
SAS code (averaging over 1000 draws):
·	data h;
do treat=0 to 1 by 1;
do subject=1 to 1000 by 1;
b=4.0164*rannor(-1) ;
do t=0 to 12 by 0.1;
if treat=0 then y=exp(-1.6308 + b -0.4043*t)
/(1+ exp(-1.6308 + b -0.4043*t));
else y=exp(-1.7454 + b -0.5657*t)
/(1+ exp(-1.7454 + b -0.5657*t));
output;
end;
end;
end;
proc sort data=h;
by t treat;
run;
proc means data=h;
var y;
by t treat;
output out=out;
run;
proc gplot data=out;
plot y*t=treat / haxis=axis1 vaxis=axis2 legend=legend1;
axis1 label=(h=2 ’Time’) value=(h=1.5)
minor=none;
order=(0 to 14 by 1)
axis2 label=(h=2 A=90 ’P(Y=1)’) value=(h=1.5)
order=(0 to 0.4 by 0.1) minor=none;
legend1 label=(h=1.5 ’Treatment: ’)
value=(h=1.5 ’A’ ’B’);
Marginal average evolutions (GLMM)’;
title h=2.5 ’
symbol1 c=black i=join w=5 l=1 mode=include;
symbol2 c=black i=join w=5 l=2 mode=include;
where stat=’MEAN’;
run;quit;run;
Introduction to Longitudinal Data Analysis
449
Average evolutions obtained from the GEE analyses:
·	P (Yij = 1)
exp(
−
1 + exp(
0.7219 −
0.7219 −
exp(
−
1 + exp(
0.6493 −
0.6493 −
0.1409tij)
0.1409tij)
−
0.2548tij)
0.2548tij)
−
=



Introduction to Longitudinal Data Analysis
450
·	In a GLMM context, rather than plotting the marginal averages, one can also plot
the proﬁle for an ‘average’ subject, i.e., a subject with random eﬀect bi = 0:
P (Yij = 1
bi = 0)
|
exp(
−
1 + exp(
1.6308 −
1.6308 −
exp(
−
1 + exp(
1.7454 −
1.7454 −
0.4043tij)
0.4043tij)
−
0.5657tij)
0.5657tij)
−
=



Introduction to Longitudinal Data Analysis
451
25.3 Example: Toenail Data Revisited
Overview of all analyses on toenail data:
·	Parameter
QUAD
PQL
MQL
GEE
Intercept group A
Intercept group B
Slope group A
Slope group B
Var. random intercepts (τ 2)
1.63 (0.44)
−
0.72 (0.24)
−
0.56 (0.17)
−
0.72 (0.17)
−
1.75 (0.45)
−
0.40 (0.05)
−
0.57 (0.06)
−
15.99 (3.02)
0.72 (0.24)
−
0.29 (0.03)
−
0.40 (0.04)
−
0.53 (0.17)
−
0.17 (0.02)
−
0.26 (0.03)
−
4.71 (0.60)
2.49 (0.29)
0.65 (0.17)
−
0.14 (0.03)
−
0.25 (0.04)
−
Conclusion:
·	GEE
|
|
<
MQL
|
|
<
PQL
|
|
<
QUAD
|
|
Introduction to Longitudinal Data Analysis
452
Model Family
·	marginal
model
↓
inference
·	&
likelihood GEE
↓
β M
↓
βM
&
random-eﬀects
model
↓
inference
·	&
marginal hierarchical
↓
β RE
↓
‘β M’
↓
(β RE, bi)
↓
‘β M’
Introduction to Longitudinal Data Analysis
453
Chapter 26
Case Study: The NTP Data
·	Research question
·	Conditional model
·	Bahadur model
·	GEE1 analyses
·	GEE2 analysis
·	Alternating logistic regressions
·	Beta-binomial model
·	Generalized linear mixed model
·	Discussion
Introduction to Longitudinal Data Analysis
454
26.1 Research Question
§	Dose-response relationship: eﬀect of dose on malformations
Regression relationship:
logit[P (Yij = 1
|
di, . . .)] = β0 + βd di
Association parameter: βd Precise meaning is model-dependent:
·	Transformed conditional odds ratio
·	Transformed correlation
·	Transformed marginal odds ratio
Introduction to Longitudinal Data Analysis
455
26.2 Conditional Model
·	Regression relationship:
logit[P (Yij = 1
|
di, Yik = 0, k
= j)] = β0 + βd di
δi = βa is conditional log odds ratio
Quadratic loglinear model
Maximum likelihood estimates (model based standard errors; empirically corrected
standard errors)
Introduction to Longitudinal Data Analysis
456
6
Outcome Par.
DEHP
EG
DYME
External
Visceral
Skeletal
β0
βd
βa
β0
βd
βa
β0
βd
βa
-2.81(0.58;0.52) -3.01(0.79;1.01) -5.78(1.13;1.23)
3.07(0.65;0.62)
2.25(0.68;0.85)
6.25(1.25;1.41)
0.18(0.04;0.04)
0.25(0.05;0.06)
0.09(0.06;0.06)
-2.39(0.50;0.52) -5.09(1.55;1.51) -3.32(0.98;0.89)
2.45(0.55;0.60)
3.76(1.34;1.20)
2.88(0.93;0.83)
0.18(0.04;0.04)
0.23(0.09;0.09)
0.29(0.05;0.05)
-2.79(0.58;0.77) -0.84(0.17;0.18) -1.62(0.35;0.48)
2.91(0.63;0.82)
0.98(0.20;0.20)
2.45(0.51;0.82)
0.17(0.04;0.05)
0.20(0.02;0.02)
0.25(0.03;0.03)
Collapsed β0
-2.04(0.35;0.42) -0.81(0.16;0.16) -2.90(0.43;0.51)
βd
βa
2.98(0.51;0.66)
0.97(0.20;0.20)
5.08(0.74;0.96)
0.16(0.03;0.03)
0.20(0.02;0.02)
0.19(0.03;0.03)
Introduction to Longitudinal Data Analysis
457
26.3 The Bahadur Model
Regression relationship:
logit[P (Yij = 1
|
di)] = β0 + βd di
βa: Fisher’s z transformed correlation
ρ: correlation
§	Introduction to Longitudinal Data Analysis
458
Outcome Parameter
DEHP
EG
DYME
External
Visceral
Skeletal
β0
βd
βa
ρ
β0
βd
βa
ρ
β0
βd
βa
ρ
-4.93(0.39)
-5.25(0.66)
-7.25(0.71)
5.15(0.56)
2.63(0.76)
7.94(0.77)
0.11(0.03)
0.12(0.03)
0.11(0.04)
0.05(0.01)
0.06(0.01)
0.05(0.02)
-4.42(0.33)
-7.38(1.30)
-6.89(0.81)
4.38(0.49)
4.25(1.39)
5.49(0.87)
0.11(0.02)
0.05(0.08)
0.08(0.04)
0.05(0.01)
0.02(0.04)
0.04(0.02)
-4.67(0.39)
-2.49(0.11)
-4.27(0.61)
4.68(0.56)
2.96(0.18)
5.79(0.80)
0.13(0.03)
0.27(0.02)
0.22(0.05)
0.06(0.01)
0.13(0.01)
0.11(0.02)
Collapsed β0
-3.83(0.27)
-2.51(0.09)
-5.31(0.40)
βd
βa
ρ
5.38(0.47)
3.05(0.17)
8.18(0.69)
0.12(0.03)
0.28(0.02)
0.12(0.03)
0.06(0.01)
0.14(0.01)
0.06(0.01)
Introduction to Longitudinal Data Analysis
459
26.4 GEE1
o	Regression relationship:
logit[P (Yij = 1
|
di)] = β0 + βd di
φ: overdispersion parameter
ρ: working correlation
Parameter estimates (model-based standard errors; empirically corrected standard
errors)
Two sets of working assumptions:
·	Independence working assumptions
·	Exchangeable working assumptions
Introduction to Longitudinal Data Analysis
460
Outcome Par.
Standard
Prentice
Linearized
External
Visceral
Skeletal
β0
βd
φ
β0
βd
φ
β0
βd
φ
-5.06(0.30;0.38) -5.06(0.33;0.38) -5.06(0.28;0.38)
5.31(0.44;0.57)
5.31(0.48;0.57)
5.31(0.42;0.57)
0.90 0.74
-4.47(0.28;0.36) -4.47(0.28;0.36) -4.47(0.28;0.36)
4.40(0.43;0.58)
4.40(0.43;0.58)
4.40(0.43;0.58)
1.00 1.00
-4.87(0.31;0.47) -4.87(0.31;0.47) -4.87(0.32;0.47)
4.89(0.46;0.65)
4.90(0.47;0.65)
4.90(0.47;0.65)
0.99 1.02
Collapsed β0
-3.98(0.22;0.30) -3.98(0.22;0.30) -3.98(0.22;0.30)
βd
φ
5.56(0.40;0.61)
5.56(0.40;0.61)
5.56(0.41;0.61)
0.99 1.04
Introduction to Longitudinal Data Analysis
461
Outcome Par.
Standard
Prentice
Linearized
External
Visceral
Skeletal
Collapsed
β0
βd
φ
ρ
β0
βd
φ
ρ
β0
βd
φ
ρ
β0
βd
φ
ρ
-4.98(0.40;0.37)
-4.99(0.46;0.37)
-5.00(0.36;0.37)
5.33(0.57;0.55)
5.32(0.65;0.55)
5.32(0.51;0.55)
0.88 0.11
0.11(0.04)
0.65 0.06
-4.50(0.37;0.37)
-4.51(0.40;0.37)
-4.50(0.36;0.37)
4.55(0.55;0.59)
4.59(0.58;0.59)
4.55(0.54;0.59)
1.00 0.08
0.11(0.05)
0.92 0.08
-4.83(0.44;0.45)
-4.82(0.47;0.44)
-4.82(0.46;0.45)
4.84(0.62;0.63)
4.84(0.67;0.63)
4.84(0.65;0.63)
0.98 0.12
0.14(0.06)
0.86 0.13
-4.05(0.32;0.31)
-4.06(0.35;0.31)
-4.04(0.33;0.31)
5.84(0.57;0.61)
5.89(0.62;0.61)
5.82(0.58;0.61)
1.00 0.11
0.15(0.05)
0.96 0.11
Introduction to Longitudinal Data Analysis
462
26.5 GEE2
Regression relationship:
logit[P (Yij = 1
|
di)] = β0 + βd di
βa: Fisher’s z transformed correlation
ρ: correlation
Working assumption: third- and fourth-order correlations are zero
Parameter estimates (empirically corrected standard errors)
o	Introduction to Longitudinal Data Analysis
463
Outcome Parameter
DEHP
EG
DYME
External
Visceral
Skeletal
Collapsed
β0
βd
βa
ρ
β0
βd
βa
ρ
β0
βd
βa
ρ
β0
βd
βa
ρ
-4.98(0.37)
-5.63(0.67)
-7.45(0.73)
5.29(0.55)
3.10(0.81)
8.15(0.83)
0.15(0.05)
0.15(0.05)
0.13(0.05)
0.07(0.02)
0.07(0.02)
0.06(0.02)
-4.49(0.36)
-7.50(1.05)
-6.89(0.75)
4.52(0.59)
4.37(1.14)
5.51(0.89)
0.15(0.06)
0.02(0.02)
0.11(0.07)
0.07(0.03)
0.01(0.01)
0.05(0.03)
-5.23(0.40)
-4.05(0.33)
5.35(0.60)
4.77(0.43)
0.18(0.02)
0.30(0.03)
0.09(0.01)
0.15(0.01)
-5.23(0.40)
-4.07(0.71)
-5.75(0.48)
5.35(0.60)
4.89(0.90)
8.82(0.91)
0.18(0.02)
0.26(0.14)
0.18(0.12)
0.09(0.01)
0.13(0.07)
0.09(0.06)
Introduction to Longitudinal Data Analysis
464
26.6 Alternating Logistic Regressions
Regression relationship:
logit[P (Yij = 1
|
di)] = β0 + βd di
Exchangeable association structure
α: log odds ratio
ψ: odds ratio
Parameter estimates (empirically corrected standard errors)
o	Introduction to Longitudinal Data Analysis
465
Outcome Parameter
DEHP
EG
DYME
External
Visceral
Skeletal
Collapsed
β0
βd
α
ψ
β0
βd
α
ψ
β0
βd
α
ψ
β0
βd
α
ψ
-5.16(0.35)
-5.72(0.64)
-7.48(0.75)
5.64(0.52)
3.28(0.72)
8.25(0.87)
0.96(0.30)
1.45(0.45)
0.79(0.31)
2.61(0.78)
4.26(1.92)
2.20(0.68)
-4.54(0.36)
-7.61(1.06)
-7.24(0.88)
4.72(0.57)
4.50(1.13)
6.05(1.04)
1.12(0.30)
0.49(0.42)
1.76(0.59)
3.06(0.92)
1.63(0.69)
5.81(3.43)
-4.87(0.49)
-3.28(0.22)
-4.92(0.34)
4.90(0.70)
3.85(0.39)
6.73(0.65)
1.05(0.40)
1.43(0.22)
1.62(0.37)
2.86(1.14)
4.18(0.92)
5.05(1.87)
-4.04(0.31)
-3.19(0.22)
-5.08(0.37)
5.93(0.63)
3.86(0.40)
7.98(0.75)
1.17(0.29)
1.40(0.22)
1.26(0.31)
3.22(0.93)
4.06(0.89)
3.53(1.09)
Introduction to Longitudinal Data Analysis
466
26.7 Beta-binomial Model
Regression relationship:
logit[P (Yij = 1
|
di, )] = β0 + βd di
βa: Fisher’s z transformed correlation
ρ: correlation
Parameter estimates (standard errors)
·	Introduction to Longitudinal Data Analysis
467
Outcome Parameter
DEHP
EG
DYME
External
Visceral
Skeletal
β0
βd
βa
ρ
β0
βd
βa
ρ
β0
βd
βa
ρ
-4.91(0.42)
-5.32(0.71)
-7.27(0.74)
5.20(0.59)
2.78(0.81)
8.01(0.82)
0.21(0.09)
0.28(0.14)
0.21(0.12)
0.10(0.04)
0.14(0.07)
0.10(0.06)
-4.38(0.36)
-7.45(1.17)
-6.21(0.83)
4.42(0.54)
4.33(1.26)
4.94(0.90)
0.22(0.09)
0.04(0.09)
0.45(0.21)
0.11(0.04)
0.02(0.04)
0.22(0.10)
-4.88(0.44)
-2.89(0.27)
-5.15(0.47)
4.92(0.63)
3.42(0.40)
6.99(0.71)
0.27(0.11)
0.54(0.09)
0.61(0.14)
0.13(0.05)
0.26(0.04)
0.30(0.06)
Collapsed β0
-3.83(0.31)
-2.51(0.09)
-5.42(0.45)
βd
βa
ρ
5.59(0.56)
3.05(0.17)
8.29(0.79)
0.32(0.10)
0.28(0.02)
0.33(0.10)
0.16(0.05)
0.14(0.01)
0.16(0.05)
Introduction to Longitudinal Data Analysis
468
26.8 Generalized Linear Mixed Model
Regression relationship:
logit[P (Yij = 1
|
di, bi)] = β0 + bi + βd di,
N (0, τ 2)
bi ∼
External malformation in DEHP study
Four ways of dealing with the integral: Laplace, adaptive Gaussian quadrature,
PQL, and MQL
Two versions of PQL and MQL: ML and REML
Parameter estimates (standard errors)
o	Introduction to Longitudinal Data Analysis
469
Eﬀect
Parameter
Laplace
QUAD
Intercept
Dose eﬀect
Intercept var.
β0
βd
τ 2
-6.02 (0.59) -5.97 (0.57)
6.50 (0.86)
6.45 (0.84)
1.42 (0.70)
1.27 (0.62)
Eﬀect
Parameter PQL (REML) PQL (ML)
Intercept
Dose eﬀect
Intercept var.
β0
βd
τ 2
-5.32 (0.40) -5.30 (0.40)
5.73 (0.65)
5.71 (0.64)
0.95 (0.40)
0.89 (0.38)
Eﬀect
Parameter MQL (REML) MQL (ML)
Intercept
Dose eﬀect
Intercept var.
β0
βd
τ 2
-5.18 (0.40) -5.17 (0.39)
5.70 (0.66)
5.67 (0.65)
1.20 (0.53)
1.10 (0.50)
Introduction to Longitudinal Data Analysis
470
26.9 Summary Table
External malformation in DEHP study
All conditional, marginal, and random-eﬀects models considered
Parameter estimates (standard errors)
For non-likelihood methods, the empirically corrected standard errors are reported
·	Introduction to Longitudinal Data Analysis
471
Family
Model
β0
βd
Association
Conditional
Quadr. loglin. (ML)
-2.81(0.58)
3.07(0.65)
LOG OR
0.18(0.04)
Quadr. loglin. (PL)
-2.85(0.53)
3.24(0.60)
LOG OR
0.18(0.04)
Marginal
Lik. Bahadur
-4.93(0.39)
5.15(0.56)
St. GEE1 (exch)
-4.98(0.37)
5.33(0.55)
St. GEE1 (ind)
-5.06(0.38)
5.31(0.57)
Prent. GEE1 (exch)
-4.99(0.37)
5.32(0.55)
Prent. GEE1 (ind)
-5.06(0.38)
5.31(0.57)
Lin. based (exch)
-5.00(0.37)
5.32(0.55)
Lin. based (ind)
-5.06(0.38)
5.31(0.57)
GEE2
ALR
-4.98(0.37)
5.29(0.55)
-.516(0.35)
5.64(0.52)
Random-eﬀects
Beta-binomial
-4.91(0.42)
5.20(0.59)
ρ
ρ
ρ
ρ
ρ
βa
ρ
0.05(0.01)
0.11 0.11 (0.04)
0.06 0.07(0.02)
0.96(0.30)
0.10(0.04)
GLLM (MQL)
-5.18(0.40)
5.70(0.66)
Int. var τ 2
1.20(0.53)
GLMM (PQL)
-5.32(0.40)
5.73(0.65)
Int. var τ 2
0.95(0.40)
GLMM (QUAD)
-5.97(0.57)
6.45(0.84)
Int. var τ 2
1.27(0.62)
Introduction to Longitudinal Data Analysis
472
26.10 Discussion
Relationship between regression model parameters:
·	conditional
|
|
<
marginal
|
|
<
random-eﬀects
|
|
·	Beta-binomial model behaves like a marginal model (similar to the linear mixed
model)
·	Marginal model parameters:
·	Mean function parameters: very similar
·	Correlation parameters:
Bahadur
|
|
<
|
GEE2
|
<
|
GEE1
<
beta-binomial
|
|
|
Introduction to Longitudinal Data Analysis
473
·	Reason: strength of constraints:
Bahadur model valid if all higher order probabilities are valid
GEE2 valid if probabilities of orders 1, 2, 3, and 4 are valid
GEE1 valid if probabilities of orders 1 and 2 are valid
beta-binomial model is unconstrained of correlations in [0, 1]
∗
∗
∗
∗
·	Correlation in Bahadur model really highly constrained:
For instance, the allowable range of βa for the external outcome in the DEHP data is
(
−
a beta-binomial model. It translates to (
0.0164; 0.1610) when β0 and βd are ﬁxed at their MLE. This range excludes the MLE under
0.0082; 0.0803) on the correlation scale.
−
Additional conditional and marginal approaches can be based on
pseudo-likelihood (Molenberghs and Verbeke 2005, Chapters 9 and 12, in
particular pages 200 and 246)
Programs: Molenberghs and Verbeke (2005, p. 219ﬀ)
o	Introduction to Longitudinal Data Analysis
474
The random eﬀects in generalized linear mixed models
o	enter linearly on the logit scale:
logit[P (Yij = 1
|
di, bi] = β0 + bi + β1 di
mean of random intercepts is 0
mean of average over litters is
mean of predicted value over litters is
−
3.8171 3.8171
−
∗
∗
∗
·	enter non-linearly on the probability scale:
P (Yij = 1
|
di, bi) =
exp(β0 + bi + β1 di)
1 + exp(β0 + bi + β1 di)
mean of random eﬀect is 0.0207
mean of average probabilities over litters is 0.0781
mean of predicted probabilities over litters is 0.0988
∗
∗
∗
Introduction to Longitudinal Data Analysis
475
Chapter 27
Case Study: Binary Analysis of Analgesic Trial
·	Research question
·	GEE
·	Alternating logistic regressions
·	Further GEE analyses
·	Generalized linear mixed model
·	Discussion
Introduction to Longitudinal Data Analysis
476
27.1 Research Question
Binary version of Global Satisfaction Assessment
GSABIN =
1 if GSA
≤
0 otherwise.



Marginal regression relationship:
3 (‘Very Good’ to ‘Moderate’),
logit[P (Yij = 1
|
tij, Xi)] = β0 + β1tij + β2t2
ij + β3Xi.
GLMM regression relationship:
logit[P (Yij = 1
|
tij, Xi, bi)] = β0 + bi + β1tij + β2t2
ij + β3Xi.
Xi: baseline pain control assessment (PCA0)
Association parameters: correlation or marginal odds ratio
o	Introduction to Longitudinal Data Analysis
477
27.2 GEE1
·	Parameter estimates (model-based standard errors; empirically corrected standard
errors)
·	Four sets of working assumptions:
·	Independence
·	Exchangeable
·	AR(1)
·	Unstructured
Introduction to Longitudinal Data Analysis
478
Eﬀect
Intercept
Time
Time2
Basel. PCA
Correlation
Eﬀect
Intercept
Time
Time2
Basel. PCA
Correlation
Correlation (1,2)
Correlation (1,3)
Correlation (1,4)
Correlation (2,3)
Correlation (2,4)
Correlation (3,4)
Parameter
IND
EXCH
β1
β2
β3
β4
ρ
Parameter
β1
β2
β3
β4
ρ
ρ12
ρ13
ρ14
ρ23
ρ24
ρ34
2.80(0.49;0.47)
2.92(0.49;0.46)
-0.79(0.39;0.34)
-0.83(0.34;0.33)
0.18(0.08;0.07)
0.18(0.07;0.07)
-0.21(0.09;0.10)
-0.23(0.10;0.10)
—
AR
0.22 UN
2.94(0.49;0.47)
2.87(0.48;0.46)
-0.90(0.35;0.33)
-0.78(0.33;0.32)
0.20(0.07;0.07)
0.17(0.07;0.07)
-0.22(0.10;0.10)
-0.23(0.10;0.10)
0.25 —
0.18 0.25
0.20 0.18
0.18 0.46
Introduction to Longitudinal Data Analysis
479
Fitted working correlation matrices:
·	REXCH =
1 0.22 0.22 0.22
1
0.22 0.22
1
0.22 1





















RAR =





















1 0.25 0.06 0.02
1
0.25 0.06
1
0.25 1










































RUN =





















1 0.18 0.25 0.20
1
0.18 0.18
1
0.46 1





















Introduction to Longitudinal Data Analysis
480
27.3 Alternating Logistic Regressions
o	Parameter estimates (empirically corrected standard errors)
Three sets of odds ratio structures:
·	Exchangeable
·	Unstructured
≡
full clustering (FULLCLUST)
·	User-deﬁned design (ZREP)
Introduction to Longitudinal Data Analysis
481
Eﬀect
Parameter
EXCH
FULLCLUST
ZREP
Intercept
Time
Time2
Basel. PCA
Log OR
Log OR(1,2)
Log OR(1,3)
Log OR(1,4)
Log OR(2,3)
Log OR(2,4)
Log OR(3,4)
Log OR par.
Log OR par.
β1
β2
β3
β4
α
α12
α13
α14
α23
α24
α34
α0
α1
2.98(0.46)
2.92(0.46)
2.92(0.46)
-0.87(0.32)
-0.80(0.32) -0.80(0.32)
0.18(0.07)
0.17(0.06)
0.17(0.07)
-0.23(0.22)
-0.24(0.10) -0.24(0.10)
1.43(0.22)
1.13(0.33)
1.56(0.39)
1.60(0.42)
1.19(0.37)
0.93(0.42)
2.44(0.48)
1.26(0.23)
1.17(0.47)
Introduction to Longitudinal Data Analysis
482
·	In the FULLCLUST structure, there is a hint that α34 is diﬀerent from the others,
with all others being equal.
·	To conﬁrm this, a Wald test can be used for the null hypothesis:
H0 : α12 = α13 = α14 = α23 = α24
Details on the test: Molenberghs and Verbeke (2005, pp. 312–313)
The reduced structure, ﬁtted with ZREP, is:
α12 = α13 = α14 = α23 = α24 = α0,
α34 = α0 + α1
At the odds ratio level, with ﬁtted values:
ψ12 =
ψ13 =
ψ14 =
ψ23 =
c
c
c
c
ψ24 =
ψ34 =
c
c
c
ψ0 = 3.53,
ψ0 · c
c
ψ1 = 11.36.
Introduction to Longitudinal Data Analysis
483
“Odds ratio matrices”:
·	ΨEXCH =



1 4.18 4.18 4.18
1
4.18 4.18
1
4.18 


ΨUN =



1 3.10 4.76 4.95
1
3.29 2.53
1
11.47 1



1



ΨZREP =
1 3.53 3.53 3.53
1
3.53 3.53
1
11.36 1



Introduction to Longitudinal Data Analysis
484
27.4 A Variety of GEE Methods
§	Methods used:
·	Ordinary logistic regression
·	Standard GEE1
·	Prentice’s GEE1
·	The linearization-based method
·	Alternating logistic regression
Exchangeably working assumption (except for logistic regression)
Parameter estimates (empirically corrected standard errors, unless for logistic
regression)
Introduction to Longitudinal Data Analysis
485
Eﬀect
Parameter Log. regr.
Standard
Prentice
Intercept
Time
Time2
Basel. PCA
Correlation
β1
β2
β3
β4
ρ
2.80(0.49)
2.92(0.46)
2.94(0.46)
-0.79(0.39) -0.83(0.33) -0.84(0.33)
0.18(0.08)
0.18(0.07)
0.18(0.07)
-0.21(0.09) -0.23(0.10) -0.23(0.10)
0.21 0.26(0.05)
Eﬀect
Parameter
Lineariz.
ALR
Intercept
Time
Time2
Basel. PCA
Corr.
Log OR
β1
β2
β3
β4
ρ
α
2.94(0.46)
2.98(0.46)
-0.84(0.33) -0.87(0.32)
0.18(0.07)
0.18(0.07)
-0.23(0.10) -0.23(0.10)
0.26(0.04)
1.43(0.22)
Introduction to Longitudinal Data Analysis
486
27.5 Generalized Linear Mixed Models
Four tools:
·	SAS procedure GLIMMIX:
·	SAS procedure NLMIXED:
MQL (= MQL1)
PQL (= PQL1)
∗
∗
·	MLwiN:
PQL1
PQL2
∗
∗
∗
∗
I: non-adaptive (Q = 10)
II: non-adaptive (Q = 10)
adaptive (Q = 10)
adaptive (Q = 20)
≡
≡
·	MIXOR
Parameter estimates (standard errors)
o	Introduction to Longitudinal Data Analysis
487
Integrand approximation
SAS GLIMMIX
MLwiN
Eﬀect
Intercept
Time
Time2
Basel. PCA
Rand. int s.d.
Rand. int var.
Par.
MQL
PQL1
PQL1
PQL2
β1
β2
β3
β4
τ
τ 2
2.91(0.53)
3.03(0.55)
3.02(0.55)
4.07(0.70)
-0.83(0.39)
-0.87(0.41)
-0.87(0.41)
-1.17(0.48)
0.18(0.08)
0.19(0.08)
0.19(0.08)
0.25(0.10)
-0.22(0.11)
-0.22(0.11)
-0.22(0.11)
-0.31(0.15)
1.06(0.25)
1.04(0.23)
1.01(0.12)
1.61(0.15)
1.12(0.53)
1.08(0.48)
1.02(0.25)
2.59(0.47)
Numerical integration
SAS NLMIXED
Eﬀect
Intercept
Time
Time2
Basel. PCA
Rand. int s.d.
Rand. int var.
Par.
I
II
MIXOR
β1
β2
β3
β4
τ
τ 2
4.07(0.71)
4.05(0.71)
4.05(0.55)
-1.16(0.47)
-1.16(0.47)
-1.16(0.45)
0.25(0.09)
0.24(0.09)
0.24(0.10)
-0.30(0.14)
-0.30(0.14)
-0.30(0.15)
1.60(0.22)
1.59(0.21)
1.59(0.21)
2.56(0.70)
2.53(0.68)
2.53(0.67)
Introduction to Longitudinal Data Analysis
488
27.6 Discussion
Results are very similar, due to a relatively weak random-eﬀects variance
PQL1 and MQL1 perform relatively poorly
The ratio between the RE and marginal parameters now is 1.37
Programs: Molenberghs and Verbeke (2005, p. 219ﬀ)
·	Introduction to Longitudinal Data Analysis
489
Chapter 28
Case Study: Ordinal Analysis of Analgesic Trial
·	Proportional odds logistic regression
·	Generalized estimating equations
·	Generalized linear mixed models
·	Analysis of the analgesic trial
Introduction to Longitudinal Data Analysis
490
28.1 Proportional Odds Logistic Regression
·	Standard logistic regression for binary data:
logit[P (Yi = 1
|
xi)] = α + βxi
An extension to ordinal data: proportional odds logistic regression
logit[P (Yi ≤
k
|
xi)] = αk + βxi,
(k = 1, . . . , c
−
A further extension poses problems with range-preserving restrictions:
logit[P (Yi ≤
|
and is usually not considered
k
xi)] = αk + βoxi,
(k = 1, . . . , c
−
An alternative model for ordinal data is the continuation-ratio model:
logit[P (Yi > k
Yi ≥
|
k, xi)] = αk + βkxi,
(k = 1, . . . , c
−
Introduction to Longitudinal Data Analysis
491
It is of use only when there is one natural directionality in the data: subjects go
from the lowest category to higher categories, without ever returning. This is
often not satisﬁed.
Proportional-odds model for the 5-point GSA outcome in the analgesic trial:
logit[P (Yij ≤
k
|
tij, Xi)] = αk + β2tij + β3t2
ij + β4Xi,
(k = 1, . . . , 4)
SAS code:
proc genmod data=m.gsa2;
title ’Analgesic, logistic regression, Ordinal’;
class patid timecls;
model gsa = time|time pca0 / dist=multinomial link=cumlogit;
run;
Note that the ‘dist’ and ‘link’ options have been adapted
§	Introduction to Longitudinal Data Analysis
492
Selected output:
·	The GENMOD Procedure
Analysis Of Parameter Estimates
Parameter
DF
Estimate
Intercept1
Intercept2
Intercept3
Intercept4
TIME
TIME*TIME
PCA0
1
1
1
1
1
1
1
-1.0048
0.5225 2.3171
4.0525 -0.2027
0.0479 -0.2141
Standard
Error
Wald 95% Confidence
Limits
0.3437 0.3407
0.3481 0.3754
0.2706 0.0545
0.0622 -1.6785
-0.1452
1.6349 3.3166
-0.7330
-0.0590
-0.3361
-0.3312
1.1903 2.9994
4.7884 0.3277
0.1547 -0.0922
Chi-
Square
8.55 2.35
44.31 116.51
0.56 0.77
11.84 Pr > ChiSq
0.0035 0.1251
<.0001
<.0001
0.4539 0.3798
0.0006 There are 5
·	−
1 = 4 intercepts, as it should.
Introduction to Longitudinal Data Analysis
493
28.2 Generalized Estimating Equations
§	The same regression model as in the PO logistic regression case is used:
tij, Xi)] = αk + β2tij + β3t2
ij + β4Xi,
k
(k = 1, . . . , 4)
logit[P (Yij ≤
|
This model is supplemented with working assumptions to obtain GEE
In the SAS procedure GENMOD, only independence working assumptions are
implemented for ordinal outcomes:
proc genmod data=m.gsa2;
title ’Analgesic, GEE, Ordinal’;
class patid timecls;
model gsa = time|time pca0 / dist=multinomial link=cumlogit;
repeated subject=patid / type=ind covb corrw within=timecls modelse;
run;
Introduction to Longitudinal Data Analysis
494
The output is structured in the same way as for PO logistic regression:
·	Analysis Of GEE Parameter Estimates
Empirical Standard Error Estimates
Parameter Estimate
Standard
Error
95% Confidence
Limits
Z Pr > |Z|
Intercept1 -1.0048
0.5225 Intercept2
2.3171 Intercept3
Intercept4
4.0525 -0.2027
TIME
TIME*TIME
0.0479 -0.2141
PCA0
0.3549 -1.7004 -0.3092
1.2218 0.3568 -0.1767
3.0363 1.5980
0.3669 4.8243
0.3938 3.2807
0.1948 0.2028 -0.6001
0.0399 -0.0304
0.1261 0.0911 -0.3927 -0.0356
-2.83
1.46 6.31
10.29 -1.00
1.20 -2.35
0.0046 0.1430
<.0001
<.0001
0.3176 0.2304
0.0187 Introduction to Longitudinal Data Analysis
495
28.3 Generalized Linear Mixed Models
A generalized linear mixed model for ordinal data:
·	logit[P (Yij ≤
k
|
Xi, Zi)] = αk + x0ijβ + z0ijbi,
(k = 1, . . . , c
−
·	This is the obvious counterpart for the PO logistic and GEE marginal models
considered above.
·	For the case of the 5-point GSA outcome in the analgesic study:
logit[P (Yij ≤
tij, Xi, bi)] = αk + bi + β2tij + β3t2
ij + β4Xi,
k
|
(k = 1, . . . , 4)
Introduction to Longitudinal Data Analysis
496
Code for the SAS procedure GLIMMIX:
proc glimmix data=m.gsa2 method=RSPL;
title ’PROC GLIMMIX analysis, ordinal, RSPL (PQL, REML)’;
class patid timecls;
nloptions maxiter=50;
model gsa = time|time pca0 / dist=multinomial link=cumlogit solution;
random intercept / subject=patid type=un;
run;
Also here, the ‘dist’ and ‘link’ functions have to be adapted to the ordinal setting.
o	Introduction to Longitudinal Data Analysis
497
Selected output:
·	Covariance Parameter Estimates
Cov
Parm
Subject
Estimate
Standard
Error
UN(1,1)
PATID
3.5348 0.4240
Solutions for Fixed Effects
Effect
GSA
Estimate
1
2
3
4
Intercept
Intercept
Intercept
Intercept
TIME
TIME*TIME
PCA0
-1.4352
0.9101 3.4720
5.6263 -0.4825
0.1009 -0.2843
Standard
Error
0.5033 0.4999
0.5084 0.5358
0.2958 0.05972
0.1249 DF
t Value
Pr > |t|
393
393
393
393
737
737
737
-2.85
1.82 6.83
10.50 -1.63
1.69 -2.28
0.0046 0.0694
<.0001
<.0001
0.1033 0.0916
0.0231 Introduction to Longitudinal Data Analysis
498
In case the procedure NLMIXED is used, more drastic changes are needed:
·	proc nlmixed data=m.gsa2 qpoints=20;
title ’Analgesic, PROC NLMIXED, ordinal, adaptive, q=20’;
parms int1=-1.5585 int2=1.0292 int3=3.8916 int4=6.2144
beta1=0.5410 beta2=-0.1123 beta3=0.3173 d=2.1082;
eta = beta1time + beta2timetime + beta3pca0 + b1;
if gsa=1 then z = 1/(1+exp(-(int1-eta)));
else if gsa=2 then z = 1/(1+exp(-(int2-eta))) - 1/(1+exp(-(int1-eta)));
else if gsa=3 then z = 1/(1+exp(-(int3-eta))) - 1/(1+exp(-(int2-eta)));
else if gsa=4 then z = 1/(1+exp(-(int4-eta))) - 1/(1+exp(-(int3-eta)));
else z = 1 - 1/(1+exp(-(int4-eta)));
if z > 1e-8 then ll = log(z);
else ll = -1e100;
model gsa  general(ll);
random b1  normal(0,dd) subject=patid;
estimate ’var(d)’ dd;
run;
Introduction to Longitudinal Data Analysis
499
Now, the general likelihood is used: a fully user-deﬁned likelihood function.
The probabilities are obtained as diﬀerences between cumulative probabilities:
P (Yij = k) = P (Yij <= k)
P (Yij <= k
1),
−
−
(k = 1, . . . 5)
with
·	P (Yij <= 0) = 0
·	P (Yij <= 5) = 1
η is the part of the linear predictor excluding the intercept
§	Introduction to Longitudinal Data Analysis
500
Selected output:
·	Parameter Estimates
Parameter Estimate
int1
int2
int3
int4
beta1
beta2
beta3
d
-1.5585
1.0292 3.8916
6.2144 0.5410
-0.1123
0.3173 2.1082
Standard
Error
0.5481 0.5442
0.5624 0.5990
0.3078 0.06187
0.1386 0.1412
DF t Value Pr > |t| Alpha
Lower
Upper Gradient
394
394
394
394
394
394
394
394
-2.84
1.89 6.92
10.37 1.76
-1.82
2.29 14.94
0.0047 0.0593
<.0001
<.0001
0.0796 0.0702
0.0226 <.0001
-2.6360
0.05 0.05 -0.04063
2.7860 0.05
5.0368 0.05
0.05 -0.06421
0.05 0.05
0.05 -0.4810 0.000235
2.0991 -0.00004
4.9973 -0.00017
7.3920 -0.00004
1.1462 -0.00008
-0.2340 0.009311 0.000019
0.5898 0.000013
0.04475 2.3858 0.000331
1.8307 Additional Estimates
Label
Estimate
Standard
Error
DF
t Value
Pr > |t|
Alpha
Lower
Upper
var(d)
4.4447 0.5952
394
7.47 <.0001
0.05 3.2746
5.6148 Introduction to Longitudinal Data Analysis
501
28.4 Analysis of the Analgesic Trial
Three approaches:
·	Logistic regression
·	GEE
·	GLMM
For GEE: (model based standard errors; empirically corrected standard errors)
MQL performs again rather poorly
§	Introduction to Longitudinal Data Analysis
502
Marginal models
Eﬀect
Parameter
OLR
GEE
Intercept 1
Intercept 2
Intercept 3
Intercept 4
Time
Time2
Basel. PCA
α1
α2
α3
α4
β2
β3
β4
-1.00(0.34)
-1.00(0.34;0.35)
0.52(0.34)
0.52(0.34;0.36)
2.32(0.35)
2.32(0.34;0.37)
4.05(0.38)
4.05(0.37;0.39)
-0.20(0.27)
-0.20(0.27;0.20)
0.05(0.05)
0.05(0.05;0.04)
-0.21(0.06)
-0.21(0.06;0.09)
Introduction to Longitudinal Data Analysis
503
Eﬀect
Parameter
MQL
PQL
N.Int.
Random-eﬀects models
Intercept 1
Intercept 2
Intercept 3
Intercept 4
Time
Time2
Basel. PCA
Rand. int s.d.
Rand. int var.
α1
α2
α3
α4
β2
β3
β4
τ
τ 2
-0.93(0.40)
-1.44(0.50)
-1.56(0.55)
0.60(0.39)
0.91(0.50)
1.03(0.54)
2.39(0.40)
3.47(0.51)
3.89(0.56)
4.13(0.42)
5.63(0.54)
6.21(0.60)
-0.30(0.28)
-0.48(0.30)
0.54(0.31)
0.06(0.06)
0.10(0.06)
-0.11(0.06)
-0.21(0.09)
-0.28(0.12)
0.32(0.14)
1.06(0.08)
1.88(0.11)
2.11(0.14)
1.13(0.16)
3.53(0.42)
4.44(0.60)
Introduction to Longitudinal Data Analysis
504
Chapter 29
Count Data: The Epilepsy Study
·	The epilepsy data
·	Poisson regression
·	Generalized estimating equations
·	Generalized linear mixed models
·	Overview of analyses of the epilepsy study
·	Marginalization of the GLMM
Introduction to Longitudinal Data Analysis
505
29.1 The Epilepsy Data
Consider the epilepsy data:
·	Introduction to Longitudinal Data Analysis
506
·	We want to test for a treatment eﬀect on number of seizures, correcting for the
average number of seizures during the 12-week baseline phase, prior to the
treatment.
·	The response considered now is the total number of seizures a patient
experienced, i.e., the sum of all weekly measurements.
·	Let Yi now be the total number of seizures for subject i:
Yi =
ni
Xi=1
Yij
where Yij was the original (longitudinally measured) weekly outcome.
Introduction to Longitudinal Data Analysis
507
Histogram:
o	As these sums are not taken over an equal number of visits for all subjects, the
above histogram is not a ‘fair’ one as it does not account for diﬀerences in ni for
this.
Introduction to Longitudinal Data Analysis
508
o	We will therefore use the following Poisson model:
Yi ∼
Poisson(λi)
ln(λi/ni) = xi0β
Note that the regression model is equivalent to
λi = ni exp(xi0β) = exp(xi0β + ln ni)
·	Since ni is the number of weeks for which the number of seizures was recorded for
subject i, exp(xi0β) is the average number of seizures per week.
o	ln ni is called an oﬀset in the above model.
In our application, the covariates in xi are the treatment as well as the baseline
seizure rate.
Introduction to Longitudinal Data Analysis
509
·	SAS statements for the calculation of outcome, oﬀset, and for ﬁtting the Poisson
model:
proc sort data=test;
by id studyweek;
run;
proc means data=test sum n nmiss;
var nseizw;
by id;
output out=result
n=n
nmiss=nmiss
sum=sum;
run;
data result;
set result;
offset=log(n-nmiss);
keep id offset sum;
run;
data first;
set test;
by id;
if first.id;
keep id bserate trt;
run;
data result;
merge result first;
by id;
run;
proc genmod data=result;
model sum=bserate trt
/ dist=poisson offset=offset;
run;
Introduction to Longitudinal Data Analysis
510
o	The treatment variable trt is coded as 0 for placebo and 1 for treated
Output from the GENMOD procedure:
Analysis Of Parameter Estimates
Parameter DF Estimate
Standard
Error
Wald 95%
Confidence Limits
Chi-
Square
Intercept
bserate
trt
Scale
1
1
1
0
0.8710 0.0172
-0.4987
1.0000 0.0218
0.0002 0.0341
0.0000 0.8283
0.0167 -0.5655
1.0000 0.9138 1596.16
0.0177 4826.14
214.18 -0.4319
1.0000 •
We obtain a highly signiﬁcant reduction in the average number of seizures in the
treated group, in comparison to the placebo group.
Introduction to Longitudinal Data Analysis
511
·	A more general model would allow the treatment eﬀect to depend on the baseline
average number of seizures:
proc genmod data=result;
model sum=bserate trt bserate*trt
/ dist=poisson offset=offset;
run;
Relevant part of the output:
·	Analysis Of Parameter Estimates
Parameter
DF Estimate
Standard
Error
Wald 95%
Confidence Limits
Chi-
Square
Intercept
bserate
trt
bserate*trt
Scale
1
1
1
1
0
0.2107 0.0450
0.2938 -0.0295
1.0000 0.0353
0.0009 0.0454
0.0010 0.0000
0.1415 0.0432
0.2047 -0.0314
1.0000 0.2799
35.60 0.0469 2286.94
41.81 0.3829
911.43 -0.0276
1.0000 Introduction to Longitudinal Data Analysis
512
§	We get a signiﬁcant interaction.
In order to explore the nature of this interaction, we estimate the treatment eﬀect
when the baseline average number of seizures equals 6, 10.5, as well as 21
(quartiles).
This is possible via inclusion of estimate statements:
proc genmod data=result;
model sum=bserate trt bserate*trt
/ dist=poisson offset=offset;
estimate ’trt, bserate=6’ trt 1 bseratetrt 6;
estimate ’trt, bserate=10.5’ trt 1 bseratetrt 10.5;
estimate ’trt, bserate=21’ trt 1 bserate*trt 21;
run;
Introduction to Longitudinal Data Analysis
513
Additional output:
·	Contrast Estimate Results
Label
Estimate
trt, bserate=6
trt, bserate=10.5
trt, bserate=21
0.1167 -0.0161
-0.3260
Standard
Error
0.0415 0.0388
0.0340 Alpha
0.05 0.05
0.05 Label
Confidence Limits
trt, bserate=6
trt, bserate=10.5
trt, bserate=21
0.0355 -0.0921
-0.3926
0.1980 0.0600
-0.2593
Chi-
Square
7.93 0.17
91.86 Pr > ChiSq
0.0049 0.6786
<.0001
·	On average, there are more seizures in the treatment group when there are few
seizures at baseline. The opposite is true for patients with many seizures at
baseline.
Introduction to Longitudinal Data Analysis
514
29.2 Generalized Estimating Equations
·	Poisson regression models will be used to describe the marginal distributions, i.e.,
the distribution of the outcome at each time point separately:
Yij = Poisson(λij)
log(λij) = β0 + β1Ti + β2tij + β3Titij
Notation:
·	Ti: treatment indicator for subject i
·	tij: time point at which jth measurement is taken for ith subject
Note that, again, the randomization would allow to set β1 equal to 0.
o	Introduction to Longitudinal Data Analysis
515
·	More complex mean models can again be considered (e.g. including polynomial
time eﬀects, or including covariates).
·	As the response is now the number of seizures during a ﬁxed period of one week,
we do not need to include an oﬀset, as was the case in the GLM ﬁtted previously
to the epilepsy data, not in the context of repeated measurements.
·	Given the long observation period, an unstructured working correlation would
require estimation of many correlation parameters.
·	Further, the long observation period makes the assumption of an exchangeable
correlation structure quite unrealistic.
·	We therefore use the AR(1) working correlation structure, which makes sense
since we have equally spaced time points at which measurements have been taken.
Introduction to Longitudinal Data Analysis
516
o	SAS code:
proc genmod data=test;
class id timeclss;
model nseizw = trt time trt*time / dist=poisson;
repeated subject=id / withinsubject=timeclss type=AR(1) corrw modelse;
run;
Relevant SAS output:
Working Correlation Matrix
Col1
Col2
Col3
.....
Col26
Col27
Row1
Row2
Row3
1.0000 0.5946
0.3535 0.5946
1.0000 0.5946
0.3535 0.5946
1.0000 ..... ......
..... ......
..... ......
......
......
......
....
......
......
......
..... ......
......
Row26
Row27
......
......
......
......
......
......
..... 1.0000
..... 0.5946
0.5946 1.0000
Introduction to Longitudinal Data Analysis
517
Analysis Of GEE Parameter Estimates
Empirical Standard Error Estimates
Parameter Estimate
Standard
Error
95% Confidence
Limits
Z Pr > |Z|
Intercept
trt
time
trt*time
1.2259 0.1681
-0.0071
-0.0183
0.1778 0.8774
0.2785 -0.3777
0.0229 -0.0519
0.0279 -0.0730
1.5743 0.7138
0.0378 0.0364
6.90 0.60
-0.31
-0.66
<.0001
0.5461 0.7574
0.5124 Analysis Of GEE Parameter Estimates
Model-Based Standard Error Estimates
Parameter Estimate
Standard
Error
95% Confidence
Limits
Z Pr > |Z|
Intercept
trt
time
trt*time
1.2259 0.1681
-0.0071
-0.0183
0.7655 0.2349
0.3197 -0.4585
0.0230 -0.0521
0.0310 -0.0790
1.6862 0.7947
0.0380 0.0425
5.22 0.53
-0.31
-0.59
<.0001
0.5991 0.7585
0.5553 Introduction to Longitudinal Data Analysis
518
o	The AR(1) correlation coeﬃcient is estimated to be equal to 0.5946
There is no diﬀerence in average evolution between both treatment groups
(p = 0.5124).
·	Note also the huge discrepancies between the results for the initial parameter
estimates and the ﬁnal results based on the GEE analysis.
Introduction to Longitudinal Data Analysis
519
29.3 Random-eﬀects Model
·	Conditionally on a random intercept bi, Poisson regression models will be used to
describe the marginal distributions, i.e., the distribution of the outcome at each
time point separately:
Yij = Poisson(λij)
log(λij) = β0 + bi + β1Ti + β2tij + β3Titij
·	Notation:
·	Ti: treatment indicator for subject i
·	tij: time point at which jth measurement is taken for ith subject
·	Similar as in our GEE analysis, we do not need to include an oﬀset, because the
response is now the number of seizures during a ﬁxed period of one week.
Introduction to Longitudinal Data Analysis
520
·	Two equivalent SAS programs:
proc nlmixed data=test;
parms int0=0.5 slope0=-0.1 int1=1 slope1=0.1 sigma=1;
if (trt = 0) then eta = int0 + b + slope0*time;
else if (trt = 1) then eta = int1 + b + slope1*time;
lambda = exp(eta);
model nseizw  poisson(lambda);
random b  normal(0,sigma**2) subject = id;
estimate ’difference in slope’ slope1-slope0;
run;
proc nlmixed data=test;
parms int0=0.5 slope0=-0.1 int1=1 slope1=0.1 sigma=1;
eta = (1-trt)int0 + trtint1 + b
·	(1-trt)slope0time + trtslope1time;
lambda = exp(eta);
model nseizw  poisson(lambda);
random b  normal(0,sigma**2) subject = id;
estimate ’difference in slope’ slope1-slope0;
run;
Introduction to Longitudinal Data Analysis
521
·	As in the MIXED procedure, CONTRAST and ESTIMATE statements can be
speciﬁed as well. However, under PROC NLMIXED, one is no longer restricted to
linear functions of the parameters in the mean structure only.
·	For example, estimation of the ratio of both slopes, as well as of the variance of
the random intercepts is achieved by adding the following ESTIMATE statements:
estimate ’ratio of slopes’ slope1/slope0;
estimate ’variance RIs’ sigma**2;
·	Inference for such functions of parameters is based on the so-called ‘delta-method’:
·	Let ψ be the vector of all parameters in the marginal model.
·	Let
ψ be the MLE of ψ
d
·	ψ is asymptotically normally distributed with mean ψ and covariance matrix
var(
ψ) (inverse Fisher information matrix).
d
d
Introduction to Longitudinal Data Analysis
522
·	The ‘delta-method’ then implies that any function F (
ψ) of
ψ is asymptotically
normally distributed with mean F (ψ) and covariance matrix equal to
d
d
var(F (
ψ)) =
d
∂F (ψ)
∂ψ0
var(
ψ)
d
∂F 0(ψ)
∂ψ
·	Hence, a Wald-type test can be constructed, replacing the parameters in
var(F (
ψ)) by their estimates
d
Relevant SAS output:
·	Parameter Estimates
Parameter Estimate
Standard
Error
DF t Value Pr > |t| Alpha
Lower
Upper Gradient
int0
slope0
int1
slope1
sigma
0.8180 0.1675
-0.01429 0.004404
0.1699 -0.01200 0.004318
0.08556 0.6478
1.0742 88
88
88
88
88
4.88 -3.24
3.81 -2.78
12.55 <.0001
0.0017 0.0003
0.0067 <.0001
0.4852 1.1509 0.006008
0.05 0.05 -0.02304 -0.00554 0.022641
0.05 0.9855 0.010749
0.05 -0.02058 -0.00342 -0.04858
1.2442 0.009566
0.05 0.3101
0.9042 Introduction to Longitudinal Data Analysis
523
Additional Estimates
Label
Estimate
Standard
Error
DF t Value Pr > |t| Alpha
Lower
Upper
difference in slope 0.002287 0.006167
0.3979 ratio of slopes
0.1838 variance RIs
0.8399 1.1539
88
88
88
0.37 2.11
6.28 0.7116
0.0376 <.0001
0.05 -0.00997
0.04923 0.05
0.7886 0.05
0.01454 1.6306
1.5192 •
The number of quadrature points was not speciﬁed, and therefore was selected
adaptively, and set equal to only one.
·	In order to check whether Q = 1 is suﬃcient, we reﬁtted the model, prespecifying
Q = 20. This produced essentially the same output.
Introduction to Longitudinal Data Analysis
524
Corresponding code for the GLIMMIX procedure is:
·	proc glimmix data=test method=RSPL;
class id trt;
model nseizw = trttime / dist=poisson solution;
random intercept time / type=UNR subject=id;
estimate ’diff slopes’ trttime 1 -1;
run;
Introduction to Longitudinal Data Analysis
525
29.4 Overview of Epilepsy Data Analyses
o	GEE analysis (empirically corrected s.e.; model based s.e.)
Eﬀect
Parameter
Estimate (s.e.)
Common intercept
Slope placebo
Slope treatment
β0
β1
β2
1.3140 (0.1435; 0.1601)
0.0142 (0.0234; 0.0185)
0.0192 (0.0178; 0.0174)
−
−
Various GLMM analyses:
·	MQL
·	PQL
·	Laplace
·	Gaussian quadrature
Introduction to Longitudinal Data Analysis
526
Eﬀect
Parameter
Estimate (s.e.)
Estimate (s.e.)
MQL
PQL
Common intercept
Slope placebo
Slope treatment
Variance of intercepts
Variance of slopes
Correlation rand.eﬀ.
β0
β1
β2
d11
d22
ρ
1.3525 (0.1492)
0.8079 (0.1261)
0.0180 (0.0144)
−
0.0151 (0.0144)
−
0.0242 (0.0094)
−
0.0191 (0.0094)
−
1.9017 (0.2986)
1.2510 (0.2155)
0.0084 (0.0014)
0.0024 (0.0006)
0.3268 (0.1039)
−
0.3394 (0.1294)
−
Laplace
QUAD
Eﬀect
Parameter
Estimate (s.e.)
Estimate (s.e.)
Common intercept
Slope placebo
Slope treatment
Variance of intercepts
Variance of slopes
Correlation rand.eﬀ.
β0
β1
β2
d11
d22
ρ
0.7740 (0.1291)
0.7739 (0.1293)
0.0244 (0.0096)
−
0.0193 (0.0096)
−
0.0245 (0.0096)
−
0.0193 (0.0097)
−
1.2814 (0.2220)
1.2859 (0.2231)
0.0024 (0.0006)
0.0024 (0.0006)
0.3347 (0.1317)
−
0.3349 (0.1318)
−
Introduction to Longitudinal Data Analysis
527
29.5 Marginalization of the Random-eﬀects Model
·	Regression coeﬃcients in GLMM need to be interpreted conditionally on the
random eﬀects bi.
Additional computations are needed for the population-averaged evolutions.
The marginal expectation of Yij measured at tij in the placebo group is
E[Yij] = E[E[Yij|
bi]]
= E [exp[(β0 + bi1) + (β1 + bi2)tij]]
= exp[β0 + β1tij]
Calculations can be done using numerical integration or numerical averaging.
SAS code and computation: Molenberghs and Verbeke (2005, pp. 343–344)
·	Introduction to Longitudinal Data Analysis
528
6
Marginal evolutions (GEE)
Marginal evolutions (integrated GLMM)
Evolutions average subjects (bi = 0)
patients & marginal evolution (bold)
Sampled predicted proﬁles for 20 placebo
Introduction to Longitudinal Data Analysis
529
o	Curvature diﬀerent in GEE and GLMM
Ordering of treatment groups diﬀerent in GEE and GLMM (although none
signiﬁcant)
·	Watch out for the eﬀects of missingness: many patients leave the study after
week 16
·	The evolution of an ‘average’ patient is completely diﬀerent from the
population-averaged evolution
Introduction to Longitudinal Data Analysis
530
Part V
Incomplete Data
Introduction to Longitudinal Data Analysis
531
Chapter 30
Setting The Scene
·	Orthodontic growth data
·	Depression trial
·	Age-related macular degeneration trial
·	Notation
·	Taxonomy
Introduction to Longitudinal Data Analysis
532
30.1 Growth Data
Taken from Potthoﬀ and Roy, Biometrika (1964)
Research question:
o	Is dental growth related to gender ?
·	The distance from the center of the pituitary to the maxillary ﬁssure was recorded
at ages 8, 10, 12, and 14, for 11 girls and 16 boys
Introduction to Longitudinal Data Analysis
533
·	Individual proﬁles:
·	Much variability between girls / boys
·	Considerable variability within girls / boys
·	Fixed number of measurements per subject
·	Measurements taken at ﬁxed time points
Introduction to Longitudinal Data Analysis
534
30.2 The Depression Trial
Clinical trial: experimental drug versus standard drug
170 patients
Response: change versus baseline in HAM D17 score
5 post-baseline visits: 4–8
·	0
2
0
1
0
0
1
0
2
e
g
n
a
h
C
Standard Drug
Experimental Drug
e
g
n
a
h
C
2
4
6
8
0
1
4
5
6
=Visit
7
8
4
5
6
Visit
7
8
Introduction to Longitudinal Data Analysis
535
30.3 Age-related Macular Degeneration Trial
Pharmacological Therapy for Macular Degeneration Study Group (1997)
An occular pressure disease which makes patients progressively lose vision
240 patients enrolled in a multi-center trial (190 completers)
Treatment: Interferon-α (6 million units) versus placebo
Visits: baseline and follow-up at 4, 12, 24, and 52 weeks
Continuous outcome: visual acuity: # letters correctly read on a vision chart
Binary outcome: visual acuity versus baseline
0 or
0
≤
≥
·	Introduction to Longitudinal Data Analysis
536
Missingness:
·	Measurement occasion
4 wks
12 wks
24 wks
52 wks Number
%
Completers
O
O
188
78.33 O
O
O
M
M
Dropouts
O
M
M
M
M
M
M
M
Non-monotone missingness
O
M
O
O
M
M
O
M
O
O
O
M
4
1
2
1
24
10.00 8
6
6
3.33 2.50
2.50 1.67
0.42 0.83
0.42 O
O
O
O
M
O
O
M
M
Introduction to Longitudinal Data Analysis
537
CRF
TRT
VISUAL0
VISUAL4
VISUAL12
VISUAL24
VISUAL52
lesion
1002
1003
1006
1007
1010
1110
1111
1112
1115
1803
1805
...
4
4
1
1
4
4
1
1
4
1
4
59
65
40
67
70
59
64
39
59
49
58
55
70
40
64
·	53
68
37
58
51
50
45
65
37
64
·	52
74
43
49
71
o	65
17
64
·	53
72
37
54
71
o	55
·	68
·	42
65
37
58
o	3
1
4
2
1
3
1
3
2
1
1
Introduction to Longitudinal Data Analysis
538
30.4 Incomplete Longitudinal Data
Introduction to Longitudinal Data Analysis
539
30.5 Scientiﬁc Question
In terms of
entire longitudinal proﬁle
In terms of
last planned measurement
In terms of
last observed measurement
§	Introduction to Longitudinal Data Analysis
540
30.6 Notation
Subject i at occasion (time) j = 1, . . . , ni
Measurement Yij
Missingness indicator Rij =
1
0
if Yij is observed,
otherwise.



Group Yij into a vector
Y i = (Yi1, . . . , Yini)0 = (Y o
i , Y m
i )
Y o
i
Y m
i
contains Yij for which Rij = 1,
contains Yij for which Rij = 0.



Group Rij into a vector Ri = (Ri1, . . . , Rini)0
Di: time of dropout: Di = 1 +
ni
j=1 Rij
P
§	Introduction to Longitudinal Data Analysis
541
30.7 Direct Likelihood/Bayesian Inference: Ignorability
MAR : f (Y o
i |
X i, θ) f (ri|
X i, Y o
i , ψ)
Mechanism is MAR
θ and ψ distinct
Interest in θ
(Use observed information matrix)



=
⇒
Lik./Bayes inference valid
Outcome type
Modeling strategy
Gaussian
Linear mixed model
Software
SAS MIXED
Non-Gaussian Gen./Non-linear mixed model SAS GLIMMIX, NLMIXED
Introduction to Longitudinal Data Analysis
542
30.8 Rubin, 1976
Ignorability: Rubin (Biometrika, 1976): 35 years ago!
Little and Rubin (1976, 2002)
Why did it take so long?
§	Introduction to Longitudinal Data Analysis
543
30.9 A Vicious Triangle
Industry
%.
&-
Academe
←−
−→
Regulatory
Academe: The R2 principle
Regulatory: Rigid procedures
scientiﬁc developments
←→
Industry: We cannot / do not want to apply new methods
§	Introduction to Longitudinal Data Analysis
544
30.10 Terminology & Confusion
o	The Ministry of Disinformation:
All directions
←−
Other directions
−→
MCAR, MAR, MNAR: “What do the terms mean?”
MAR, random dropout, informative missingness, ignorable,
censoring,. . .
Dropout from the study, dropout from treatment, lost to follow up,. . .
“Under MAR patients dropping out and patients not dropping out are
similar.”
Introduction to Longitudinal Data Analysis
545
30.11 A Virtuous Triangle
Industry
%.
&-
Academe
←−
−→
Regulatory
FDA/Industry Workshops
DIA/EMA Meetings
The NAS Experience
§	Introduction to Longitudinal Data Analysis
546
30.12 The NAS Experience: A Wholesome Product
§	FDA
NAS
−→
−→
the working group
Composition
Encompassing:
·	terminology/taxonomy/concepts
·	prevention
·	treatment
Introduction to Longitudinal Data Analysis
547
30.13 Taxonomy
Missingness pattern: complete — monotone — non-monotone
Dropout pattern: complete — dropout — intermittent
Model framework: SEM — PMM — SPM
Missingness mechanism: MCAR — MAR — MNAR
Ignorability: ignorable — non-ignorable
Inference paradigm: frequentist — likelihood — Bayes
§	Introduction to Longitudinal Data Analysis
548
30.14 The NAS Panel
Name
Rod Little
Ralph D’Agostino
Kay Dickerson
Scott Emerson
John Farrar
Constantine Frangakis
Joseph Hogan
Geert Molenberghs
Susan Murphy
James Neaton
Andrea Rotnitzky
Dan Scharfstein
Joseph Shih
Jay Siegel
Hal Stern
Specialty
biostat
biostat
epi
biostat
epi
biostat
biostat
biostat
stat
biostat
stat
biostat
biostat
biostast
stat
Aﬃliation
U Michigan
Boston U
Johns Hopkins
U Washington
U Penn
Johns Hopkins
Brown U
U Hasselt & K.U.Leuven
U Michigan
U Minnesota
Buenos Aires & Harvard
Johns Hopkins
New Jersey SPH
J&J
UC at Irvine
Introduction to Longitudinal Data Analysis
549
30.15 Modeling Frameworks & Missing Data Mechanisms
f (yi, ri|
Xi, θ, ψ)
Selection Models: f (yi|
Xi, θ) f (ri|
Xi, yo
i , ym
i , ψ)
MCAR
−→
MAR
−→
MNAR
f (ri|
Xi, ψ)
f (ri|
Xi, yo
i , ψ)
f (ri|
Xi, yo
i , ym
i , ψ)
Pattern-mixture Models: f (yi|
Shared-parameter Models: f (yi|
Xi, ri, θ) f (ri|
Xi, bi, θ) f (ri|
Xi, ψ)
Xi, bi, ψ)
Introduction to Longitudinal Data Analysis
550
30.16 Frameworks and Their Methods
f (yi, ri|
Xi, θ, ψ)
Selection Models: f (yi|
Xi, θ) f (ri|
Xi, yo
i , ym
i , ψ)
MCAR/simple
−→
MAR
−→
MNAR
CC?
LOCF?
single imputation?
...
direct likelihood!
direct Bayesian!
multiple imputation (MI)!
IPW
⊃
W-GEE!
joint model!?
sensitivity analysis?!
d.l. + IPW = double robustness! (consensus)
Introduction to Longitudinal Data Analysis
551
30.17 Frameworks and Their Methods: Start
f (yi, ri|
Xi, θ, ψ)
Selection Models: f (yi|
Xi, θ) f (ri|
Xi, yo
i , ym
i , ψ)
MCAR/simple
−→
MAR
−→
MNAR
direct likelihood!
direct Bayesian!
multiple imputation (MI)!
IPW
⊃
W-GEE!
d.l. + IPW = double robustness!
Introduction to Longitudinal Data Analysis
552
30.18 Frameworks and Their Methods: Next
f (yi, ri|
Selection Models: f (yi|
Xi, θ, ψ)
Xi, θ) f (ri|
Xi, yo
i , ym
i , ψ)
MCAR/simple
−→
MAR
−→
MNAR
joint model!?
sensitivity analysis!
PMM
MI (MGK, J&J)
local inﬂuence
interval ignorance
IPW based
Introduction to Longitudinal Data Analysis
553
30.19 Overview
MCAR/simple CC
LOCF
biased
ineﬃcient
not simpler than MAR methods
MAR
direct likelihood easy to conduct
direct Bayes
Gaussian & non-Gaussian
weighted GEE
MI
MNAR
variety of methods
strong, untestable assumptions
most useful in sensitivity analysis
Introduction to Longitudinal Data Analysis
554
30.20 Selection Models versus Pattern-mixture Models:
A Paradox!?
Glynn, Laird and Rubin (1986)
Two measurements (Y1, Y2)
Y1 always observed.
Y2 observed (R = 1) or missing (R = 0).
·	Introduction to Longitudinal Data Analysis
555
Selection model versus pattern-mixture model
f (y1, y2)g(r = 1
|
f (y1, y2)g(r = 0
|
y1, y2) = f1(y1, y2)p(r = 1)
y1, y2) = f0(y1, y2)p(r = 0)
or
f (y1, y2)g(y1, y2) = f1(y1, y2)p
f (y1, y2)[1
−
g(y1, y2)] = f0(y1, y2)[1
p]
−
of which the ratio yields:
f0(y1, y2) =
1
g(y1, y2)
−
g(y1, y2)
·	1
p
−
p
.f1(y1, y2)
The right hand side is identiﬁable
the left hand side is not. . .
←→
o	Introduction to Longitudinal Data Analysis
556
30.21 Illustration
Data:
20 30
10 40
75
25
LOCF:
20 30
75
0
10 40
0 25
CC:
20 30
10 40
0
0
0
0
MAR:
20 30
30 45
10 40
5 20
=
⇒
=
⇒
=
⇒
95 30
10 65
20 30
10 40
50 75
15 60
=
⇒
=
⇒
=
⇒
θ = 95
200 = 0.475 [0.406; 0.544] (biased & too narrow)
b
θ = 20
100 = 0.200 [0.122; 0.278] (biased & too wide)
b
θ = 50
200 = 0.250 [0.163; 0.337]
b
Introduction to Longitudinal Data Analysis
557
Chapter 31
Proper Analysis of Incomplete Data
·	Simple methods
·	Bias for LOCF and CC
·	Direct likelihood inference
·	Weighted generalized estimating equations
Introduction to Longitudinal Data Analysis
558
31.1 Incomplete Longitudinal Data
Introduction to Longitudinal Data Analysis
559
Data and Modeling Strategies
Introduction to Longitudinal Data Analysis
560
Modeling Strategies
Introduction to Longitudinal Data Analysis
561
31.2 Simple Methods
MCAR
Complete case analysis:
Last observation carried forward:
delete incomplete subjects
Standard statistical software
Loss of information
Impact on precision and power
Missingness
= MCAR
bias
⇒
⇒
·	impute missing values
Standard statistical software
Increase of information
Constant proﬁle after dropout:
unrealistic
Usually bias
⇒
·	Introduction to Longitudinal Data Analysis
562
6
Quantifying the Bias
Dropouts
tij = 0
Probability
p0
Treatment indicator
Ti = 0, 1
Completers
tij = 0, 1
Probability
1
−
Treatment indicator
p0 = p1
Ti = 0, 1
E(Yij) = β0 + β1Ti + β2tij + β3Titij
E(Yij) = γ0+γ1Ti+γ2tij+γ3Titij
CC
MCAR
MAR
σ[(1
−
−
0
p1)(β0 + β1 −
p0)(β0 −
(1
−
−
γ1)
γ0 −
γ0)]
LOCF
p0)β2 −
(p1 −
−
p1(γ0 + γ1 + γ2 + γ3) + (1
(1
p1)β3
p0(γ0 + γ2)
−
−
p1)(β0 + β1 −
−
(1
−
γ0 −
σ[(1
−
−
p0)β0 −
(1
γ1)
−
p1)(β0 + β1)
γ1 −
−
γ3
p0)(β0 −
γ0)]
Introduction to Longitudinal Data Analysis
563
31.3 Ignorability
§	Let us decide to use likelihood based estimation.
The full data likelihood contribution for subject i:
L∗(θ, ψ
Y i, Di)
|
f (Y i, Di|
∝
θ, ψ).
Base inference on the observed data:
with
L(θ, ψ
Y i, Di)
|
f (Y o
i , Di|
∝
θ, ψ)
f (Y o
i , Di|
θ, ψ) =
=
Z
Z
f (Y i, Di|
i , Y m
f (Y o
i |
θ, ψ)dY m
i
Y o
θ)f (Di|
i , Y m
i , ψ)dY m
i .
Introduction to Longitudinal Data Analysis
564
Under a MAR process:
f (Y o
i , Di|
θ, ψ) =
Z
f (Y o
i , Y m
i |
θ)f (Di|
θ)f (Di|
Y o
i , ψ),
= f (Y o
i |
Y o
i , ψ)dY m
i
The likelihood factorizes into two components.
o	Introduction to Longitudinal Data Analysis
565
31.4 Ignorability: Summary
Likelihood/Bayesian + MAR
&
Frequentist + MCAR
Introduction to Longitudinal Data Analysis
566
31.5 Direct Likelihood/Bayesian Inference: Ignorability
MAR : f (Y o
i |
X i, θ) f (ri|
X i, Y o
i , ψ)
Mechanism is MAR
θ and ψ distinct
Interest in θ
(Use observed information matrix)



=
⇒
Lik./Bayes inference valid
Outcome type
Modeling strategy
Software
Gaussian
Non-Gaussian
Linear mixed model
Gen./Non-linear mixed model
SAS MIXED
SAS GLIMMIX, NLMIXED
Introduction to Longitudinal Data Analysis
567
31.6 Original, Complete Orthodontic Growth Data
Mean
unstructured
= slopes
= slopes
= slopes
Covar
unstructured
unstructured
unstructured
CS
par
18
14
13
6
1
2
3
7
Introduction to Longitudinal Data Analysis
568
6
6
31.7 Trimmed Growth Data: Simple Methods
Method
Complete case
LOCF
Unconditional mean
Conditional mean
Model
7a
2a
7a
1
Mean
= slopes
quadratic
= slopes
unstructured
Covar
CS
unstructured
CS
unstructured
par
5
16
5
18
distorting
Introduction to Longitudinal Data Analysis
569
31.8 Trimmed Growth Data: Direct Likelihood
Mean Covar # par
7
= slopes CS
6
Introduction to Longitudinal Data Analysis
570
6
31.9 Growth Data: Comparison of Analyses
·	Model
·	Unstructured group by time mean
·	Unstructured covariance matrix
o	Data
·	Complete cases
·	LOCF imputed data
·	All available data
Analysis methods
·	Direct likelihood
ML
REML
∗
∗
·	MANOVA
·	ANOVA per time point
Introduction to Longitudinal Data Analysis
571
Principle Method
Boys at Age 8 Boys at Age 10
Original
Direct likelihood, ML
22.88 (0.56)
23.81 (0.49)
1964
Direct likelihood, REML
ANOVA per time point
≡
MANOVA
22.88 (0.58)
23.81 (0.51)
22.88 (0.61)
23.81 (0.53)
Direct Lik. Direct likelihood, ML
22.88 (0.56)
23.17 (0.68)
1987
Direct likelihood, REML
22.88 (0.58)
23.17 (0.71)
CC
1987
LOCF
1987
MANOVA
ANOVA per time point
Direct likelihood, ML
Direct likelihood, REML
ANOVA per time point
Direct likelihood, ML
Direct likelihood, REML
ANOVA per time point
≡
≡
24.00 (0.48)
24.14 (0.66)
22.88 (0.61)
24.14 (0.74)
24.00 (0.45)
24.14 (0.62)
MANOVA
24.00 (0.48)
24.14 (0.66)
24.00 (0.51)
24.14 (0.74)
22.88 (0.56)
22.97 (0.65)
MANOVA
22.88 (0.58)
22.97 (0.68)
22.88 (0.61)
22.97 (0.72)
Introduction to Longitudinal Data Analysis
572
31.9.1 Behind the Scenes
R completers
N
−
↔
R “incompleters”
Yi1
Yi2
















N
∼
µ1
µ2
























,








σ11 σ12
σ22
















o	Conditional density
Yi2|
µ1
µ2
µ2
freq. & lik.
frequentist
likelihood
yi1 ∼
1
N
µ1 =
d
1
R
1
N
µ2 =
g
µ2 =
d
N (β0 + β1yi1, σ22.1)
N
Xi=1
yi1
R
Xi=1
yi2
R
Xi=1



yi2 +
N
Xi=R+1 (cid:20)
y2 +
β1(yi1 −
d
y1)
(cid:21)



Introduction to Longitudinal Data Analysis
573
31.9.2 Growth Data: Further Comparison of Analyses
Principle
Original
Method
Boys at Age 8 Boys at Age 10
Direct likelihood, ML
22.88 (0.56)
23.81 (0.49)
Direct Lik. Direct likelihood, ML
22.88 (0.56)
23.17 (0.68)
CC
LOCF
Direct likelihood, ML
24.00 (0.45)
24.14 (0.62)
Direct likelihood, ML
22.88 (0.56)
22.97 (0.65)
Data
Mean
Covariance Boys at Age 8 Boys at Age 10
Complete
Unstructured Unstructured
Unstructured CS
Unstructured Independence
Incomplete Unstructured Unstructured
Unstructured CS
Unstructured Independence
22.88 22.88
22.88 22.88
22.88 22.88
23.81 23.81
23.81 23.17
23.52 24.14
Introduction to Longitudinal Data Analysis
574
Introduction to Longitudinal Data Analysis
575
31.9.3 Growth Data: SAS Code for Model 1
IDNR AGE SEX MEASURE
SAS code:
·	8
10
12
14
8
12
14
2
2
2
2
2
2
2
1
1
1
1
·	. .
3
3
3
·	. .
21.0 20.0
21.5 23.0
20.5 24.5
26.0 proc mixed data = growth method = ml;
class sex idnr age;
model measure = sex age*sex / s;
repeated age / type = un
subject = idnr;
run;
·	Subjects in terms of IDNR blocks
·	age ensures proper ordering of observations
within subjects!
Introduction to Longitudinal Data Analysis
576
31.9.4 Growth Data: SAS Code for Model 2
IDNR AGE SEX MEASURE
·	SAS code:
8
10
12
14
8
12
14
2
2
2
2
2
2
2
1
1
1
1
·	. .
3
3
3
·	. .
21.0 20.0
21.5 23.0
20.5 24.5
26.0 data help;
set growth;
agecat = age;
run;
proc mixed data = growth method = ml;
class sex idnr agecat;
model measure = sex age*sex / s;
repeated agecat / type = un
subject = idnr;
run;
·	Time ordering variable needs to be cate-
gorical
Introduction to Longitudinal Data Analysis
577
31.10 Analysis of the Depression Trial
Complete case analysis:
%cc(data=depression, id=patient, time=visit, response=change, out={cc});
performs analysis on CC data set
⇒
LOCF analysis:
%locf(data=depression, id=patient, time=visit, response=change, out={locf});
performs analysis on LOCF data
⇒
Direct-likelihood analysis:
⇒
ﬁt linear mixed model to incomplete data
§	Introduction to Longitudinal Data Analysis
578
Treatment eﬀect at visit 8 (last follow-up measurement):
·	Method Estimate (s.e.) p-value
CC
LOCF
MAR
-1.94
-1.63
-2.38
(1.17)
0.0995 (1.08)
0.1322 (1.16)
0.0419 Observe the slightly signiﬁcant p-value under the MAR model
Introduction to Longitudinal Data Analysis
579
Chapter 32
Analysis of the ARMD Trial
Model for continuous outcomes:
·	Yij = βj1 + βj2Ti + εij
with:
·	Ti = 0 for placebo and Ti = 1 for interferon-α
·	tj (j = 1, . . . , 4) refers to the four follow-up measurements
·	unstructured variance-covariance matrix
Introduction to Longitudinal Data Analysis
580
o	Marginal mean for GEE:
logit[P (Yij = 1
|
Ti, tj)] = βj1 + βj2Ti
Model for GLMM with random interecept:
logit[P (Yij = 1
|
Ti, tj, bi)] = βj1 + bi + βj2Ti
with
·	bi ∼
N (0, τ 2)
Introduction to Longitudinal Data Analysis
581
·	Complete case analysis preparation (continuous outcome):
%cc(data=armd155,id=subject,time=time,response=diff,
out=armdcc2);
Complete case analysis preparation (discrete outcome):
%cc(data=armd111,id=subject,time=time,response=bindif,
out=armdcc);
Preparing for LOCF analysis (continuous outcome):
%locf(data=armd155,id=subject,time=time,response=diff,
out=armdlocf2);
Preparing for LOCF analysis (discrete outcome):
%locf(data=armd111,id=subject,time=time,response=bindif,
out=armdlocf);
Introduction to Longitudinal Data Analysis
582
Program for linear mixed model (PROC MIXED):
·	proc mixed data=armdcc2 method=ml;
class time treat subject;
model diff = time treat*time / noint solution ddfm=kr;
repeated time / subject=subject type=un;
run;
Introduction to Longitudinal Data Analysis
583
o	Program for GEE (PROC GENMOD):
proc genmod data=armdcc;
class time treat subject;
model bindif = time treat*time / noint dist=binomial;
repeated subject=subject / withinsubject=time type=exch modelse;
run;
Program for GEE (PROC GEE):
proc gee data=armdcc;
class time treat subject;
model bindif = time treat*time / noint dist=binomial;
repeated subject=subject / withinsubject=time type=exch modelse;
run;
Introduction to Longitudinal Data Analysis
584
Program for GLMM (PROC GLIMMIX):
·	proc glimmix data=armdcc method=gauss(q=20);
nloptions maxiter=50 technique=newrap;
class time treat subject;
model bindif = time treat*time / noint solution dist=binary;
random intercept
run;
/ subject=subject type=un g gcorr;
Introduction to Longitudinal Data Analysis
585
Program for GLMM (PROC NLMIXED):
·	data help; set armdcc;
time1=0; if time=1 then time1=1;
time2=0; if time=2 then time2=1;
time3=0; if time=3 then time3=1;
time4=0; if time=4 then time4=1;
run;
proc nlmixed data=help qpoints=20 maxiter=100 technique=newrap;
title ’CC - mixed - numerical integration’;
eta = beta11time1+beta12time2+beta13time3+beta14time4
+b
+(beta21time1+beta22time2+beta23time3+beta24time4)
*(2-treat);
p = exp(eta)/(1+exp(eta));
model bindif  binary(p);
random b  normal(0,tautau) subject=subject;
estimate ’tau^2’ tautau;
run;
Introduction to Longitudinal Data Analysis
586
Eﬀect
Parameter
CC
LOCF
direct lik.
Parameter estimates (standard errors) for linear mixed model
Intercept 4
Intercept 12
Intercept 24
Intercept 52
Treatm. eﬀ. 4
Treatm. eﬀ. 12
Treatm. eﬀ. 24
Treatm. eﬀ. 52
Treatm. eﬀ. 4
Treatm. eﬀ. 12
Treatm. eﬀ. 24
Treatm. eﬀ. 52
Treatm. eﬀ. (overall)
β11
β21
β31
β41
β12
β22
β32
β42
β12
β22
β32
β42
-3.24(0.77)
-3.48(0.77)
-3.48(0.77)
-4.66(1.14)
-5.72(1.09)
-5.85(1.11)
-8.33(1.39)
-8.34(1.30)
-9.05(1.36)
-15.13(1.73)
-14.16(1.53)
-16.21(1.67)
2.32(1.05)
2.20(1.08)
2.20(1.08)
2.35(1.55)
3.38(1.53)
3.51(1.55)
2.73(1.88)
2.41(1.83)
3.03(1.89)
4.17(2.35)
3.43(2.15)
4.86(2.31)
p-values
0.0282 0.1312
0.1491 0.0772
0.1914 0.0432
0.0287 0.1891
0.1119 0.1699
0.0435 0.0246
0.1096 0.0366
0.1234 Introduction to Longitudinal Data Analysis
587
Eﬀect
Parameter
CC
PQL
LOCF
direct lik.
Int.4
Int.12
Int.24
Int.52
Trt.4
Trt.12
Trt.24
Trt.52
R.I. s.d.
R.I. var.
Int.4
Int.12
Int.24
Int.52
Trt.4
Trt.12
Trt.24
Trt.52
R.I. s.d.
R.I. var.
β11
β21
β31
β41
β12
β22
β32
β42
τ
τ 2
β11
β21
β31
β41
β12
β22
β32
β42
τ
τ 2
-1.19(0.31)
-1.05(0.31)
-1.35(0.32)
-1.97(0.36)
0.45(0.42)
0.58(0.41)
0.55(0.42)
0.44(0.47)
1.42(0.14)
2.03(0.39)
Numerical integration
-1.05(0.28)
-1.18(0.28)
-1.30(0.28)
-1.89(0.31)
0.24(0.39)
0.68(0.38)
0.50(0.39)
0.39(0.42)
1.53(0.13)
2.34(0.39)
-1.73(0.42)
-1.53(0.41)
-1.93(0.43)
-2.74(0.48)
0.64(0.54)
0.81(0.53)
0.77(0.55)
0.60(0.59)
2.19(0.27)
4.80(1.17)
-1.63(0.39)
-1.80(0.39)
-1.96(0.40)
-2.76(0.44)
0.38(0.52)
0.98(0.52)
0.74(0.52)
0.57(0.56)
2.47(0.27)
6.08(1.32)
-1.00(0.26)
-1.19(0.28)
-1.26(0.29)
-2.02(0.35)
0.22(0.37)
0.71(0.37)
0.49(0.39)
0.46(0.46)
1.40(0.13)
1.95(0.35)
-1.50(0.36)
-1.73(0.37)
-1.83(0.39)
-2.85(0.47)
0.34(0.48)
1.00(0.49)
0.69(0.50)
0.64(0.58)
2.20(0.25)
4.83(1.11)
Introduction to Longitudinal Data Analysis
588
Chapter 33
Weighted Generalized Estimating Equations
·	General Principle
·	Analysis of the analgesic trial
·	Analysis of the ARMD trial
·	Analysis of the depression trial
Introduction to Longitudinal Data Analysis
589
33.1 General Principle
MAR and non-ignorable !
Standard GEE inference correct only under MCAR
Under MAR: weighted GEE
Robins, Rotnitzky & Zhao (JASA, 1995)
Fitzmaurice, Molenberghs & Lipsitz (JRSSB, 1995)
Decompose dropout time
Di = (Ri1, . . . , Rin) = (1, . . . , 1, 0, . . . , 0)
§	Introduction to Longitudinal Data Analysis
590
o	Weigh a contribution by inverse dropout probability
νidi ≡
P [Di = di] =
di−
1
Yk=2
(1
P [Rik = 0
|
−
Ri2 = . . . = Ri,k
1 = 1])
−
×
P [Ridi = 0
|
Ri2 = . . . = Ri,di−
1 = 1]I
T
di≤
{
}
Adjust estimating equations
N
Xi=1
1
νidi ·
∂µi
∂β0
1
V −
i
(yi
−
µi) = 0
Introduction to Longitudinal Data Analysis
591
33.2 Computing the Weights
o	Predicted values from (PROC GENMOD) output
The weights are now deﬁned at the individual measurement level:
·	At the ﬁrst occasion, the weight is w = 1
·	At other than the last ocassion, the weight is the already accumulated weight,
multiplied by 1
−
the predicted probability
·	At the last occasion within a sequence where dropout occurs the weight is
multiplied by the predicted probability
·	At the end of the process, the weight is inverted
Introduction to Longitudinal Data Analysis
592
33.3 The Analgesic Trial
single-arm trial with 530 patients recruited (491 selected for analysis)
analgesic treatment for pain caused by chronic nonmalignant disease
treatment was to be administered for 12 months
we will focus on Global Satisfaction Assessment (GSA)
GSA scale goes from 1=very good to 5=very bad
GSA was rated by each subject 4 times during the trial, at months 3, 6, 9, and 12.
§	Introduction to Longitudinal Data Analysis
593
Research questions:
·	Evolution over time
·	Relation with baseline covariates: age, sex, duration of the pain, type of pain,
disease progression, Pain Control Assessment (PCA), . . .
·	Investigation of dropout
Frequencies:
o	GSA
Month 3
Month 6
Month 9
Month 12
1
2
3
4
5
55
112
151
52
15
14.3%
29.1%
39.2%
13.5%
3.9%
Tot
385
38
84
115
51
14
302
12.6%
27.8%
38.1%
16.9%
4.6%
40
67
76
33
11
227
17.6%
29.5%
33.5%
14.5%
4.9%
30
66
97
27
3
223
13.5%
29.6%
43.5%
12.1%
1.4%
Introduction to Longitudinal Data Analysis
594
Missingness:
·	Measurement occasion
Month 3
Month 6
Month 9
Month 12
Number
%
O
O
O
O
O
O
O
O
M
M
M
M
Completers’ pattern
O
O
163
41.2 Dropout patterns
O
M
M
M
M
M
Non-monotone patterns
M
O
O
M
O
O
M
M
O
O
M
O
O
M
O
M
51
51
63
30
7
2
18
2
1
1
3
12.91 12.91
15.95 7.59
1.77 0.51
4.56 0.51
0.25 0.25
0.76 O
O
O
M
O
M
M
M
O
O
O
O
Introduction to Longitudinal Data Analysis
595
33.4 Analysis of the Analgesic Trial
A logistic regression for the dropout indicator:
·	logit[P (Di = j
Di ≥
|
j,
·
)] = ψ0 + ψ11I(GSAi,j
−
+ψ13I(GSAi,j
1 = 3) + ψ14I(GSAi,j
−
+ψ2PCA0i + ψ3PFi + ψ4GDi
with
·	GSAi,j
1 the 5-point outcome at the previous time
−
) is an indicator function
·	I(
·
·	PCA0i is pain control assessment at baseline
·	PFi is physical functioning at baseline
·	GDi is genetic disorder at baseline are used)
1 = 1) + ψ12I(GSAi,j
1 = 2)
−
1 = 4)
−
Introduction to Longitudinal Data Analysis
596
Eﬀect
Intercept
Previous GSA= 1
Previous GSA= 2
Previous GSA= 3
Previous GSA= 4
Basel. PCA
Phys. func.
Genetic disfunc.
Par. Estimate (s.e.)
ψ0
ψ11
ψ12
ψ13
ψ14
ψ2
ψ3
ψ4
-1.80 (0.49)
-1.02 (0.41)
-1.04 (0.38)
-1.34 (0.37)
-0.26 (0.38)
0.25 (0.10)
0.009 (0.004)
0.59 (0.24)
There is some evidence for MAR: P (Di = j
Di ≥
|
j) depends on previous GSA.
Furthermore: baseline PCA, physical functioning and genetic/congenital disorder.
o	Introduction to Longitudinal Data Analysis
597
GEE and WGEE:
·	logit[P (Yij = 1
|
tj, PCA0i)] = β1 + β2tj + β3t2
j + β4PCA0i
Eﬀect
Parameter
GEE
WGEE
Intercept
Time
Time2
Basel. PCA
β1
β2
β3
β4
2.95 (0.47)
2.17 (0.69)
-0.84 (0.33)
-0.44 (0.44)
0.18 (0.07)
0.12 (0.09)
-0.24 (0.10)
-0.16 (0.13)
A hint of potentially important diﬀerences between both
·	Introduction to Longitudinal Data Analysis
598
33.5 Analgesic Trial: Steps for WGEE in SAS
1.	Preparatory data manipulation:
%dropout(...)
2.	Logistic regression for weight model:
proc genmod data=gsac;
class prevgsa;
model dropout = prevgsa pca0 physfct gendis / pred dist=b;
ods output obstats=pred;
run;
3.	Conversion of predicted values to weights:
...
%dropwgt(...)
Introduction to Longitudinal Data Analysis
599
4. Weighted GEE analysis:
proc genmod data=repbin.gsaw;
scwgt wi;
class patid timecls;
model gsabin = time|time pca0 / dist=b;
repeated subject=patid / type=un corrw within=timecls;
run;
Introduction to Longitudinal Data Analysis
600
33.6 Analgesic Trial: Steps for WGEE in SAS,
Using PROC GEE
o	Available since SAS 9.4 (SAS/STAT 13.2)
Preparation:
data gsaw;
set gsaw;
by patid;
prevgsa = lag(gsa);
if first.id then prevgsa = 1;
time = time-1;
timeclss = time;
run;
Introduction to Longitudinal Data Analysis
601
Weighted GEE analysis:
·	ods graphics on;
proc gee data=gsaw plots=histogram;
class patid timecls prevgsa;
model gsabin = time|time pca0 / dist=bin;
repeated subject=patid / within=timecls corr=un;
missmodel prevgsa pca0 physfunt gendist / type=obslevel;
run;
Introduction to Longitudinal Data Analysis
602
33.7 Analysis of the ARMD Trial
·	Model for the weights:
Di ≥
logit[P (Di = j
|
j)] = ψ0 + ψ1yi,j
1 + ψ2Ti + ψ31L1i + ψ32L2i + ψ34L3i
−
+ψ41I(tj = 2) + ψ42I(tj = 3)
with
·	yi,j
1 the binary outcome at the previous time ti,j
−
1 = tj
−
1 (since time is
−
common to all subjects)
·	Ti = 1 for interferon-α and Ti = 0 for placebo
·	Lki = 1 if the patient’s eye lesion is of level k = 1, . . . , 4 (since one dummy
variable is redundant, only three are used)
) is an indicator function
·	I(
·
Introduction to Longitudinal Data Analysis
603
Results for the weights model:
·	Eﬀect
Intercept
Previous outcome
Treatment
Lesion level 1
Lesion level 2
Lesion level 3
Time 2
Time 3
Parameter Estimate (s.e.)
ψ0
ψ1
ψ2
ψ31
ψ32
ψ33
ψ41
ψ42
0.14 (0.49)
0.04 (0.38)
-0.86 (0.37)
-1.85 (0.49)
-1.91 (0.52)
-2.80 (0.72)
-1.75 (0.49)
-1.38 (0.44)
Introduction to Longitudinal Data Analysis
604
GEE:
·	with
logit[P (Yij = 1
|
Ti, tj)] = βj1 + βj2Ti
·	Ti = 0 for placebo and Ti = 1 for interferon-α
·	tj (j = 1, . . . , 4) refers to the four follow-up measurements
·	Comparison between CC, LOCF, and GEE analyses
SAS code: Molenberghs and Verbeke (2005, Section 32.5)
Results:
o	Introduction to Longitudinal Data Analysis
605
Standard GEE
Observed data
Eﬀect Par.
CC
LOCF
Unweighted
WGEE
Int.4
Int.12
Int.24
Int.52
Tr.4
Tr.12
Tr.24
Tr.52
Corr.
β11
β21
β31
β41
β12
β22
β32
β42
ρ
-1.01(0.24;0.24) -0.87(0.20;0.21) -0.87(0.21;0.21) -0.98(0.10;0.44)
-0.89(0.24;0.24) -0.97(0.21;0.21) -1.01(0.21;0.21) -1.78(0.15;0.38)
-1.13(0.25;0.25) -1.05(0.21;0.21) -1.07(0.22;0.22) -1.11(0.15;0.33)
-1.64(0.29;0.29) -1.51(0.24;0.24) -1.71(0.29;0.29) -1.72(0.25;0.39)
0.40(0.32;0.32) 0.22(0.28;0.28) 0.22(0.28;0.28) 0.80(0.15;0.67)
0.49(0.31;0.31) 0.55(0.28;0.28) 0.61(0.29;0.29) 1.87(0.19;0.61)
0.48(0.33;0.33) 0.42(0.29;0.29) 0.44(0.30;0.30) 0.73(0.20;0.52)
0.40(0.38;0.38) 0.34(0.32;0.32) 0.44(0.37;0.37) 0.74(0.31;0.52)
0.39 0.44
0.39 0.33
Introduction to Longitudinal Data Analysis
606
Linearization-based GEE
Observed data
Eﬀect Par.
CC
LOCF
Unweighted
WGEE
Int.4
Int.12
Int.24
Int.52
Tr.4
Tr.12
Tr.24
Tr.52
β11
β21
β31
β41
β12
β22
β32
β42
σ2
τ 2
-1.01(0.24;0.24)
-0.87(0.21;0.21)
-0.87(0.21;0.21)
-0.98(0.18;0.44)
-0.89(0.24;0.24)
-0.97(0.21;0.21)
-1.01(0.22;0.21)
-1.78(0.26;0.42)
-1.13(0.25;0.25)
-1.05(0.21;0.21)
-1.07(0.23;0.22)
-1.19(0.25;0.38)
-1.64(0.29;0.29)
-1.51(0.24;0.24)
-1.71(0.29;0.29)
-1.81(0.39;0.48)
0.40(0.32;0.32)
0.22(0.28;0.28)
0.22(0.29;0.29)
0.80(0.26;0.67)
0.49(0.31;0.31)
0.55(0.28;0.28)
0.61(0.28;0.29)
1.85(0.32;0.64)
0.48(0.33;0.33)
0.42(0.29;0.29)
0.44(0.30;0.30)
0.98(0.33;0.60)
0.40(0.38;0.38)
0.34(0.32;0.32)
0.44(0.37;0.37)
0.97(0.49;0.65)
ρ
Corr.
Introduction to Longitudinal Data Analysis
0.62 0.39
0.39 0.57
0.44 0.44
0.62 0.39
0.39 1.29
1.85 0.59
607
33.8 WGEE for the ARMD Trial in SAS: PROC GENMOD
1.	Data manipulation prior to estimating the parameters in the weight
model:
%dropout(data=armd111,id=subject,time=time,response=bindif,
out=armdhlp);
2.	Code for the weight model:
proc genmod data=armdhlp descending;
class trt prev lesion time;
model dropout = prev trt lesion time / pred dist=binomial;
ods output obstats=pred;
run;
Introduction to Longitudinal Data Analysis
608
3. Data manipulation to prepare for WGEE:
data pred;
set pred;
keep observation pred;
run;
data armdhlp;
merge pred armdhlp;
run;
%dropwgt(data=armdhlp,id=subject,time=time,pred=pred,
dropout=dropout,out=armdwgee);
Introduction to Longitudinal Data Analysis
609
4. WGEE using PROC GENMOD:
proc genmod data=armdcc;
weight wi;
class time treat subject;
model bindif = time treat*time / noint dist=binomial;
repeated subject=subject / withinsubject=time type=exch modelse;
run;
Introduction to Longitudinal Data Analysis
610
33.9 ARMD Trial Analyzed With PROC GEE
Eﬀect Parameter
observation
subject
Weights
Int.4
Int.12
Int.24
Int.52
Tr.4
Tr.12
Tr.24
Tr.52
Corr.
β11
β21
β31
β41
β12
β22
β32
β42
ρ
-0.95 (0.20)
-0.98 (0.35)
-1.03 (0.22)
-1.77 (0.30)
-1.03 (0.23)
-1.11 (0.29)
-1.52 (0.30)
-1.72 (0.37)
0.32 (0.28)
0.78 (0.56)
0.65 (0.29)
1.83 (0.47)
0.39 (0.30)
0.71 (0.49)
0.30 (0.39)
0.72 (0.47)
0.38 0.33
Introduction to Longitudinal Data Analysis
611
o	Subject-level weights:
·	A single weight is calculated for the entire subject
·	The original proposal of the method
·	Slightly easier to calculate
·	Less precise
Observation-level weights:
·	A separate weight is calculated for every measurement within a subject
·	A later modiﬁcation
·	Increased precision
·	Default in SAS
Introduction to Longitudinal Data Analysis
612
33.10 WGEE for the ARMD Trial in SAS: PROC GEE
1.	Preparatory steps (creation of the variables needed in the
MISSMODEL statement):
data help;
set armdwgee;
by subject;
prevbindif=lag(bindif);
if first.id then prevbindif=1;
time2=0;
if time=2 then time2=1;
time3=0;
if time=3 then time3=1;
run;
Introduction to Longitudinal Data Analysis
613
2. PROC GEE code:
proc gee data=help;
class time treat subject lesion;
model bindif = time treat*time / noint dist=binomial;
repeated subject=subject / withinsubject=time type=exch corrw
modelse;
missmodel prevbindif treat lesion time2 time3 / type=obslevel;
run;
Introduction to Longitudinal Data Analysis
614
33.11 Analysis of the Depression Trial
Response: create binary indicator ybin for HAM D17 > 7
Model for dropout:
o	logit[P (Di = j
Di ≥
|
j)] = ψ0 + ψ1yi,j
1 + γTi
−
with
·	yi,j
1: the binary indicator at the previous occasion
−
·	Ti: treatment indicator for patient i
Introduction to Longitudinal Data Analysis
615
·	SAS code:
·	Preparing the dataset:
%dropout(data=depression,id=patient,time=visit,response=ybin,out=dropout);
producing:
·	dropout indicates whether missingness at a given time occurs
·	prev contains outcome at the previous occasion
·	The logistic model for dropout:
proc genmod data=dropout descending;
class trt;
model dropout = prev trt / pred dist=b;
output out=pred p=pred;
run;
·	The weights can now be included in the GENMOD program which speciﬁes the
GEE, through the WEIGHT or SCWGT statements:
Introduction to Longitudinal Data Analysis
616
proc genmod data=study descending;
weight wi;
class patient visitclass trt;
model ybin = trt visit trtvisit basval basvalvisit / dist=bin;
repeated subject=patient / withinsubject=visitclass type=cs corrw;
run;
Results:
·	WGEE
GEE
Eﬀect
est.
(s.e.)
p-value
est.
(s.e.)
p-value
Treatment at visit 4
-1.57
(0.99)
Treatment at visit 5
-0.67
(0.65)
Treatment at visit 6
0.62 (0.56)
Treatment at visit 7
-0.57
(0.37)
Treatment at visit 8
-0.84
(0.39)
0.11 0.30
0.27 0.12
0.03 -0.24
(0.57)
0.09 (0.40)
0.17 (0.34)
-0.43
(0.35)
-0.71
(0.38)
0.67 0.82
0.62 0.22
0.06 Introduction to Longitudinal Data Analysis
617
Chapter 34
Multiple Imputation
·	General idea
·	Estimation
·	Hypothesis testing
·	Use of MI in practice
·	Analysis of the growth data
·	Analysis of the ARMD trial
·	Creating monotone missingness
Introduction to Longitudinal Data Analysis
618
34.1 General Principles
Valid under MAR
An alternative to direct likelihood and WGEE
Three steps:
1.	The missing values are ﬁlled in M times =
⇒
M complete data sets
2.	The M complete data sets are analyzed by using standard procedures
3.	The results from the M analyses are combined into a single inference
Rubin (1987), Rubin and Schenker (1986), Little and Rubin (1987)
·	Introduction to Longitudinal Data Analysis
619
·	Multiple imputation (M = 5 imputations):
Imputation
Combination
........................................
........................................
.........................................................................................................................................................................................................................................................
.........................................................................................................................................................................................................................................................................................................
..........................................................................................................................................................................................................................
....................................................................................................................................................................................... ........................................
.......................................................................................................................................................................................................................... ........................................
......................................................................................................................................................................................................................................................................................................... ........................................
Analysis
....................................
.............................................................................................................................................................
....................................................................................................................................................................................... ........................................
Imputed 1
Results 1
Imputed 2
....................................................................................................................................................................................... ........................................
Results 2
Imputed 3
....................................................................................................................................................................................... ........................................
Results 3
Imputed 4
....................................................................................................................................................................................... ........................................
Results 4
Imputed 5
....................................................................................................................................................................................... ........................................
Results 5
.........................................................................................................................................................................................................................................................
......................................................................................................................................................................................................................................................................................................... ........................................
.......................................................................................................................................................................................................................... ........................................
....................................................................................................................................................................................... ........................................
..........................................................................................................................................................................................................................
........................................
.........................................................................................................................................................................................................................................................................................................
........................................
Observed
data
Final
results
Introduction to Longitudinal Data Analysis
620
34.1.1
Informal Justiﬁcation
o	We need to estimate θ from the data (e.g., from the complete cases)
Plug in the estimated ˆθ and use
to impute the missing data.
f (ym
i |
i , ˆθ)
yo
We need to acknowledge that ˆθ is a random variable; its uncertainty needs to be
included in the imputation process
o	Given this distribution we:
·	draw a random θ∗ from the distribution of ˆθ
·	put this θ∗ in to draw a random Y m
i
f (ym
i |
from
yo
i , θ∗).
Introduction to Longitudinal Data Analysis
621
34.1.2 The Algorithm
1.	Draw θ∗ from its posterior distribution
2.	Draw Y m
∗i
from f (ym
i |
yo
i , θ∗).
3.	To estimate β, then calculate the estimate of the parameter of interest, and its
estimated variance, using the completed data, (Y o, Y m
∗):
The within imputation variance is
ˆβ = ˆβ(Y ) = ˆβ(Y o, Y m
∗)
U =
Var( ˆβ)
d
4.	Repeat steps 1, 2 and 3 a number of M times
m
ˆβ
& U m
(m = 1, . . . , M )
⇒
Introduction to Longitudinal Data Analysis
622
34.1.3 Pooling Information
o	With M imputations, the estimate of β is
ˆβ∗ = P
m
ˆβ
·	M
m=1
M
Further, one can make normally based inferences for β with
where
ˆβ∗)
(β
−
∼
N (0, V ),
total:
V = W + 


M + 1
M 


B
within:
W = P
between:
B = P
m=1 U m
M
M
m=1( ˆβ
M
m
m
ˆβ∗)0
−
ˆβ∗)( ˆβ
1
−
M
−
Introduction to Longitudinal Data Analysis
623
34.1.4 Hypothesis Testing
Two “sample sizes”:
·	N : The sample size of the data set
·	M : The number of imputations
Both play a role in the asymptotic distribution (Li, Raghunathan, and Rubin 1991)
o	H0 : θ = θ0
↓
p = P (Fk,w > F )
Introduction to Longitudinal Data Analysis
624
where
k :
length of the parameter vector θ
Fk,w ∼
F
F =
(θ∗
−
θ0)0W −
1(θ∗
k(1 + r)
θ0)
−
4.	
w = 4 + (τ
r =
1 +
1
k 


−
1
M 


1 +
(1
−
2τ −
r
2







tr(BW −
τ = k(M
−
Limiting behavior:
·	F M
→∞
−→
Fk,
∞
= χ2/k
Introduction to Longitudinal Data Analysis
625
34.2 Use of MI in Practice
Many analyses of the same incomplete set of data
A combination of missing outcomes and missing covariates
As an alternative to WGEE: MI can be combined with classical GEE
MI in SAS:
·	Imputation Task:
PROC MI
Analysis Task:
↓
PROC “MYFAVORITE”
Inference Task:
↓
PROC MIANALYZE
Introduction to Longitudinal Data Analysis
626
34.3 MI Analysis of the Orthodontic Growth Data
The same Model 1 as before
Focus on boys at ages 8 and 10
Results
Method
Boys at Age 8 Boys at Age 10
Original Data
22.88 (0.58)
23.81 (0.51)
Multiple Imputation
22.88 (0.66)
22.69 (0.81)
Between-imputation variability for age 10 measurement
Conﬁdence interval for Boys at age 10: [21.08,24.29]
o	Introduction to Longitudinal Data Analysis
627
34.4 MI Analysis of the ARMD Trial
M = 10 imputations
GEE:
GLMM:
logit[P (Yij = 1
|
Ti, tj)] = βj1 + βj2Ti
logit[P (Yij = 1
|
Ti, tj, bi)] = βj1 + bi + βj2Ti,
N (0, τ 2)
bi ∼
Ti = 0 for placebo and Ti = 1 for interferon-α
tj (j = 1, . . . , 4) refers to the four follow-up measurements
Imputation based on the continuous outcome
§	Introduction to Longitudinal Data Analysis
628
Results:
·	Eﬀect
Par.
GEE
GLMM
Int.4
Int.12
Int.24
Int.52
Trt.4
Trt.12
Trt.24
Trt.52
R.I. s.d.
R.I. var.
β11
β21
β31
β41
β12
β22
β32
β42
τ
τ 2
-0.84(0.20)
-1.46(0.36)
-1.02(0.22)
-1.75(0.38)
-1.07(0.23)
-1.83(0.38)
-1.61(0.27)
-2.69(0.45)
0.21(0.28)
0.32(0.48)
0.60(0.29)
0.99(0.49)
0.43(0.30)
0.67(0.51)
0.37(0.35)
0.52(0.56)
2.20(0.26)
4.85(1.13)
Introduction to Longitudinal Data Analysis
629
34.5 SAS Tools for MI
MONOTONE: Input data must be monotone only!
There are several options:
reg: every ‘later’ variable is regressed on the earlier ones (normal model)
logistic: every ‘later’ variable is regressed on the earlier ones (logistic model)
discrim: version for categorical data where discriminant functions are used
regpmm: predictive mean matching (a choice is made from observations with
value similar to the predicted mean)
propensity: Some detail:
∗
∗
The propensity score is the conditional probability of assignment to a
particular treatment given a vector of observed covariates.
A propensity score is generated for each variable with missing values to
indicate the probability of that observations being missing.
Introduction to Longitudinal Data Analysis
630
The observations are then grouped based on these propensity scores, and an
approximate Bayesian bootstrap imputation (Rubin 1987) is applied to each
group.
Less suitable when also relations betwen variables (association, regression) is
of interest.
Hence of limited use.
∗
∗
∗
MCMC: can be used with monotone and non-monotone missingness:
·	a multivariate normal model is considered
·	Monte Carlo draws are taken from it
·	Default in SAS
Introduction to Longitudinal Data Analysis
631
FCS: a extension of MONOTONE that also allows for non-monotone missingness:
·	Fill in ‘starting values’ for the missing values (e.g., from a predictive model
based on completers)
·	Cycle repeatedly through the following two steps:
Build a predictive model for the jth variable given all the others (recall that
missing values are ‘temporarily’ ﬁlled in)
Based on the updated model, draw values for the missing ones
∗
∗
·	Same options as MONOTONE, except for propensity scores
Introduction to Longitudinal Data Analysis
632
34.6 SAS Code for MI
1.	Preparatory data analysis so that there is one line per subject
2.	The imputation task (default MCMC method):
proc mi data=armd13 seed=486048 out=armd13a simple nimpute=10 round=0.1;
var lesion diff4 diff12 diff24 diff52;
by treat;
run;
Note that the imputation task is conducted on the continuous outcome ‘diﬀ
·
indicating the diﬀerence in number of letters versus baseline
’,
3.	Then, data manipulation takes place to deﬁne the binary indicators and to create
a longitudinal version of the dataset
Introduction to Longitudinal Data Analysis
633
4. The imputation task using FCS:
proc mi data=m.armd13 seed=486048 simple out=m.armd13fcs nimpute=30 round=0.01;
fcs reg(diff4=lesion);
fcs reg(diff12=lesion diff4);
fcs reg(diff24=lesion diff4 diff12);
fcs reg(diff52=lesion diff4 diff12 diff24);
var lesion diff4 diff12 diff24 diff52;
by treat;
run;
Introduction to Longitudinal Data Analysis
634
5. Details on intermediate steps:
·	Dichotomization of imputed data:
proc sort data=m.armd13a;
by imputation subject;
run;
data m.armd13a;
set m.armd13a;
bindif4=0; if diff4 <= 0
bindif12=0; if diff12 <= 0
bindif24=0; if diff24 <= 0
bindif52=0; if diff52 <= 0
if diff4=. then bindif4=.;
if diff12=. then bindif12=.;
if diff24=. then bindif24=.;
if diff52=. then bindif52=.;
run;
then bindif4=1;
then bindif12=1;
then bindif24=1;
then bindif52=1;
Introduction to Longitudinal Data Analysis
635
·	Transforming horizontal dataset in vertical dataset:
data m.armd13b;
set m.armd13a;
array x (4) bindif4 bindif12 bindif24 bindif52;
array y (4) diff4 diff12 diff24 diff52;
do j=1 to 4;
bindif=x(j);
diff=y(j);
time=j;
output;
end;
run;
·	Creating dummies:
data m.armd13c;
set m.armd13b;
time1=0;
time2=0;
time3=0;
time4=0;
Introduction to Longitudinal Data Analysis
636
trttime1=0;
trttime2=0;
trttime3=0;
trttime4=0;
if time=1 then time1=1;
if time=2 then time2=1;
if time=3 then time3=1;
if time=4 then time4=1;
if (time=1 & treat=1) then trttime1=1;
if (time=2 & treat=1) then trttime2=1;
if (time=3 & treat=1) then trttime3=1;
if (time=4 & treat=1) then trttime4=1;
run;
proc sort data=m.armd13cs;
by imputation subject time;
run;
Introduction to Longitudinal Data Analysis
637
6. The analysis task (GEE):
proc gee data=armd13c;
class time subject;
by imputation;
model bindif = time1 time2 time3 time4 trttime1 trttime2 trttime3 trttime4
/ noint dist=binomial covb;
repeated subject=subject / withinsubject=time type=exch modelse;
ods output GEEEmpPEst=gmparms parminfo=gmpinfo CovB=gmcovb;
run;
Deletion of redundant parameter information:
data gmpinfo;
set gmpinfo;
if parameter=’Prm1’ then delete;
run;
Introduction to Longitudinal Data Analysis
638
7. The analysis task (GLMM):
proc nlmixed data=armd13c qpoints=20 maxiter=100 technique=newrap cov ecov;
by imputation;
eta = beta11time1+beta12time2+beta13time3+beta14time4+b
+beta21trttime1+beta22trttime2+beta23trttime3+beta24trttime4;
p = exp(eta)/(1+exp(eta));
model bindif  binary(p);
random b  normal(0,tautau) subject=subject;
estimate ’tau2’ tautau;
ods output ParameterEstimates=nlparms
CovMatParmEst=nlcovb
AdditionalEstimates=nlparmsa
CovMatAddEst=nlcovba;
run;
Introduction to Longitudinal Data Analysis
639
8. The inference task (GEE):
proc mianalyze parms=gmparms covb=gmcovb parminfo=gmpinfo wcov bcov tcov;
modeleffects time1 time2 time3 time4 trttime1 trttime2 trttime3 trttime4;
run;
9.	The inference task (GLMM):
proc mianalyze parms=nlparms covb=nlcovb wcov bcov tcov;
modeleffects beta11 beta12 beta13 beta14 beta21 beta22 beta23 beta24;
run;
Introduction to Longitudinal Data Analysis
640
Chapter 35
Creating Monotone Missingness
·	When missingness is non-monotone, one might think of several mechanisms
operating simultaneously:
·	A simple (MCAR or MAR) mechanism for the intermittent missing values
·	A more complex (MNAR) mechanism for the missing data past the moment of
dropout
·	Analyzing such data are complicated, especially with methods that apply to
dropout only
Introduction to Longitudinal Data Analysis
641
·	Solution:
·	Generate multiple imputations that render the datasets monotone missing, by
including into the MI procedure:
mcmc impute=monotone;
·	Apply method of choice to the so-completed multiple sets of data
·	Note: this is diﬀerent from the monotone method in PROC MI, intended to fully
complete already monotone sets of data:
o	The MONOTONE statement takes monotone patterns as input and returns
completed data.
The MCMC statement with impute=monotone option takes data with also
non-monotone patterns as input and returns monotonized data.
Introduction to Longitudinal Data Analysis
642
35.1 Example: Creating Monotone Missingness to Then
Apply Weighted GEE
o	Consider again the analgesic trial
Multiple imputation to create monotone missingness:
proc mi data=m.gsa4 seed=459864 simple nimpute=10
round=0.1 out=m.gsaimput;
title ’Monotone multiple imputation’;
mcmc impute = monotone;
var pca0 physfct gsa1 gsa2 gsa3 gsa4;
run;
Introduction to Longitudinal Data Analysis
643
Preparation of the data in vertical format, so that the data can be used in
ordinary GEE:
·	data m.gsaw;
set m.gsa4;
array y (4) gsa1 gsa2 gsa3 gsa4;
do j=1 to 4;
gsa=y(j);
time=j;
timecls=time;
gsabin=.;
if gsa=1 then gsabin=1;
if gsa=2 then gsabin=1;
if gsa=3 then gsabin=1;
if gsa=4 then gsabin=0;
if gsa=5 then gsabin=0;
output;
end;
run;
Introduction to Longitudinal Data Analysis
644
Standard GEE:
·	proc gee data=m.gsaw plots=histogram;
title ’Standard GEE for GSA data’;
class patid timecls;
model gsabin = time|time pca0 / dist=bin;
repeated subject=patid / within=timecls corr=un;
run;
Steps to prepare the data for weighted GEE, including deﬁnition of the ‘previous’
outcome:
·	data m.gsaimput02;
set m.gsaimput;
array y (4) gsa1 gsa2 gsa3 gsa4;
do j=1 to 4;
gsa=y(j);
time=j;
timecls=time;
Introduction to Longitudinal Data Analysis
645
patid2=1000*imputation+patid;
output;
end;
run;
proc sort data=m.gsaimput02;
by imputation patid2;
run;
data m.gsaimput03;
set m.gsaimput02;
by patid2;
prevgsa = lag(gsa);
if time=1 then prevgsa = 1;
timeclss = time;
run;
Introduction to Longitudinal Data Analysis
646
data m.gsaimput03;
set m.gsaimput03;
if gsa<=3.5 then gsabin=1;
if gsa>3.5 then gsabin=0;
gsabin=gsabin+gsa-gsa;
run;
Introduction to Longitudinal Data Analysis
647
Weighted GEE, where weights are created at observation level:
·	ods graphics on;
proc gee data=m.gsaimput03 plots=histogram;
title ’Weighted GEE for GSA Data Based on Multiple
Imputation to Monotonize - OBSLEVEL’;
by imputation;
class patid timecls;
model gsabin = time|time pca0 / dist=bin covb;
repeated subject=patid / within=timecls corr=un ecovb;
missmodel prevgsa pca0 physfct / type=obslevel;
ods output GEEEmpPEst=gmparms parminfo=gmpinfo
modelinfo=modelinfo GEERCov=gmcovb;
run;
proc mianalyze parms=gmparms parminfo=gmpinfo covb=gmcovb;
title ’Multiple Imputation Analysis After Weighted GEE for GSA Data’;
modeleffects intercept time time*time pca0;
run;
Introduction to Longitudinal Data Analysis
648
To use weights at subject rather than observation level:
missmodel prevgsa pca0 physfct / type=sublevel;
Evidently, using these monotonized data, also standard GEE can be used:
o	ods graphics on;
proc gee data=m.gsaimput03 plots=histogram;
title ’Standard GEE for GSA Data Based on
Multiple Imputation to Monotonize’;
by imputation;
class patid timecls;
model gsabin = time|time pca0 / dist=bin covb;
repeated subject=patid / within=timecls corr=un ecovb;
ods output GEEEmpPEst=gmparms parminfo=gmpinfo
modelinfo=modelinfo GEERCov=gmcovb;
run;
Introduction to Longitudinal Data Analysis
649
Files:
·	analg11(met-proc-gee).sas
·	analg11(met-proc-gee).lst
Overview of results:
o	Introduction to Longitudinal Data Analysis
650
Eﬀect
Par.
Est.(s.e.)
p-value
Est.(s.e.)
p-value
Standard GEE
Without MI
After MI
2.90(0.46)
-0.81(0.32)
0.17(0.07)
-0.23(0.10)
2.87(0.45)
-0.83(0.32)
0.18(0.06)
-0.21(0.10)
0.0124 0.0083
0.0178 0.0087
0.0058 0.0253
Weighted GEE (after MI)
Observation level
Subject level
2.74(0.46)
-0.76(0.33)
0.17(0.07)
-0.19(0.10)
2.62(0.60)
-0.71(0.40)
0.16(0.08)
-0.21(0.12)
0.0231 0.0155
0.0384 0.0747
0.0444 0.0853
Intercept
Time
Time2
PCA0
Intercept
Time
Time2
PCA0
β0
β1
β2
β3
β0
β1
β2
β3
Introduction to Longitudinal Data Analysis
651
·	The dropout model is similar but slightly diﬀerent than the one used with PROC
GENMOD.
o	Weighted GEE leads to increased standard errors, as observed before.
This eﬀect is less pronounced when weigths are constructed at observation level,
rather than at subject level.
·	A typical output for one of the imputed datasets takes the form (ﬁrst imputation
out of ten; with weights at observation level):
Introduction to Longitudinal Data Analysis
652
Parameter Estimates for Response Model
with Empirical Standard Error
Parameter Estimate
Standard
Error
95% Confidence
Limits
Z Pr > |Z|
Intercept
TIME
TIME*TIME
PCA0
2.4299 -0.5881
0.1392 -0.1797
0.5890 1.2755
0.3912 -1.3548
0.0794 -0.0165
0.1173 -0.4096
3.5843 0.1787
0.2949 0.0501
4.13 -1.50
1.75 -1.53
<.0001
0.1328 0.0796
0.1254 Parameter Estimates for Missingness Model
Parameter
Estimate
Standard
Error
95% Confidence
Limits
Z
Pr > |Z|
Intercept
prevgsa
PCA0
PHYSFCT
3.1335 -0.1974
-0.2495
-0.0079
0.4060 0.0822
0.0956 0.0037
2.3377 -0.3585
-0.4370
-0.0151
3.9293 -0.0363
-0.0621
-0.0007
7.72 -2.40
-2.61
-2.16
<.0001
0.0163 0.0091
0.0311 Introduction to Longitudinal Data Analysis
653
Part VI
Topics in Methods and Sensitivity Analysis for Incomplete
Data
Introduction to Longitudinal Data Analysis
654
Chapter 36
An MNAR Selection Model and Local Inﬂuence
·	The Diggle and Kenward selection model
·	Mastitis in dairy cattle
·	An informal sensitivity analysis
·	Local inﬂuence to conduct sensitivity analysis
Introduction to Longitudinal Data Analysis
655
36.1 A Full Selection Model
MNAR :
f (Y i|
θ)f (Di|
Z
Y i, ψ)dY m
i
f (Y i|
θ)
f (Di|
Y i, ψ)
Linear mixed model
Logistic regressions for dropout
Y i = Xiβ + Zibi + εi
logit [P (Di = j
Di ≥
|
j, Yi,j
1, Yij)]
−
= ψ0 + ψ1Yi,j
1 + ψ2Yij
−
Diggle and Kenward (JRSSC 1994)
Introduction to Longitudinal Data Analysis
656
36.2 Mastitis in Dairy Cattle
·	Infectious disease of the udder
Leads to a reduction in milk yield
High yielding cows more susceptible?
But this cannot be measured directly be-
cause of the eﬀect of the disease: ev-
idence is missing since infected cause
have no reported milk yield
Introduction to Longitudinal Data Analysis
657
Model for milk yield:
Yi1
Yi2
















N
∼
















µ
µ + ∆








,








σ2
1
ρσ1σ2
ρσ1σ2
σ2
2
















Model for mastitis:
logit [P (Ri = 1
|
Yi1, Yi2)] = ψ0 + ψ1Yi1 + ψ2Yi2
= 0.37 + 2.25Yi1 −
0.29Yi1 −
= 0.37
−
2.54Yi2
2.54(Yi2 −
Yi1)
LR test for H0 : ψ2 = 0 : G2 = 5.11
§	Introduction to Longitudinal Data Analysis
658
36.3 Criticism
−→
Sensitivity Analysis
“. . . , estimating the ‘unestimable’ can be accomplished only by making
modelling assumptions,. . . . The consequences of model misspeci-
ﬁcation will (. . . ) be more severe in the non-random case.”
(Laird
1994)
Change distributional assumptions
(Kenward 1998)
Local and global inﬂuence methods
Pattern-mixture models
Several plausible models
or
ranges of inferences
Semi-parametric framework
(Scharfstein et al 1999)
o	Introduction to Longitudinal Data Analysis
659
36.4 Kenward’s Sensitivity Analysis
·	Deletion of #4 and #5
0.08 5.11
−→
G2 for ψ2:
⇒
·	Cows #4 and #5 have unusually large
increments
·	Kenward conjectures: #4 and #5 ill dur-
ing the ﬁrst year
Kenward (SiM 1998)
·	Introduction to Longitudinal Data Analysis
660
36.5 Local Inﬂuence
§	Verbeke, Thijs, Lesaﬀre, Kenward (Bcs 2001)
Perturbed MAR dropout model:
logit [P (Di = 1
|
Yi1, Yi2)]
= ψ0 + ψ1Yi1 + ωiYi2
Likelihood displacement:
LD(ω) = 2
Lω=0
(cid:20)
θ,
ψ
c
d
(cid:18)
(cid:19) −
Lω=0
θω,
c
ψω(cid:19) (cid:21) ≥
d
(cid:18)
0
Introduction to Longitudinal Data Analysis
661
36.5 Local Inﬂuence
§	Verbeke, Thijs, Lesaﬀre, Kenward (Bcs 2001)
Perturbed MAR dropout model:
logit [P (Di = 1
|
Yi1, Yi2)]
= ψ0 + ψ1Yi1 + ωiYi2
or ψ0 + ψ1Yi1 + ωi(Yi2 −
Yi1)
Likelihood displacement:
LD(ω) = 2
Lω=0
(cid:20)
θ,
ψ
c
d
(cid:18)
(cid:19) −
Lω=0
θω,
c
ψω(cid:19) (cid:21) ≥
d
(cid:18)
0
Introduction to Longitudinal Data Analysis
662
36.5.1 Likelihood Displacement
LD(ω)
ωj
PPPPPPPPPPPPPPPP
ω = 0
s
h
Local inﬂuence direction h
↑
normal curvature Ch
T0
Local inﬂuence for θ and ψ:
·	ωi
Ch = Ch(θ) + Ch(ψ)
Introduction to Longitudinal Data Analysis
663
36.5.2 Computational Approaches
Measuring local inﬂuence:
Fit for continuous outcomes:
o	Expression for Ch:
Ch = 2
|
h0 ∆0 ¨L−
1 ∆ h
|
·	Fit MAR model:
·	linear mixed model for outcomes
·	logistic regression for dropout
Choices for h
·	Direction of the ith subject
·	evaluate closed-form expressions for
local inﬂuence
Ci
·	Direction hmax of maximal curva-
⇒
ture Cmax
Introduction to Longitudinal Data Analysis
664
36.6 Application to Mastitis Data
Removing #4, #5 and #66
G2 = 0.005
⇒
hmax: diﬀerent signs for (#4,#5) and #66
o	Introduction to Longitudinal Data Analysis
665
36.6.1
Interpretable Components of Ci(ψ)
P (Ri = 1) [1
P (Ri = 1)]
−
Yi1
Yi2 −
or
Yi1)
E(Yi2 |
Yi1
−
Vi ( 1 Yi1 )



× (cid:18)
j
X
1
Yi1
1
Yj1
Vj
(cid:18)
(cid:19)
( 1 Yj1 )
(cid:19)
−1



666
Introduction to Longitudinal Data Analysis
36.7 Global Inﬂuence Analysis
o	MAR versus MNAR model
For a variety of subsets:
·	All data
·	Removal of:
(53,54,66,69): from local inﬂuence on Yi2
(4,5): from Kenward’s informal analysis
(66): additional one identiﬁed from local inﬂuence on Yi2 −
(4,5,66): frol local inﬂuence on Yi2 −
Yi1
∗
∗
∗
∗
Yi1
Introduction to Longitudinal Data Analysis
667
Eﬀect
Measurement model:
Intercept
Time eﬀect
First variance
Second variance
Correlation
Dropout model:
Intercept
First measurement
Second measurement
-2 loglikelihood
Eﬀect
Measurement model:
Intercept
Time eﬀect
First variance
Second variance
Correlation
Dropout model:
Intercept
First measurement
Second measurement
-2loglikelihood
G2 for MNAR
Parameter
all
(53,54,66,69)
µ
∆
σ2
1
σ2
2
ρ
ψ0
ψ1
ω = ψ2
5.77(0.09)
0.72(0.11)
0.87(0.12)
1.30(0.20)
0.58(0.07)
-2.65(1.45)
0.27(0.25)
0
280.02 5.69(0.09)
0.70(0.11)
0.76(0.11)
1.08(0.17)
0.45(0.08)
-3.69(1.63)
0.46(0.28)
0
246.64 Parameter
all
(53,54,66,69)
µ
∆
σ2
1
σ2
2
ρ
ψ0
ψ1
ω = ψ2
5.77(0.09)
0.33(0.14)
0.87(0.12)
1.61(0.29)
0.48(0.09)
0.37(2.33)
2.25(0.77)
-2.54(0.83)
274.91 5.11
5.69(0.09)
0.35(0.14)
0.76(0.11)
1.29(0.25)
0.42(0.10)
-0.37(2.65)
2.11(0.76)
-2.22(0.86)
243.21 3.43
MAR
(4,5)
5.81(0.08)
0.64(0.09)
0.77(0.11)
1.30(0.20)
0.72(0.05)
-2.34(1.51)
0.22(0.25)
0
237.94 MNAR
(4,5)
5.81(0.08)
0.40(0.18)
0.77(0.11)
1.39(0.25)
0.67(0.06)
-0.77(2.04)
1.61(1.13)
-1.66(1.29)
237.86 0.08
(66)
(4,5,66)
5.75(0.09)
0.68(0.10)
0.86(0.12)
1.10(0.17)
0.57(0.07)
5.80(0.09)
0.60(0.08)
0.76(0.11)
1.09(0.17)
0.73(0.05)
-2.77(1.47)
0.29(0.24)
0
264.73 -2.48(1.54)
0.24(0.26)
0
220.23 (66)
(4,5,66)
5.75(0.09)
0.34(0.14)
0.86(0.12)
1.34(0.25)
0.48(0.09)
0.45(2.35)
2.06(0.76)
-2.33(0.86)
261.15 3.57
5.80(0.09)
0.63(0.29)
0.76(0.11)
1.10(0.20)
0.73(0.05)
-2.77(3.52)
0.07(1.82)
0.20(2.09)
220.23 0.005
Introduction to Longitudinal Data Analysis
668
Chapter 37
Local Inﬂuence for the ARMD Trial
Eﬀect
Parameter
Ignorable
MCAR
MAR
MNAR
Int. 4
Int.12
Int.24
Int.52
Trt. 4
Trt. 12
Trt. 24
Trt. 52
Int.
Previous
Current
β11
β21
β31
β41
β12
β22
β32
β42
ψ0
ψ1
ψ2
Measurement model
54.00 (1.47)
53.01 (1.60)
49.20 (1.74)
43.99 (1.79)
-3.11 (2.10)
-4.54 (2.29)
-3.60 (2.49)
-5.18 (2.59)
54.00 (1.46)
53.01 (1.59)
49.20 (1.73)
43.99 (1.78)
-3.11 (2.07)
-4.54 (2.25)
-3.60 (2.46)
-5.18 (2.57)
Dropout model
-2.79 (0.17)
54.00 (1.47)
53.01 (1.60)
49.19 (1.74)
43.99 (1.79)
-3.11 (2.09)
-4.54 (2.29)
-3.60 (2.50)
-5.18 (2.62)
54.00 (1.47)
52.98 (1.60)
49.06 (1.74)
43.52 (1.82)
-3.11 (2.10)
-4.67 (2.29)
-3.80 (2.50)
-5.71 (2.63)
-1.86 (0.46)
-0.020 (0.009)
-1.81 (0.47)
0.016 (0.022)
-0.042 (0.023)
-2 log-likelihood
6488.7 6782.7
6778.4 6775.9
Treatment eﬀect at 1 year (p-value)
0.046 0.044
0.048 0.030
Introduction to Longitudinal Data Analysis
669
# 27
28
10
0
0
0
0
0
4
0
0
0
0
0
3
0
0
0
0
0
2
0
0
0
0
0
1
0
C i
139
154
114
0
50
100
150
200
Introduction to Longitudinal Data Analysis
670
C i(θ)
68
185
0
0
5
0
0
4
0
0
3
0
0
2
0
0
1
0
0
50
100
150
200
Introduction to Longitudinal Data Analysis
671
C i(β)
6
·	0
5
·	0
4
·	0
3
·	0
2
·	0
1
·	0
0
·	0
0
50
100
150
200
Introduction to Longitudinal Data Analysis
672
C i(α)
68
185
0
0
5
0
0
4
0
0
3
0
0
2
0
0
1
0
0
50
100
150
200
Introduction to Longitudinal Data Analysis
673
# 27
28
10
0
0
0
0
0
4
0
0
0
0
0
3
0
0
0
0
0
2
0
0
0
0
0
1
0
C i(ψ)
139
154
114
0
50
100
150
200
Introduction to Longitudinal Data Analysis
674
# 27
10
28
hmax,i
139
114
154
0
·	1
5
·	0
0
·	0
5
·	0
0
·	1
0
50
100
150
200
Introduction to Longitudinal Data Analysis
675
Active
y
t
i
u
c
A
l
a
u
s
V
i
0
0
1
0
8
0
6
0
4
0
2
0
68
139
28
10
154
114
4
12
24
Introduction to Longitudinal Data Analysis
Weeks
52
676
Placebo
y
t
i
u
c
A
l
a
u
s
V
i
0
0
1
0
8
0
6
0
4
0
2
0
27
185
4
12
24
Introduction to Longitudinal Data Analysis
Weeks
52
677
Eﬀect
Int. 4
Int.12
Int.24
Int.52
Trt. 4
Trt.12
Trt.24
Trt.52
Parameter
Set 1
MAR
Set 2
MAR
Set 3
MAR
Measurement model
β11
β21
β31
β41
β12
β22
β32
β42
54.14(1.51)
53.09(1.64)
49.56(1.77)
44.40(1.82)
-3.13(2.17)
-4.48(2.36)
-3.80(2.56)
-5.45(2.66)
Dropout model
-1.90(0.47)
-0.019(0.010)
6535.3 0.040
54.30(1.47)
53.16(1.59)
49.31(1.74)
44.00(1.79)
-3.28(2.08)
-4.55(2.26)
-3.55(2.48)
-5.06(2.59)
53.84(1.48)
52.94(1.60)
49.44(1.73)
44.38(1.78)
-2.95(2.07)
-4.47(2.26)
-3.85(2.44)
-5.56(2.55)
-1.90(0.47)
-0.019(0.010)
6606.9 0.051
-1.85(0.46)
-0.020(0.009)
6706.4 0.029
Intercept
Previous
-2 log-likelihood
Treatm. eﬀ. at 1 year (p-value)
ψ0
ψ1
Introduction to Longitudinal Data Analysis
678
Eﬀect
Int. 4
Int.12
Int.24
Int.52
Trt. 4
Trt.12
Trt.24
Trt.52
Parameter
Set 1
MNAR
Set 2
MNAR
Set 3
MNAR
Measurement model
β11
β21
β31
β41
β12
β22
β32
β42
54.15(1.49)
53.06(1.62)
49.46(1.75)
43.97(1.84)
-3.13(2.11)
-4.63(2.29)
-4.04(2.49)
-6.12(2.66)
Dropout model
-1.85(0.49)
0.018(0.022)
-0.044(0.024)
6532.7 0.021
54.30(1.46)
53.13(1.59)
49.20(1.72)
43.58(1.82)
-3.28(2.06)
-4.69(2.24)
-3.79(2.44)
-5.72(2.61)
53.84(1.47)
52.91(1.59)
49.31(1.72)
43.90(1.82)
-2.95(2.05)
-4.60(2.23)
-4.04(2.42)
-6.09(2.58)
-1.85(0.49)
0.017(0.022)
-0.043(0.024)
6604.4 0.028
-1.81(0.47)
0.017(0.022)
-0.043(0.024)
6703.8 0.018
Intercept
Previous
Current
-2 log-likelihood
Treatm. eﬀ. at 1 year (p-value)
ψ0
ψ1
ψ2
Introduction to Longitudinal Data Analysis
679
Chapter 38
Mechanism for Growth Data
By the way, how did Little and Rubin
delete data from the growth data set ?
Introduction to Longitudinal Data Analysis
680
38.1 Modeling Missingness
Candidate model for missingness:
logit[P (Ri = 0
|
yi)] = ψ0 + ψ1yij,
with j = 1, 2, 3, or 4
When j = 2, then MNAR, else MAR.
Results:
Mechanism Eﬀects Deviance
Yi1
Yi3
Yi4
Yi2
MAR
MAR
MAR
MNAR
19.51 7.43
2.51 2.55
p-value
<0.0001
0.0064 0.1131
0.1105 Including covariates:
·	Introduction to Longitudinal Data Analysis
681
Boys : logit[P (Ri = 0
|
Girls : logit[P (Ri = 0
|
yi1, xi = 0)] =
yi1, xi = 1)] =
These models are interpreted as follows:
·	(22
−
(20.75
∞
∞
yi1)
yi1)
−
Boys : P (Ri = 0
|
yi1, xi = 0) =
1
0.5 0



if yi1 < 22,
if yi1 = 22,
if yi1 > 22.
Girls : P (Ri = 0
|
yi1, xi = 1) = 


1
0
if yi1 < 20.75,
if yi1 > 20.75.
Introduction to Longitudinal Data Analysis
682
Chapter 39
Interval of Ignorance / Bodyguard
·	The Toenail Data
·	The Fluvoxamine Trial
·	The Slovenian Public Opinion Survey
·	MAR and MNAR analyses
·	Informal sensitivity analysis
·	Interval of ignorance & interval of uncertainty
Introduction to Longitudinal Data Analysis
683
39.1 Fluvoxamine Trial: Side Eﬀects
89 13
57 65
26
49
2 0
14
Post-marketing study of ﬂuvoxamine in psychiatric patients
Absence versus presence of side eﬀects
Two measurement occasions
315 subjects:
·	224 completers, 75 drop out after ﬁrst, 2 non-monotone, 14 without follow up
Questions:
·	Do side eﬀects evolve over time ?
·	Are both measurements dependent ?
o	Introduction to Longitudinal Data Analysis
684
39.2 The Slovenian Plebiscite
Rubin, Stern, and Vehovar (1995)
Slovenian Public Opinion (SPO) Survey
Four weeks prior to decisive plebiscite
Three questions:
1.	Are you in favor of Slovenian independence ?
2.	Are you in favor of Slovenia’s secession from Yugoslavia ?
3.	Will you attend the plebiscite ?
Political decision: ABSENCE≡NO
Primary Estimand: θ: Proportion in favor of independence
§	Introduction to Longitudinal Data Analysis
685
Slovenian Public Opinion Survey Data:
·	Independence
Secession
Attendance
Yes
No
Yes
No
∗
Yes
No
∗
Yes
No
∗
Yes
No
∗
1191
8
107
158
7
18
90
1
19
8
0
3
68
14
43
2
2
8
∗
21
4
9
29
3
31
109
25
96
Introduction to Longitudinal Data Analysis
686
39.3 Slovenian Public Opinion: 1st Analysis
§	Pessimistic: All who can say NO will say NO
1439
2074
θ =
c
= 0.694
Optimistic: All who can say YES will say YES
1439 + 159 + 144 + 136
2074
=
1878
2076
θ =
c
= 0.904
Resulting Interval:
[0.694; 0.904]
θ
∈
Introduction to Longitudinal Data Analysis
687
§	Resulting Interval:
[0.694; 0.904]
θ
∈
Complete cases: All who answered on 3 questions
1191 + 158
1454
θ =
c
= 0.928 ?
Available cases: All who answered on both questions
1191 + 158 + 90
1549
θ =
c
= 0.929 ?
Introduction to Longitudinal Data Analysis
688
39.4 Slovenian Public Opinion: 2nd Analysis
o	Missing at Random:
Non-response is allowed to depend on observed, but not on unobserved outcomes:
·	Based on two questions:
·	Based on three questions:
θ = 0.892
c
θ = 0.883
c
Missing Not at Random (NI):
Non-response is allowed to depend on unobserved measurements:
θ = 0.782
c
Introduction to Longitudinal Data Analysis
689
39.5 Slovenian Public Opinion Survey
Estimator
Pessimistic bound
Optimistic bound
Complete cases
Available cases
MAR (2 questions)
MAR (3 questions)
MNAR
θ
c
0.694 0.904
0.928 ?
0.929 ?
0.892 0.883
0.782 Introduction to Longitudinal Data Analysis
690
39.6 Slovenian Plebiscite: The Truth ?
θ =0.885
Estimator
Pessimistic bound
Optimistic bound
Complete cases
Available cases
MAR (2 questions)
MAR (3 questions)
MNAR
θ
c
0.694 0.904
0.928 ?
0.929 ?
0.892 0.883
0.782 Introduction to Longitudinal Data Analysis
691
39.7 MAR in 3 Frameworks
Selection models
f (ri|
yi, ψ) = f (ri|
yo
i , ψ)
Pattern-mixture models
f (ym
i |
i , ri, θ) = f (ym
yo
i |
yo
i , θ)
Shared-parameter models
?
Introduction to Longitudinal Data Analysis
692
39.8 MAR in Selection Models
Diggle and Kenward (ApStat 1994)
f (ri|
yi, ψ) = f (ri|
yo
i , ψ)
Longitudinal data:
logit [P (Di = j
j, yij, yi,j
Di ≥
|
= 0
ψ2 6
ψ2 = 0
ψ1 = ψ2 = 0
1)] = ψ0 + ψ1yi,j
−
1 + ψ2yij
−
←→
←→
←→
MNAR
MAR
MCAR
No dependence on the future (NFD): built in
o	Introduction to Longitudinal Data Analysis
693
39.9 MAR in Pattern-mixture Models
Molenberghs, Michiels, Kenward, and Diggle (Statistica Neerlandica 1998)
Thijs, Molenberghs, Michiels, Verbeke, and Curran (Biostatistics 2002)
f (ym
i |
i , ri, θ) = f (ym
yo
i |
yo
i , θ)
For longitudinal data: ACMV: available case missing value restrictions:
t
2,
∀
≥
· · ·
Practical implementation: doable!
∀
s < t : f (yit|
yi1,
, yi,t
1, di = s) = f (yit|
yi1,
· · ·
, yi,t
1, di ≥
−
t)
−
o	f (yit|
yi1,
· · ·
, yi,t
1, di = s) =
−
n
Xd=s




P
αdfd(yi1, . . . , yi,s
−
ni
d=s αdfd(yi1, . . . , yi,s
−
fd(ys|




yi1, . . . , yi,s
−
Introduction to Longitudinal Data Analysis
694
39.10 Non-future Dependence in Pattern-Mixture Models
Kenward, Molenberghs, and Thijs (Biometrika 2003)
Within every pattern:
·	Past: Build a model for the observed data
·	Present, given past: For the ﬁrst unobserved time, given the past: Free
choice!
·	Future, given past and present: Use ACMV-type restrictions
Named NFMV: non-future missing values
Equivalence:
§	SeM: NFD
⇐⇒
PMM: NFMV
Introduction to Longitudinal Data Analysis
695
39.11 MAR in Shared-parameter Models
Creemers, Hens, Aerts, Molenberghs, Verbeke, and Kenward (2008)
Conventional
f (yi|
Xi, bi, θ) f (ri|
Xi, bi, ψ)
Extended
f (yo
i |
gi, hi, ji, `i)f (ym
i |
yo
i , gi, hi, ki, mi) f (ri|
gi, ji, ki ni)
∩
f (yo
i |
R
MAR
gi, ji, ki)f (bi) dbi
∪
R
f (yo
i |
gi, hi, ji)f (ym
yo
i , gi, hi, ki)f (ri|
i |
gi, ji)f (ri|
yo
gi, hi)f (ym
i , gi, hi)f (bi) dbi
i |
f (yo
i )
gi, ji)f (bi) dbi
f (yo
i |
R
Introduction to Longitudinal Data Analysis
696
39.12 MAR in Shared-parameter Models
Creemers, Hens, Aerts, Molenberghs, Verbeke, and Kenward (2008)
Extended
f (yo
i |
gi, hi, ji, `i)f (ym
i |
yo
i , gi, hi, ki, mi) f (ri|
gi, ji, ki ni)
f (yo
i |
R
MAR
gi, ji, ki)f (bi) dbi
∪
R
f (yo
i |
gi, hi, ji)f (ym
yo
i , gi, hi, ki)f (ri|
i |
gi, ji)f (ri|
yo
gi, hi)f (ym
i , gi, hi)f (bi) dbi
i |
f (yo
i )
gi, ji)f (bi) dbi
f (yo
i |
R
Sub-class MAR
f (yo
i |
ji, `i)f (ym
i |
∪
yo
i , mi)f (ri|
ji, ni)
Introduction to Longitudinal Data Analysis
697
39.13 Slovenian Public Opinion Survey: An MNAR Model
Family
Baker, Rosenberger, and DerSimonian (1992)
Counts: Yr1r2jk
Questions: j, k = 1, 2
Non-response: r1, r2 = 0, 1
§	E(Y11jk) = mjk
E(Y10jk) = mjkβjk
E(Y01jk) = mjkαjk
E(Y00jk) = mjkαjkβjkγjk
·	αjk: non-response on independence question
·	βjk: non-response on attendance question
·	γjk: interaction between both non-response indicators
Introduction to Longitudinal Data Analysis
698
39.14 A More Formal Look
Statistical Uncertainty
@
(cid:0)
(cid:0)
@
(cid:0)
(cid:0)
@
@
@
(cid:0)
(cid:0)
(cid:0)
(cid:0)(cid:9)
@
@
@R
Statistical Imprecision
Statistical Ignorance
Introduction to Longitudinal Data Analysis
699
Statistical Imprecision: Due to ﬁnite sampling
Fundamental concept of mathematical statistics
Consistency, eﬃciency, precision, testing,. . .
Disappears as sample size increases
§	Statistical Ignorance: Due to incomplete observations
Received less attention
Can invalidate conclusions
Does not disappear with increasing sample size
§	Kenward, Goetghebeur, and Molenberghs (StatMod 2001)
Introduction to Longitudinal Data Analysis
700
39.14.1 Monotone Patterns
R = 1
Y1,11 Y1,12
Y1,21 Y1,22
↑
R = 1
Y1,11 Y1,12
Y1,21 Y1,22
R = 0
Y0,1
Y0,2
↑
R = 0
Y0,11 Y0,12
Y0,21 Y0,22
Introduction to Longitudinal Data Analysis
701
39.14.2 Models for Monotone Patterns
R = 1
Y1,11 Y1,12
Y1,21 Y1,22
↑
R = 1
Y1,11 Y1,12
Y1,21 Y1,22
R = 0
Y0,1
Y0,2
↑
R = 0
Y0,11 Y0,12
Y0,21 Y0,22
µr,ij = pijqr
ij,
|
(i,j=1,2;r=0,1)
Model
1.	MCAR
2.	MAR
3.	MNAR(0)
4.	MNAR(1)
5.	MNAR(2)
qr
ij
|
qr
qr
i
|
qr
j
|
logit(qr
|
ij) = α + βi + γj
qr
ij
|
Par.
Observed d.f.
Complete d.f.
4
5
5
6
7
Non-saturated
Non-saturated
Saturated
Saturated
Non-saturated
Non-saturated
Overspeciﬁed
Non-saturated
Overspeciﬁed
Saturated
Introduction to Longitudinal Data Analysis
702
39.14.3 Sensitivity Parameter Method
Sensitivity Parameter: A minimal set η
Estimable Parameter: µ, estimable, given η
Procedure:
·	Given η, calculate parameter and C.I. for µ
·	Set of parameter estimates: region of ignorance
·	Set of interval estimates: region of uncertainty
·	Single parameter case: ‘region’ becomes ‘interval’
Introduction to Longitudinal Data Analysis
703
39.15 Side Eﬀects: Monotone Patterns
Parameter
Model 1/2 Model 3
Model 4
Model 5
1st Margin
2nd Margin
Log O.R.
O.R.
II
IU
II
IU
II
IU
II
IU
0.43 0.43
0.43 0.43
[0.37;0.48]
[0.37;0.48]
[0.37;0.48]
[0.37;0.48]
0.64 0.59
[0.49;0.74]
[0.49;0.74]
[0.58;0.70]
[0.53;0.65]
[0.43;0.79]
[0.43;0.79]
2.06 2.06
[1.52;2.08]
[0.41;2.84]
[1.37;2.74]
[1.39;2.72]
[1.03;2.76]
[0.0013;2.84]
7.81 7.81
[4.57;7.98]
[1.50;17.04]
[3.95;15.44]
[4.00;15.24]
[2.79;15.74]
[1.0013;32.89]
Introduction to Longitudinal Data Analysis
704
39.16 Side Eﬀects: Non-Monotone Patterns
Marg. Prob.
Odds Ratio
Model
d.f. G2
P
First
Second
Orig.
Log
BRD1
BRD2
BRD3
BRD4
BRD7
BRD9
Model 10:II
Model 10:IU
6
7
7
7
8
8
9
9
4.5 0.104 0.43[0.37;0.49] 0.64[0.58;0.71] 7.80[3.94;15.42] 2.06[1.37;2.74]
1.7 0.192 0.43[0.37;0.48] 0.64[0.58;0.70] 7.81[3.95;15.44] 2.06[1.37;2.74]
2.8 0.097 0.44[0.38;0.49] 0.66[0.60;0.72] 7.81[3.95;15.44] 2.06[1.37;2.74]
1.7 0.192 0.43[0.37;0.48] 0.58[0.49;0.68] 7.81[3.95;15.44] 2.06[1.37;2.74]
0.0 0.0
0.0 0.0
·	0.44[0.38;0.49] 0.61[0.53;0.69] 7.81[3.95;15.44] 2.06[1.37;2.74]
·	0.43[0.38;0.49] 0.66[0.60;0.72] 7.63[3.86;15.10] 2.03[1.35;2.71]
[0.425;0.429]
[0.47;0.75]
[4.40;7.96]
[1.48;2.07]
[0.37;0.49]
[0.41;0.80]
[2.69;15.69]
[0.99;2.75]
Introduction to Longitudinal Data Analysis
705
39.17 Slovenian Public Opinion Survey: Identiﬁable Models
Model
BRD1
BRD2
BRD3
BRD4
BRD5
BRD6
BRD7
BRD8
BRD9
Structure
d.f.
(α, β)
(α, βj)
(αk, β)
(α, βk)
(αj, β)
(αj, βj)
(αk, βk)
(αj, βk)
(αk, βj)
6
7
7
7
7
8
8
8
8
loglik
-2495.29
-2467.43
-2463.10
-2467.43
-2463.10
-2431.06
-2431.06
-2431.06
-2431.06
θ
0.892 0.884
0.881 0.765
0.844 0.819
0.764 0.741
0.867 C.I.
[0.878;0.906]
[0.869;0.900]
[0.866;0.897]
[0.674;0.856]
[0.806;0.882]
[0.788;0.849]
[0.697;0.832]
[0.657;0.826]
[0.851;0.884]
Introduction to Longitudinal Data Analysis
706
39.18 SPO: An MNAR “Interval”
θ =0.885
Estimator
[Pessimistic; optimistic]
Complete cases
Available cases
MAR (2 questions)
MAR (3 questions)
MNAR
MNAR “interval”
c
θ
[0.694;0.904]
0.928 0.929
0.892 0.883
0.782 [0.741;0.892]
Introduction to Longitudinal Data Analysis
707
39.19 SPO: Interval of Ignorance
Model
BRD1
BRD2
BRD3
BRD4
BRD5
BRD6
BRD7
BRD8
BRD9
Model 10
Model 11
Model 12
Structure
(α, β)
d.f.
6
(α, βj)
(αk, β)
(α, βk)
(αj, β)
(αj, βj)
(αk, βk)
(αj, βk)
(αk, βj)
(αk, βjk)
(αjk, βj)
7
7
7
7
8
8
8
8
9
9
loglik
-2495.29
-2467.43
-2463.10
-2467.43
-2463.10
-2431.06
-2431.06
-2431.06
-2431.06
-2431.06
θ
0.892 0.884
0.881 0.765
0.844 0.819
0.764 0.741
0.867 [0.762;0.893]
C.I.
[0.878;0.906]
[0.869;0.900]
[0.866;0.897]
[0.674;0.856]
[0.806;0.882]
[0.788;0.849]
[0.697;0.832]
[0.657;0.826]
[0.851;0.884]
[0.744;0.907]
-2431.06
[0.766;0.883]
[0.715;0.920]
(αjk, βjk)
10
-2431.06
[0.694;0.904]
Introduction to Longitudinal Data Analysis
708
39.20 Every MNAR Model Has Got an MAR Counterpart
Molenberghs, Beunckens, Sotto, and Kenward (JRSSB 2008)
Creemers, Hens, Aerts, Molenberghs, Verbeke, and Kenward (2008)
Fit an MNAR model to a set of incomplete data
Change the conditional distribution of the unobserved outcomes, given the
observed ones, to comply with MAR
Resulting new model has exactly the same ﬁt as the original MNAR model
The missing data mechanism has changed
This implies that deﬁnitively testing for MAR versus MNAR is not possible
o	Introduction to Longitudinal Data Analysis
709
39.21 MAR Counterpart to Pattern-mixture Models
f (yi
o, yi
m, ri| c
θ,
ψ) = f (yi
d
ri,
o
|
↓
h(yi
o, yi
m, ri| c
θ,
ψ) = f (yi
d
ri,
o
|
θ) f (ri|d
c
ψ) f (yi
θ) f (ri|d
c
ψ) f (yi
yi
o, ri,
m
|
θ)
c
yi
o,
m
|
θ,
ψ)
c
d
·	Starting from PMM is “natural”: clear separation into:
·	fully observable components
·	entirely unobserved component
Introduction to Longitudinal Data Analysis
710
39.22 MAR Counterpart to Selection Models
f (yi
o, yi
m, ri| c
θ,
ψ) = f (yi
m
o, yi
d
θ) f (ri|
| c
yi
o, yi
m,
ψ)
d
↓
f (yi
o, yi
m, ri| c
θ,
ψ) = f (yi
d
o
|
ri,
θ,
c
θ,
ψ) f (ri| c
d
ψ) f (yi
m
d
yi
o, ri,
|
θ,
ψ)
c
d
↓
h(yi
o, yi
m, ri| c
θ,
ψ) = f (yi
d
o
|
ri,
θ,
c
θ,
ψ) f (ri| c
d
ψ) f (yi
m
d
yi
o,
|
θ,
ψ)
c
d
Introduction to Longitudinal Data Analysis
711
39.23 MAR Counterpart to Shared-parameter Models
f (yi
o, yi
m, ri|
bi) = f (yo
i |
gi, hi, ji, `i) f (ym
i |
yo
i , gi, hi, ki, mi) f (ri|
gi, ji, ki ni)
↓
h(yi
o, yi
m, ri|
bi) = f (yo
i |
gi, hi, ji, `i) h(ym
i |
yo
i , mi) f (ri|
gi, ji, ki ni)
with
h(ym
i |
yo
i , mi) =
Zgi Zhi Zki
f (ym
i |
yo
i , gi, hi, ki, mi)dgidhidki
Introduction to Longitudinal Data Analysis
712
39.24 Slovenian Public Opinion Survey: Counterpart Added
Structure
d.f.
loglik
Model
BRD1
BRD2
BRD3
BRD4
BRD5
BRD6
BRD7
BRD8
BRD9
Model 10
Model 11
(α, β)
(α, βj)
(αk, β)
(α, βk)
(αj, β)
(αj, βj)
(αk, βk)
(αj, βk)
(αk, βj)
(αk, βjk)
(αjk, βj)
6
7
7
7
7
8
8
8
8
9
9
-2495.29
-2467.43
-2463.10
-2467.43
-2463.10
-2431.06
-2431.06
-2431.06
-2431.06
θ
0.892 0.884
0.881 0.765
0.844 0.819
0.764 0.741
0.867 C.I.
[0.878;0.906]
θMAR
b
0.8920 [0.869;0.900]
0.8915 [0.866;0.897]
0.8915 [0.674;0.856]
0.8915 [0.806;0.882]
0.8915 [0.788;0.849]
0.8919 [0.697;0.832]
0.8919 [0.657;0.826]
0.8919 [0.851;0.884]
0.8919 -2431.06
[0.762;0.893]
[0.744;0.907]
0.8919 -2431.06
[0.766;0.883]
[0.715;0.920]
0.8919 Model 12
(αjk, βjk)
10
-2431.06
[0.694;0.904]
0.8919 Introduction to Longitudinal Data Analysis
713
39.25 Slovenian Public Opinion Survey: Incomplete Data
Observed
BRD7
BRD9
≡
≡
≡
BRD7(MAR)
BRD9(MAR):
≡
1439
16
78
16
159
32
144
54
136
BRD1
≡
BRD1(MAR):
1381.6 101.7
24.2 41.4
182.9 8.1
179.7 18.3
136.0 BRD2
≡
BRD2(MAR):
1402.2 108.9
15.6 22.3
159.0 32.0
181.2 16.8
136.0 Introduction to Longitudinal Data Analysis
714
39.26 Slovenian Public Opinion Survey: Complete-data
Prediction
BRD1
≡
BRD1(MAR):
1381.6 101.7
41.4 24.2
170.4 12.5
5.1 3.0
176.6 13.0
5.3 3.1
121.3 9.0
2.1 3.6
BRD2:
BRD2(MAR):
BRD7:
BRD9:
BRD7(MAR)
≡
BRD9(MAR):
1402.2 108.9
22.3 15.6
147.5 11.5
13.2 18.8
179.2 13.9
2.9 2.0
105.0 8.2
9.4 13.4
1402.2 108.9
22.3 15.6
147.7 11.3
13.3 18.7
177.9 12.5
4.3 3.3
121.2 9.3
2.3 3.2
1439
16
1439
16
1439
16
78
16
78
16
78
18
3.2 155.8
32.0 0.0
142.4 44.8
9.2 1.6
0.4 112.5
23.1 0.0
150.8 8.2
16.0 16.0
142.4 44.8
9.2 1.6
66.8 21.0
7.1 41.1
148.1 10.9
11.8 20.2
141.5 38.4
2.5 15.6
121.3 9.0
2.1 3.6
Introduction to Longitudinal Data Analysis
715
39.27 Slovenian Public Opinion Survey: Collapsed
(Marginalized) Predictions
BRD1
≡
BRD1(MAR):
BRD2:
BRD2(MAR):
BRD7:
BRD9:
BRD7(MAR)
≡
BRD9(MAR):
1849.9 136.2
55.4 32.4
1833.9 142.5
57.5 40.2
1849.0 142.0
48.5 34.5
1585.0 391.1
80.3 17.6
1799.7 152.0
82.3 40.7
1849.9 136.3
57.4 30.4
=
⇒
=
⇒
=
⇒
=
⇒
=
⇒
=
⇒
θ = 89.2%
b
θ = 88.4%
b
θ = 89.2%
b
θ = 76.4%
b
θ = 86.7%
b
θ = 89.2%
b
Introduction to Longitudinal Data Analysis
716
39.28 Toenail Data: Unaﬀected Nail Length
We opt for the following SPM:
·	logit [P (Rij = 1
|
Ri,j
gi, Ti, tj, β) = β0 + gi + β1Ti + β2tj + β3Titj
E(Yij|
1 = 0, gi, Ti, tj, γ)] = γ0 + γ01gi + γ1Ti + γ2tj + γ3Titj
−
with
o	Yij: unaﬀected nail length for subject i at occasion j
·	tj: time at which the jth measurement is made
·	Ti: treatment indicator for subject i
·	gi: normal random eﬀect
Introduction to Longitudinal Data Analysis
717
Parameter estimates (standard errors):
·	Unaﬀected nail length
Dropout
Eﬀect
Par. Estimate (s.e.)
Par. Estimate (s.e.)
Mean structure parameters
Intercept
Treatment
Time
Treatment-by-time
β0
β1
β2
β3
2.510 (0.247)
0.255 (0.347)
0.558 (0.023)
0.048 (0.031)
γ0
γ1
γ2
γ3
-3.127 (0.282)
-0.538 (0.436)
0.035 (0.041)
0.040 (0.061)
Variance-covariance structure parameters
Residual variance
Scale factor
Rand. int. variance
σ2
τ 2
6.937(0.248)
6.507 (0.630)
γ01
01τ 2
γ2
-0.076 (0.057)
0.038 (0.056)
Introduction to Longitudinal Data Analysis
718
Graphical representation of predictions for incomplete portions:
o	MNAR model: Y m
i |
yo
i , gi ∼
N (Xiβ + Z m
i gi, σ2Ii)
·	MAR counterpart: Y m
i |
yo
i ∼
N (Xiβ, dJi + σ2Ii)
(dashed lines)
(solid lines)
Introduction to Longitudinal Data Analysis
719
39.29 Toenail Data: Severity of Infection
πgi1i2rt = πg
πi1|
g ·
πi2|
·
i1gt ·
πr
g
|
Variable
Index
0
Complete ﬁrst measurement
Incomplete last measurement
Dropout indicator
Treatment arm
Latent class
i1
i2
r
t
g
non-severe
non-severe
1
severe
severe
dropout
completer
standard
experimental
class 0
class 1
Introduction to Longitudinal Data Analysis
720
39.30 Toenail Data: Severity of Infection
πg =
eαg
1 + eα
g =
πi1|
e(β0+β1g)i1
1 + eβ0+β1g
πi2|
i1gt =
e(γ0+γ1i1+γ2g+γ3i1g+γ4t)i2
1 + eγ0+γ1i1+γ2g+γ3i1g+γ4t
g =
πr
|
e(δ0+δ1g)r
1 + eδ0+δ1g
Model
Restriction
Mechanism
Implication
Bin1
Bin2
β1 = 0
MNAR
Bin1
= Bin1(MAR)
γ2 = γ3 = 0
MAR
Bin2 = Bin2(MAR)
Introduction to Longitudinal Data Analysis
721
6
Standard treatment
Experimental treatment
Completers
Dropouts
Completers
Dropouts
77
42
5
9
Observed data
10
3
79
42
3
3
Fit of Model ‘Bin1’
76.85 40.60
5.66 7.99
9.04 4.62
0.34 0.90
9.38 5.52
81.21 45.62
2.43 3.63
9.36 5.19
0.15 0.41
Fit of Model ‘Bin1(MAR)’
77.12 40.61
5.39 7.98
8.77 4.62
0.61 0.91
9.38 5.52
81.32 45.63
2.32 3.63
9.24 5.18
0.26 0.41
11
6
9.51 5.60
9.51 5.59
Fit of Model ‘Bin2’
‘Bin2(MAR)’
75.86 41.50
5.58 8.15
9.72 3.74
0.72 0.73
10.44 4.47
2.40 3.72
10.27 4.20
0.31 0.34
10.58 4.53
≡
80.16 46.61
Introduction to Longitudinal Data Analysis
722
39.31 Conclusion: Correspondence Between Model Families
Molenberghs, Michiels, Kenward, and Diggle (Statistica Neerlandica 1998)
Kenward, Molenberghs, and Thijs (Biometrika 2003)
Creemers, Hens, Aerts, Molenberghs, Verbeke, and Kenward (2008)
SeM : MCAR
PMM : MCAR
l
SPM : MCAR
l
⊂
⊂
⊂
MAR
l
ACMV
l
Theorem 1
∪
Subfamily 1
NFD
l
NFMV
interior
l
Theorem 2
∪
Subfamily 2
⊂
⊂
⊃
⊂
⊂
⊂
⊂
⊂
⊂
general MNAR
l
general MNAR
l
general MNAR
Introduction to Longitudinal Data Analysis
723
6
39.32 Conclusion: Counterparts to Models
Molenberghs, Beunckens, Sotto, and Kenward (JRSSB 2008)
Creemers, Hens, Aerts, Molenberghs, Verbeke, and Kenward (2008)
Verbeke and Molenberghs (2008)
MNAR model =
⇒
MAR model:
·	Observed data: same ﬁt
·	Unobserved data given observed data: MAR prediction
Holds more generally:
o	Introduction to Longitudinal Data Analysis
724
39.33 Conclusion: Counterparts to Models
Enriched data
Coarse data
Augmented data
Incomplete data
Censored data
Grouped data
Random eﬀects
Latent classes
Latent variables
Mixtures
Introduction to Longitudinal Data Analysis
725
Chapter 40
Pattern-mixture Models
·	A selection model for the vorozole study
·	Initial pattern-mixture models for the vorozole study
·	Principles of pattern-mixture models
·	Connection between selection models and pattern-mixture models
Introduction to Longitudinal Data Analysis
726
40.1 The Vorozole Study
§	open-label study in 67 North American centers
postmenopausal women with metastatic breast cancer
452 patients, followed until disease progression/death
two groups:
vorozole 2.5 mg
1
×
←→
megestrol acetate 40 mg
4
×
several outcomes: response rate, survival, safety,. . .
focus: quality of life: total Function Living Index: Cancer (FLIC)
a higher score is more desirable
Introduction to Longitudinal Data Analysis
727
40.2 A Selection Model for the Vorozole Study
Eﬀect
Parameter Estimate (s.e.)
Fixed-Eﬀect Parameters:
Treatment eﬀect: p = 0.5822
time
time
∗
time
∗
time2
baseline
treatment
time2
∗
baseline
Variance Parameters:
Random intercept
Serial variance
Serial association
Measurement error
β0
β1
β2
β3
β4
d
τ 2
λ
σ2
7.78 (1.05)
-0.065 (0.009)
0.086 (0.157)
-0.30 (0.06)
0.0024 (0.0005)
105.42 77.96
7.22 77.83
Introduction to Longitudinal Data Analysis
728
40.2.1 The Dropout Model
MAR
MNAR
First
0.080(0.341)
Extended 0.033(0.401)
−
0.014(0.003)basei
−
0.013(0.003)basei
−
0.047(0.010) yi,j−1
−
0.033(0.004)yi,j
−
0.023(0.005) yi,j−2+yi,j−1
yi,j−2
1
−
2
−
2
0.53 −
1.38 0.015basei
0.076yi,j
−
0.021basei
0.064yi,j
−
−
−
1 + 0.035yij
−
1 + 0.057yij
−
0.0027yi,j
2
−
Dropout increases with:
low score
negative trend
lower baseline value
§	Introduction to Longitudinal Data Analysis
729
40.3 Pattern-mixture Analysis of the Vorozole Study:
Proﬁles
Introduction to Longitudinal Data Analysis
730
Introduction to Longitudinal Data Analysis
731
40.3.1 Two Pattern-Mixture Models
treatment
Includes: time
∗
pattern
treatment
Includes: time
∗
∗
Assessment of Treatment Eﬀect
Selection model: p = 0.5822 (1 df; output)
PMM1: p = 0.6868 (1 df; output)
PMM2:
p = 0.2403 (13 df; output)
p = 0.3206 (1 df; delta method)
Introduction to Longitudinal Data Analysis
732
40.3.2 Estimating Marginal Eﬀects From PMM
§	Pattern-membership probabilities:
π1, . . . , πt, . . . , πT .
The marginal eﬀects:
Their variance:
where
and
β` =
n
Xt=1
β`tπt,
` = 1, . . . , g
Var(β1, . . . , βg) = AV A0
V = 





Var(β`t)
0
0
Var(πt)






A =
∂(β1, . . . , βg)
∂(β11, . . . , βng, π1, . . . , πn)
Introduction to Longitudinal Data Analysis
733
40.3.3 Considerations
Models ﬁtted over the observation period within a certain pattern
How do we extrapolate beyond dropout time ?
Making the model simple enough ?
Formal identifying restrictions ?
·	. . ?
o	Introduction to Longitudinal Data Analysis
734
40.4 PMM: Three Strategies
(1a)
Simple model per pattern:
Yi = Xiβ (di) + Zibi + εi
bi ∼
εi ∼
N (0, D (di))
N (0, Σi (di))
(1b)
Pattern as covariate:
Yi = Xiβ + Zibi + diθ + εi
(2)
Identifying restrictions:
CCMV: Complete Case Missing Values
ACMV: Available Case Missing Values
NCMV: Neighbouring Case Missing Values
Introduction to Longitudinal Data Analysis
735
40.4.1
Identifying Restrictions
Pattern 3
Pattern 2
Pattern 1
Introduction to Longitudinal Data Analysis
736
Eﬀect
Initial
CCMV
NCMV
ACMV
Pattern 1:
Time
Time∗base
Time∗treat
Time2
Time2∗base
σ11
σ12
σ22
σ13
σ23
σ33
Pattern 2:
Time
Time∗base
Time∗treat
Time2
Time2∗base
σ11
σ12
σ22
σ13
σ23
σ33
Pattern 3:
Time
Time∗base
Time∗treat
Time2
Time2∗base
σ11
σ12
σ22
σ13
σ23
σ33
3.40(13.94)
-0.11(0.13)
0.33(3.91)
131.09(31.34)
53.85(14.12)
-0.46(0.12)
-0.95(1.86)
-18.91(6.36)
0.15(0.05)
170.77(26.14)
151.84(29.19)
292.32(44.61)
29.91(9.08)
-0.26(0.08)
0.82(0.95)
-6.42(2.23)
0.05(0.02)
206.73(35.86)
96.97(26.57)
174.12(31.10)
87.38(30.66)
91.66(28.86)
262.16(44.70)
13.21(15.91)
-0.16(0.16)
-2.09(2.19)
-0.84(4.21)
0.01(0.04)
151.91(42.34)
59.84(40.46)
201.54(65.38)
55.12(58.03)
84.99(48.54)
245.06(75.56)
29.78(10.43)
-0.29(0.09)
-1.68(1.21)
-4.45(2.87)
0.04(0.02)
175.59(27.53)
147.14(29.39)
297.38(46.04)
57.22(37.96)
71.58(36.73)
212.68(101.31)
29.91(9.08)
-0.26(0.08)
0.82(0.95)
-6.42(2.23)
0.05(0.02)
206.73(35.86)
96.97(26.57)
174.12(31.10)
87.38(30.66)
91.66(28.86)
262.16(44.70)
7.56(16.45)
-0.14(0.16)
-1.20(1.93)
-2.12(4.24)
0.03(0.04)
134.54(32.85)
119.76(40.38)
257.07(86.05)
49.88(44.16)
99.97(57.47)
241.99(79.79)
33.74(11.11)
-0.33(0.10)
-1.56(2.47)
-7.00(3.80)
0.07(0.03)
176.49(27.65)
149.05(29.77)
299.40(47.22)
89.10(34.07)
107.62(47.59)
264.57(76.73)
29.91(9.08)
-0.26(0.08)
0.82(0.95)
-6.42(2.23)
0.05(0.02)
206.73(35.86)
96.97(26.57)
174.12(31.10)
87.38(30.66)
91.66(28.86)
262.16(44.70)
4.43(18.78)
-0.11(0.17)
-0.41(2.52)
-0.70(4.22)
0.02(0.04)
137.33(34.18)
97.86(38.65)
201.87(80.02)
61.87(43.22)
110.42(87.95)
286.16(117.90)
28.69(11.37)
-0.29(0.10)
-2.12(1.36)
-4.22(4.20)
0.05(0.04)
177.86(28.19)
146.98(29.63)
297.39(46.04)
99.18(35.07)
166.64(66.45)
300.78(77.97)
29.91(9.08)
-0.26(0.08)
0.82(0.95)
-6.42(2.23)
0.05(0.02)
206.73(35.86)
96.97(26.57)
174.12(31.10)
87.38(30.66)
91.66(28.86)
262.16(44.70)
Introduction to Longitudinal Data Analysis
737
40.4.2 Pattern As Covariate
treat
treat
treat
base
base
base
treat
∗
base
Pattern
1
2
3
1
2
3
1
2
3
1
2
3
treat
base
∗
∗
Eﬀect
Time
Time
Time
Time
∗
Time
∗
Time
∗
Time
∗
Time
∗
Time
∗
Time
∗
Time2
Time2
Time2
Time2
Time2
σ11
σ12
σ22
σ13
σ23
σ33
Estimate (s.e.)
7.29(15.69)
37.05(7.67)
39.40(9.97)
5.25(6.41)
3.48(5.46)
3.44(6.04)
-0.21(0.15)
-0.34(0.06)
-0.36(0.08)
-0.06(0.04)
-9.18(2.47)
-7.70(2.29)
1.10(0.74)
0.07(0.02)
173.63(18.01)
117.88(17.80)
233.86(26.61)
89.59(24.56)
116.12(34.27)
273.98(48.15)
Introduction to Longitudinal Data Analysis
738
40.4.3 Plot for Three Diﬀerent Strategies
Introduction to Longitudinal Data Analysis
739
40.5 Connection SEM–PMM
Molenberghs, Michiels, Kenward, and Diggle (Stat Neerl 1998)
Kenward, Molenberghs, and Thijs (Bka 2002)
Selection Models: f (Di|
Y i, ψ)
←→
f (Y i|
Di, θ) : Pattern-Mixture Models
f(Di)
f(Di
SeM : MCAR
PMM : MCAR
l
⊂
⊂
|
Yi1, . . . , Yi,j−1)
MAR
l
ACMV
f(Di
|
Yi1, . . . , Yi,j−1, Yij)
future
¬
l
future
¬
f(Di
|
⊂
⊂
Yi1, . . . , Yi,j−1, Yij, . . . , Yin)
MNAR
l
general
⊂
⊂
Introduction to Longitudinal Data Analysis
740
Chapter 41
PMM Analysis of the ARMD Trial
Eﬀect
Intercept 4
Intercept 12
Intercept 24
Intercept 52
Treatment 4
Treatment 12
Treatment 24
Treatment 52
Intercept 4
Intercept 12
Intercept 24
Intercept 52
Treatment 4
Treatment 12
Treatment 24
Treatment 52
Parameter
CCMV
ACMV
Estimate (standard error)
54.00(1.47)
52.87(1.68)
48.65(2.00)
44.19(2.14)
-3.11(2.10)
-4.18(2.48)
-4.36(3.83)
-5.04(3.86)
p-values
β11
β21
β31
β41
β12
β22
β32
β42
54.00(1.47)
52.92(1.61)
49.16(1.87)
44.69(2.54)
-3.11(2.10)
-4.07(2.30)
-5.14(3.61)
-2.33(4.93)
β11
β21
β31
β41
β12
β22
β32
β42
—
< .0001
< .0001
< .0001
—
0.092 0.271
0.211 —
< .0001
< .0001
< .0001
—
0.077 0.173
0.647 NCMV
54.00(1.47)
52.86(1.63)
48.77(1.78)
44.00(1.80)
-3.11(2.10)
-4.40(2.42)
-4.19(2.62)
-4.89(2.70)
—
< .0001
< .0001
< .0001
—
0.069 0.110
0.071 Introduction to Longitudinal Data Analysis
741
Chapter 42
Sensitivity Analysis Based on Multiple Imputation
Multiple imputation in its basic form: MAR
Various ways to deviate from this:
·	Apply shift and/or inﬂation factor to imputed data (in some groups)
·	Apply inﬂation factor to imputed data (in some groups)
·	Use a ‘placebo’ rather than an ‘active’ predictive distributions
·	More generally, base predictive distribution on any subset of your choice
·	Use pattern-mixture models with NCMV, CCMV,. . .
Implemented in PROC MI using the MNAR statement.
§	Introduction to Longitudinal Data Analysis
742
SAS code for a shift to the treatment group in the ARMD data:
·	proc mi data=m.armd13 seed=486048 simple out=m.armd13as1
nimpute=10 round=0.1;
title ’Shift multiple imputation’;
class treat;
var lesion diff4 diff12 diff24 diff52;
fcs reg;
mnar adjust (diff12 / shift=10 adjustobs=(treat=’2’));
mnar adjust (diff24 / shift=15 adjustobs=(treat=’2’));
mnar adjust (diff52 / shift=20 adjustobs=(treat=’2’));
by treat;
run;
Introduction to Longitudinal Data Analysis
743
o	SAS code for a subgroup adjustment in the ARMD data:
proc mi data=m.armd13 seed=486048 simple out=m.armd13as2 nimpute=10;
title ’Model multiple imputation’;
class treat;
var lesion diff4 diff12 diff24 diff52;
fcs reg;
mnar model (diff4 / modelobs= (treat=’1’));
mnar model (diff12 / modelobs= (treat=’1’));
mnar model (diff24 / modelobs= (treat=’1’));
mnar model (diff52 / modelobs= (treat=’1’));
run;
Suppose we want to undertake NCMV adjustment:
·	The method can be applied to monotone data only
·	First MI (10 imputations): Start by making the data monotone (under MAR)
·	Second MI (1 imputation): Then apply NCMV to the monotonized data
·	The end result is 10 imputations, as we want
Introduction to Longitudinal Data Analysis
744
SAS code for NCMV in the ARMD data:
·	proc mi data=m.armd13 seed=486048 simple out=m.armd13as3 nimpute=10;
title ’Montone imputation’;
var lesion diff4 diff12 diff24 diff52;
mcmc impute=monotone;
by treat;
run;
proc mi data=m.armd13as3 seed=486048 simple out=m.armd13as4 nimpute=1;
title ’Model multiple imputation’;
var lesion diff4 diff12 diff24 diff52;
monotone reg;
mnar model (diff4 diff12 diff24 diff52 / modelobs=ncmv);
by treat;
run;
Introduction to Longitudinal Data Analysis
745
In the latter case, SAS prints what predictive distributions have been used:
·	Observations Used for Imputation Models Under MNAR Assumption
Observations
Nonmissing lesion, diff4; Missing diff12, ..., diff52
Nonmissing lesion, ..., diff12; Missing diff24, diff52
Nonmissing lesion, ..., diff24; Missing diff52
Complete Cases
Imputed
Variable
diff4
diff12
diff24
diff52
·	Results:
·	GEE
·	GLMM
Introduction to Longitudinal Data Analysis
746
Eﬀect
Par.
MAR
shift
placebo
NCMV
Int.4
Int.12
Int.24
Int.52
Trt.4
Trt.12
Trt.24
Trt.52
Int.4
Int.12
Int.24
Int.52
Trt.4
Trt.12
Trt.24
Trt.52
R.I. s.d.
R.I. var.
β11
β21
β31
β41
β12
β22
β32
β42
β11
β21
β31
β41
β12
β22
β32
β42
τ
τ 2
-0.73(0.20)
-0.71(0.19)
-0.56(0.19)
-0.82(0.20)
0.07(0.28)
0.29(0.27)
-0.10(0.27)
-0.43(0.30)
Generalized estimating equations
-0.82(0.20)
-0.97(0.22)
-1.07(0.23)
-1.66(0.27)
0.17(0.29)
0.56(0.29)
0.41(0.30)
0.41(0.35)
Generalized linear mixed models
-1.46(0.36)
-1.75(0.38)
-1.83(0.38)
-2.71(0.45)
0.32(0.48)
0.99(0.49)
0.67(0.51)
0.53(0.57)
2.21(0.26)
4.90(1.14)
-1.32(0.36)
-1.27(0.35)
-1.01(0.34)
-1.47(0.36)
0.12(0.50)
0.50(0.48)
-0.19(0.48)
-0.74(0.51)
2.28(0.25)
5.21(1.15)
-0.81(0.20)
-0.98(0.22)
-1.05(0.22)
-1.58(0.29)
0.17(0.28)
0.56(0.29)
0.39(0.29)
0.32(0.35)
-1.39(0.35)
-1.67(0.37)
-1.78(0.38)
-2.62(0.46)
0.25(0.48)
0.91(0.48)
0.62(0.48)
0.45(0.56)
2.17(0.25)
4.72(1.09)
-0.83(0.21)
-1.06(0.21)
-1.00(0.22)
-1.59(0.27)
0.17(0.29)
0.67(0.28)
0.34(0.29)
0.32(0.35)
-1.42(0.35)
-1.80(0.36)
-1.70(0.38)
-2.64(0.44)
0.24(0.48)
1.09(0.47)
0.53(0.49)
0.45(0.55)
2.16(0.24)
4.66(1.05)
Introduction to Longitudinal Data Analysis
747
Chapter 43
Overview
MCAR/simple
CC
LOCF
biased
ineﬃcient
MAR
direct likelihood
easy to conduct
direct Bayes
Gaussian & non-Gaussian
not simpler than MAR methods
weighted GEE
MI
variety of methods
strong, untestable assumptions
most useful in sensitivity analysis
Introduction to Longitudinal Data Analysis
748
