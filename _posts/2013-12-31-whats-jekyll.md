---
layout: post
category: spatial
title: What and why spatial statistics?
---

## Spatial Epidemiology
Concerns the analysis of the spatial/geographical distribution of the incidence or prevalence of disease. Also called:

- geographical epidemiology
- environmental epidemiology
- medical geography
- small area health study
- spatial biostatistics

> How to analyse geographical data when we have geographical information available?

### Non-Spatial Analysis
- Using conventional statistical methods
- The results are independent of the spatial arrangement of the geographical entities 
- Observations or entities are assumed to be independent and identically distributed, or in some occasions temporal dependence are also explored.


|     asas    | sasa |       asas       |
|-------------|------|------------------|
| adasdsaasas | sasa | asasadasdasdsads |



|          |   Variable 1  |   Variable 2  |   Variable n  |
|----------|---------------|---------------|---------------|
| Entity 1 | attribute~1~  | attribute~12~ | attribute~1n~ |
| Entity 2 | attribute~21~ | attribute~22~ | attribute~2n~ |
| Entity m | attribute~m1~ | attribute~m2~ | attribute~mn~ |

### Spatial analysis
- When geographical information available, data are called geo-referenced
- Use of spatial statistical methods
- The results depend on the spatial arrangement of the geographical entities 
- It can also include temporal dependence

|          | Geographical |    |   attribute   |               |               |
|----------|--------------|----|---------------|---------------|---------------|
|          | X            | Y  | Variable 1    | Variable 2    | Variable n    |
| Entity 1 | X~1~           | Y1 | attribute~1~  | attribute~12~ | attribute~1n~ |
| Entity 2 | X~2~           | Yg | attribute~21~ | attribute~22~ | attribute~2n~ |
| Entity m | X~m~           | Ym | attribute~m1~ | attribute~m2~ | attribute~mn~ |

### Areas of Applications
Earliest example of spatial epidemiology: Dr. John Snow’s study of spatial distribution of cholera victims. Epidemic around the Broad Street pump in London (1854)

Other example: In an epidemiological investigation, we might wish to analyze lung, breast,
colorectal, and cervical cancer rates by county and year in a particular state, with smoking, mammography, and other important screening and staging information also available at some level.

![Source: Limburg cancer registry data (Belgium)](a1.png)

Researchers in diverse areas such as climatology, ecology, environmental health, and real estate marketing are increasingly faced with the task of analyzing data that are highly multivariate, with many important predictors and response variables, geographically referenced, and often presented as maps and temporally correlated, as in longitudinal or other time series structures.

### Statistical Inference
Public health professionals who collect such data are charged not only with surveillance, but also statistical inference tasks, such as __modeling__ of trends and correlation structures, __estimation__ of underlying model parameters, __hypothesis testing__ (or comparison of competing models), __prediction__ of observations at unobserved times or locations. All naturally accomplished through hierarchical modeling (which can be implemented via Markov chain Monte Carlo (MCMC) methods or other frequentist approaches)

### Spatial Statistics Books
Pioneer Book: Paelinck and Klaassen (1979): which focused the attention of regional scientists on the need for specialized econometric methods to deal with estimation and specification problems caused by spatial data.

Books in the field:  Cressie (1990, 1993): the legendary “bible” of spatial statistics, but rather high mathematical level and lacks modern hierarchical modeling/computing. Wackernagel (1998): terse; only geostatistics. Chiles and Delfiner (1999): only geostatistics. Stein (1999a): theoretical treatise on kriging.

More descriptive presentations: Bailey and Gatrell (1996) focuses on description of pattern, tests of hypotheses, and interpolation, the authors continually stress the role of visualization in understanding spatial phenomena. Fotheringham and Rogerson (1994): deal with the integration of GIS and spatial analysis. Haining (1990):include data description, map interpolation, exploratory and explanatory analyses.

More recent books: Banerjee, S., Carlin, B.P. and Gelfand, A.E. (2004) Hierarchical modeling of analysis for spatial data, CRC Press. Waller, LA, and Gotway, CA. (2004) Applied spatial statistics for public health, Wiley & Sons. Bivand, R.S., Pebesma, E.J., Gomez-Rubio, V. (2008) Applied Spatial Data Analysis with R, Springer. Lawson, A. (2009) Bayesian disease mapping. Hierarchical modeling in spatial epidemiology, Chapman & Hall. Gelfand, A.E., Diggle, P., Fuentes, M. and Guttorp, P. (2011) Handbook of Spatial Statistics, CRC press. Cressie, N. and Wikle, C. (2011) Statistics for spatio-temporal data.