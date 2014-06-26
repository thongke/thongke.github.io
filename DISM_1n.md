# What and why spatial statistics?
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

## Properties and Nature of Spatial Data and Spatial Process
What is Special about Spatial Data?

- Location, Location, Location: where matters 
- Spatial Dependence is the rule
    - spatial interaction, contagion, externalities, spill-overs, copycatting
    - The first law of geography (Tobler): everything depends on everything else, but closer things more.
- Spatial Heterogeneity (or Non-stationarity)
    - The second law of geography (a law of spatial heterogeneity): conditions vary (“smoothly”) over the Earth's surface
- Pertains to the spatial or regional differentiation observed in the value of a variable: Spatial drift (e.g., a trend surface), spatial association.
- The properties of geographical data present a fundamental challenge to classical (conventional) statistics.
- They violate classical assumptions of independence and homogeneity (stationarity) and render classical methods inefficient or inappropriate.

### Nature of Spatial Data
Spatially referenced data “georeferenced”: “attribute” data associated with location: __where matters__. Example, spatial Objects

- points: x, y coordinates => cities, stores, crimes, accidents
- lines: arcs, from node, to node => road network, transmission lines
- polygons: series of connected arcs => provinces, cities, census tracts

![](a2.png)

![](a3.png)

### Spatial Process: Basic Definitions
#### Regionalized Variable
Any variable distributed in space is said to be “regionalized” or “spatial”. Some examples of regionalized variables are: 

1. The price of gold on the NYSE in time (one dimension - longitudinal data) 
2. Total snowfall during January in the province Luxembourg (two dimensions)
3. Magnesium content in the soil (in ppm), measured at varying depths in an agricultural field (three dimensions)

We view a regionalized variable as a function $f(\mathbf{s})$ which adopts a value at every point *s* in the appropriately defined space. The “site” or “location” **s~i~** is bold to indicate that it may be multidimensional. For example, over some field, we might view **s~i~** as ($x_{1i},x_{2i}$), where $x_1, x_2$ refer to northing and easting coordinates. In general no more than 4 dimensions are used.

#### Random Function
The regionalized variables in general will possess both random and structured spatial characteristics. Let $y(\mathbf{s_i})$ denote the observed value of the variable of interest (ex: magnesium content) at location $\mathbf{s_i}$. This can be regarded as a particular realization of the random variable $y(\mathbf{s_i})$ at the point $\mathbf{s_i}$. The set of random variables: $\{Y(\mathbf{s}) : \mathbf{s}\in \mathcal{R}\}$ where $\mathcal{R}$ is the region of interest (ex: the agricultural field), is called a random function where

- $Y(\mathbf{s})$ : the data response at location $\mathbf{s} = (x_1; x_2)$ or $\mathbf{s} = (x_1; x_2; x_3)$
- $\mathcal{R} = $ the set of all spatial locations in the study area.

## Classes of Spatial Data
As outlined in Cressie's book, spatial data generally fall into one of three categories:

- Spatially Continuous (Geostatistical or point-referenced) Data
    - $\mathcal{R}$ is a fixed subset of the plane of positive area (2-D) or volume (3-D).
    - $Y(\mathbf{s})$ is a random variable at each of the __infinite continuous__ locations $s \in \mathcal{R}$
- Area (Lattice) Data 
    - $\mathcal{R} = {\mathbf{s_1}; \mathbf{s_2};...;\mathbf{s_n}}$ is a fixed regular or irregular lattice on the plane.
    - $Y(sidam)$ is a random variable at each location $sidam, i =  1;...;n$
- Spatial Point Process Data
    - $\mathcal{R} = {\mathbf{s_1}; \mathbf{s_2};...;\mathbf{s_n}}$  is a random collection of points on the plane.
    - $Y(\mathbf{s})$ is not specified or is a random variable at a location $\mathbf{s}\in \mathcal{R}$ (marked point process).


### Geostatistical (Point-Referenced) Data
The term “geostatistics” was coined by Matheron (1962, 1963) to describe the statistical methodology for examining ore reserves from spatially distributed data in an ore body.

![](Selection_001.png)

More generally, geostatistics refers to data from a random process $\{Y(\mathbf{s}): \mathbf{s}\in \mathcal{R}\}$ where $\mathcal{R}\subset \mathfrak{R}^n $ fixed, and $\mathbf{s}$ is allowed to vary continuously over $\mathcal{R}$. For example, if $\mathcal{R}\subset \mathfrak{R}^n $ is an agricultural field, and we are measuring the magnesium content in the soil(*Y*), then we can measure Y at any point $\mathbf{s}$ within the field $\mathcal{R}$.

- Important concepts: Stationarity, lsotropy and Variogram.

#### Research Questions
- interest focuses on modeling continuous spatial variation across space
- spatial interpolation 
- estimating how spatial dependence varies with distance

#### Methods of Analysis (for geostatistical data)
- __Linear Interpolation__: Based in general on inverse-distance weighting, this method does not account for the variability in the data, measurements errors, distributional assumptions, etc. Although basic linear interpolation methods do __not attempt to account for spatial correlation__ in the data, they can produce __predictions of equivalent quality__ to methods which do account for spatial variability (such as kriging, described below).

- __Variograms__: The variogram __is a function of the distance and relative orientation__ of pairs of points which __describes the degree of correlation__ between such points. Choice of a variogram model effectively assigns a spatial correlation structure to the data. This model can then be fit to the data and used to make predictions for the response variable at any point in the spatial domain $\mathcal{R}$ using kriging, as described below.

- __Kriging__: The method of kriging attempts to __model the variability in the data as a function of distance, through the variogram and make predictions__. “Kriging” is named after D.G. Krige, a South African mining engineer, who developed empirical methods for predicting the spatial distribution of ore grades based on a sample.

### Area (Lattice) Data
Under the assumption that we have data from the spatial process $\{Y(\mathbf{s}): \mathbf{s}\in \mathcal{R}$, lattice data refers to the case where $\mathcal{R}$ is some __countable collection of spatial sites__.

![](a4.png)

In other words, data can only be observed at the sites in $\mathcal{R}$, and all subsequent __inference applied only to those sites__. 

- Important concepts: Spatial Weights and Spatial Autocorrelation

#### Research Questions
- interest focuses on statistical inference
- estimation, specification tests
- hypothesis test on spatial randomness of attributes : value and location

#### Methods of Analysis:
- __Tests for Spatial Autocorrelation__: In these discrete-space data, a common test performed is to assess __whether there is any spatial autocorrelation__ present in the data. Some common tools for studying autocorrelation include the Geary’s C and Moran's I statistics, and their corresponding randomization tests.
- __Median Polish__: Iterative algorithm for removing large-scale row and column trend effects, before analysis on the small-scale spatial dependence is performed.
- __Nearest-Neighbor Models__: These models make use of a spatial Markovian assumption to model the data (i.e.: The value of the random variable at a given site only depends on the values at a specified set of neighboring sites). These types of models lead to a number of so-called “auto”-models, (Autonormal, AutoPoisson, Autologistic, etc.) __where the distribution of the random variable at a given site depends on itself through its dependence with the nearest neighbors__.

### Spatial Point Processes
For this third type of spatial process $\{Y(\mathbf{s}) : \mathbf{s}\in \mathcal{R}\}, \mathcal{R}$ is again considered to be a random index set for the locations of the process, but here, __the data do not consist of realizations of some random variable at a given site__. The __data are the locations of the sites themselves__, and the __collection of all the sites__ is the event of interest. In this way, a measure such as the count of the number of items over any subset of $\mathcal{R}$ might be the key variable.

![](a5.png)

If a spatial point process $\{Y(\mathbf{s}) : \mathbf{s}\in \mathcal{R}\}$ also consists of measurements $\{Z(\mathbf{s}) : \mathbf{s}\in \mathcal{R}\}$ taken at the locations indicated by *Y*, the process is __known as a marked point process__, where the measurements in *Z* are known as the “marks”.

#### Research Questions
- Research Question
- interest focuses on detecting absence of __spatial randomness (cluster statistics)__ 
- clustered points vs dispersed points

#### Methods of Analysis
- __Quadrat Methods__: These methods involve counting the number of events (trees) in subsets of the region of study $\mathcal{R}$, and __comparing the observed frequencies to the expected frequencies under a Poisson process__. The quadrats themselves are __usually taken to be rectangular, but can be any shape__ desired.
- __Kernel Estimation__: Letting $\lambda(.)$ represent the intensity parameter function for the Poisson process, kernel estimation can be used to estimate $\lambda$.
- __Distance Methods__: These methods take advantage of the __precise distances between points__, and often use nearest-neighbor methods to examine the distribution of sites (trees) in a given area. Such methods include the Ripley's K- and L-functions, as well as the F- and G-functions.

# Exploring Areal Unit Data
## Mapping Count Data
A disease map provides instant visual information on variation of that disease throughout space. Naive use of mapping of health indicators can be __misleading__. The choices of shadings, the scaling of the mapped index quantity, the number of risk classes and their delimitation have to be determined with care. This choice depends on the range of variation, the precision of the estimates and the need for comparability over multiple maps. When working with observational data, regions are not necessarily comparable.

### LIKAR data
The Limburg cancer registry contains all records in Limburg from 1996-2005 which have been classified as cancers[^LIKAR]

[^LIKAR]:http://likas.edm.uhasselt.be/

![Cervix cancer](a6.png)

![Disease Number (Urinary bladder cancer among Males)](a7.png)

### Issues with Crude Counts
Interpretation of crude counts is limited: (1) the population density should be accounted for and (2) areas of high density of ‘at risk’ population would tend to yield high incidence of case-events.

Crude counts of mortality/disease cannot be used to express an increased risk in certain regions. The incidence rate (IR, also called ‘crude rate’) expresses the number of new cases of disease occurring in a population of disease-free individuals in a specific period of time

$$IR = \dfrac{\text{Total number of cases observed in study period}}{\text{Total number of people at risk}}$$

![Cervix Cancer: Crude Counts](a8.png)

![Bladder Cancer: Crude Counts](a9.png)

### Age-specific disease rates
Are regions comparable? Most diseases affect peop.le of certain age disproportionally. In general, there is an increasing incidence of cancer with age.

![Cervix Cancer Incidence](a10.png)

![Urinary Bladder Cancer Incidence](a11.png)

The age-distribution in different regions is not the same. We expect more cases for the region with more residents in the higher-age (and higher-risk) categories.


![Age Distribution Males Age Distribution Males](q1.png)
![Age Distribution Males Age Distribution Females](q2.png)

### Making Rates Comparable: Standardization
Incidence rates reflect estimated average risks for a study population. Thus, populations containing more people in higher age ranges will have higher summary incidence rates than those of younger populations. As a result, the incidence rates for two regions may appear different, but this difference may be due to the different age distributions within the regions rather than to a difference in the underlying age-specific risk of disease. Conclusion: the spatial variation of background population should also be accounted for.

Rate standardization is a mechanism to adjust summary rates to remove the effect of known risk factors (such as age) and make rates from different populations comparable. The idea of rate standardisation is to select a standard population, and adjust observed rates to reflect the age distribution within the standard population. Possible standard populations includes a superpopulation containing the study population e.g. the Belgian, European or World population or the total subpopulation (if we are interested in standardising some subset of individuals, e.g., from a particular region within the study area).

Assume we have $G$ age groups $g (g =1,...,G)$, $y_{gi}$ is number of cases in age group $g$ for study population $i$. $n_{gi}$ is number of people at risk in age group $g$ for study population $i$. $r_{gi} = y_{gi}/n_{gi}$  is the observed incidence proportion in age group $g$ for the study population $i$. $y_{i} = \sum_g y_{gi}$, $n_{i} = \sum_g n_{gi}$. $y_g^S, n_g^S, r_g^S, y^S, n^S$ denote the same quantities for the standard population.

#### Direct Standardisation
How many cases would we __observe in the standard population__ if the observed age-speciﬁc rates of disease applied? {Mausner and Kamer, 1985).

- The expected number of cases in age group $g$ for the standard population is 
$$E^S_{gi} = r_{gi} n^S_g = \dfrac{y_{gi}}{n_{gi}} n^S_g$$
The overall expected rate for the standard population is
$$\dfrac{E^S_i}{n^S} = \dfrac{\sum_g E^S_{gi}}{\sum_g n^S_g}$$
This is the *directly standardised rate*: a weighted average of the observed age-specific rates in the study population, with weights corresponding to the numbers at risk in each group within the standard population.
- _Comparative mortality ﬁgure_ (CMF) is defined as
$$CMF_i = \dfrac{E_i^S}{y^S} = \dfrac{\sum_g E_{gi}^S}{\sum_g y^S_g}$$
This is related to the directly standardised rate:
$$\dfrac{E^S_i}{n^S} = CMF_i \times \dfrac{\sum_g y_{g}^S}{\sum_g n^S_g} = CMF_i \times r^S$$
Thus, the directly standardised rate is equal to the CMF multiplied by the crude
incidence proportion in the standard population.
- Direct standardisation requires the following data (Inskip 1988)
    - Age-specific rates for the study population observed, $r_{gi}$
    - Number of people at risk in the standard population $n_g^S$
    - Total number of cases observed in the standard population $\sum_g y^S_g$

Direct standardisation requires accurate assessment of the age-specific incidence proportions for the observed study population $r_{gi}$. For rare diseases, $r_{gi}$ may be statistically unstable, especially when $n_g$ is small. If $n_g$ is small, addition of a single case could drastically change the value of $r_{gi}$. Also, age-specific incidence counts $y_g$, may not be as readily available as the total observed incidence counts and the age-specific population counts within the study population observed, due to confidentiality issues.

#### Indirect Standardisation
What would be the __number of cases expected in the study population__ if people in the study population contracted the disease at the same rate as people in the standard population? (Mausner and Kramer, 1985).

- The expected number of cases in age group $g$ for the study population is
$$E_{gi} = r_g^S n_{gi} = \dfrac{y_{g}^S}{n_{g}^S} n_{gi}$$
The overall _expected number of cases_ for the study population is
$$E_i  = \sum_g E_{gi} $$
The _standardised mortality ratio_ (SMR) is defined as
$$SMR_i = \dfrac{y_i}{E_i}$$
When referring to incidence rather than mortality, this is called standardised morbidity ratio (SMR) or standardised incidence ratio (SIR). $SMR > 1$ indicate more cases observed in the study population than expected based on the age-specific incidence proportions from the standard population. These SMRs are often the basis for atlases of disease risk (Pickle et al. 1999)
- The _indirectly standardized rated_ is the product of the SMR and the crude rate in the standard population:
$$SMR_i \times \dfrac{\sum_g y^S_g}{\sum_g n^S_g}$$
- Indirect standardisation requires the following data (Inskip 1988)
    - Age-specific rates for the standard population $r_g^S$
    - Number of people at risk in the study population $n_{gi}$
    - Total number of cases observed in the study population $y_{i}$

More stable, requires least information. Often, aggregate (marginal) standards are used in indirect standardisation. The SMR for region $i$ is:
$$SMR_i = \dfrac{y_i}{e_i} = \dfrac{\sum_g y_{ig}}{\sum_g E_{ig}} = \dfrac{\sum_g y_{ig}}{\sum_g n_{ig} (\frac{\sum_i y_{ig}}{\sum_i n_{ig}})}$$
This is called _internal standardisation_. Internal standardisation is in some sense "cheating", since we are "losing a degree of freedom" by estimating the rate $\frac{\sum_i y_{ig}}{\sum_i n_{ig}}$ from our current data. Alternatively, one can refer to an existing standard table of age-adjusted rates for the disease (which are available for many types of cancer). This is called external standardisation.

![Cervix Cancer: indirect standardisation](q3.png)

![Bladder Cancer: indirect standardisation](q4.png)

### Summary
__Direct Age Standardisation__: Adjust observed rates to reflect rates that we would observe if the observed age-specific rates applied to the standard population. Apply the observed age-specific rates directly to the standard population.

__Indirect Age Standardisation__: Adjust observed rates to reflect rates that we would observe if the population standard’s age-specific rates applied to the study population. Apply the age-specific rates from the standard population to estimate indirectly the numbers of cases expected in each age group in the observed study population.

The choice between direct and indirect standardisation often reduces to the type of data available. Often, the standard population is much larger than the study population, and the values $r_j^S$  may provide more stable estimates than $r_j = y_j/n_j$. For indirect standardisation, an internal standard population (superpopulation containing the regions of interest) is preferred as compared to an external standard population (entirely seperate population), for comparability of the regions (Breslow and Day 1975).

In addition to age standardisation, we could also standardise rates to compensate for other risk strata (e.g. gender, race, etc).
>Always be cautious with interpreting and comparing
maps!

## Making Chloropleth maps
To visualize attribute data of spatial tessellations, a map called “Choropleth map"is used. Choropleth map is a map showing attribute data of a spatial tessellation by colors and textures. Categorical variable is directly visualized by colors and textures: different colors indicate different categories. Numerical variable is first categorized into several classes, and then visualized by a progression of colors and textures. To make a Choropleth map of a numerical variable, we determine (1) classification scheme (2) class number (3) class boundaries and (4) colors and textures.

### Classification of Schemes
There are four schemes for categorizing numerical variables: (1) Equal interval scheme(2) Quantile scheme(3) Nonuniform scheme(4) Irregular interval scheme.

#### Equal interval scheme
Equal interval scheme categorizes a numerical variable by an equal interval value. If we specify interval value, GIS calculates boundary values, classifies spatial units into categories, and visualizes the categories by different colors or textures.

![](q5.png)

#### Quantile scheme
Quantile scheme categorizes a numerical variable so that every category has the same number of spatial units. In this scheme we give the number of categories instead of interval value. GIS then calculates boundary values, classifies spatial units into categories, and visualizes the categories by different colors or textures.

![](q6.png)

#### Nonequal interval scheme
When we are interested in a certain range of attribute value, we want to use finer categories in the range. In this case we use nonequal interval scheme. Typical examples include monotonically increasing (decreasing) interval scheme, in which the interval monotonically increases (decreases) with attribute value.

![](q7.png)

#### Irregular interval scheme
Frequency distribution of attribute value often shows “breakpoints”, where the frequency suddenly drops. In such a case, we can obtain a natural classification of attribute variable if we take the breakpoints as boundaries of intervals.

![](q8.png)

### Class Number and Colors
In theory, we can use any number of categories in classification of numerical variables. You may think that the more categories you use, the better map you obtain. In practice, however, we can discriminate only a limited number of colors used for visualizing categories. It is not always useful to increase class number.

Class boundary greatly affects the appearance of Choropleth maps that represent numerical variables. We should be careful when classification scheme involves subjective choice of class boundaries, as seen in nonequal and irregular interval schemes.

Black/white, colors, red-to-green, rainbow, light-to-dark colors, 

#### Effect of the size of spatial units
In Choropleth map, map readers tend to pay attention to large polygons while they often overlook small polygons. Attribute data of large polygons are more influential than those of small polygons in the perception of map readers. Large polygons are emphasized as a result of visualization. This often leads to misunderstanding of the spatial distribution of attribute values.

Example: effect of categorization

![Equal Intervals: for fairly uniformly distributed data](q9.png)

![Quantile map: 4 classes, lowest quartile, second quartile, third and highest](q10.png)

![Standard Deviates: divide data into units of standard deviation around the mean](w1.png)

![Log-Scale: interpretation on log-scale](w2.png)

Alternative: independent categorization

- Bi-chromatic range from red to green
- Based on a uniform log-scale division (Knorr-Held and Raser, 2000)
- 7 categories with a flexion zone in yellow centered around the median
- Dark red regions indicate a high risk of disease 
- Dark green indicate a small risk

![](w3.png)

    2.2.3.R

- Visualizing data in the attribute space and geographical space simultaneously
- Useful for exploring spatial stationary (homogeneity) of spatial patterns and process 

![](w4.png)

Several maps can be produced depending on the style used when we define the R function `classlntervals`.

- The `"fixed"` style permits a `"classlnterval"` object to be specified with given breaks
- The `"equal"` style divides the range of the variable into n parts.
- The `"quantile"` style provides quantile breaks
- The `"sd"` style chooses breaks based on pretty of the centred and scaled variables
- The `"hclust"` style uses hclust to generate the breaks using hierarchical clustering
- The `“bclust” `style uses bclust to generate the breaks using bagged clustering;

## Defining the Neighbourhood Structure
 Defining the relationship between regions. Spatial weights define the spatial relationships among spatial objects (e.g., polygons, rasters, points). Spatial weights are used to identify spatial contiguity or neighbourhood of a given object. Spatial matrix (for $n$ objects, there will be $n \times n$ pairs of relationships).

### Spatial Weights
Depending the case we choose (Rook’s or Queen’s case) the weight matrix can be constructed and give a weight equals 1 to those cells in the matrix reflecting contiguity in the space.

![](Selection_002.png)

Example

![](w5.png)

![](Selection_003.png)

Other weights matrices can be used as well. The distance between spatial locations, can be used to construct the weights, in general the weights are $1/h_{ij}^{\alpha}$1 and $\alpha$ is taken to be 1. Other approaches consider a threshold distance $h_t$, and if the distance is larger than $h_t$, then the weight is taken equal to 0.  In previous example, the weight matrix is constructed based on the distances between the centroids of each spatial region.

Example

![](Selection_004.png)

![](Selection_005.png)

    2.3.2 Code in R

In `spdep` package, this is implemented in the `dnearneigh` function, which takes as input a matrix of coordinates and a lower and upper distance bound (as well as, optionally, a variable containing region ids). All points that are within the defined distance band from each other are categorized as neighbors. It is important to note that the distance calculation __applies only to Euclidean__ distance, which requires that the coordinates are projected (not simply latitude, longitude). A slight complication is that the value for the lower and upper distance bound must be specified beforehand. Typically the lower bound is zero, but the upper distance bound should be such that no observations become islands.

## Spatial Autocorrelation
Measuring Dependency: The objective is to measure __how strong__ the tendency is for observations from __nearby regions to be more (or less) alike than observations from regions farther apart__, to judge whether any apparent tendency is sufficiently strong that it is unlikely to be due to chance alone.

### Definition Autocorrelation
Complicated name, simple concept... expresses the amount of spatial dependence (how much proximity matters in spatial data). Correlation is the key notion, it indicates how much two properties vary together. Correlation of a variable with itself through space, is a variable in a location correlated with its values in nearby places? Spatial + auto + correlation.

If there is any systematic pattern in the spatial distribution of a variable, it is said to be spatially autocorrelated. If nearby or neighboring areas are more alike, this is positive spatial autocorrelation. Negative autocorrelation describes patterns in which neighboring areas are unlike. Random patterns exhibit no spatial autocorrelation. 

It measures the extent to which the occurrence of an event in an areal unit __constrains, or makes more probable__, the occurrence of an event in a neighboring areal unit.

#### Why important
Most statistics are based on the assumption that the values of observations in each sample are independent of one another. Positive spatial autocorrelation may violate this, if the samples were taken from nearby areas. Assumption that observations have been selected randomly, not valid. Estimates obtained will be biased and overly precise.

Biased because areas with higher concentration of events will have a greater impact on the model estimate. Overestimate precision because, since events tend to be concentrated, there are actually fewer number of independent observations than are being assumed.

Goals of spatial autocorrelation:

- Measure the strength of spatial autocorrelation in a map
- Test the assumption of independence or randomness

### Global Indicators of Spatial Autocorrelation
Let $Y_i$ denote the response at the i~th~ location  $i= (1;...;n)$. Let $S_{ij}$ be a measure of how similar or dissimilar the responses are at locations $i$ and $j$. Let $W_{ij}$ be a measure of the spatial proximity of locations $i$ and $j$. For future reference, define matrices $Sdam : (S_{ij})$ and $Wdam =(W_{ij})$. Most global indexes of spatial autocorrelation are of the form
$$\dfrac{\sum_{i=1}^n \sum_{j=1}^n W_{ij} S_{ij}}{\sum_{i=1}^n \sum_{j=1}^n W_{ij}}$$
This is also called general cross-product statistics (Mantel, 1967)

### Moran’s I
One of the oldest indicators of spatial autocorrelation (Moran, 1950). Still a defacto standard for determining spatial autocorrelation. Applied to __zones or points__ with continuous variables associated with them. Compares the value of the variable at any one location with the value at all other locations.
$$I = \dfrac{n\sum_{i} \sum_{j} W_{ij} (Y_i-\bar{Y})(Y_j-\bar{Y})}{(\sum_{i} \sum_{j} W_{ij})\sum_i(Y_i-\bar{Y})^2 }$$

Note that $I$ resembles ordinary correlation coefficient (but is not restricted to $|I| < 1$) (covariation divided by variance).

#### Properties Moran’s I
Expected value of $I$ under assumption of independence $(E(I))$ is $-\frac{1}{n-1}$

- If $I > E(I)$ : positive autocorrelation (clustered pattern)
- If $I < E(I)$ : negative autocorrelation (regular pattern)

Testing for Spatial Autocorrelation, can be done using randomisation or Normal approximation. Randomization distribution obtained by reassign data values among the $n$ fixed locations. If $I$ lies in the tails of this distribute, reject assumption of independence. Normal approximation under independence ($n > 25$):
$$Z(I) = \dfrac{I - E(I)}{\sqrt{Var(I)}}$$
$$Var(I) = \dfrac{n^2 \sum_{ij} W_{ij}^2 + 3(\sum_{ij} W_{ij})^2 - n \sum_{i} (\sum_j W_{ij})^2}{(n^2-1) (\sum_{ij} W_{ij})^2}$$
Compare observed $Z$-score to a standard normal distribution.

#### Notes on Moran’s I
Spatial similarity is assessed by deviations of each regional count with overall mean. Does this assess clustering? Regional at-risk population sizes may result in variations. __Replace counts with incidence proportions__, to remove impact of population heterogeneity. This removes heterogeneity in the expected counts, __but local incidence proportions remain heterogeneous because of differences in sample size__. An alternative approach using constant risk hypothesis (Walter, 1992):
$$I_{cr} = \dfrac{\sum_j\sum_j W_{ij} \left(\dfrac{Y_i-rn_i}{\sqrt{rn_i}}\right) \left(\dfrac{y_j-rn_i}{\sqrt{rn_i}}\right)}{\left(\sum_i \sum_j W_{ij}\right)}$$
which assume risk $r$ is constant in every area and compare observed and expected values.

### Geary’s c
A similar statistic is proposed by Geary (Geary, 1954). Interaction is not the cross-product of the deviations from the mean, but the __deviations in intensities__ of each observation location with one another
$$c = \dfrac{(n-1)\sum_{i} \sum_{j} W_{ij} (Y_i-Y_j)^2}{2(\sum_{i} \sum_{j} W_{ij})\sum_i (Y_i-\bar{Y})^2}$$

If regions $i$ and $j$ have similar values $(Y_i-Y_j)^2$ will be small. Scaled by overall variation around mean regional observation Y. Also called, Geary's contiguity ratio.

#### Properties of Geary’s c
Value ranges between 0 and 2

- 0 indicates perfect positive spatial correlation
- 2 indicates perfect negative spatial autocorrelation

If value of any one zone are spatially unrelated to any other zone, the expected
value of $c$ will be 1. Does not provide identical inference because it emphasizes the differences in values between pairs of observations, rather than the covariation between the pairs. Geary's $c$ corresponds to the Durbin-Watson statistic (used to test for serial autocorrelation)

Testing can be done in a similar way as before: using normality assumption or via randomisation. __Moran’s I gives a more global indicator__, whereas the Geary's coefficient $c$ __is more sensitive to differences in small__ neighborhoods. Geary's $c$ can be adjusted in similar way as $I_{cr}$, by replacing $Y_i$, by standardised value.

    2.4.5 Code in R

Test values may contradict, should use as a exploratory tool and checking further by modelling.

### Local Indicators of Spatial Autocorrelation (LISA)
Global indicators of spatial association assess patterns of spatial similarity summarised over the entire study area to detect clustering. 

What if we want to find out local pockets of mutually similar deviations from the overall mean? Use a local indicator, compares local value to that of its neighbours and detect individual clusters.

#### Local Moran’s I index
Basic form of LISA for region $i$ is:
$$\sum_i W_{ij} S_{ij}$$
with $W_{ij}$ row-standardised weights (rows sum to 1)
Local Moran’s I index can be calculated as:

$$I_i = n  \dfrac{(Y_i-\bar{Y}) \sum_j W_{ij} (Y_j - \bar{Y})}{\sum_i(Y_i - \bar{Y})}$$

#### Moran Scatter Plot
A Moran Scatter Plot can be usefull to detect local pattern of spatial association, if the plot is divided in quadrants, then

- In the first quadrant are the locations in which the value is high and the neighbors are also having high values (positive autocorrelation) 
- The second contain locations with low values surrounded by locations with high values (negative autocorrelation)
- The third contain locations with low values and the neighbors have also low values (positive autocorrelation) 
- The last one contain those locations with high values surrounded by locations with low values.

```
2.4.7 Code in R
```

#### Use of LISA
Identify Hot Spots which usually cannot identify with global index.

- significant local clusters __in the absence of global autocorrelation__, some complications in the presence of global autocorrelation (extra heterogeneity)
- significant local outliers: high surrounded by low and vice versa

Indicate Local Instability: 

- local deviations from global pattern of spatial autocorrelation

# Modeling Areal Unit Data
## Aggregate Count Data
In spatial epidemiology, focus is often on tract count data. Common to use the standardised mortality (or morbidity) ratio (SMR). Notation:

- $Y_i$ observed number of cases of disease in county $i$
- $E_i$ expected number of cases of disease in country $i$

$Y_i$ are random, but the $E_i$, are thought of as fixed and known.

![Cervix Cancer in Limburg, Belgium](w6.png)

### Traditional Models and Methods
- The usual model for the $Y_i$ is the Poisson model:
$Y_i \sim Poisson(E_i\theta_i)$
where $\theta$ is the _relative risk_ of disease in region $i$. This corresponds to the log-Likelihood:
$$l = \sum_{i=1}^m Y_i  \ln(E_i\theta_i) - \sum_{i=1}^m E_i\theta_i$$
Maximum likelihood estimator of $\theta_i$ is the SMR $Y_i/E_i$
- The estimate variance is
$$Var(SMR_i) = Var(Y_i)/E_i^2 = \hat{\theta}/E_i = Y_i/E_i^2 $$
Thus, the estimated standard error is $s_i  \sqrt{Y_i}/E_i$
- Wald-based confidence intervals for the SMRs can then be calculated
$$[SMR-1.96 s_i, SMR + 1.96 s_i]$$
This assumes that the SMRs are normally distributed. This is a bit awkward since the data are discrete and SMR can be positive only. 
- Alternatively, assume $\log(SMR_i)$ is roughly normally distributed. Using the delta method, one can find that 
$$Var[\log(SMR_i)] = \dfrac{1}{SMR_i^2} Var(SMR_i) = \dfrac{E_i^2}{Y_i^2}\dfrac{Y_i}{E_i^2} = \dfrac{1}{Y_i}$$
An approximate 95% CI for $\log(SMR_i)$ is thus
$$[\log(SMR_i) - 1.96/\sqrt{Y_i}, \log(SMR_i) + 1.96/\sqrt{Y_i}]$$
Back-transforming gives an approximate 95% CI for $SMR_i$
$$[SMR_i \exp(-1.96/\sqrt{Y_i}), SMR_i \exp(1.96/\sqrt{Y_i})]$$
This is equal to the CI [$SMR_i$/errfac, $SMR_i \times$ errfac] with error factor (Clayton and Hills, 1995)
$$\text{errfac} = \exp\left(z_{1-\alpha/2} \sqrt{\dfrac{1}{Y_i}}\right)$$

#### Excess Risk?
Table: Cervix Cancer in Limburg, Belgium

| Community  | Y  |  E   |  SIR  | lower |   upper    | Community | Y  |  E   | SIR  | lower | upper |
|------------|----|------|-------|-------|------------|-----------|----|------|------|-------|-------|
| ALKEN      |  7 |    5 |  1.41 |  0.67 | 2.96       | KINROOl   |  3 |  5.1 | 0.59 |  0.19 |  1.82 |
| AS         |  5 |  3.2 |  1.54 |  0.64 | 3.7        | KORTES    |  4 |  3.6 | 1.11 |  0.42 |  2.95 |
| BERlNGEN   | 14 | 17.3 |  0.81 |  0.48 | 1.37       | LANAKEN   |  8 | 11.1 | 0.72 |  0.36 |  1.44 |
| BILZEN     |  7 | 13.4 |  0.52 |  0.25 | 1.09       | LEOPOL-   | 10 |  6.4 | 1.56 |  0.84 |   2.9 |
| BOCHOLT    |  6 |  5.2 |  1.15 |  0.52 | 2.56       | LOMMEL    |  9 | 13.6 | 0.66 |  0.34 |  1.27 |
| BORGLOON   |  7 |  5.1 |   138 |  0.66 | 2.89       | LUMMEN    |  6 |  6.3 | 0.96 |  0.43 |  2.13 |
| BREE       |  1 |  6.5 |  0.15 |  0.02 | 1.09       | MAASEIK   | 12 | 10.5 | 1.14 |  0.65 |  2.01 |
| DIEPENBEEK |  7 |  7.7 |  0.91 |  0.43 | 1.91       | MAASM     | 18 | 15.8 | 1.14 |  0.72 |  1.81 |
| DlLSEN     | 18 |  8.2 |   2.2 |  1.39 | 3.49       | MEEUWEN   |  8 |  5.2 | 1.52 |  0.76 |  3.05 |
| GENK       | 28 | 27.6 |  1.01 |   0.7 | 1.47       | NEERPELT  |  6 |  6.9 | 0.87 |  0.39 |  1.93 |
| GINGELOM   |  7 |  3.8 |  1.83 |  0.87 | 3.83       | NIEUWERK  |  2 |  3.1 | 0.65 |  0.16 |   2.6 |
| HALEN      |  1 |    4 |  0.25 |  0.03 | 1.76       | OPGLAB    |  2 |  3.9 | 0.52 |  0.13 |  2.07 |
| HAM        |  6 |  4.2 |  1.41 |  0.64 | 3.15       | OVERPELT  |  5 |  5.7 | 0.88 |  0.37 |  2.12 |
| HAMONT     |  3 |  6.2 |  0.49 |  0.16 | 1.5        | PEER      |  7 |  6.4 | 1.09 |  0.52 |  2.28 |
| HASSELT    | 35 | 33.9 |  1.03 |  0.74 | 1.44       | RIEMST    |  2 |  7.4 | 0.27 |  0.07 |  1.08 |
| HECHTEL    |  4 |  4.9 |  0.82 |  0.31 | 2.19       | SINTVTRUI | 38 | 18.8 | 2.02 |  1.47 |  2.77 |
| HEERS      |  3 |  3.3 |  0.91 |  0.29 | 2.81       | TESSEN    |  4 |  7.2 | 0.55 |  0.21 |  1.47 |
| HERK       |  7 |  5.4 |   1.3 |  0.62 | 2.73       | TONGER    | 19 | 14.9 | 1.28 |  0.81 |     2 |
| HERSTAPPE  |  1 |    0 | 22.18 |  3.12 | **157.48** | VOEREN    |  0 |    2 | O    |       |       |
| HEUSDEN    |  8 | 12.9 |  0.62 |  0.31 | 1.24       | WELLEN    |  2 |  3.1 | 0.64 |  0.16 |  2.57 |
| HOESELT    |  7 |  4.2 |  1.67 |   0.8 | 3.5        | ZONHOV    |  6 |  8.6 | 0.7  |  0.31 |  1.56 |
| HOUTH      |  8 | 12.3 |  0.65 |  0.33 | 1.31       | ZUTEN-    |  2 |    3 | 0.68 |  0.17 |   2.7 |

Suppose we wish to test whether the true relative risk in region $i$ is elevated or not
$$H_0:\theta_i = 1 \text{ versus } H_A:\theta_i > 1$$

Under the null hypothesis, $Y\sim Poisson(E_i)$. Thus the (one-sided) p-value for this test is
$$p = P(X \geq Y_i|E_i) = 1- P(X < Y_i|E_i) = 1 - \sum_{z=0}^{Y_i-1} \dfrac{\exp(-E_i)E_i^z}{z!}$$
If $p < 0.05$ we would typically reject $H_0$ and conclude that there is statistically significant excess risk in region $i$.

#### Drawbacks of the traditional Poisson model:
- Can __yield large changes in estimate with relatively small changes__ in expected value (since they are based on ratio estimators)
- when a (close to) zero expectation is found, the SMR will be very large for any
positive count
- the zero SMRs do not distinguish variation in expected count
- variance of SMR is proportional to $1/e_i$
- it is a saturated estimate of relative risk, and hence not parsimonious

Thus, naive use of disease mapping on rare diseases or small areas can be __very misleading__ (Molli, 1999). Possibly, the most extreme SMRs are those based on only a few cases. On the contrary, the most extreme p-values of tests comparing SMRs to unity or confidence intervals excluding unity may simply identify areas with large populations.

### Spatial Smoothers
So we want to smooth (take out noise)the data? When interested in the expected values we might expect to do some smoothing. To smooth the observed quantity $Y_i$ replace with
$$\widehat{Y}_i = \dfrac{\sum_i w_{ij} Y_{j}}{w_{i+}}$$
More generally, we could include the value actually observed for unit $i$, and revise our smoother to
$$(1 - \alpha)Y_i + \alpha Y_i$$
For $0 < \alpha < 1$, this is a linear (convex) combination of 'shrinkage‘ form.

Finally, we could try model-based smoothing, i.e. based on the mean of the predictive distribution $E(Y_i|\text{data})$. This leads to hierarchical spatial modeling. The basic idea is to “borrow” information from neighboring regions to produce a better estimate of the rate associated with each region. A “Better” estimate means more stable and less noisy. In this way we separate out the “signal” (i.e., the spatial pattern) from the noise.

## Hierarchical Bayesian Methods
Think of the true underlying relative risks $\theta_i$, as random effects, to allow ‘borrowing of strength’ across regions. Appropriate if we want to estimate and map the underlying risk surface. The random effects can be high dimensional, and are couched in a nonnormal (Poisson) likelihood. Thus, hierarchical Bayesian modeling seems natural.

Assume the following model
$$Y_i \sim Poisson(E_i\theta_i)$$
To circumvent shortcomings of the SMR as a relative risk estimator (extra-poisson variation), we want to control the behavior of the $\theta$. This can be done by assuming that the parameters (the true risks) come from a common underlying distribution. This also allow the procedure to 'borrow strength’ across the various regions in order to come up with an improved estimate for the relative risk in each region. This is the idea of random effects and hierarchical modelling.

### Poisson-Gamma Model
Assume that the number of deaths is each area follow a Poisson distribution (likelihood) $y_i \sim Poisson(e_i\theta_i)$. Let the relative risks follow a gamma distribution (random effect, prior) 
$$\theta_i \sim Gamma(a, b)$$
with mean $m_{\theta_i} = a/b$ and variance $v_{\theta_i} = a/b^2$. As a result, the relative risk has the following distribution (posterior) 
$$\theta_i \sim Gamma(a + y_i , b + e_i )$$
This emerges in closed form thanks to the conjugacy of the gamma prior with the Poisson likelihood.

#### Posterior Mean
The posterior mean of the relative risk 67; is given by 
$$E[\theta_i |y_i ] $$
$$= \dfrac{a + y_i}{b + e_i} $$
$$= \dfrac{\dfrac{m^2_{\theta_i}}{v_{\theta_i}} + y_i}{\dfrac{m_{\theta_i}}{v_{\theta_i}} +  e_i} $$
$$= \left(\dfrac{e_i}{\dfrac{m_{\theta_i}}{v_{\theta_i}} +e_i}\right)\dfrac{y_i}{e_i} + \left(\dfrac{\dfrac{m_{\theta_i}}{v_{\theta_i}}}{\dfrac{m_{\theta_i}}{v_{\theta_i}}+e_i}\right)m_{\theta_i} $$
$$= C_i SMR_i + (1 - C_i ) \dfrac{a}{b}$$
The posterior mean of $\theta_i$ is a weighted average of the data-based SMR for the i^th^ area, and the relative risk in the overall map (the prior mean $m_{\theta_i}$).

#### Shrinkage
$$E[\theta_i |y_i ] = C_i SMR_i + (1 - C_i ) \dfrac{a}{b}$$

For __rare diseases and small areas__ $C_i$ is small, posterior mean tends towards a global mean $a/b$ and __prior is highly informative__.

For areas with a lot of data, posterior mean is close to $y_i/e_i$, prior is weakly informative.

__Example__: Suppose we have a $Gamma(4,4)$ prior. In this case, we assume that $E[\theta_i] = 1$, and $Var[\theta_i] = 1/4$. Suppose in area $i$ we observe $y_i$ = 27 cases, when $e_i$ = 21 were expected, or risk $y_i /e_i$ = 1.29. As a result, we obtain as posterior distribution
$$Gamma(y_i + a, e_i + b) = Gamma(31, 25)$$
This distribution has mean 31/25 = 1.24 indicates slightly elevated risk (24%).

![](Selection_006.png)

The probability that the true risk is bigger than 1, is $P (\theta = 1|y_i ) = 0.863$. This is derived exactly from the gamma distributionll. Can also be derived empirically as the proportion of samples that are greater than 1.

A $100(1 - a)\%$ confidence interval for $\theta$ can be derived by taking the upper and lower $\alpha/2$ points of the $Gamma(31, 25)$ posterior. This is called the equal-tail credibility interval. We obtain (0.84127 1.1713).

#### Estimation of the model parameters
$a$ and $b$ in the Poisson-Gamma model are so-called hyperparameters. In Empirical Bayes approach hyperparameters are estimated from the data and yield acceptable point estimates of the rates but underestimates their uncertainty.

In Full Bayesian approach, the prior distributions for $a$ and $b$ are specified as well. This expresses our ignorance or prior knowledge about $a$ and $b$. For example, $a \sim exponential(0.1)$ and $b \sim exponential(0.1)$. Sensitivity analysis to investigate the influence of the choice of hyperprior on estimates of relative risk can be done. If __data are scarce, the choice of a suitable combination of hyperparameters is important__.

    3.2.2 WinBugs Code: Poisson-Gamma model

![Cervix Cancer in Limburg, Belgium](w7.png)
![Cervix Cancer in Limburg, Belgium](w8.png)

Table: Cervix Cancer in Limburg, Belgium, Posterior Parameter Estimates:

| node |  mean  |   sd   | MC error |  2.5%  | median | 97.5%  |
|------|--------|--------|----------|--------|--------|--------|
| a    |   3.57 | 0.9012 |  0.03898 |  2.062 |  3.489 |   5.59 |
| b    |  3.495 | 0.9429 |  0.04048 |  1.903 |  3.408 |  5.595 |
| mean |  1.032 | 0.1062 |  0.00142 | 0.8409 |  1.025 |   1.26 |
| var  | 0.3213 |  0.112 | 0.004055 | 0.1676 | 0.2998 | 0.6033 |


Table: Cervix Cancer in Limburg, Belgium, Posterior Estimates of Relative Risk:

|    node   |  mean  |   sd   | MC error |  2.5%  | median | 97.5%  |
|-----------|--------|--------|----------|--------|--------|--------|
| theta[1]  |  1.254 | 0.3851 |  0.00438 | 0.6168 |  1.214 |  2.118 |
| theta[2]  |  1.279 | 0.4519 |   0.0052 | 0.5546 |  1.221 |  2.318 |
| theta[3]  | 0.8426 | 0.1992 | 0.002209 | 0.4992 | 0.8267 |  1.286 |
| theta[4]  | 0.6217 | 0.1914 | 0.002385 | 0.3043 | 0.6029 |  1.048 |
| theta[5]  |  1.099 | 0.3589 | 0.004094 |  0.515 |  1.053 |  1.921 |
| theta[6]  |  1.233 | 0.3786 | 0.004608 | 0.6102 |  1.192 |  2.076 |
| theta[7]  | 0.4518 |  0.217 | 0.003479 | 0.1262 | 0.4199 | 0.9512 |
| theta[8]  | 0.9475 | 0.2941 | 0.003245 | 0.4549 | 0.9187 |   1.61 |
| theta[9]  |  1.852 | 0.4053 | 0.004673 |  1.146 |  1.828 |  2.747 |
| theta[10] |  1.015 | 0.1811 | 0.002093 | 0.6924 |  1.003 |  1.398 |

![caterpillar plot: Cervix Cancer in Limburg, Belgium](w9.png)

![probability of theta greater than 1.0 ](w10.png)

#### Remarks
Gamma prior is mathematically convenient and leads to robust estimates but, __covariate adjustment is difficult__ and __does not cope with spatial correlation between risks in nearby areas__. Extending this model to allow for such spatial correlations among the $\theta_i$ also difficult (we would need a multivariate version of the gamma distribution). As an alternative, __use a multivariate normal distribution for prior__.

### Poisson-Lognormal Model
he number of deaths is each area follow a Poisson distribution (likelihood) $y_i \sim Poisson(e_i \theta_i )$. Use normal prior distribution on the log-relative risks 
$$log(\theta_i ) = \alpha + x_i \beta + \nu_i,\qquad \nu_i \sim N(0,\sigma^2_{\nu})$$
$\nu_i$ is the heterogeneity random effect, capturing extra-Poisson variability in the log-relative risks. $x_i$ are explanatory spatial covariates (at region-level), having parameter coefficients $\beta$. Lognormal model for the relative risk is more flexible.

Parameters $\nu_i$ are called random effects, this represents the residual (unexplained) (log) relative risk in area $i$ after adjusting for known covariates ($x_i\beta$) and overall mean risk ($\alpha$). $\nu_i$ captures the effects of unknown or unmeasured area level covariates. The variance of the random effects ($\sigma^2_{\nu}$) reflects the amount of extra-Poisson variation in the data.

    3.2.4 WinBugs Code: Poisson-Lognormal model

![Example 1: Cervix Cancer in Limburg, Belgium](e1.png)

![Example 1: Cervix Cancer in Limburg, Belgium](e2.png)

Table: Cervix Cancer in Limburg, Belgium, Posterior Parameter Estimates:

|  node |   mean  |    sd   | MC error |   2.5%  |  median  |  97.5%  |
|-------|---------|---------|----------|---------|----------|---------|
| alpha | -0.0576 | 0.07584 |  0.00141 | -0.2107 | -0.05575 | 0.08602 |
| tau.u |   15.07 |   13.79 |   0.5335 |   4.491 |     11.8 |    46.5 |
| var   | 0.09471 | 0.05194 | 0.001591 |  0.0215 |  0.08475 |  0.2228 |

![probability of theta greater than 1.0](e3.png)

### Conditional autoregressive (CAR) model
Many statistical methods assume that observations are independent. However, data that occur at locations close together in space are likely to be correlated. Dependence between observations is a realistic assumption. Poisson-Gamma and Poisson-Lognormal model do not account for possible __unknown spatial structures__ (e.g. environmental effects). These are called __independent prior models__. Prior distributions should allow for spatial correlation (spatially structured priors).

#### Gaussian (autonormal) model
$$p(y_i |y_j , j \neq i) = N\left(\sum_j b_{ij} y_j , \tau_i^2\right)$$
where $b =1$ if neighbors and o otherwise, the mean is the sum average response in all neighbors. Using Brook's Lemma we can obtain
$$p(y_1 , y_2 , ... , y_n ) \propto \exp\left(-\dfrac{1}{2} \mathbf{y'} D^{-1} (I - B)\mathbf{y}\right)$$
where $B = \{b_{ij}\}$ and $D$ is diagonal with $D_{ii} = \tau^2_i$
- Suggests a multivariate normal distribution with $\mu_Y = 0$ and
$$\sum_Y = (I - B)^{-1}D$$
$D^{-1} (I - B)$ symmetric requires $\dfrac{b_{ij}}{\tau_i^2}  = \dfrac{b_{ji}}{tau_j^2}$ for all $i, j$
- Returning to $W$ , let $b_{ij} = w_{ij} /w_{i}+$ and $\tau_i^2 = \tau^2/w_{i+}$, so
$$p(y_1 , y_2 , ... , y_n ) \propto \exp\left(-\dfrac{1}{2\tau^2} \mathbf{y'} (D_w - W)\mathbf{y}\right)$$
where $D_w$ is diagonal with $(D_w )_{ii} = w_{i+}$ and thus
$$p(y_1 , y_2 , ... , y_n ) \propto \exp\left(-\dfrac{1}{2\tau^2} w_{ij}(y_i-y_j)^2\right)$$
This is the intrinsic autoregressive (__IAR__) model
- __Improper__ distribution since $(D_w - W)^{-1} = 0$, so requires a constraint, say $\sum_i y_i = 0$
- Not a data model, __a random effects model__

#### Proper CAR Model
Proper version: replace $D_w - W$ by $D_w - \rho W$ , and choose $\rho$ such that
$$\sum_y = (D_w - \rho W )^{-1}$$
exists. This in turn implies
$$Y_{i} |Y_{j\neq i} \sim N (\rho \sum_j w_{ij} Y_{j} , \tau^2 /m_i )$$

Or the mean is proportional to sum weighted average. Using $\rho$ makes distribution proper and adds parametric flexibility. When $\rho$ = 0 interpretable as independence.

#### Issues
Why should we expect $Y_i$ to be a proportion of average of neighbors - a sensible spatial interpretation? Calibration of $\rho$ as a correlation, e.g.

- $\rho$ = 0.80 yields 0.1 $\leq$ Moran’s I $\leq$ 0.15
- $\rho$ = 0.90 yields 0.2 $\leq$ Moran’s I $\leq$ 0.25
- $\rho$ = 0.99 yields Moran’s I $\leq$ 0.5

So, __used with random effects, scope of spatial pattern__ may be limited.

#### Bayesian Specification of CAR Model
Introduced by Clayton and Kaldor (1987) in an empirical Bayes setting. Developed by Besag et al. (1991) in a fully Bayes implementation. Also called the Besag, York and Mollie model. __Widely used__ in the area of disease mapping!

#### Convolution Model
$$y_i \sim Poisson(e_i \theta_i )$$
$$\log(\theta_i ) = \alpha + \beta x_i + u_i + v_i$$

- $\alpha$ is an overall level of the relative risk
- $v_i$ is the __uncorrelated__ heterogeneity
$$v_i \sim N (0, \sigma_v^2 )$$
- $u_i$ is the __correlated__ heterogeneity, a spatial correlation structure is used and estimation of the risk in any area depends on neighboring areas

The conditional autregressive model uses
$$[u_i |u_j , i \neq j, \tau_u^2 ] \sim N ( \bar{u_i} , \sigma_i^2 )$$
where

- $\bar{u_i} = \dfrac{1}{\sum_j w_{ij}} \sum_j u_j w_{ij}$
- $\sigma_i^2 = \dfrac{\sigma_u^2}{\sum_j w_{ij}} $
- $w_{ij} = 1$ is adjacent (and 0 otherwise) 

$u_i$ is __smoothed towards the mean risk in the set of neighboring areas__, with __variance inversely proportional to the number of neighbors__.

### Interpretation
Area-specific random effects are decomposed into clustering or correlated heterogeneity component and uncorrelated heterogeneity component. $\sigma_u^2$ and $\sigma_v^2$ control the variability of $u$ and $v$. 

The model allows the data to decide __how much of the residual disease risk is due to spatially structured variation__, and __how much is unstructured over-dispersion__, $\sigma_u^2 /\sigma_v^2$ small reflects an unstructured heterogeneity, $\sigma_u^2 /\sigma_v^2$ large indicates that a spatially structured variation dominates.

$\sigma_u^2$ controls the strength of spatial similarity induced by the CAR prior. Small values of $\sigma_u^2$ indicate stronger spatial correlation between neighboring regions. Setting $\sigma_u^2 = 0$ produces complete shrinkage to a common value.

    3.2.6 WinBUGS Code: Convolution model

`u[] ∼ car.normal(adj[], weight[], num[], tau.u)`

`adj[]` is a vector listing the adjacent areas for each area. The adjacency matrix can be generated using the Adjacency Tool from the Map menu. `weight[]` is a vector of the same length as `adj[]`, giving weights associated with each pair or areas; these weights can be generated using a loop in the WinBUGS code. `num[]` is a vector giving the number of neighbors for each area.

#### Running the MCMC
Two separate chains starting from different initial values. Convergence was checked by visual examination of time series plots of samples for each chain. Gelman and Rubin diagnostics were calculated. The first 20,000 samples were discarded as burn-in. Each chain was run for a further 10,000 iterations. __Estimates of the relative risk show less variation than the observed SMR__.

![Example 1: Cervix Cancer in Limburg, Belgium](e4.png)

![Example 1: Cervix Cancer in Limburg, Belgium](e5.png)

![Example 1: Cervix Cancer in Limburg, Belgium, with spacial structured](e6.png)

![Example 1: Cervix Cancer in Limburg, Belgium, heterogeneity](e7.png)

Table: Cervix Cancer in Limburg, Belgium, Posterior Parameter Estimates:

|   node  |    mean    |    sd   | MC error |   2.5%   |  median  |  97.5%  |
|---------|------------|---------|----------|----------|----------|---------|
| alpha   | -0.03222   | 0.06181 | 0.001151 |  -0.1577 | -0.03046 | 0.08741 |
| sigma.u | **0.1844** |  0.1465 | 0.008371 | 0.001734 |   0.1595 |   0.537 |
| Sigma.v | 0.01105    | 0.02497 | 0.001945 | 1.785E-4 | 0.001787 | 0.09169 |
| tau.u   | 50.8       |   182.0 |    14.65 |    1.862 |    6.268 |   584.2 |
| tau.v   | 1122.0     |  1529.0 |    97.37 |    10.94 |    559.8 |  5608.0 |

The variability of the relative risk is attributed less to the uncorrelated heterogeneity than to the spatially structured effects.

![Cervix Cancer in Limburg, Belgium, probability of RR greater than 1](e8.png)

#### Which model is the best one

Table: Model comparison with DIC:

| MODEL |   Dbar  |   Dhat  |   pD   |   DIC   |
|-------|---------|---------|--------|---------|
| PG    | 198.745 | 171.833 | 26.912 | 225.656 |
| PL    | 206.355 | 189.152 | 17.203 | 223.558 |
| CAR   | 208.222 | 194.186 | 14.036 | 222.258 |
| CONV  | 206.613 | 191.743 | 14.870 | 221.483 |

Convolution model is preferred

### Other Correlation Models
Other correlation models available in WinBugs: Multivariate normal model (similar to Bayesian kriging): uses distances between municipalities $d_{ij}$, e.g.
$$u \sim N_n (0, \tau \Sigma)$$
with $\Sigma_{ij} = cov(u_i , u_j ) = \exp (−\alpha d_{ij} )$

More flexible spatial form (2 parameters instead of 1). Stationary model. But requires inversion of large $n \times n$ covariance matrix. Proper CAR model with other way to specify correlation structured.

### Importing maps in WinBUGS

    3.2.8 code

## Ecological analysis
Example 2: South Africa Unemployment Data. Census Data in South Africa: 2001,2011. Number of households with access to piped water, % of 20+ year-olds without any schooling. Around 5% of population sampled.
>Is the number of households with access to piped water related to the education?

![South Africa Neighborhood Structure](e9.png)

![Proportion of households with access to piped water](e10.png)

![Proportion of 20+ year-olds without schooling](r1.png)

### Binomial Model
This can be modelled using a binomial likelihood:
$$y_i ∼ Binomial(n_i , p_i )$$
$$\logit(p_i ) = \alpha + \beta x_i + u_i + v_i$$

- $\alpha$ is an overall level of the relative risk, $\beta$ is the covariate effect, $v_i$ is the __uncorrelated__ heterogeneity with $v_i \sim N (0, \sigma_v^2 )$, $u_i$ is the __correlated__ heterogeneity a spatial correlation structure is used estimation of the risk in any area depends on neighboring areas

$$[u_i |u_j , i\neq j, \tau u_2 ] \sim N ( \bar{u}_i , \sigma_i^2 ) $$

    3.3.2 Code in WinBugs

#### Predicted Proportions

![Predicted Proportion Piped Waterd](r2.png)

![Standard Deviation of Prediction](r3.png)

#### CAR (CH) component

![Unstructured Heterogeneity Terms](r4.png)

#### UH component
![Spatially Structured Heterogeneity Terms](r5.png)

#### Is % no-schooling significant?

|  node |    mean    |    sd    |     LL    |   median  |     UL    |
|-------|------------|----------|-----------|-----------|-----------|
| alpha |  0.2169000 | 0.021360 |  0.172900 |  0.216900 |  0.260800 |
| beta  | -7.1140000 | 0.589400 | -8.305000 | -7.121000 | -5.889000 |

__An effect at the aggregate level cannot be interpreted on the individual level: ecological fallacy!__

#### Residuals

![Residuals](r6.png)
