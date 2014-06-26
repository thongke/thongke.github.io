# Point- Referenced Data

## Basics
Basic tool is a spatial process, $\{Y(\mathbf{s}), \mathbf{s} \in D\}$, where $D \subset Rcal^r$. Note that time series follows this approach with $r = 1$; we will
usually have $r$ = 2 or 3. We begin with essenials of point-level data modeling,
including stationarity, isotropy, and variograms - key ellements of the ‘Matheron school’. Then we add the spatial (typically Gaussian) process modeling that enables likelihood (and Bayesian) inference in these settings.

### Scallops Data
![](Selection_007.png)

![](Selection_008.png)

![](Selection_009.png)

![](Selection_010.png)

### Bivariate Interpolation Methods (EDA)
- Bivariate linear interpolation
- Bivariate spline interpolation
- Generalized additive models
- ...

```
3.1.R code
```
### Three standard statistic
Three ways of providing a description of how the data are related (correlated) with distance (and direction):

- __Covariogram__: Describes the way in which the __deviations of observations from their mean__ value __at different__ locations on the map __co-vary__. The covariance of a spatial stochastic process at any two locations (sites) $s$ and $s + h$ is defined as
$$COV [X(s), X(s + h)] = [X(s) - \overline{X(s)}] · [X(s) - \overline{X(s + h)}]$$
- Correlogram: 
$$COR[X(s), X(s + h)] = \dfrac{COV [X(s), X(s + h)]}{SD[X(s)]·SD[X(s+h)]}$$ 
where SD refers the standard devuations.
- Variogram (Semi-variogram): The variance is defined as a __half the expected squared difference__ between paired data values
$$VAR[X(s + h) - X(s)] = 2 · \gamma(h)$$
The $\gamma(d)$ function is known as variogram (or semi-variogram)

### Stationarity
Suppose our spatial process has a mean $\mu(\mathbf{s}) = E(Y (\mathbf{s}))$ and
that the variance of $Y(\mathbf{s})$ exists for adamll $\mathbf{s} \in D$. 

- The process is said to be __strictly stationary__ (also called __strong__ stationary) if, for any given $n \leq 1$, any set of $n$ sites $\{\mathbf{s}_1 , ... , \mathbf{s}_n \}$ and any $h \in IR r$ , the distribution of $\{Y (\mathbf{s}_1 ), ... , Y (\mathbf{s}_n )\}$ is the same as that of $\{Y (\mathbf{s}_1 + \mathbf{h}), ... , Y (\mathbf{s}_n + \mathbf{h})\}$. Thus, a process is strictly stationary if the __joint distribution__ of the process __only depends on the relative locations__, not on the actual locations themselves.
- A less restrictive condition is given by __weak stationarity__ (also called __second-order stationarity__). A process is weakly stationary if $\mu(s) \equiv \mu$ and $$Cov(Y(\mathbf{s}), Y(\mathbf{s} + \mathbf{h})) = C(\mathbf{h})$$ for all $\mathbf{h} \in IR r$ such that $s$ and $s + h$ both lie within $D$. Thus, a process is weakly stationary if the __expected value__ of the process __is constant__ and the __covariance__ of the process at any two sites __depends only on the relative locations__ of the two sites, not on the actual locations themselves.  Weak stationarity implies that the covariance relationship between the values of the process at any two locations can be summarized by a covariance function $C(\mathbf{h})$ (sometimes called a __covariogram__), and this function depends only on the separation vector $\mathbf{h}$. Note that with all variances assumed to exist, __strong stationarity implies weak stationarity. The converse is not true in general, but it holds for Gaussian processes__.

###Variograms
Suppose we assume $E[Y (\mathbf{s} + \mathbf{h}) - Y (\mathbf{s})] = 0$ and define
$$E[Y (\mathbf{s} + \mathbf{h}) - Y (\mathbf{s})]^2 = Var[Y (\mathbf{s} + \mathbf{h}) - Y (\mathbf{s})] = 2\gamma(\mathbf{h})$$
This expression makes sense when the __left-hand side depends only on__ $\mathbf{h}$, and __not the particular choice of__ $\mathbf{s}h$. If this is the case, we say the process is __intrinsically stationary__. The function $2\gamma(\mathbf{h})$ is then called the variogram, and $\gamma(\mathbf{h})$ is called the semivariogram.

Note that intrinsic stationarity __requires only the first and second moments__ of the differences $Y(\mathbf{s} + \mathbf{h}) - Y(\mathbf{s})$. __It says noting about the joint distribution__ of a collection of variables $Y (\mathbf{s}_1 ), ... , Y (\mathbf{s}_n )$, and __thus provides no likelihood__.

### Relationship between $C(h)$ and $\gamma(h)$
$2\gamma(\mathbf{h}) = Var[Y (\mathbf{s} + \mathbf{h}) - Y (\mathbf{s})]$

$= Var[Y (\mathbf{s} + \mathbf{h})] + Var[Y (\mathbf{s})] - 2Cov[Y (\mathbf{s} + \mathbf{h}), Y (\mathbf{s})]$

$= C(\mathbf{0}) + C(\mathbf{0}) - 2C(\mathbf{h})$

$= 2[C(\mathbf{0}) - C(\mathbf{h})]$

Thus,

$$\gamma(\mathbf{h}) = C(\mathbf{0}) - C(\mathbf{h})$$

So given $C$, we are able to determine $\gamma$, But in general, we can not always recover $C$ from $\gamma$.

In the relationship $\gamma(\mathbf{h}) = C(\mathbf{0}) - C(\mathbf{h})$ we can add a constant on the right side so $C(\mathbf{h})$ is not identified. If the spatial process is ergodic[^ergodic], then $C(\mathbf{h}) \rightarrow 0$ as $||\mathbf{h}|| \rightarrow \infty$, with $||\mathbf{h}||$ the length of the vector $\mathbf{h}$.

[^ergodic]: In mathematics, the term ergodic is used to describe a dynamical system which, broadly speaking, __has the same behavior averaged over time as averaged over the space of all the system's states__ (phase space). In physics the term is used to imply that a system satisfies the ergodic hypothesis of thermodynamics. In statistics, the term describes a __random process__ for __which the time average of one sequence of events is the same as the ensemble average__. In other words, for a Markov chain, as one increases the steps, there exists a positive probability measure at step  n  that is __independent__ of probability distribution at __initial step 0__ (Feller, 1971, p. 271).

Taking the limit of both sides of $\gamma(\mathbf{h}) = C(\mathbf{0}) - C(\mathbf{h})$ as $\mathbf{h} \rightarrow \infty$, we then have that $\lim_{||\mathbf{h}||\rightarrow\infty} \gamma(\mathbf{h}) = C(\mathbf{0})$. Thus, we have:
$$C(\mathbf{h}) = C(\mathbf{0}) - \gamma(\mathbf{h}) = \lim_{||\mathbf{h}||\rightarrow\infty} \gamma(\mathbf{h}) - \gamma(\mathbf{h})$$
So $C(\mathbf{h})$ is well defined if $\lim_{||\mathbf{h}||\rightarrow\infty} \gamma(\mathbf{h})$ __exists__. In this case, __intrinsic stationarity implies weak__ (second-order) stationarity. Thus, __weak stationarity implies intrinsic stationarity, but the converse is not true in GENERAL__.

### Isotropy
If the __semivariogram__ $\gamma(\mathbf{h})$ __depends__ upon the separation vector **only through its length** $\mathbf{h}$ , then we say that the process is isotropic; if not, we say it is anisotropic. For an isotropic process, $\gamma(\mathbf{h})$ is a real-valued function of a univariate argument, and can be written as $\gamma(||\mathbf{h}||)$.

If the process **is intrinsically stationary and isotropic, it is also called homogeneous**. Isotrophic processes are popular because of their simplicity, interpretability, and because a number of relatively simple parametric forms are available as candidates for $C$ (and $\gamma$).

## Some common isotropic covariograms
![Some common isotropic covariograms](Selection_011.png)
![Some common isotropic covariograms](Selection_012.png)

![](Selection_013.png)
![](Selection_014.png)
![](Selection_015.png)

### Spherical semivariogram
Selection_016.png

- While $\gamma(0) = 0$ by definition, $\gamma(0^{+}) = \lim_{t\rightarrow 0} \gamma(t) = \tau^2 $, this quantity is called the **nugget**.
- $\lim_{t\rightarrow \infty} \gamma(t) = \tau^2 + \sigma^2$ ; this asymptotic value of the semivariogram is called the sill.
- The sill minus the nugget, $\sigma^2$ , is called the partial sill 
- The value $t = 1/\phi$ at which $\gamma(t)$ reaches its ultimate level (the sill) is called the range, $R = 1/\phi$. The parameter $\phi$ is the decay parameter.

![Spherical semivariogram](Selection_018.png)

### Linear semivariogram
![](Selection_017.png)

Note that $\gamma(t) \rightarrow  \infty$ as $t \rightarrow  \infty$, and so the semivariogram does **not correspond to weakly stationary** process. The process in intrinsically stationary. The nugget is $\tau^2$ , but the sill and range ar both infinite.

### Exponential semivariogram
![](Selection_019.png)

The sill is only reached asymptotically, meaning that strictly speaking, the range $R = 1/\phi$ is infinite. To define the notion of an *effective range*, for $t > 0$,

$C(t) = \lim_{u\rightarrow \infty} \gamma(u) - \gamma(t)$

$= \tau^2 + \sigma^2 - \tau^2 + \sigma^2 (1 - \exp(-\phi t))$

$= \sigma^2 \exp(-\phi t)$

However, with $\gamma(h) = C(0) - C(h)$ we set $C(0) = \tau^2 + \sigma^2$.

![](Selection_020.png)

Then the **correlation between two points** distance $t$ apart is $\exp(-\phi t)$; note that $\exp(-\phi t) = 1$ for $t = 0^+$ and $\exp(-\phi t) = 0$ for $t = \infty$, as expected. We define the **effective range**, $t_0$ , as the distance **at which this correlation has dropped to only 0.05**. Setting $\exp(-\phi t_0 )$ equal to this value, we obtain $t_0 \approx 3/\phi$ (since $\log(0.05) \approx -3$). 

Finally, the form of $C(t)$ shows why the **nugget** $\tau^2$ is often viewed as a ‘**nonspatial** effect variance’ and the **partial sill** $\sigma^2$ as a ‘**spatial** effect variance’.

### Mat`ern Correlation function

![](Selection_021.png)

Versatile family. $K_\nu$ is bessel function of order $\nu$, with $\nu$ being a smoothness parameter:

- $\nu$ = 1/2: exponential
- $\nu$ = 3/2: results in convenient closed form for $C(t), \gamma(t)$
- $\nu \rightarrow \infty$ : Gaussian

### Variogram or Covariogram
Cressie (1991) argues in favor of variograms: **Variogram is defined in cases when covariogram is not**. Classical estimation of **variograms is more robust**.

In practice, modeling often done via **covariograms (more intuitive) even though weak stationarity is not usually verified**. Rough verification: Look at plot of $\widehat{C(h)}$ (method of moments estimate) and see if $\widehat{C(h)}$ appears to **go to 0 for large** $h$.

If process **is not weakly stationary**, try to **remove trends to produce weak** stationarity in the residual process. In general it is best to first try to satisfy the stationarity assumption (by removing trends, etc.) before trying to fit more complex (e.g. nonstationary) models.

### Empirical variogram
Empirical variogram is a method of moments estimate:

![](Selection_022.png)

where $N(h)$ is the number of pairs of samples that are at a distance of $h$ apart from each other. Usually need to ‘grid up’ the space in intervals $0 < h_1 < h_2 < ... < h_K$. The method of moments empirical variogram is **sensitive to outliers and squared differences may not be well behaved**. Robust estimator due to Hawkins and Cressie (1984):

![](Selection_023.png)

This estimator is **approximately unbiased**. Another alternative (replacing above mean with median): 

![](Selection_024.png)

![Spherical semivariogram](Selection_025.png)


### Variogram model fitting
$2\overline{\gamma(h)} , \widetilde{\gamma(h)}$ are **invalid** variograms, so spatial predictions using them may have negative variances. An ad-hoc solution (produces valid variogram) is plot empirical variogram and try different parametric forms or can **visually** estimate sill, nugget, range and thus obtain corresponding parameters of parametric variogram. More formal (**old**) approach: Use least squares, weighted least squares, generalized least squares etc. You can pursue full likelihood or Bayesian inference (have a model for the data).

After building an experimental variogram, we need to fit a theoretical function in order to model the spatial variation. The theoretical model to be selected should be that model which best fits the data. Some useful models: Gaussian, Exponential, Spherical models.

# Point-Referenced Data - Modeling
## Interpolation and Spatial Prediction
In exposure assessment, interest is in **predicting** the exposure at a location **where we have not recorded an observation**, say at location $s_o$. Interpolation is the process of obtaining a value for a variable of interest $(Y(s_o)$ at an unsampled location, **based on the surrounding** measurements (linear bivariate interpolation, inverse distance interpolation). 

When probabilistic models are used for interpolation, they are referred to as methods for spatial prediction (**kriging**). The **spatial predictions have standard errors that quantify the uncertainty associated with the interpolated values**.

## Inverse-Distance Interpolation
A deterministic approach, not a probabilistic model. A weighted average of neighboring values: 

![](Selection_026.png)

The weight given to each observation is a function of the distance $d_{0,i}$ between the observed locations $s_i$ and the unobserved location $s_0$. The weighting power $p$ is selected to control **how fast the weights tend to zero as the distance increases**. Distance powers between 1 and 3 are typically chosen. Taking $p = 2$ refers to the inverse-distance-squared interpolator.

![](Selection_027.png)
![](Selection_028.png)
![](Selection_029.png)

![](r8.png)

    5.Interpolating.R

Advantage: Simple conceptually, computationally fast, **exact** interpolator (interpolated surface passes through the original observations) 

Disadvantage: The mapped surfaces have flat-topped peaks and flat-bottomed valleys (concentric contours around data points). No underlying statistical model and no easy measure of the uncertainty associated with the value interpolated.

## Kriging
Geostatistical technique for optimal spatial prediction. Named by Matheron (1963) in honor of D.G. Krige, a South African mining engineer whose seminal work on empirical methods for geostatistical data inspired the general approach. 

Many different types of kriging, differing by underlying assumptions and analytical goals: Simple kriging uses linear prediction (i.e. predictor is linear combination of observed data values) assuming a **known mean**. Ordinary kriging uses linear prediction with a **constant unknown mean**. Universal kriging uses linear prediction with a **nonstationary mean** structure. Filtered kriging is kriging with measurement error. Lognormal kriging is optimal spatial prediction based on the lognormal distribution.

Trans-Gaussian kriging is spatial prediction based on transformations of the data. Co-kriging is multivariate spatial prediction. Indicator kriging is probabilitymapping under indicator functions (binary data). Probability kriging is probability mapping based on both indicator functions of the data and the original data. Disjunctive kriging is nonlinear prediction based on univariate functions of the data. Bayesian kriging is incorporates prior information about the mean and/or covariance functions into spatial prediction. Block kriging is optimal linear prediction of areal data from point data,... 

### Ordinary Kriging
Assume that we have data $\mathbf{Y} = (Y (s_1 ), ... , Y (s_n ))$ and want to predict the value of the process $Y (.)$ at an unobserved location, $Z(s_0 ) (s_0 \in D)$. Assume $Y (.)$ **is intrinsic stationary** (i.e. **constant unknown mean** $\mu$ and known semivariogram $\gamma(h))$. A linear predictor is takes the form
$$Y (s_0 ) = \lambda_i Y (s_i )$$
**Instead** of specifying the weights as an **arbitrary function** of distance (inverse-distance interpolation), determine the weights **based on the data using the semivariogram**. Determine weights based on two criteria: unbiasedness and
minimum mean-squared prediction error. The best linear prediction is found using two criteria:

- Unbiasedness requires $E( Y (s_0 )) = \mu = E(Y (s_0 ))$, which requires $\sum_{i=1}^n \lambda_i = 1$ 
- Minimize the mean-squared prediction error (MSPE) defined as 
$$E[( Y (s_0 ) - Y (s_0 ))^2 ]$$
This results in solving a constrained optimization problem, e.g. using Lagrange multipliers. This gives a system of equations, referred to as the kriging equations: 
$$\sum_{j=1}^N \lambda_j \gamma(s_i - s_j ) + m = \gamma(s_0 - s_i )$$
$$\sum_{i=1}^N \lambda_i = 1$$

Note that no distributional assumptions are required for the $Y (s_i )$. The kriging weights are given by:

![](Selection_030.png)

Note that we need to calculate these kriging weights $(\lambda_1 , ... , \lambda_n )$ for each prediction location $s_0$. The matrix $\Gamma_0$ **depends only on the data** locations, and **not on the prediction locations**, so we need only to invert $\Gamma_0$ once, and then multiply by the associated $\gamma_0$ vector to obtain a prediction for any $s_0 \in D$,

Difficulties: **Limitation of constant mean, Mean surface unknown, Variogram unknown**. If we put estimates of both into the kriging equations we **fail to  take into account the uncertainty in these estimates**. Therefor, we turn to Gaussian processes and likelihood based methods, to work with the covariance function.

    5.kriging.R

![Simple kriging](r9.png)
![Simple ordinary kriging](r91.png)
![Simple universal kriging](e92.png)

### Kriging with Gaussian processes
Consider the case where we have **no covariates**, but only responses $Y (s_i)$ (*ordinary kriging with Gaussian process*)
$$\mathbf{Y} = \mathbf{\mu} + \mathbf{\varepsilon} \text{ with } \mathbf{\varepsilon} \sim \mathbf{N (0, \Sigma)}$$
Consider the case where we have covariate values $\mathbf{x}(s_i )$
(*universal kriging with Gaussian process*)
$$\mathbf{Y} = \mathbf{X}\beta + \mathbf{\varepsilon} \text{ with } \mathbf{\varepsilon} \sim \mathbf{N (0,}\Sigma)$$
For a spatial covariance structure without nugget effect, we specify
$$\Sigma = \sigma^2 H(\phi)\text{ with }(H(\phi))_{ij} = \rho(\phi; d_{ij} )$$
with $\rho$ a **valid correlation function** depending on the distance between $s_i$ and $s_j$ (covariance functions). For a model having a nugget effect, we instead set
$$\Sigma = \sigma^2 H(\phi) + \tau^2 I$$

We seek the function $f (\mathbf{y})$ that minimizes the mean-squared prediction error $E[Y (s_0 ) - f (\mathbf{y}))^2 \mathbf{y}]$. This expression equals
$$E [(Y (s_0 ) - E[Y (s_0 )|\mathbf{y}])^2 |\mathbf{y}] + [E(Y (s_0 )|\mathbf{y}) - f (\mathbf{y})]^2$$
since the expectation of the cross-product equals 0.
Since the second term is nonnegative, we have
$$E[Y (s_0 ) - f (\mathbf{y}))^2 |\mathbf{y}] \geq E (Y (s_0 ) - E[Y (s_0 )|\mathbf{y}])^2 |\mathbf{y}$$
for any function $f (\mathbf{y})$

Equality holds if and only if $f (\mathbf{y}) = E(Y (s_0 )|\mathbf{y})$. **Thus, the conditional expectation** of $Y (s_0 )$ given the data (i.e. posterior mean of $Y (s_0 )$) is the predictor that minimizes the error.

We now turn to the estimation of this predictor. Assume (the unrealistic situation) that the parameters ($\beta, \sigma, \phi, \tau^2 $) are known. Suppose

![](Selection_031.png)

with $\Omega_{21} = \Omega^T_{12}$. Then $p(Y_1 |Y_2 )$ is normal with mean and variance
$$E(\mathbf{Y}_1 |\mathbf{Y}_2 ) = \mu_1 + \Omega_{12} \Omega^{-1}_{22} (\mathbf{Y}_2 - \mu_2 )$$
$$Var(\mathbf{Y}_1 |\mathbf{Y}_2 ) = \Omega_{11} - \Omega_{12} \Omega^{-1}_{22} \Omega_{21}$$

In our setting, $Y_1 = Y (s 0 )$ and $Y_2 = \mathbf{y}$ (all observed data), meaning that

$\Omega_{11} = \sigma^2 + \tau^2 ;$

$\Omega_{12} = \gamma^T ;$

$\Omega_{22} = \Sigma = \sigma^2 H(\phi) + \tau^2 I$

with $\gamma^T = \sigma^2 \rho(\phi; d_{01} ), ... , \sigma^2 \rho(\phi; d_{0n} )$. This results in
$$E(Y (s 0 )|\mathbf{y}) = X_0 \beta + \gamma^T \Sigma^{-1} (y - X\beta)$$
$$V ar(Y (s 0 )|\mathbf{y}) = \sigma^2 + \tau^2 - \gamma^T \Sigma^{-1} \gamma$$
The weights $\gamma^T \Sigma^{-1}$ are **known as the kriging weights**. Note that this solution **assumes we have observed the covariate value** $X_0 = X(s_0 )$ at the ‘new’ location $s_0$.

Next, consider the **more realistic scenario where the model parameters are unknown**, then 
$$\widehat{f (\mathbf{y})} = \mathbf{X}_0 \hat{\beta} + \hat{\gamma}^T \Sigma^{-1} (\mathbf{y} - \mathbf{X} \hat{\beta})$$
with $\hat{\gamma}^T = (\hat{\sigma}^2 \rho( \hat{\phi}; d_{01} ), ... , \hat{\sigma}^2 \rho( \hat{\phi}; d_{0n} ))$, $\hat{\Sigma} = \hat{\sigma}^2 H( \hat{\phi})$ and $\hat{\beta} = (X^T \hat{\Sigma}^{-1} X)^{-1} X^T \hat{\Sigma}^{-1} \mathbf{y}$ the usual weighted least squares estimator of $\beta$.

The result $\widehat{f (\mathbf{y})}$ can be written as $\mathbf{\lambda^T} \mathbf{y}$ where 
$$\lambda = \hat{\Sigma}^{-1} \hat{\gamma} + \hat{\Sigma}^{-1} X (X^T \hat{\Sigma}^{-1}X)^{-1} (X_0 - X^T \hat{\Sigma}^{-1}\hat{\gamma})$$

#### Remarks

- if the number of predictors $p > 0$, the above best linear unbiased prediction method is referred to as **universal kriging**
- if $p = 1$ (and X does not include coordinates) the term kriging with external drift is used
- if $p = 0$ and $X_0 \equiv 1$ this is called ordinary kriging
- if $\beta$ is **assumed a priori known** this results corresponds to simple kriging.
- If $x_0$ is **unobserved**, we can estimate $x_0$ and $Y (s_0 )$ jointly by iterating between the previous formula and a corresponding one for $x_0$ , namely $x_0 = X^T \lambda$, which arises simply by multiplying both sides of the previous equation by $X^T$ and simplifying. 
    - This is essentially an EM (expectation-maximization) algorithm, with the calculations of $x_0$ being the $E$ step and the updating of $\lambda$ being the $M$ step. 
    - In the classical framework, restricted maximum likelihood (REML) estimates are often plugged in above, and shown to have certain optimal properties.\

```
Maximum.likelihood.R
```

### Bayesian Kriging
The likelihood is given by
$$Y |\theta \sim N (X\beta, \sigma^2 H(\phi) + \tau^2 I)$$
Typically, independent priors are chosen for the parameters
$$p(\theta) = p(\beta)p(\sigma^2 )p(\tau^2 )p(\phi)$$
Useful candidates are multivariate normal for $\beta$ and inverse gamma for $\sigma^2$ and $\tau^2$. Specification of $p(\phi)$ depends upon choice of $\rho$ function; a uniform or gamma prior is usually selected. Informativeness: $p(\beta)$ can be ‘flat’ (improper), but $\phi$ and at least **one of** $\sigma^2$ and $\tau^2$ **require informative** priors.

We can recast the foregoing in a hierarchical setup by considering conditional likelihood on the spatial random effects $w = (w(s_1 ), ... , w(s_n ))$. 

- First stage:
$$Y |\theta, w \sim N (X\beta + w, \tau^2 I)$$
The $Y (s_i )$ are conditionally independent given the random effects $w$
- Second stage:
$$w|\sigma^2 , \phi \sim N (0, \sigma^2 H(\phi))$$
- Third stage: priors on $(\beta, \sigma^2 , \tau^2 , \phi)$
- Often we need to predict the response $Y$ at a new site $s_0$ with associated covariates $x_0 ≡ x(s_0 )$. Predictive distribution:
$$p(y_0 |\mathbf{y}, X, x_0 ) = \int p(y_0 , \theta|\mathbf{y}, X, x_0 )d\theta$$
$$= p(y_0 |\mathbf{y}, X, x_0 , \theta)p(\theta|\mathbf{y}, X)d\theta$$
$p(y_0 |\mathbf{y}, X, x_0 , \theta)$ is normal since $p(y_0 , y|\mathbf{y}, X, x_0 , \theta)$ is. Easy Monte Carlo estimate using composition with Gibbs draws of $\theta$.
- Suppose we want to predict at a set of $m$ new sites, say $S_0 = \{s_{01} , ... , s_{0m} \}$. We could individually predict each site ‘independently’ using method of the previous slide. But, joint prediction is also of interest, e.g. bivariate predictive distributions to reveal pairwise dependence, to reflect posterior associations in the realized surface. Form the unobserved vector $\mathbf{Y}_0 = (Y (s_{01} ), ... , Y (s_{0m} ))$, set $X_0$ as covariate matrix at $S_0$ and compute 
$$p(\mathbf{y}_0 |\mathbf{y}, X, X_0 ) = p(\mathbf{y}_0 , \theta|\mathbf{y}, X, X_0 , \theta)p(\theta|y, X)d\theta$$. Estimation can be performed again using **methods of composition** sampling.
```
5.bayes.R
```

### Spatial GLM
Some data sets preclude Gaussian modeling: $Y (s)$ might not even be continuous Example $Y (s)$ is a binary or count variable such as precipitation or deposition was measurable or not, number of insurance claims by residents of a single family home at s price is high or low for home at location $s$. **Replace Gaussian likelihood by an appropriate exponential family member** (Diggle, Tawn and Moyeed (1988))

- First stage: $Y (s_i )$ are conditionally independent given $\beta$ and $w(s_i)$. The likelihood $f (y(s_i )|\beta, w(s_i ), \gamma)$ is such that
$$g(E(Y (s_i ))) = \eta(s_i ) = X^T (s_i )\beta + w(s_i )$$
where $\eta$ is a canonical link function (such as a logit) and $\gamma$ is a dispersion parameter
- Second stage: model $w(s)$ as a Gaussian process: 
$$w|\sigma^2 , \phi \sim N (\mathbf{0}, \sigma^2 H(\phi))$$
- Third stage: priors and hyperpriors

# Point Pattern Data
For a **specified, bounded region** $D$, a set of locations $s_i , i = 1, 2, ... , n$. The **locations are random** (the event locations), **only locations, no variables at locations**. Primary goal is detect patterns, e.g. complete randomness, i.e. no common causative factor between event locations, clustering, i.e. events occurring more frequently in close proximity to one another and regular/systematic pattern. Marks or marked patterns is used to comparison of patterns.

Examples: pattern of trees in a forest, pattern disease cases (possible cases and controls), breast cancer cases (treatment option), bovine tuberculosis (also over time).

## Patterns
**Complete Spatial Randomness (CSR)**
: Event is equally likely to occur at any location within the study area, regardless of the locations of other events. Events follow a uniform distribution across the study areas and independent of one another. CSR is boundary condition between spatial processes that are more clustered than random and that are more regular than random.

**Clustering**
: Reflecting areas with increases occurrences of events, more cases then under CSR (attraction, contagion)

**Regular**
: Distances between events is **larger than under CSR** and occurs when there is inhibition/competition among points.

![Six realizations (data sets) based on complete spatial randomness, with 30 event locations distributed indepen- dently (no interactions between events) and uniformly (events equally likely at any locations) across the unit square.](Selection_034.png)

![Three examples each of realizations (data sets) arising from a spatial point process more clustered that complete spatial randomness (top row), and a spatial point process more regular than complete spatial randomness (bottom row). Each data set contains 30 event locations within the unit square.](Selection_033.png)

![Realizations of two hypothetical spatial point processes containing both clustering and regularity at different spatial scales.](Selection_035.png)

The realizations from 3 CSR model: contains collections of nearby events (apparent clusters), contains locations with **large gaps between areas**.

**Clustering** occurs by chance, making visual assessment of overall pattern difficult. Observed patterns do not always fall neatly into one of the three classes. The level of spatial scale is important to describe the observed pattern.

## F and G functions
To measure the degree of accomplishment of the CSR, several functions can be computed on the data. These measures can be used to **test for CSR**. Distance based methods:

- $G$-function measures the distribution of the **distances from an arbitrary event to its nearest event**. $G(d)$ is the nearest neighbor distance, **event to event**, i.e. 
$$G(d) = Pr(\text{nearest event} \leq d)$$
- $F$-function measures the distribution of all **distances from an arbitrary point of the plane to its nearest event**. $F(d)$ is the nearest neighbor distance, **point to event**, i.e. 
$$F (d) = Pr(\text{nearest event} \leq d)$$
- The $G$-function can be estimated as 
$$\hat{G} (d) = \dfrac{\text{#} \{d_i : d_i \leq d, \forall i\}}{n} $$ 
where $d_i = \min_j \{d_{ij} , \forall j \neq i\}$ and $n$ is the total number of points. $\hat{G}$ is the empirical cdf of the $n$ nearest neighbor distances (nearest neighbor distance for $s_1$ , for $s_2$ , etc.). **Under CSR**: 
$$G(d) = 1 - \exp(-\lambda \pi d^2 )$$
Compatibility with CSR of the point pattern can be assessed by plotting the empirical function $\hat{G}$ against the theoretical expectation.
- The $F$-function can be estimated as
$$\hat{F} (d) = \dfrac{\text{#} \{d_i : d_i \leq d, \forall i\}}{m} $$
where $d_i = \min_j \{d_{ij} , \forall j \neq i\}$ are the minimal distance from
$m$ **random locations****** $s_j$ to the event locations $s_i$. $F$ is the empirical cdf arising from the $m$ nearest neighbor distances associated with a **randomly selected set** of $m$ points in $D$. The $F$-function is often called the **empty space function, because it measures the average space left between events**. **Under CSR**: 
$$F (d) = 1 - \exp(-\lambda \pi d^2 )$$
Hence, compare the estimated value of $F$-function to its theoretical value.
- Edge correction if $d > b_i$ , distance form $s_i$ to edge of $D$.
- Compare $\hat{G}$ with $G$, if $\hat{G} > G$ clustered pattern, if $\hat{G} < G$ more **regular** pattern. Compare $\hat{F}$ with $F$ if $\hat{F} > F$ more regular pattern, if $\hat{F} < F$ **clustered** pattern.

Monte Carlo (simulation-based) methods of inference to assess significance. No need to rely on asymptotic properties. Procedure:

- Calculate the test statistic value based on the data observed
- Calculate the same statistic for a large number (N~sim~ ) of data sets simulated independently under the null hypothesis
- A histogram of the statistic values associated with simulated data sets provides description of the null-distribution. The estimated p-value (one-sided test) is
$$P (T \geq T_{obs} |H_0 ) = \dfrac{l}{N_{sim} + 1} $$
with $T_{(1)} \geq T_{(2)} \geq ... \geq T_{(l)} ≥ T_{(obs)} ≥ T_{(l+1})≥...≥T_{(N_{sim}}$ is the ordering test statistics. For functions, point-wise envelopes under CSR can be obtained in a similar way.

![Example of three point patterns re-scaled to ﬁt in the unit square. On the left, spatial distribution of the location of cell centres (Ripley, 1977); in the middle, Japanese black pine saplings (Numata, 1961); and on the right, saplings of California redwood trees (Strauss, 1975)](Selection_036.png)

![Envelopes and observed values of the G function for three point patterns](Selection_037.png)

![Envelopes and observed values of the F function for three point patterns](Selection_038.png)

## Spatial Point Processes
A spatial point process describes a stochastic process where **each random variable represents the location of an event in space**. A realization of the process is a collection of locations generated under the spatial point process model and represents a data set resulting from a particular model (either observed or simulated). 

Two important underlying concepts for modeling spatial point patterns are *stationarity* when the process is **invariant to translation** within the space, *isotropic*  when the process is **invariant to rotation** about the origin. *Stationarity + isotropy* events **depend only on the distance** separating their locations an **not on their orientation** to each other.  These properties offer a notion of replication within the data.

### Homogeneous spatial Poisson point process
Defined by the following criteria:

1. The number of events occurring within a **finite region** $A$ is a random variable **following a Poisson** distribution with mean $\lambda |A|$ for some positive constant $\lambda$  and $|A|$ denoting the area of $A$.
2. Given the total number of events $N$ occurring within an area $A$, the location of the $N$ events represent an **independent** random sample of $N$ locations, where each point is **equally likely to be chosen as an event**

Criterion 1 introduces idea of an intensity $\lambda$, representing the **number of events expected per unit area**. Criterion 2 represents general concept of CSR (events uniformly distributed across the study area).

Estimation of $\lambda$: divide total number of events observed by the total area: 
$$\hat{\lambda}  = \dfrac{N}{|A|}$$
Since the **intensity of events is constant at all locations**, the process is **called homogeneous**. Note that a homogeneous spatial Poisson point process is **stationary and isotropic**.

The definition of homogeneous spatial Poisson point process provides a two-stage approach for generating realization from CSR in study area $D$:

- Generate the total number of points $N (D)$ from a $Poi(\lambda |D|)$ 
- Place events within $D$ according to a uniform distribution
- If $D$ rectangular: generate $u$ and $v$ coordinates using uniform random number generators on the intervals corresponding on the width and height of $D$ 
- If $D$ nonrectangular: embed $D$ within a larger rectable $R$, and generate event locations within $R$ until $N (D)$ events occur within $D$.

### Heterogeneous spatial Poisson point process
The poisson process is homogeneous when the intensity $\lambda$ is constant across the study area. CSR may not be an appropriate model for the lack of clustering, since a the population at risk is not necessarily uniformly across space (people tend to live in towns and cities), or environmental factors such as humidity, quality of soil, affect the chance of an event. As an alternative to assumption of CSR, consider the **hypothesis of constant risk** as a model of ‘no clustering’. Under hypothesis of constant risk **each person has the same risk** of disease during the observation period, **regardless of location**, we expect more cases in areas with more people at risk.

Clusters of cases in high population areas could Violate CSR but not the constant risk hypothesis.

![Example of a process that appears clustered with respect to CSR but not clustered with respect to a hypothesis of constant risk. The box represents an area of high population density. The set of event locations is the same in both plots. In the plot on the left, we observe a higher intensity of events in the high-population area, consistent with a constant risk hypothesis but inconsistent with CSR. In the plot on the right, the cluster occurs outside the area of high population density, reﬂecting an area of increased local risk and a clustered process with respect to both CSR and constant risk](Selection_039.png)

Heterogeneous spatial Poisson point process defined by the following criteria:

1. The number of events occurring within a finite region $A$ is a random variable following a Poisson distribution with mean $\int_A \lambda (s)ds$ 
2. Given the total number of events $N$ occurring within an area $A$, the **location** of the $N$ events **represent** an **independent random sample** of $N$ locations, with the **probability of sampling** a particular point $s$  **proportional** to $\lambda (s)$

The constant risk hypothesis requires a spatially varying intensity function $\lambda (s)$ and the events are distributed according to a spatial density function proportional to $\lambda (s)$, Heterogeneity implies **nonstationarity** (no longer translation invariant). Heterogeneity in the intensity function $\lambda (s)$ results in an anisotropic process only if $\lambda (s)$ is anisotropic (i.e. not symmetric around the origin).

The intensity function $\lambda (s)$ is a first-order (mean) property, describing the **expected density of events in any location of the region**. Under a heterogeneous spatial poisson point process, **clusters occur solely due to heterogeneities in the intensity function** and **individual event locations remain independent** of one another.

![Example intensity function, \lambda (s), for a heterogeneous Poisson point process defined for s = (u, v) and u, v \in (0, 20).](Selection_040.png) 

![Six simulated realizations of 100 events each from an inhomogeneous Poisson process with intensity shown in Figure 5.5.](Selection_041.png)

Crude approach: imagine a refined grid over $D$. Then
$$\lambda (\partial s) = \partial s \lambda (s)ds ≈ \lambda (s)|\partial s|$$
So, for grid cell $A_l$ , assume $\lambda$  is constant over $A_l$ and estimate with $N (A_l )/|A_l |$. Two dimensional step function, tile function: like a two-dimensional histogram.

More sophisticated: a kernel intensity estimate (like a kernel density estimate). This is a non-parametric estimator of the intensity. For every $s \in D$: 
$$\lambda_\tau (s) = \dfrac{1}{h^2} \sum_i\kappa \left(\dfrac{||s - s_i ||}{h} \right)/q(||s||)$$
$\kappa$ is a bivariate and symmetrical kernel function (usually a bivariate normal pdf). $h$ is the bandwidth (controls smoothness of $\lambda$ ). $q(||s||)$ is a border correction to compensate for the missing observations that occur when $s$ is close to the border of the region.

Example kernel smoothing in one dimension, and effect of bandwidth:

![Kernel intensity estimates based on the same set of 20 event locations and four different kernel bandwidths [kernel variances reﬂect the square of the bandwidth (standard deviation)]](Selection_042.png) 

Different kernel functions can be used (overview given by Silverman, 1986) Example: Product kernel based on two univariate Gaussian kernels The two-dimensional quartic kernel (also known as biweight) 
$$\kappa(u) = \dfrac{3}{\pi} (1 - ||u||^2)^2 $$
if $u \in (-1, 1)$, and $\kappa(u) = 0$ otherwise. **A possible approach to select optimal bandwidth involves e.g. the mean squared error** of the kernel smoothing estimator.

![Early medieval grave site locations. Circled locations denote those grave sites where the person shows evidence of a particular tooth defect (“affected individuals“). The polygon surrounding the points represents the edge of the study area Arrows indicate two locations, each containing two grave sites occurring at such small distances that visual distinction of event locations is difﬁcult at the scale of the ﬁgure. The question of interest: Do the burial sites of affected individuals tend to cluster?](Selection_043.png) 

![Kernel smoothed density estimate (proportional to the intensity function) for affected loca- tions (grave sites with tooth defect) in an early medieval grave site data set. Bandwidth set to 872.24 in the u direction and 997.45 in the 1; direction based on Scott’s rule [equation (5.3)] for Gaussian kernels (see text).](Selection_044.png)

![Kernel smoothed spatial density estimate (proportional to the intensity function) for non- affected locations (grave sites without tooth defect) in early medieval gave site data set. Bandwidth set to 695.35 in the u direction and 734.82 in the 1; direction based on Scott’s rule [equation (5.3)] for Gaussian kernels (see text).](Selection_045.png) 

Alternatively, a parametric or semi-parametric form for the intensity may be of interest (e.g. to include covariates). The log-likelihood of a realisation of $n$ independent events with intensity $\lambda (s)$ is
$$L(\lambda ) = \sum_{i=1}^n \log \lambda (s_i ) -\int_A \lambda (s)ds$$
where $\int_A \lambda (s)ds$ is the expected number of cases in region $A$. Diggle (2003) suggests a log-linear model 
$$\log \lambda (s) = \sum_{j=1}^p \beta_j z_j (s)$$
using covariates $z_j (s)$ measured at location $s$. Standard statistical methods, such as MLE, can be used to estimate the parameters.

![Location of maple trees from the Lansing data set and their estimated parametric intensity using model (72)](Selection_046.png)

Log-intensity specified by $log(\lambda (s)) =  \alpha + \beta_1 x + \beta_2 x_2 + \beta_3  x^2_1 + \beta_4 x^2_2 + \beta_5 x_1 x_2$

Generating a realization under the constant risk hypothesis can be done by generating from a heterogeneous spatial Poisson point process given a known intensity function $\lambda (s)$: 

- Compute $\lambda_{max} = \max_{s\in D} \lambda (s)$
- Sample $n$ from $Poi(\lambda_{max} |D|)$
- Given $n$, sample $n$ locations uniformly over $D$
- ‘Thin’ by retaining $s_i$ with probability $\lambda (s_i )/\lambda_{max}$

In summary: Let $N (A)$ be the number of points in area $A$. Distribution of $N (A)$ is driven by the intensity surface $\lambda (s)$, based on Poisson process
$$N (A) \sim Po(\lambda (A)) \text{ where } \lambda (A) = \int_A \lambda(s) ds$$
Special processes:

- $\lambda (s) = \lambda$ : homogeneous Poisson process, spatial homogeneity, complete spatial randomness (CSR), $\lambda (A) = \lambda |A|$ 
- $\lambda (s)$ nonconstant, fixed: nonhomogeneous (heterogeneous) Poisson process 
- $\lambda (s)$ random: Cox process

K functions

- We may be interested in how often events occur within a given distance to other events. 
- Intensity function = first-order property of a point process: informs on the mean 
- Relative position of events = second-order property: informs on the interrelationship between events. Second-order properties **measure the strength and type of the interactions between events of the point process**.

They are interesting when **studying clustering or competition between** events
- The $K$ function considers the number of points within distance $d$ of an arbitrary point. Formally,
$$K(d) = \dfrac{E(\text{# of events within d of an arbitraty event})}{\lambda}$$
Also called **reduced second moment measure**. To compute this function, Ripley (1976) proposed the estimate
$$\hat{K}(d) = (n(n - 1))^{-1} |A| \sum_{i=1}^n \sum_{j\neq i} w_{ij}^{-1} |\{s_j : d(s_i , s_j ) \leq d\}|$$
with $w_{ij}$ weights equal to **the proportion of the area inside the region** $A$ **of the circle centered at** $s_i$ and radius $d(s_i , s_j )$. Under CSR, the value of $K(d$) is $\pi d^2$ (the area of a circle of radius $d$). For **processes more regular, we would expect fewer events within distance** $d$ of a randomly chosen event than under CSR. For processes **more clustered, we would expect more events** within a given distance than under CSR.

Under CSR, $K(d) = \lambda \pi d^2 /\lambda  = \pi d^2$
$$\hat{K}(d) = (n \hat{\lambda} )^{-1} \sum_i \sum_j (w_{ij} )^{-1} I(||s_i - s_j || \leq d) $$
with $\hat{\lambda}  = n/|D|$. $w_{ij}$ is an edge correction, the proportion of the circumference of the circle centered at $s_i$ with radius $||s_i - s_j ||$ within $D$.

Compare $\hat{K}(d)$ with $K(d) = \pi d^2$, regularity implies $K(d) < \pi d^2$, clustering implies $K(d) > \pi d^2$. Alternatively, compare $\hat{L}(d) = \left(\dfrac{\hat{K}(d)}{\pi} \right)^{1/2}$  with $d$.

![Envelopes and observed values of Ripley’s K-function for three point patterns](Selection_047.png)

The $K$-function **does not suggest where clusters occur**, but rather, at what distances events tend to occur from other events with respect to distances expected under CSR. If event location occur within an irregularly shaped, nonconvex polygon, CSR over a rectangle does assign event locations uniformly throughout the rectangle (including locations outside polygon). In this case the K function might **entirely describe the occurrence of events within the enclosing polygon**. The above methodology always makes the comparison with CSR, but it can be adapted to compare directly patterns (e.g. affected versus nonaffected sites) by using random labeling hypothesis, random labeling the set of locations of all events observed. This offers one approach for assessing the similarity in the underlying processes. 

![early medieval grave site locations. Circled locations denote those grave sites where the person shows evidence of a particular tooth defect (“affected individuals“). The polygon surrounding the points represents the edge of the study area Arrows indicate two locations, each containing two grave sites occurring at such small distances that visual distinction of event locations is difﬁcult at the scale of the ﬁgure. The question of interest: Do the burial sites of affected individuals tend to cluster?](Selection_048.png)


![](Selection_049.png)

![](Selection_050.png)

![L plot for the early medieval grave site data shown in Figure 5.8, compared to a random labeling hypothesis (see text). The solid line illustrates L(h) - h for the case data (comparable to the middle plot of Figure 5.13). The dashed, dotted, and dash-dotted “envelopes” represent the minimum, maximum, 2.5th percentile, 97.5th percentile, and median values, respectively, based on 499 random samples of 30 sites each from the set of 143 sites.](Selection_051.png)

### First and Second-order properties
- First order (intensity function) and second order (K function) analysis provide different but complementary insight into the analysis of spatial point patterns. 
- Second-order analysis gives insight into global aspects: Are there general patterns of clustering and/or regularity with respect to CSR or another pattern
- First-order analysis provide local insights Where do the patterns appear to differ? 

### Other Point Processes
**Poisson cluster processes**: defines spatial point processes wherein each event belongs to a particular (latent) cluster.

**Contagion/inhibition processes**: focuses on direct modeling of interevent interactions wherein the occurrence of an event raises or lowers the probability of observing subsequent events nearby (e.g. for modeling of infectious diseases or territory of an animal) 

**Cox processes**: considers the intensity function as random quantity drawn from some probability distribution (e.g. heterogeneity changes from year to year)

Distinguishing between possible underlying processes based on observed data can be problematic. With replicate realizations, one may be able to distinguish the patterns.