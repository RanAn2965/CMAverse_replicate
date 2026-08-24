# Code Review: `cmsens.R`

This guide provides an in-depth code review, architectural blueprint, and mathematical breakdown of the **`cmsens`** function in the `CMAverse` package. `cmsens` serves as the primary gateway for assessing the robustness of causal mediation analyses against two major threats to validity: **unmeasured confounding** and **measurement error (misclassification)**.

---

## 1. Executive Summary & Design Philosophy

The core design philosophy of `cmsens.R` is **post-estimation robustness evaluation**. Rather than integrating sensitivity analysis directly into the primary model-fitting stage, `cmsens()` acts as an analytical wrapper around a completed class `cmest` object.

```
                   ┌────────────────────────────────────────┐
                   │  Primary Mediation Object (class cmest)│
                   └───────────────────┬────────────────────┘
                                       │
                                       ▼
                       ┌──────────────────────────────┐
                       │      cmsens(sens = ...)      │
                       └──────┬────────────────┬──────┘
                              │                │
            sens = "uc"       ▼                ▼       sens = "me"
     ┌─────────────────────────────┐    ┌─────────────────────────────┐
     │   Unmeasured Confounding    │    │      Measurement Error      │
     ├─────────────────────────────┤    ├─────────────────────────────┤
     │  - Standardizes variables   │    │  - Error models specified   │
     │  - Transforms difference -> │    │  - Wraps fits in rcreg()    │
     │    Risk Ratio scale         │    │    or simexreg()            │
     │  - Calculates E-Values via  │    │  - Re-runs estinf() loop    │
     │    EValue::evalues.RR       │    │  - Out: Corrected PE / SEs  │
     └─────────────────────────────┘    └─────────────────────────────┘
```

By decoupling estimation from sensitivity checks, `cmsens()` allows researchers to run multiple sensitivity scenarios (e.g., varying measurement error standard deviations or testing different confounding limits) without having to re-specify their underlying baseline mediation models.

---

## 2. Core Architectural Scaffolding

### A. Parameter Inheritance and Extraction
The function begins by extracting the active data, variables, regressions, point estimates (`effect.pe`), and standard errors (`effect.se`) directly from the parsed `cmest` object:
```r
data <- object$data
model <- object$methods$model
reg.input <- object$reg.input
effect.pe <- object$effect.pe
effect.se <- object$effect.se
outcome <- object$variables$outcome
```
This ensures perfect alignment between the sensitivity checks and the original study's specifications (e.g., exposure-mediator interactions, covariates, and sample weights).

### B. Environment Sharing via Lexical Scoping
Just like `cmest.R`, `cmsens.R` relies on R's dynamic lexical scoping overrides to run its internal pipelines without passing complex lists of state variables:
```r
environment(regrun) <- environment()
regs <- regrun()
...
environment(estinf) <- environment()
estinf_res <- estinf()
```
By assigning the local environment of `cmsens()` to the helper functions `regrun` and `estinf` (housed in external package files), these helpers gain direct read/write access to all locally assigned variables (e.g., modified outcome and mediator regressions, measurement error structures, and bootstrapping parameters).

---

## 3. Section-by-Section Code Logic

### Branch A: Unmeasured Confounding (`sens = "uc"`)
This branch calculates **E-values**, which quantify the minimum strength of association that an unmeasured mediator-outcome or exposure-outcome confounder would need to have with both the exposure and the outcome to fully explain away the observed direct or indirect effects.

#### 1. Variable Mapping and Exclusions
The code first filters out unsupported configurations. E-values for unmeasured confounding are currently not supported for survival outcomes (fitted by `coxph` or `survreg`), so it throws an error if survival terms are detected:
```r
if (inherits(reg.input$yreg, "coxph") | inherits(reg.input$yreg, "survreg") | ...)
    stop("sensitivity analysis for unmeasured confounding currently doesn't support survival outcomes")
```
It then maps the target causal indices based on the primary model selected:
* For standard models (`rb`, `wb`, `msm`, `gformula`, `ne`), it reports E-values for the six major effects: `cde`, `pnde`, `tnde`, `pnie`, `tnie`, and `te`.
* For `iorw` (Inverse Odds Ratio Weighting), it restricts estimation to the three identifiable effects: `te`, `pnde`, and `tnie`.

#### 2. Difference-Scale to Ratio-Scale Transformation
E-values are mathematically defined on the **Risk Ratio (RR)** or **Rate Ratio** scale. If the primary analysis fitted a continuous outcome on a difference scale, `cmsens` implements the standard conversion framework developed by VanderWeele and Ding (2017).

It first standardizes the point estimate and standard error by dividing them by the standard deviation of the outcome variable:
$$d = \frac{\text{effect.pe}}{\text{sd}(Y)}$$
$$sd = \frac{\text{effect.se}}{\text{sd}(Y)}$$

The standardized mean difference $d$ is then transformed into an approximate Risk Ratio ($RR$) using the logistic approximation:
$$\text{estRR} = \exp(0.91 \times d)$$
$$\text{lowerRR} = \exp(0.91 \times d - 1.78 \times sd)$$
$$\text{upperRR} = \exp(0.91 \times d + 1.78 \times sd)$$

This transformation is executed in code inside a loop over the effect indices:
```r
outcome.se <- sd(data[, outcome], na.rm = TRUE)
d <- unname(effect.pe / outcome.se)
sd <- unname(effect.se / outcome.se)
for (i in out.index) {
  evalues_mid <- evalues.RR(est = exp(0.91 * d[i]), 
                            lo = exp(0.91 * d[i] - 1.78 * sd[i]), 
                            hi = exp(0.91 * d[i] + 1.78 * sd[i]))
  ...
}
```

#### 3. Standard Ratio-Scale E-Values
If the outcome is binary or count, the effects are already on the ratio scale. `cmsens` extracts the pre-calculated 95% confidence intervals (`ci_lo`, `ci_hi`) directly from the class `cmest` object and passes them straight to the E-value calculator:
```r
evalues_mid <- evalues.RR(est = effect.pe[i], lo = ci_lo[i], hi = ci_hi[i])
```

---

### Branch B: Sensitivity Analysis for Measurement Error (`sens = "me"`)
This branch corrects for systematic measurement error or misclassification in a single covariate (such as a mismeasured mediator or baseline confounder). It provides two core validation methods: **Regression Calibration (RC)** and **Simulation-Extrapolation (SIMEX)**.

#### 1. Input Validation and Configuration
The code enforces that the error model is correctly specified depending on the variable type:
* **Continuous variables** require a numeric vector of measurement error standard deviations (`MEerror`):
  $$\text{Reliability Ratio} = 1 - \frac{\sigma^2_{\text{error}}}{\sigma^2_{\text{observed}}}$$
* **Categorical variables** require a list of misclassification transition probability matrices (`MEerror`), representing the probability of observing category $j$ given that the true category is $i$.

```r
if (MEvartype == "continuous") {
  if (!is.vector(MEerror)) stop("MEerror should be a vector of error standard deviations...")
  out$ME$reliabilityRatio <- 1 - MEerror/sd(data[, MEvariable], na.rm = TRUE)
} else if (MEvartype == "categorical") {
  if (!is.list(MEerror)) stop("MEerror should be a list of misclassification matrices...")
}
```

#### 2. Dynamic Model Wrapping and Swapping
`cmsens` iterates through each specified level of measurement error in `MEerror`. For each level, it clones the original regression parameters and dynamically replaces standard modeling calls with specialized measurement-error-corrected regressions:

##### A. Regression Calibration (`MEmethod = "rc"`)
Regression calibration replaces the mismeasured continuous covariate $W$ with its predicted true value $E[X \mid W, Z]$ estimated from a calibration model, then fits the primary model using this imputed value. In code, this wraps the target regression objects in the internal helper function `rcreg()`:
```r
regs_mid[[r]] <- eval(bquote(rcreg(reg = .(regs_mid[[r]]), 
                                   formula = .(formula(regs_mid[[r]])), 
                                   data = .(data), 
                                   MEvariable = .(MEvariable), 
                                   MEerror = .(MEerror[[i]]), 
                                   variance = .(variance), 
                                   nboot = .(nboot.rc))))
```

##### B. SIMEX (`MEmethod = "simex"`)
The SIMEX algorithm adds simulated pseudo-errors of increasing variance $(1 + \lambda)\sigma^2_u$ to the observed variable, fits the regression models over multiple simulation rounds ($B$), and performs a quadratic extrapolation back to a theoretical zero-error state ($\lambda = -1$). In code, this wraps the regressions in the internal helper function `simexreg()`:
```r
regs_mid[[r]] <- eval(bquote(simexreg(reg = .(regs_mid[[r]]), 
                                      formula = .(formula(regs_mid[[r]])), 
                                      data = .(data), 
                                      MEvariable = .(MEvariable), 
                                      MEvartype = .(MEvartype), 
                                      MEerror = .(MEerror[[i]]), 
                                      variance = .(variance), 
                                      lambda = .(lambda), 
                                      B = .(B))))
```

#### 3. Re-Running the Estimation Loop
Once the standard regressions are replaced with their error-corrected equivalents (`rcreg` or `simexreg`), `cmsens()` binds its local environment to the master estimation engine `estinf()`:
```r
environment(estinf) <- environment()
estinf_res <- estinf()
```
This forces `estinf()` to run the counterfactual simulations or parameter integrations over the adjusted models, automatically producing error-corrected causal point estimates and bootstrapped standard errors.

---

## 4. Helper Dependencies

The `cmsens` module depends heavily on external R libraries and internal package functions to execute calculations:

| Dependency | Origin | Function & Purpose |
| :--- | :--- | :--- |
| **`EValue`** | External Package | `evalues.RR()`: Calculates the E-values for point estimates and confidence intervals on the relative risk scale. |
| **`simex`** | External Package | `check.mc.matrix()`: Validates that user-defined misclassification matrices for categorical SIMEX are row-stochastic. |
| **`rcreg`** | Internal Helper | Performs regression calibration on linear, generalized linear, or survival model components. |
| **`simexreg`** | Internal Helper | Implements SIMEX simulation and quadratic extrapolation back to the error-free state ($\lambda = -1$). |
| **`estinf`** | Internal Helper | Re-calculates point estimates, bootstrapping distributions, and p-values for all causal mediation metrics. |
