# Deep-Dive Code Review & Technical Notes: `cmest.R` in `CMAverse`

This document provides a comprehensive breakdown of the internal logic, architectural patterns, helper dependencies, and bootstrap loops of the master function `cmest()` in the `CMAverse` package.

---

## 1. Executive Structural Overview

The `cmest()` function is the central "air traffic controller" for causal mediation analysis in `CMAverse`. Instead of housing thousands of lines of regression and calculation logic in a single file, the package is highly modular. `cmest.R` coordinates input validation, runs regressions, performs causal estimation, and manages standard error inferences using two external scripts: `regrun.R` (for model fitting) and `estinf.R` (for causal effect estimation and statistical inference).

```
[ User Input (Data & Spec) ]
             │
             ▼
      ┌──────────────┐
      │   cmest()    │  ◄─── Handles validation, imputation, and progress bars
      └──────┬───────┘
             │
      ┌──────┴───────┐
      │  regrun()    │  ◄─── Translates text specifications and fits models
      └──────┬───────┘
             │
      ┌──────┴───────┐
      │   estinf()   │  ◄─── Computes effects (Parametric/Imputation)
      └──────────────┘
```

---

## 2. The Core Architectural Secret: Lexical Environment Sharing

One of the most notable programming patterns in `CMAverse` is how `cmest.R` communicates with `regrun.R` and `estinf.R`. Rather than passing dozens of arguments through function calls, `cmest.R` dynamically overrides the lexical environment of these functions at runtime.

### The Code Snippet
```r
environment(regrun) <- environment()
regs <- regrun()

# ... Subsequent calculations ...

environment(estinf) <- environment()
out <- c(out, estinf())
```

### How This Works
* **Environment Binding**: `environment(regrun) <- environment()` resets the parent/enclosing environment of `regrun` to the **active execution frame** of `cmest()`.
* **Shared Memory Access**: When `regrun()` and `estinf()` execute, they can read and write all local variables inside `cmest()` (such as `data`, `exposure`, `mediator`, `yreg`, `mreg`, etc.) as if they were local variables.
* **Benefits & Trade-offs**: This avoids passing long, fragile parameter lists. However, it tightly couples the scripts together, meaning any modification to variable names in `cmest.R` will immediately break `regrun.R` or `estinf.R`.

---

## 3. Step-by-Step Code Logic Flow

The execution of `cmest()` can be broken down into four distinct phases:

### Phase 1: Argument Validation and Input Guardrails
The first ~300 lines of `cmest()` are dedicated to rigorous sanity checks and error throwing:
* **Missing Data Handlers**: If missing values are detected in `data` and multiple imputation is inactive (`multimp = FALSE`), the function stops with an error. If `multimp = TRUE`, it prepares arguments for the `mice` package.
* **Weighting Restrictions**: It enforces theoretical rules from causal literature. For example, because weights computed from continuous variables are highly unstable, it throws an error if weighting-based models (`wb`, `msm`, or `iorw`) are combined with non-categorical mediators or exposures.
* **Case-Control Checks**: If `casecontrol = TRUE`, it verifies that the outcome is binary. It also enforces that the user either specifies the outcome prevalence (`yprevalence`) or flags the case as rare (`yrare = TRUE`), which is necessary to apply correction factors to odds ratios.

### Phase 2: Translation and Model Fitting via `regrun()`
The local variables are handed to `regrun()`, which is responsible for fitting the statistical models:
* **Character-to-Model Translation**: Users can specify regressions as character strings (e.g., `yreg = "logistic"`, `mreg = list("linear")`). `regrun()` translates these into R formulas (e.g., `Y ~ A + M + C`) and calls standard regression functions.
* **Custom Models**: If users provide pre-fitted R model objects (e.g., from `lm` or `coxph`), `regrun()` bypasses formula compilation and extracts the coefficients directly.
* **Weighting Computations**: For weighted estimation approaches (like MSM or MSMs with post-exposure confounding), `regrun()` automatically calculates stable or unstabilized inverse-probability weights.

### Phase 3: Estimation via `estinf()`
Once the regressions are fit, `estinf()` is invoked to calculate the point estimates:
* **Closed-Form Parameter Functions (`paramfunc`)**:
  * Limited to a single mediator and no post-exposure confounding.
  * It extracts the regression coefficients (the $\beta$ and $\theta$ parameters) and plugs them directly into analytic closed-form formulas (like those derived in Valeri & VanderWeele, 2013).
* **Direct Counterfactual Imputation (`imputation`)**:
  * The default and most flexible pathway. It works by simulating individual-level counterfactuals.
  * For each individual $i$ in the dataset, it uses the fitted mediator models to simulate what their mediator value *would have been* under the active exposure level ($M_{a,i}$) and the control exposure level ($M_{a^*,i}$) by drawing randomly from the mediator's estimated conditional distribution.
  * It then plugs these simulated counterfactual mediator values back into the outcome model (`yreg`) to predict each individual’s counterfactual outcome.
  * Taking the sample averages of these individual predictions imputes the population-level counterfactual expectations (e.g., $E[Y_aM_{a^*}]$), which are then subtracted or divided to compute direct, indirect, and interactive effects.

---

## 4. Anatomy of the Bootstrap Loop

When `inference = "bootstrap"` is selected, `cmest()` bypasses standard asymptotic variance calculations (like the Delta method) and executes a non-parametric bootstrap routine.

### Step-by-Step Bootstrap Mechanics

1. **Progress Bar Setup**:
   ```r
   if (inference == "bootstrap") {
     env <- environment()
     counter <- 0
     progbar <- txtProgressBar(min = 0, max = nboot, style = 3)
   }
   ```
   An interactive progress bar is displayed to the user.

2. **Bootstrapping over Multiple Imputations**:
   If multiple imputation is turned on (`multimp = TRUE`), `cmest()` performs bootstrapping in a nested fashion:
   * It draws bootstrap samples from each of the imputed datasets generated by `mice()`. This accounts for both **imputation uncertainty** and **sampling variance** simultaneously.

3. **Resampling with Replacement**:
   In each bootstrap step:
   * It draws a random sample of size $N$ with replacement from the active dataset.
   * It re-fits all regressions using `regrun()` and re-estimates all causal effects using `estinf()` on this resampled dataset.
   * The calculated point estimates are appended as a new row in a bootstrap matrix.
   * The progress bar is incremented: `setTxtProgressBar(progbar, counter)`.

4. **Confidence Interval Construction**:
   Once the bootstrap loop completes, `cmest()` extracts the standard deviation of each column to serve as the standard error (`effect.se`). It then constructs $95\%$ confidence intervals:
   * **Percentile Method (`per`)**: It takes the raw $2.5^{\text{th}}$ and $97.5^{\text{th}}$ quantiles of the bootstrap distribution.
   * **BCa Method (`bca`)**: Calls the `boot` package to adjust the endpoints based on median bias and acceleration (skewness), which is highly recommended when the bootstrap distribution is asymmetric.

---

## 5. Key Helper & Methodological Dependencies

To function correctly, `cmest.R` imports and relies heavily on the following external package dependencies:

* **`msm`**: Used exclusively for `deltamethod()` to calculate standard errors when `inference = "delta"` is selected (only available under `paramfunc` estimation).
* **`nnet`**: Used to fit multinomial logistics (`multinom`) for nominal categorical mediators.
* **`MASS`**: Used to fit ordered categorical variables (`polr`) and overdispersed count data (`glm.nb`).
* **`survival`**: Used to fit Cox proportional hazards (`coxph`) and accelerated failure time (`survreg`) survival outcomes.
* **`survey`**: Used to fit weighted generalized linear models (`svyglm`) to accommodate complex survey designs or weights calculated in-house.
* **`medflex`**: Specifically imported to execute the *Natural Effect Model* pathway using the `neImpute` function under the hood.
* **`boot`**: Used to compute bias-corrected and accelerated (BCa) bootstrap confidence intervals.
* **`mice`**: Handles missing data through multiple imputation prior to model estimation.

---

## 6. Guide to Reading and Navigating the R Code

When opening `cmest.R` to modify or study it, use these line-range landmarks as your map:
* **Lines 1–135**: Extensive Roxygen2 documentation detailing parameters, methodological references, and user examples.
* **Lines 136–160**: Imports from core R packages and helper namespaces.
* **Lines 161–330**: Argument restriction and safety checks (checking dataset formatting, rare-disease options, categorical requirements).
* **Lines 331–345**: Running the regressions by binding the local lexical environment to `regrun()`.
* **Lines 346–360**: Initializing and managing the Bootstrap progress bar if `bootstrap` is selected.
* **Lines 361–370**: Invoking `estinf()` to perform the final causal estimations and return the compiled `cmest` S3 object.
