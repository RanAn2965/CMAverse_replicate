# =================================================================================
# CAVERSE: FULLY ANNOTATED MASTER ESTIMATION FUNCTION (cmest)
# =================================================================================
# This file contains the complete, mathematically unmodified source code of the 
# central cmest() function from the CMAverse package. 
#
# No active code has been altered. Educational, step-by-step technical comments and 
# annotations have been added inline to clarify the validation logic, environment-sharing 
# mechanics, and execution pathways.
# =================================================================================

#' Causal Mediation Analysis
#'
#' @export
cmest <- function(
  # Core Data and Model Selection
  data = NULL,           # Dataset (dataframe)
  model = "rb",          # Causal approach: "rb" (regression-based), "wb" (weighting-based),
                         # "iorw" (inverse odds-ratio weighting), "ne" (natural effect model),
                         # "msm" (marginal structural model), or "gformula" (g-formula)
  full = TRUE,           # TRUE = output full list of causal effects; FALSE = reduced list
  casecontrol = FALSE,   # TRUE = case-control study design (binary outcome)
  yrare = NULL,          # TRUE = rare outcome assumption for case-control studies
  yprevalence = NULL,    # Population prevalence of the outcome for case-control corrections
  
  # Methodological Configurations
  estimation = "imputation", # Point estimate method: "paramfunc" (analytical formulas)
                             # or "imputation" (direct counterfactual imputation)
  inference = "bootstrap",   # Variance method: "delta" (Delta method) or "bootstrap"
  
  # Structural Variables
  outcome = NULL,        # Variable name of the outcome (Y)
  event = NULL,          # Event indicator for survival outcomes (coxph or AFT models)
  exposure = NULL,       # Variable name of the exposure (A)
  mediator = NULL,       # Vector of variable name(s) of the mediator(s) (M)
  EMint = NULL,          # TRUE = allow exposure-mediator interaction in outcome regression
  basec = NULL,          # Vector of baseline confounders (C) not affected by exposure
  postc = NULL,          # Vector of post-exposure/mediator-outcome confounders (L) affected by exposure
  
  # Model Objects / Formulas
  yreg = NULL,           # Regression model for the outcome Y
  mreg = NULL,           # List of regressions for each mediator variable in M
  wmnomreg = NULL,       # List of nominator regressions for MSM weights
  wmdenomreg = NULL,     # List of denominator regressions for MSM weights
  ereg = NULL,           # Exposure regression model (used in weighting / iorw)
  postcreg = NULL,       # List of regressions for post-exposure confounders L
  
  # Causal Contrasts and Reference Values
  astar = 0,             # Reference/control value of exposure (a*)
  a = 1,                 # Active value of exposure (a)
  mval = NULL,           # List of controlled values for mediators (m)
  yval = NULL,           # Category of Y at which risk/odds ratio effects are computed
  basecval = NULL,       # Specific values of C to condition on if using "paramfunc"
  
  # Bootstrap and Simulation Parameters
  nboot = 200,           # Number of bootstrap resamples (if inference = "bootstrap")
  boot.ci.type = "per",  # Bootstrap interval type: "per" (percentile) or "bca" (bias-corrected)
  nRep = 5,              # Samples drawn per unit in Natural Effect models
  multimp = FALSE,       # TRUE = run Multiple Imputation chained equations (mice) first
  args_mice = NULL       # Optional list of auxiliary arguments passed to mice()
) {

  # -------------------------------------------------------------------------------
  # Step 1: Initialize Call and Output Records
  # -------------------------------------------------------------------------------
  # match.call() captures the exact function call arguments entered by the user
  # for later extraction in the summary() and print() methods.
  cl <- match.call()
  out <- list(call = cl)

  # -------------------------------------------------------------------------------
  # Step 2: Input Validation - Data Integrity Guardrails
  # -------------------------------------------------------------------------------
  if (is.null(data)) stop("Unspecified data")
  data <- as.data.frame(data)
  
  # Missing data is highly problematic for causal counterfactual estimation.
  # Unless the user explicitly enables Multiple Imputation (multimp = TRUE),
  # CMAverse halts execution if any NAs are present.
  if (sum(is.na(data)) > 0 && !multimp) {
    stop("NAs in the data; delete rows with NAs from the data or set multimp = TRUE")
  }
  
  out$data <- data
  n <- nrow(data)  # Sample size, used extensively downstream for simulation/boot loops

  # -------------------------------------------------------------------------------
  # Step 3: Methodological Scope Constraints
  # -------------------------------------------------------------------------------
  # Confirm the user selected a supported causal paradigm.
  if (!model %in% c("rb", "wb", "iorw", "ne", "gformula", "msm")) {
    stop("Select model from 'rb', 'wb', 'iorw', 'ne', 'gformula', 'msm'")
  }
  out$methods$model <- model

  if (!is.logical(full)) stop("full should be TRUE or FALSE")
  out$methods$full <- full

  # -------------------------------------------------------------------------------
  # Step 4: Case-Control and Selection Bias Corrections
  # -------------------------------------------------------------------------------
  if (!is.logical(casecontrol)) stop("casecontrol should be TRUE or FALSE")
  out$methods$casecontrol <- casecontrol
  
  if (casecontrol) {
    # Case-control designs require a binary outcome
    if (length(unique(data[, outcome])) != 2) {
      stop("When casecontrol is TRUE, the outcome must be binary")
    }
    # To compute valid population-level risks/odds from a retrospectively sampled
    # case-control cohort, we must either assume the disease is extremely rare (yrare = TRUE)
    # or adjust weights using the true population prevalence (yprevalence).
    if (is.null(yprevalence) && yrare != TRUE) {
      stop("When casecontrol is TRUE, specify yprevalence or set yrare to be TRUE")
    }
    # Imputation methods require sampling from continuous conditional distributions.
    # Without true prevalence, direct counterfactual simulation is biased, forcing
    # the user to fall back onto parametric closed-form models.
    if (is.null(yprevalence) && !(estimation %in% c("para", "paramfunc"))) {
      stop("When casecontrol is TRUE, specify yprevalence or use estimation = 'paramfunc'")
    }
    if (!is.null(yprevalence)) {
      if (!is.numeric(yprevalence)) stop("yprevalence should be numeric")
      out$methods$yprevalence <- yprevalence
    } else {
      out$methods$yrare <- yrare
    }
  }

  # -------------------------------------------------------------------------------
  # Step 5: Parser/Abbreviation Standardizer
  # -------------------------------------------------------------------------------
  # Standardize four-letter parameter abbreviations to complete words for internal logic.
  if (estimation == "para") estimation <- "paramfunc"
  if (estimation == "impu") estimation <- "imputation"
  if (inference == "delt")  inference <- "delta"
  if (inference == "boot")  inference <- "bootstrap"

  # -------------------------------------------------------------------------------
  # Step 6: Verify Model-Estimation-Inference Alignments
  # -------------------------------------------------------------------------------
  # Regression-based (rb) models can use analytical formulas (paramfunc) or simulations.
  # Weighting-based, G-formula, MSMs, and IORW methods can ONLY be estimated via Imputation.
  if (model == "rb" && !estimation %in% c("paramfunc", "imputation")) {
    stop("When model = 'rb', select estimation from 'paramfunc', 'imputation'")
  }
  if (model != "rb" && !estimation == "imputation") {
    stop("Use estimation = 'imputation'")
  }

  # Validate parametric configuration constraints
  if (estimation == "paramfunc") {
    # Closed-form solutions require string names to fetch corresponding coefficients
    if (!is.character(yreg) | length(yreg) != 1) {
      stop("When estimation = 'paramfunc', select yreg from 'linear', 'logistic', 'loglinear', 'poisson', 'quasipoisson', 'negbin', 'coxph', 'aft_exp', 'aft_weibull'")
    }
    if (!yreg %in% c("linear", "logistic", "loglinear", "poisson", "quasipoisson", "negbin", "coxph", "aft_exp", "aft_weibull")) {
      stop("When estimation = 'paramfunc', select yreg from 'linear', 'logistic', 'loglinear', 'poisson', 'quasipoisson', 'negbin', 'coxph', 'aft_exp', 'aft_weibull'")
    }
    # Closed-form analytical integration is mathematically intractable for multiple mediators
    if (length(mediator) > 1) {
      stop("'paramfunc' only supports a single mediator")
    }
    if (!is.character(mreg[[1]]) | length(mreg[[1]]) != 1) {
      stop("When estimation = 'paramfunc', select mreg[[1]] from 'linear', 'logistic', 'multinomial'")
    }
    if(!mreg[[1]] %in% c("linear", "logistic", "multinomial")) {
      stop("When estimation = 'paramfunc', select mreg[[1]] from 'linear', 'logistic', 'multinomial'")
    }
  }
  out$methods$estimation <- estimation

  # Delta method standard errors are calculated analytically and can only be used 
  # with closed-form parameter function estimation. Imputations must use bootstrapping.
  if (!inference %in% c("delta", "bootstrap")) {
    stop("Select inference from 'delta', 'bootstrap'")
  }
  if (estimation == "imputation" && inference == "delta") {
    stop("Use inference = 'bootstrap' when estimation = 'imputation'")
  }
  out$methods$inference <- inference

  if (inference == "bootstrap") {
    if (!is.numeric(nboot)) stop("nboot should be numeric")
    if (!boot.ci.type %in% c("per", "bca")) stop("Select boot.ci.type from 'per', 'bca'")
    out$methods$nboot <- nboot
    out$methods$boot.ci.type <- boot.ci.type
  }

  # -------------------------------------------------------------------------------
  # Step 7: Variables & Covariates Role Mapping
  # -------------------------------------------------------------------------------
  # Outcome Variable
  if (length(outcome) == 0) stop("Unspecified outcome")
  if (length(outcome) > 1)  stop("length(outcome) > 1")
  out$variables$outcome <- outcome

  # Survival Outcome Event Checks
  if (is.character(yreg)) {
    if (yreg %in% c("coxph", "aft_exp", "aft_weibull") && !is.null(event)) {
      out$variables$event <- event
    }
  } else {
    if (!is.null(event)) warning("event is ignored when yreg is not character")
  }

  # Exposure Variable
  if (length(exposure) == 0) stop("Unspecified exposure")
  if (length(exposure) > 1)  stop("length(exposure) > 1")
  out$variables$exposure <- exposure

  # Mediator Variables
  if (length(mediator) == 0) stop("Unspecified mediator")
  out$variables$mediator <- mediator

  # Exposure-Mediator Interaction
  if (model != "iorw") {
    if (!is.logical(EMint)) stop("EMint should be TRUE or FALSE")
    out$variables$EMint <- EMint
  }

  # Confounders (Baseline and Post-exposure)
  if (!is.null(basec)) out$variables$basec <- basec
  
  if (length(postc) != 0) {
    # Post-exposure confounders (L) introduce intermediate pathways affected by exposure
    # that confound mediator-outcome relations. Only MSMs and G-formulas can identify these.
    if (!model %in% c("msm", "gformula")) {
      stop("When postc is not empty, select model from 'msm' and 'gformula'")
    }
    out$variables$postc <- postc
  }

  # -------------------------------------------------------------------------------
  # Step 8: Assemble Registration Inputs
  # -------------------------------------------------------------------------------
  # Bundle the appropriate models together based on the active causal method.
  if (model == "rb")        out$reg.input <- list(yreg = yreg, mreg = mreg)
  if (model == "wb")        out$reg.input <- list(yreg = yreg)
  if (model == "gformula")  out$reg.input <- list(yreg = yreg, mreg = mreg)
  if (model == "msm")       out$reg.input <- list(yreg = yreg, mreg = mreg, wmnomreg = wmnomreg, wmdenomreg = wmdenomreg)
  if (model == "iorw")      out$reg.input <- list(yreg = yreg, ereg = ereg)
  if (model == "ne")        out$reg.input <- list(yreg = yreg)
  
  # Exposure models (ereg) are fitted when confounders are present in weighting models
  if (length(basec) != 0 && model %in% c("wb", "msm")) {
    out$reg.input$ereg <- ereg
  }
  # Post-exposure confounder models (postcreg) are fitted in G-formula approaches
  if (length(postc) != 0 && model == "gformula") {
    out$reg.input$postcreg <- postcreg
  }

  # -------------------------------------------------------------------------------
  # Step 9: Establish Reference Value Contours
  # -------------------------------------------------------------------------------
  if (is.null(a) | is.null(astar)) stop("Unspecified a or astar")

  if (model != "iorw") {
    if (!is.list(mval)) stop("mval should be a list")
    if (length(mval) != length(mediator)) stop("length(mval) != length(mediator)")
    for (p in 1:length(mval)) {
      if (is.null(mval[[p]])) stop(paste0("Unspecified mval[[", p, "]]"))
    }
  }

  if (!(model == "rb" && estimation == "paramfunc" && length(basec) != 0) && !is.null(basecval)) {
    warning("basecval is ignored")
  }

  if (model == "ne") {
    if (!is.numeric(nRep)) stop("nRep should be numeric")
    out$methods$nRep <- nRep
  }

  # -------------------------------------------------------------------------------
  # Step 10: Multiple Imputation Validation & Setup
  # -------------------------------------------------------------------------------
  out$multimp <- list(multimp = multimp)
  if (!is.logical(multimp)) stop("multimp should be TRUE or FALSE")

  if (multimp) {
    if (!is.null(args_mice) && !is.list(args_mice)) stop("args_mice should be a list")
    if (is.null(args_mice)) args_mice <- list()
    args_mice$print <- FALSE  # Mute mice console prints inside the main function
    if (!is.null(args_mice$data)) warning("args_mice$data is overwritten by data")
    args_mice$data <- data
    out$multimp$args_mice <- args_mice
  }

  # Weighting-based methods, natural effect models, and MSMs do not support active NA values.
  # Missing data must be resolved prior to weight calculation to avoid calculation failures.
  if (!multimp && model %in% c("wb", "iorw", "ne", "msm")) {
    if (sum(is.na(data)) > 0) {
      stop("Selected model doesn't support missing values; use multimp = TRUE")
    }
  }

  # -------------------------------------------------------------------------------
  # Step 11: Lexical Scoping Magic (Regressions via regrun)
  # -------------------------------------------------------------------------------
  # R's scoping environment manipulation:
  # environment(regrun) <- environment() dynamically changes the parent environment of
  # the external regrun function to be the active execution workspace of cmest().
  # This lets regrun directly read and write cmest's local variables (yreg, data, etc.)
  # without copying or passing massive parameter lists.
  environment(regrun) <- environment()
  regs <- regrun()  # Executes model fitting

  # Map fitted models back to local objects
  yreg <- regs$yreg
  ereg <- regs$ereg
  mreg <- regs$mreg
  wmnomreg <- regs$wmnomreg
  wmdenomreg <- regs$wmdenomreg
  postcreg <- regs$postcreg

  # -------------------------------------------------------------------------------
  # Step 12: Progress Bar Initialization (Inference via Bootstrap)
  # -------------------------------------------------------------------------------
  if (inference == "bootstrap") {
    env <- environment()
    counter <- 0
    # Creates an ASCII progress bar (1 to nboot) to show boots updates on console
    progbar <- txtProgressBar(min = 0, max = nboot, style = 3)
  }

  # -------------------------------------------------------------------------------
  # Step 13: Estimation and Inference (via estinf)
  # -------------------------------------------------------------------------------
  # Bind the execution environment again for the effect calculator helper 'estinf'.
  # This allows estinf to execute closed-form or direct counterfactual simulation
  # and bootstrap loops while modifying our final 'out' list object directly.
  environment(estinf) <- environment()
  out <- c(out, estinf())

  class(out) <- "cmest"  # Assign S3 class for print/summary dispatcher
  return(out)
}
# =================================================================================
# END OF cmest FUNCTION CODE
# =================================================================================
