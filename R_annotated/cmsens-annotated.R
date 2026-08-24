#' Sensitivity Analysis For Unmeasured Confounding and Measurement Error
#' 
#' \code{cmsens} is used to conduct sensitivity analysis for unmeasured confounding via 
#' the \emph{E-value} approach by Vanderweele et al. (2017) and Smith et al. (2019), and 
#' sensitivity analysis for measurement error via \emph{regression calibration} by Carroll 
#' et al. (1995) and \emph{SIMEX} by Cook et al. (1994) and Küchenhoff et al. (2006).
#' 
#' @param object an object of class \code{cmest}.
#' @param sens sensitivity analysis for unmeasured confounding or measurement error. 
#' \code{uc} represents unmeasured confounding and \code{me} represents measurement error. 
#' See \code{Details}.
#' @param MEmethod method for measurement error correction. \code{rc} represents regression 
#' calibration and \code{simex} represents SIMEX. See \code{Details}.
#' @param MEvariable variable measured with error.
#' @param MEvartype type of the variable measured with error. Can be \code{continuous} or 
#' \code{categorical} (first 3 letters are enough).
#' @param MEerror a vector of standard deviations of the measurement error (when \code{MEvartype} 
#' is \code{continuous}) or a list of misclassification matrices (when \code{MEvartype} 
#' is \code{categorical}).
#' @param lambda a vector of lambdas for SIMEX. Default is \code{c(0.5, 1, 1.5, 2)}.
#' @param B number of simulations for SIMEX. Default is \code{200}.
#' @param nboot.rc number of boots for correcting the var-cov matrix of coefficients 
#' with regression calibration. Default is \code{400}.
#' @param x an object of class \code{cmsens}
#' @param digits minimal number of significant digits. See \link{print.default}.
#' @param ... other arguments.
#' 
#' @details 
#' 
#' \strong{Sensitivity Analysis for Unmeasured Confounding}
#' 
#' Currently, sensitivity analysis for unmeasured confounding are available when the outcome 
#' regression model is fitted by \link{lm}, \link{glm}, \link{glm.nb}, \link[mgcv]{gam}, 
#' \link[nnet]{multinom}, \link[MASS]{polr}.
#' 
#' All E-values are reported on the risk or rate ratio scale. If the causal effects are estimated 
#' on the difference scale (i.e., the outcome is continuous), they are transformed into risk 
#' ratios using the transformation described by Vanderweele et al. (2017).
#' 
#' \strong{Sensitivity Analysis for Measurement Error}
#' 
#' Currently, sensitivity analysis for measurement error are available:
#' 
#' 1) when the regression 
#' model involving the variable measured with error is fitted by \link{lm}, 
#' \link{glm} (with family \code{gaussian}, \code{binomial} or \code{poisson}), \link[nnet]{multinom}, 
#' \link[MASS]{polr}, \link[survival]{coxph} or \link[survival]{survreg} and \code{model} is 
#' \code{rb} or \code{gformula};
#' 2) when \code{estimation} is \code{paramfunc}.
#' 
#' Sensitivity analysis for measurement error only supports a single variable measured with error. 
#' Regression calibration requires that the variable measured with error be an independent 
#' continuous variable in the regression it's involved in. SIMEX supports a continuous or 
#' categorical variable measured with error. Quadratic extrapolation method is implemented for 
#' SIMEX.
#' 
#' @return 
#' If \code{sens} is \code{uc}, an object of class \code{cmsens.uc} is returned:
#' \item{call}{the function call,}
#' \item{evalues}{a data frame in which the first three columns are point estimates, lower limits 
#' of 95\% confidence intervals and upper limits of 95\% confidence intervals of causal effects on 
#' the risk or rate ratio scale and the last three columns are E-values on the risk or rate ratio scale,}
#' 
#' If \code{sens} is \code{me}, an object of class \code{cmsens.me} is returned:
#' \item{call}{the function call,}
#' \item{ME}{a list which might contain \code{MEmethod}, \code{MEvariable}, \code{MEvartype}, 
#' \code{MEerror}, \code{lambda}, \code{B}, \code{nboot.rc} and reliability ratio (which is calculated 
#' by \code{1 - MEerror[i]/sd(data[, MEvariable])} for \code{i=1,...,length(MEerror)} when 
#' \code{MEvartype} is \code{continuous}),}
#' \item{naive}{naive causal mediation analysis results,}
#' \item{sens}{a list of causal mediation analysis results after correcting errors in \code{MEerror},}
#' ...
#' 
#' @seealso \link{ggcmsens}, \link{cmdag}, \link{cmest}, \link{rcreg}, \link{simexreg}.
#' 
#' @references 
#' VanderWeele TJ, Ding P (2017). Sensitivity analysis in observational research: introducing the 
#' E-Value. Annals of Internal Medicine. 167(4): 268 - 274.
#' 
#' Smith LH, VanderWeele TJ (2019). Mediational E-values: Approximate sensitivity analysis for 
#' unmeasured mediator-outcome confounding. Epidemiology. 30(6): 835 - 837.
#' 
#' Carrol RJ, Ruppert D, Stefanski LA, Crainiceanu C (2006). Measurement Error in Nonlinear Models: 
#' A Modern Perspective, Second Edition. London: Chapman & Hall.
#' 
#' Cook JR, Stefanski LA (1994). Simulation-extrapolation estimation in parametric measurement error 
#' models. Journal of the American Statistical Association, 89(428): 1314 - 1328.
#' 
#' Küchenhoff H, Mwalili SM, Lesaffre E (2006). A general method for dealing with misclassification 
#' in regression: the misclassification SIMEX. Biometrics. 62(1): 85 - 96.
#' 
#' Stefanski LA, Cook JR (1995). Simulation-extrapolation: the measurement error jackknife. 
#' Journal of the American Statistical Association. 90(432): 1247 - 56.
#' 
#' Valeri L, Lin X, VanderWeele TJ (2014). Mediation analysis when a continuous mediator is measured 
#' with error and the outcome follows a generalized linear model. Statistics in 
#' medicine, 33(28): 4875 – 4890.
#' 
#' @examples 
#' 
#' \dontrun{
#' library(CMAverse)
#' 
#' # 10 boots are used for illustration
#' naive <- cmest(data = cma2020, model = "rb", outcome = "contY", 
#' exposure = "A", mediator = c("M1", "M2"), 
#' basec = c("C1", "C2"), EMint = TRUE, 
#' mreg = list("logistic", "multinomial"), yreg = "linear", 
#' astar = 0, a = 1, mval = list(0, "M2_0"), 
#' estimation = "imputation", inference = "bootstrap", nboot = 10)
#' 
#' exp1 <- cmsens(object = naive, sens = "uc")
#' exp2 <- cmsens(object = naive, sens = "me", MEmethod = "rc", 
#' MEvariable = "C1", MEvartype = "con", 
#' MEerror = c(0.1, 0.2))
#' summary(exp2)
#' 
#' # B = 10 is used for illustration
#' exp3 <- cmsens(object = naive, sens = "me", MEmethod = "simex", 
#' MEvariable = "M1", MEvartype = "cat", 
#' MEerror = list(matrix(c(0.95,0.05,0.05,0.95), nrow = 2), 
#' matrix(c(0.9,0.1,0.1,0.9), nrow = 2)), B = 10)
#' summary(exp3)
#' }
#' 
#' @importFrom stats model.frame printCoefmat family sd getCall
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom EValue evalues.RR
#' @importFrom simex check.mc.matrix
#' @importFrom ggplot2 ggplot geom_errorbar aes geom_point ylab geom_hline position_dodge2 
#' scale_colour_hue theme element_blank facet_grid
#' 
#' @export
cmsens <- function(object = NULL, sens = "uc", MEmethod = "simex", MEvariable = NULL, 
                   MEvartype = NULL, MEerror = NULL, lambda = c(0.5, 1, 1.5, 2), 
                   B = 200, nboot.rc = 400) {
  
  # Ensure the master causal mediation object is provided
  if (is.null(object)) stop("Unspecified object")
  
  cl <- match.call()
  out <- list(call = cl)
  
  # --- EXTRACT PARAMETERS FROM THE PRIMARY CMEST OBJECT ---
  data <- object$data
  model <- object$methods$model
  reg.input <- object$reg.input
  effect.pe <- object$effect.pe
  effect.se <- object$effect.se
  outcome <- object$variables$outcome
  
  #################################################################################################
  ############################ Sensitivity Analysis for Unmeasured Confounding #######################
  #################################################################################################
  
  if (sens == "uc") {
    
    # Extract the original fitted regressions. Handle multiple imputation results if present.
    if (object$multimp$multimp) {
      reg.output <- object$reg.output[[1]]
    } else {
      reg.output <- object$reg.output
    }
    
    # GUARDRAIL: Survival outcomes are currently unsupported for E-value sensitivity analyses
    if (inherits(reg.input$yreg, "coxph") | inherits(reg.input$yreg, "survreg") | 
        (is.character(reg.input$yreg) && reg.input$yreg %in% c("coxph", "aft_exp", "aft_weibull"))
    ) stop("sensitivity analysis for unmeasured confounding currently doesn't support survival outcomes")
    
    # STEP 1: Determine which causal effect indices to report
    # Standard models estimate 6 effects (cde, pnde, tnde, pnie, tnie, te)
    if (model %in% c("rb", "wb", "msm", "gformula", "ne")) out.index <- 1:6
    # Inverse odds ratio weighting (iorw) only estimates 3 effects (te, pnde, tnie)
    if (model == "iorw") out.index <- 1:3
    
    # Subset point estimates and standard errors to the active indices
    effect.pe <- effect.pe[out.index]
    effect.se <- effect.se[out.index]
    
    # Retrieve the appropriate outcome regression model (yregTot for iorw, yreg for others)
    yreg <- switch((model == "iorw") + 1, "1" = reg.output$yreg, "2" = reg.output$yregTot)
    
    # Unwrap measurement error regressions (rc or simex class objects) back to their naive form
    yreg4class <- switch(("rcreg" %in% class(yreg) | "simexreg" %in% class(yreg)) + 1, "1" = yreg, "2" = yreg$NAIVEreg)
    
    if (inherits(yreg4class, "lm")) yreg_family <- family(yreg)$family
    
    evalues <- c()
    
    # STEP 2: Calculate E-values
    # Scenario A: Continuous outcomes on the difference scale (lm with linear, Gamma, etc. families)
    if (inherits(yreg4class, "lm") && (yreg_family %in% c("gaussian","Gamma","inverse.gaussian","quasi"))) {
      
      # Standardize the point estimates and SEs by the standard deviation of the outcome
      # to obtain Cohen's d values.
      outcome.se <- sd(data[, outcome], na.rm = TRUE)
      d <- unname(effect.pe / outcome.se)
      sd <- unname(effect.se / outcome.se)
      
      # Transform standardized effects to Risk Ratio (RR) equivalents using the VanderWeele (2017) approx:
      # RR ≈ exp(0.91 * d), with lower/upper limits calculated on the 1.78 * sd scaling factor.
      for (i in out.index) {
        evalues_mid <- evalues.RR(
          est = exp(0.91 * d[i]), 
          lo = exp(0.91 * d[i] - 1.78 * sd[i]), 
          hi = exp(0.91 * d[i] + 1.78 * sd[i])
        )
        evalues_mid <- c(evalues_mid[1, ], evalues_mid[2, ])
        evalues <- rbind(evalues, evalues_mid)
      }
      
    } else {
      # Scenario B: Categorical or count outcomes already estimated on the Risk/Rate/Odds Ratio scale
      ci_lo <- object$effect.ci.low
      ci_hi <- object$effect.ci.high
      
      # Map estimates and original confidence intervals directly to the EValue calculator
      for (i in out.index) {
        evalues_mid <- evalues.RR(est = effect.pe[i], lo = ci_lo[i], hi = ci_hi[i])
        evalues_mid <- c(evalues_mid[1, ], evalues_mid[2, ])
        evalues <- rbind(evalues, evalues_mid)
      }
    }
    
    # Format and label the output E-value matrix
    colnames(evalues) <- c("estRR", "lowerRR", "upperRR", "Evalue.estRR", "Evalue.lowerRR", "Evalue.upperRR")
    rownames(evalues) <- names(effect.pe)
    out$evalues <- evalues
    class(out) <- "cmsens.uc"
    
  #################################################################################################
  ############################ Sensitivity Analysis for Measurement Error #########################
  #################################################################################################
    
  } else if (sens == "me") {
    
    # GUARDRAILS AND ARGUMENT VALIDATION:
    if(!model %in% c("rb", "gformula")) stop("Measurement error correction currently only supports model = 'rb' or 'gformula'")
    if (length(MEvariable) == 0) stop("Unspecified MEvariable")
    if (length(MEvariable) != 1) stop("Currently only supports one variable measured with error")
    if (length(MEvartype) != length(MEvariable)) stop("length(MEvartype) != length(MEvariable)")
    
    # Standardize types and variable class checks
    if (MEvartype == "con") MEvartype <- "continuous"
    if (MEvartype == "cat") MEvartype <- "categorical"
    
    if (MEvartype == "continuous") {
      if (!is.vector(MEerror)) stop("MEerror should be a vector of error standard deviations for a continuous variable measured with error")
    } else if (MEvartype == "categorical") {
      if (!is.list(MEerror)) stop("MEerror should be a list of misclassification matrices for a categorical variable measured with error")
    }
    
    # Save the measurement error metadata
    if (MEmethod == "rc") out$ME <- list(MEmethod = MEmethod, MEvariable = MEvariable, MEvartype = MEvartype, MEerror = MEerror)
    if (MEmethod == "simex") out$ME <- list(MEmethod = MEmethod, MEvariable = MEvariable, MEvartype = MEvartype, MEerror = MEerror, lambda = lambda, B = B)
    
    # Calculate reliability ratios for continuous mismeasured variables: 1 - (error_variance / total_variance)
    if (MEvartype == "continuous") out$ME$reliabilityRatio <- 1 - MEerror/sd(data[, MEvariable], na.rm = TRUE)
    
    out$naive <- object
    
    # Re-extract execution parameters from the baseline cmest object
    n <- nrow(data)
    full <- object$methods$full
    casecontrol <- object$methods$casecontrol
    yrare <- object$methods$yrare
    yprevalence <- object$methods$yprevalence
    estimation <- object$methods$estimation
    inference <- object$methods$inference
    outcome <- object$variables$outcome
    event <- object$variables$event
    exposure <- object$variables$exposure
    mediator <- object$variables$mediator
    EMint <- object$variables$EMint
    basec <- object$variables$basec
    postc <- object$variables$postc
    yreg <- reg.input$yreg
    ereg <- reg.input$ereg
    mreg <- reg.input$mreg
    wmnomreg <- reg.input$wmnomreg
    wmdenomreg <- reg.input$wmdenomreg
    postcreg <- reg.input$postcreg
    a <- object$ref$a
    astar <- object$ref$astar
    mval <- object$ref$mval
    yval <- object$ref$yval
    basecval <- object$ref$basecval
    nboot <- object$methods$nboot
    boot.ci.type <- object$methods$boot.ci.type
    nRep <- object$methods$nRep
    multimp <- object$multimp$multimp
    args_mice <- object$multimp$args_mice
    
    # Bind environment and fit initial regressions
    environment(regrun) <- environment()
    regs <- regrun()
    
    # Delta method requires analytical variance estimation during the measurement error correction
    if (inference == "delta") {
      variance <- TRUE
      if (MEmethod == "rc") out$ME$nboot.rc <- nboot.rc
    } else {
      variance <- FALSE
    }
    
    sens <- list()
    
    # LOOP OVER EACH SPECIFIED LEVEL OF MEASUREMENT ERROR / MISCLASSIFICATION
    for (i in 1:length(MEerror)) {
      regs_mid <- regs
      
      # --- METHOD 1: REGRESSION CALIBRATION ---
      if (MEmethod == "rc") {
        if (MEvartype != "continuous") stop("Regression calibration only supports a continuous variable measured with error")
        
        # Replace traditional regression fits with rcreg() modeling
        for (r in 1:length(regs_mid)) {
          if (!is.null(regs_mid[[r]])) {
            if (inherits(regs_mid[[r]], "list")) {
              # Iterate across lists of models (such as multiple mediator regressions)
              regs_mid[[r]] <- lapply(1:length(regs_mid[[r]]), function(x) {
                eval(bquote(rcreg(
                  reg = .(regs_mid[[r]][[x]]), 
                  formula = .(formula(regs_mid[[r]][[x]])), 
                  data = .(data), 
                  weights = .(model.frame(regs_mid[[r]][[x]])$'(weights)'), 
                  MEvariable = .(MEvariable), 
                  MEerror = .(MEerror[[i]]), 
                  variance = .(variance), 
                  nboot = .(nboot.rc)
                )))
              })
            } else {
              # Handle standalone single regression models
              regs_mid[[r]] <- eval(bquote(rcreg(
                reg = .(regs_mid[[r]]), 
                formula = .(formula(regs_mid[[r]])), 
                data = .(data), 
                weights = .(model.frame(regs_mid[[r]])$'(weights)'), 
                MEvariable = .(MEvariable), 
                MEerror = .(MEerror[[i]]), 
                variance = .(variance), 
                nboot = .(nboot.rc)
              )))
            }
          }
        }
        
      # --- METHOD 2: SIMULATION-EXTRAPOLATION (SIMEX) ---
      } else if (MEmethod == "simex") {
        
        # Replace traditional regression fits with simexreg() modeling
        for (r in 1:length(regs_mid)) {
          if (!is.null(regs_mid[[r]])) {
            if (inherits(regs_mid[[r]], "list")) {
              regs_mid[[r]] <- lapply(1:length(regs_mid[[r]]), function(x) {
                eval(bquote(simexreg(
                  reg = .(regs_mid[[r]][[x]]), 
                  formula = .(formula(regs_mid[[r]][[x]])), 
                  data = .(data), 
                  weights = .(model.frame(regs_mid[[r]][[x]])$'(weights)'), 
                  MEvariable = .(MEvariable), 
                  MEvartype = .(MEvartype), 
                  MEerror = .(MEerror[[i]]), 
                  variance = .(variance), 
                  lambda = .(lambda), 
                  B = .(B)
                )))
              })
            } else {
              regs_mid[[r]] <- eval(bquote(simexreg(
                reg = .(regs_mid[[r]]), 
                formula = .(formula(regs_mid[[r]])), 
                data = .(data), 
                weights = .(model.frame(regs_mid[[r]])$'(weights)'), 
                MEvariable = .(MEvariable), 
                MEvartype = .(MEvartype), 
                MEerror = .(MEerror[[i]]), 
                variance = .(variance), 
                lambda = .(lambda), 
                B = .(B)
              )))
            }
          }
        }
      } else {
        stop("Unsupported MEmethod; use 'rc' or 'simex'")
      }
      
      # Re-map the newly corrected measurement-error models to the execution variables
      yreg <- regs_mid$yreg
      mreg <- regs_mid$mreg
      ereg <- regs_mid$ereg
      wmnomreg <- regs_mid$wmnomreg
      wmdenomreg <- regs_mid$wmdenomreg
      postcreg <- regs_mid$postcreg
      
      # Setup progress bar if bootstrapping is the chosen inference method
      if (inference == "bootstrap") {
        env <- environment()
        counter <- 0
        progbar <- txtProgressBar(min = 0, max = nboot, style = 3)
      }
      
      # Execute causal effect estimation using the error-corrected regressions
      environment(estinf) <- environment()
      estinf_res <- estinf()
      
      # Store point estimates, SEs, intervals, p-values, and model structures for this error level
      sens[[i]] <- estinf_res[c("effect.pe", "effect.se", "effect.ci.low", "effect.ci.high", "effect.pval")]
      sens[[i]]$reg.output <- estinf_res$reg.output
    }
    
    out$sens <- sens
    class(out) <- "cmsens.me"
    
  } else {
    stop("Unsupported sens; use 'uc' or 'me'")
  }
  
  return(out)
}

#' @describeIn cmsens Print results of \code{cmsens.uc} nicely
#' @export
print.cmsens.uc <- function(x, ...) {
  cat("Sensitivity Analysis For Unmeasured Confounding \n")
  cat("\nEvalues on the risk or rate ratio scale: \n")
  print(x$evalues)
}

#' @describeIn cmsens Print results of \code{cmsens.me} nicely
#' @export
print.cmsens.me <- function(x, ...) {
  cat("Sensitivity Analysis For Measurement Error \n \n")
  cat("The variable measured with error: ")
  cat(x$ME$MEvariable)
  cat("\nType of the variable measured with error: ")
  cat(x$ME$MEvartype)
  cat("\n\n")
  
  for (i in 1:length(x$sens)) {
    cat(paste0("# Measurement error ", i, ": \n"))
    if (x$ME$MEvartype == "continuous") print(x$ME$MEerror[i])
    if (x$ME$MEvartype == "categorical") print(x$ME$MEerror[[i]])
    cat(paste0("\n## Error-corrected regressions for measurement error ", i, ": \n\n"))
    
    # PRINT ERROR CORRECTED REGRESSIONS: Handle structural printing based on multiple imputation status
    if (!x$naive$multimp$multimp) {
      regnames <- names(x$sens[[i]]$reg.output)
      for (name in regnames) {
        if (name == "yreg") {
          cat("### Outcome regression:\n")
          if (inherits(x$sens[[i]]$reg.output$yreg, "svyglm")) {
            x$sens[[i]]$reg.output$yreg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$yreg,design = getCall(x$sens[[.(i)]]$reg.output$yreg)$design, family = getCall(x$sens[[.(i)]]$reg.output$yreg)$family, evaluate = FALSE)))
            x$sens[[i]]$reg.output$yreg$survey.design$call <- eval(bquote(as.call(update(summary(x$sens[[.(i)]]$reg.output$yreg)$survey.design, data = getCall(summary(x$sens[[.(i)]]$reg.output$yreg)$survey.design)$data, weights = getCall(summary(x$sens[[.(i)]]$reg.output$yreg)$survey.design)$weights, evaluate = FALSE))))
            print(x$sens[[i]]$reg.output$yreg)
          } else if (inherits(x$sens[[i]]$reg.output$yreg, "rcreg")|inherits(x$sens[[i]]$reg.output$yreg, "simexreg")) {
            x$sens[[i]]$reg.output$yreg$call <- eval(bquote(update(x$sens[[i]]$reg.output$yreg, reg = getCall(x$sens[[.(i)]]$reg.output$yreg)$reg, data=getCall(x$sens[[.(i)]]$reg.output$yreg)$data, weights=getCall(x$sens[[.(i)]]$reg.output$yreg)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output$yreg)
          } else {
            x$sens[[i]]$reg.output$yreg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$yreg, data=getCall(x$sens[[.(i)]]$reg.output$yreg)$data, weights=getCall(x$sens[[.(i)]]$reg.output$yreg)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output$yreg)
          }
        }
        if (name == "yregTot") {
          cat("### Outcome regression for the total effect: \n")
          if (inherits(x$sens[[i]]$reg.output$yregTot, "rcreg")|inherits(x$sens[[i]]$reg.output$yregTot, "simexreg")) {
            x$sens[[i]]$reg.output$yregTot$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$yregTot, reg = getCall(x$sens[[.(i)]]$reg.output$yregTot)$reg, data=getCall(x$sens[[.(i)]]$reg.output$yregTot)$data, weights=getCall(x$sens[[.(i)]]$reg.output$yregTot)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output$yregTot)
          } else {
            x$sens[[i]]$reg.output$yregTot$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$yregTot,data=getCall(x$sens[[.(i)]]$reg.output$yregTot)$data, weights=getCall(x$sens[[.(i)]]$reg.output$yregTot)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output$yregTot)
          }
        }
        if (name == "yregDir") {
          cat("### Outcome regression for the direct effect: \n")
          if (inherits(x$sens[[i]]$reg.output$yregDir, "rcreg")|inherits(x$sens[[i]]$reg.output$yregDir, "simexreg")) {
            x$sens[[i]]$reg.output$yregDir$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$yregDir, reg = getCall(x$sens[[.(i)]]$reg.output$yregDir)$reg, data=getCall(x$sens[[.(i)]]$reg.output$yregDir)$data, weights=getCall(x$sens[[.(i)]]$reg.output$yregDir)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output$yregDir)
          } else {
            x$sens[[i]]$reg.output$yregDir$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$yregDir,data=getCall(x$sens[[.(i)]]$reg.output$yregDir)$data, weights=getCall(x$sens[[.(i)]]$reg.output$yregDir)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output$yregDir)
          }
        }
        if (name == "ereg") {
          cat("### Exposure regression for weighting: \n")
          if (inherits(x$sens[[i]]$reg.output$ereg, "rcreg")|inherits(x$sens[[i]]$reg.output$ereg, "simexreg")) {
            x$sens[[i]]$reg.output$ereg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$ereg, reg = getCall(x$sens[[.(i)]]$reg.output$ereg)$reg, data=getCall(x$sens[[.(i)]]$reg.output$ereg)$data, weights=getCall(x$sens[[.(i)]]$reg.output$ereg)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output$ereg)
          } else {
            x$sens[[i]]$reg.output$ereg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$ereg,data=getCall(x$sens[[.(i)]]$reg.output$ereg)$data, weights=getCall(x$sens[[.(i)]]$reg.output$ereg)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output$ereg)
          }
        }
        if (name == "mreg") {
          cat("### Mediator regressions: \n")
          for (j in 1:length(x$sens[[i]]$reg.output$mreg)) {
            if (inherits(x$sens[[i]]$reg.output$mreg[[j]], "svyglm")) {
              x$sens[[i]]$reg.output$mreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$mreg[[.(j)]], design = getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$design, family = getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$family, evaluate = FALSE)))
              x$sens[[i]]$reg.output$mreg[[j]]$survey.design$call <- eval(bquote(as.call(update(summary(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$survey.design, data = getCall(summary(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$survey.design)$data, weights = getCall(summary(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$survey.design)$weights, evaluate = FALSE))))
              print(x$sens[[i]]$reg.output$mreg[[j]])
            } else if (inherits(x$sens[[i]]$reg.output$mreg[[j]], "rcreg") | inherits(x$sens[[i]]$reg.output$mreg[[j]], "simexreg")){
              x$sens[[i]]$reg.output$mreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$mreg[[.(j)]], reg=getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$reg, data=getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$data, weights=getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$weights, evaluate = FALSE)))
              print(x$sens[[i]]$reg.output$mreg[[j]])
            } else {
              x$sens[[i]]$reg.output$mreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$mreg[[.(j)]], data=getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$data, weights=getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$weights, evaluate = FALSE)))
              print(x$sens[[i]]$reg.output$mreg[[j]])
            }
            if (j < length(x$sens[[i]]$reg.output$mreg)) cat("\n")
          }
        }
        if (name == "wmdenomreg") {
          cat("### Mediator regressions for weighting (denominator): \n")
          for (j in 1:length(x$sens[[i]]$reg.output$wmdenomreg)) {
            if (inherits(x$sens[[i]]$reg.output$wmdenomreg[[j]], "rcreg") | inherits(x$sens[[i]]$reg.output$wmdenomreg[[j]], "simexreg")) {
              x$sens[[i]]$reg.output$wmdenomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]], reg=getCall(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]])$reg, data=getCall(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]])$data, weights=getCall(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]])$weights, evaluate = FALSE)))
              print(x$sens[[i]]$reg.output$wmdenomreg[[j]])
            } else {
              x$sens[[i]]$reg.output$wmdenomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]], data=getCall(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]])$data, weights=getCall(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]])$weights, evaluate = FALSE)))
              print(x$sens[[i]]$reg.output$wmdenomreg[[j]])
            }
            if (j < length(x$sens[[i]]$reg.output$wmdenomreg)) cat("\n")
          }
        }
        if (name == "wmnomreg") {
          cat("### Mediator regressions for weighting (nominator): \n")
          for (j in 1:length(x$sens[[i]]$reg.output$wmnomreg)) {
            if (inherits(x$sens[[i]]$reg.output$wmnomreg[[j]], "rcreg") | inherits(x$sens[[i]]$reg.output$wmnomreg[[j]], "simexreg")) {
              x$sens[[i]]$reg.output$wmnomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]], reg=getCall(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]])$reg, data=getCall(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]])$data, weights=getCall(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]])$weights, evaluate = FALSE)))
              print(x$sens[[i]]$reg.output$wmnomreg[[j]])
            } else {
              x$sens[[i]]$reg.output$wmnomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]], data=getCall(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]])$data, weights=getCall(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]])$weights, evaluate = FALSE)))
              print(x$sens[[i]]$reg.output$wmnomreg[[j]])
            }
            if (j < length(x$sens[[i]]$reg.output$wmnomreg)) cat("\n")
          }
        }
        if (name == "postcreg") {
          cat("### Regressions for mediator-outcome confounders affected by the exposure: \n")
          for (j in 1:length(x$sens[[i]]$reg.output$postcreg)) {
            if (inherits(x$sens[[i]]$reg.output$postcreg[[j]], "rcreg") | inherits(x$sens[[i]]$reg.output$postcreg[[j]], "simexreg")) {
              x$sens[[i]]$reg.output$postcreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]], reg=getCall(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]])$reg, data=getCall(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]])$data, weights=getCall(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]])$weights, evaluate = FALSE)))
              print(x$sens[[i]]$reg.output$postcreg[[j]])
            } else {
              x$sens[[i]]$reg.output$postcreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]], data=getCall(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]])$data, weights=getCall(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]])$weights, evaluate = FALSE)))
              print(x$sens[[i]]$reg.output$postcreg[[j]])
            }
            if (j < length(x$sens[[i]]$reg.output$postcreg)) cat("\n")
          }
        }
        cat("\n")
      }
    } else {
      # Handle printing for multiple imputation fits (looping across imputed datasets)
      for (m in 1:length(x$sens[[i]]$reg.output)) {
        cat(paste("### Regressions with imputed dataset", m, "\n\n"))
        regnames <- names(x$sens[[i]]$reg.output[[m]])
        for (name in regnames) {
          if (name == "yreg") {
            cat("#### Outcome regression: \n")
            if (inherits(x$sens[[i]]$reg.output[[m]]$yreg, "svyglm")) {
              x$sens[[i]]$reg.output[[m]]$yreg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg, design = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$design, family = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$family, evaluate = FALSE)))
              x$sens[[i]]$reg.output[[m]]$yreg$survey.design$call <- eval(bquote(as.call(update(summary(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$survey.design, data = getCall(summary(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$survey.design)$data, weights = getCall(summary(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$survey.design)$weights, evaluate = FALSE))))
              print(x$sens[[i]]$reg.output[[m]]$yreg)
            } else if (inherits(x$sens[[i]]$reg.output[[m]]$yreg, "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$yreg, "simexreg")){
              x$sens[[i]]$reg.output[[m]]$yreg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg, reg = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$reg, data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$data, weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$weights, evaluate = FALSE)))
              print(x$sens[[i]]$reg.output[[m]]$yreg)
            } else {
              x$sens[[i]]$reg.output[[m]]$yreg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg, data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$data, weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$weights, evaluate = FALSE)))
              print(x$sens[[i]]$reg.output[[m]]$yreg)
            }
          }
          if (name == "yregTot") {
            cat("#### Outcome regression for the total effect: \n")
            if (inherits(x$sens[[i]]$reg.output[[m]]$yregTot, "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$yregTot, "simexreg")){
              x$sens[[i]]$reg.output[[m]]$yregTot$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot, reg = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot)$reg, data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot)$data, weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot)$weights, evaluate = FALSE)))
              print(x$sens[[i]]$reg.output[[m]]$yregTot)
            } else {
              x$sens[[i]]$reg.output[[m]]$yregTot$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot, data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot)$data, weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot)$weights, evaluate = FALSE)))
              print(x$sens[[i]]$reg.output[[m]]$yregTot)
            }
          }
          if (name == "yregDir") {
            cat("#### Outcome regression for the direct effect: \n")
            if (inherits(x$sens[[i]]$reg.output[[m]]$yregDir, "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$yregDir, "simexreg")){
              x$sens[[i]]$reg.output[[m]]$yregDir$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir, reg = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir)$reg, data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir)$data, weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir)$weights, evaluate = FALSE)))
              print(x$sens[[i]]$reg.output[[m]]$yregDir)
            } else {
              x$sens[[i]]$reg.output[[m]]$yregDir$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir, data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir)$data, weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir)$weights, evaluate = FALSE)))
              print(x$sens[[i]]$reg.output[[m]]$yregDir)
            }
          }
          if (name == "ereg") {
            cat("#### Exposure regression for weighting: \n")
            if (inherits(x$sens[[i]]$reg.output[[m]]$ereg, "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$ereg, "simexreg")){
              x$sens[[i]]$reg.output[[m]]$ereg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg, reg = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg)$reg, data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg)$data, weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg)$weights, evaluate = FALSE)))
              print(x$sens[[i]]$reg.output[[m]]$ereg)
            } else {
              x$sens[[i]]$reg.output[[m]]$ereg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg, data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg)$data, weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg)$weights, evaluate = FALSE)))
              print(x$sens[[i]]$reg.output[[m]]$ereg)
            }
          }
          if (name == "mreg") {
            cat("#### Mediator regressions: \n")
            for (j in 1:length(x$sens[[i]]$reg.output[[m]]$mreg)) {
              if (inherits(x$sens[[i]]$reg.output[[m]]$mreg[[j]], "svyglm")) {
                x$sens[[i]]$reg.output[[m]]$mreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]], design = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$design, family = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$family, evaluate = FALSE)))
                x$sens[[i]]$reg.output[[m]]$mreg[[j]]$survey.design$call <- eval(bquote(as.call(update(summary(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$survey.design, data = getCall(summary(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$survey.design)$data, weights = getCall(summary(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$survey.design)$weights, evaluate = FALSE))))
                print(x$sens[[i]]$reg.output[[m]]$mreg[[j]])
              } else if (inherits(x$sens[[i]]$reg.output[[m]]$mreg[[j]], "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$mreg[[j]], "simexreg")){
                x$sens[[i]]$reg.output[[m]]$mreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]], reg=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$reg, data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$data, weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$weights, evaluate = FALSE)))
                print(x$sens[[i]]$reg.output[[m]]$mreg[[j]])
              } else {
                x$sens[[i]]$reg.output[[m]]$mreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]], data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$data, weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$weights, evaluate = FALSE)))
                print(x$sens[[i]]$reg.output[[m]]$mreg[[j]])
              }
              if (j < length(x$sens[[i]]$reg.output[[m]]$mreg)) cat("\n")
            }
          }
          if (name == "wmdenomreg") {
            if (!is.null(x$sens[[i]]$reg.output[[m]]$wmdenomreg)) {
              cat("#### Mediator regressions for weighting (denominator): \n")
              for (j in 1:length(x$sens[[i]]$reg.output[[m]]$wmdenomreg)) {
                if (inherits(x$sens[[i]]$reg.output[[m]]$wmdenomreg[[j]], "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$wmdenomreg[[j]], "simexreg")){
                  x$sens[[i]]$reg.output[[m]]$wmdenomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]], reg=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]])$reg, data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]])$data, weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]])$weights, evaluate = FALSE)))
                  print(x$sens[[i]]$reg.output[[m]]$wmdenomreg[[j]])
                } else {
                  x$sens[[i]]$reg.output[[m]]$wmdenomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]], data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]])$data, weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]])$weights, evaluate = FALSE)))
                  print(x$sens[[i]]$reg.output[[m]]$wmdenomreg[[j]])
                }
                if (j < length(x$sens[[i]]$reg.output[[m]]$wmdenomreg)) cat("\n")
              }
            }
          }
          if (name == "wmnomreg") {
            if (!is.null(x$sens[[i]]$reg.output[[m]]$wmnomreg)) {
              cat("#### Mediator regressions for weighting (nominator): \n")
              for (j in 1:length(x$sens[[i]]$reg.output[[m]]$wmnomreg)) {
                if (inherits(x$sens[[i]]$reg.output[[m]]$wmnomreg[[j]], "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$wmnomreg[[j]], "simexreg")){
                  x$sens[[i]]$reg.output[[m]]$wmnomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]], reg=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]])$reg, data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]])$data, weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]])$weights, evaluate = FALSE)))
                  print(x$sens[[i]]$reg.output[[m]]$wmnomreg[[j]])
                } else {
                  x$sens[[i]]$reg.output[[m]]$wmnomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]], data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]])$data, weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]])$weights, evaluate = FALSE)))
                  print(x$sens[[i]]$reg.output[[m]]$wmnomreg[[j]])
                }
                if (j < length(x$sens[[i]]$reg.output[[m]]$wmnomreg)) cat("\n")
              }
            }
          }
          if (name == "postcreg") {
            if (!is.null(x$sens[[i]]$reg.output[[m]]$postcreg)) {
              cat("#### Regressions for mediator-outcome confounders affected by the exposure: \n")
              for (j in 1:length(x$sens[[i]]$reg.output[[m]]$postcreg)) {
                if (inherits(x$sens[[i]]$reg.output[[m]]$postcreg[[j]], "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$postcreg[[j]], "simexreg")){
                  x$sens[[i]]$reg.output[[m]]$postcreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]], reg=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]])$reg, data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]])$data, weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]])$weights, evaluate = FALSE)))
                  print(x$sens[[i]]$reg.output[[m]]$postcreg[[j]])
                } else {
                  x$sens[[i]]$reg.output[[m]]$postcreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]], data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]])$data, weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]])$weights, evaluate = FALSE)))
                  print(x$sens[[i]]$reg.output[[m]]$postcreg[[j]])
                }
                if (j < length(x$sens[[i]]$reg.output[[m]]$postcreg)) cat("\n")
              }
            }
          }
          cat("\n")
        }
      }
    }
  }
}
