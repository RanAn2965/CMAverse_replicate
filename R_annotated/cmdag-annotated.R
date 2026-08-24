#' Directed Acyclic Graph (DAG) Visualization for Causal Mediation Analysis
#'
#' This script contains the annotated version of the `cmdag` function from the `CMAverse` R package.
#' It visualizes the causal pathways between an exposure, mediator(s), outcome, and confounders.
#' 
#' @param outcome Name of the outcome variable (Y)
#' @param exposure Name of the exposure variable (A)
#' @param mediator Vector of variable name(s) for the mediator(s) (M)
#' @param basec (Optional) Vector of baseline confounders not affected by the exposure (C)
#' @param postc (Optional) Vector of post-exposure confounders affected by the exposure (L)
#' @param x.outcome,y.outcome 2D coordinates for plotting the outcome node
#' @param x.exposure,y.exposure 2D coordinates for plotting the exposure node
#' @param x.mediator,y.mediator 2D coordinates for plotting the mediator node(s)
#' @param x.basec,y.basec 2D coordinates for plotting baseline confounders
#' @param x.postc,y.postc 2D coordinates for plotting post-exposure confounders
#' @param caption.width Max character width before wrapping the caption text
#' @param caption.size Text size (pts) for the caption
#' @param ... Additional graphical arguments passed down to ggdag()
#'
#' @importFrom ggplot2 ggproto theme labs element_text
#' @importFrom ggdag dagify ggdag theme_dag_blank
#' @importFrom stringr str_wrap
#'
#' @export
cmdag <- function(outcome = NULL, 
                  exposure = NULL, 
                  mediator = NULL, 
                  basec = NULL, 
                  postc = NULL, 
                  x.outcome = 4, 
                  x.exposure = 0, 
                  x.mediator = 2, 
                  x.basec = 2, 
                  x.postc = 2, 
                  y.outcome = 0, 
                  y.exposure = 0, 
                  y.mediator = 1, 
                  y.basec = 2, 
                  y.postc = -0.5, 
                  caption.width = 50, 
                  caption.size = 10, 
                  ...) {
  
  # -------------------------------------------------------------------------
  # 1. INPUT VALIDATION
  # -------------------------------------------------------------------------
  # Causal mediation analysis strictly requires an exposure, mediator, and outcome.
  # If any of these core nodes are missing, the function stops immediately.
  if (is.null(outcome) | is.null(exposure) | is.null(mediator)) {
    stop("Unspecified outcome, exposure or mediator")
  }
  
  # -------------------------------------------------------------------------
  # 2. GLOBAL ENVIRONMENT WORKAROUND
  # -------------------------------------------------------------------------
  # Standard R packages sometimes run into scope issues with custom ggplot2 layers
  # in sister packages (like ggdag). To ensure that 'ggproto' methods resolve
  # properly during visualization, this helper dynamically copies 'ggproto' 
  # from the ggplot2 namespace directly into the global execution environment.
  assign2glob <- function(key, val, pos) {
    assign(key, val, envir = as.environment(pos))
  }
  assign2glob("ggproto", ggplot2::ggproto, 1L)
  
  # -------------------------------------------------------------------------
  # 3. CAUSAL MODEL BRANCHING (4 STRUCTURAL SCENARIOS)
  # -------------------------------------------------------------------------
  # The function determines which DAG to draw based on the presence of
  # baseline confounders (basec / C) and post-exposure confounders (postc / L).
  
  # === SCENARIO A: Simple Mediation Model (No Confounders) ===
  if (is.null(basec) && is.null(postc)) {
    # Define structural equations:
    # - Outcome (Y) is affected by Exposure (A) and Mediator (M)
    # - Mediator (M) is affected by Exposure (A)
    dag <- dagify(Y ~ A + M, 
                  M ~ A, 
                  coords = list(
                    x = c(A = x.exposure, Y = x.outcome, M = x.mediator), 
                    y = c(A = y.exposure, Y = y.outcome, M = y.mediator)
                  ))
    
    # Generate plot, remove coordinate axes, wrap caption labels, and style
    p <- ggdag(dag, ...) + 
      theme_dag_blank() + 
      labs(caption = paste0(
        "A (exposure): ", str_wrap(exposure, width = caption.width), "\n",
        "M (mediator): ", str_wrap(paste(mediator, collapse = ", "), width = caption.width), "\n",
        "Y (outcome): ", str_wrap(outcome, width = caption.width)
      )) + 
      theme(plot.caption = element_text(hjust = 0, face = "italic", size = caption.size))
    
    print(p)
    
  # === SCENARIO B: Baseline Confounding Only (C) ===
  } else if (!is.null(basec) && is.null(postc)) {
    # Define structural equations:
    # - Outcome (Y) is affected by Exposure (A), Mediator (M), and Baseline Confounders (C)
    # - Mediator (M) is affected by Exposure (A) and Baseline Confounders (C)
    # - Exposure (A) is affected by Baseline Confounders (C) [representing selection/confounding]
    dag <- dagify(Y ~ A + M + C, 
                  M ~ A + C, 
                  A ~ C, 
                  coords = list(
                    x = c(A = x.exposure, Y = x.outcome, M = x.mediator, C = x.basec), 
                    y = c(A = y.exposure, Y = y.outcome, M = y.mediator, C = y.basec)
                  ))
    
    # Generate plot, strip backgrounds, wrap metadata, and style
    p <- ggdag(dag, ...) + 
      theme_dag_blank() + 
      labs(caption = paste0(
        "A (exposure): ", str_wrap(exposure, width = caption.width), "\n",
        "M (mediator): ", str_wrap(paste(mediator, collapse = ", "), width = caption.width), "\n",
        "Y (outcome): ", str_wrap(outcome, width = caption.width), "\n",
        "C (confounders not affected by the exposure): ", str_wrap(paste(basec, collapse = ", "), width = caption.width)
      )) + 
      theme(plot.caption = element_text(hjust = 0, face = "italic", size = caption.size))
    
    print(p)
    
  # === SCENARIO C: Post-Exposure Confounding Only (L) ===
  } else if (is.null(basec) && !is.null(postc)) {
    # Define structural equations:
    # - Outcome (Y) is affected by Exposure (A), Mediator (M), and Post-Exposure Confounders (L)
    # - Mediator (M) is affected by Exposure (A) and Post-Exposure Confounders (L)
    # - Post-Exposure Confounder (L) is causally downstream of Exposure (A) [A -> L]
    dag <- dagify(Y ~ A + M + L, 
                  M ~ A + L, 
                  L ~ A, 
                  coords = list(
                    x = c(A = x.exposure, Y = x.outcome, M = x.mediator, L = x.postc), 
                    y = c(A = y.exposure, Y = y.outcome, M = y.mediator, L = y.postc)
                  ))
    
    # Generate plot, strip backgrounds, wrap metadata, and style
    p <- ggdag(dag, ...) + 
      theme_dag_blank() + 
      labs(caption = paste0(
        "A (exposure): ", str_wrap(exposure, width = caption.width), "\n",
        "M (mediator): ", str_wrap(paste(mediator, collapse = ", "), width = caption.width), "\n",
        "Y (outcome): ", str_wrap(outcome, width = caption.width), "\n",
        "L (confounders affected by the exposure): ", str_wrap(paste(postc, collapse = ", "), width = caption.width)
      )) + 
      theme(plot.caption = element_text(hjust = 0, face = "italic", size = caption.size))
    
    print(p)
    
  # === SCENARIO D: Baseline & Post-Exposure Confounding (Full Model) ===
  } else if (!is.null(basec) && !is.null(postc)) {
    # Define structural equations representing the most complex causal mediation scenario:
    # - Outcome (Y) is affected by A, M, C, and L
    # - Mediator (M) is affected by A, C, and L
    # - Baseline Confounders (C) affect Exposure (A)
    # - Exposure (A) affects Post-Exposure Confounders (L)
    dag <- dagify(Y ~ A + M + C + L, 
                  M ~ A + C + L, 
                  A ~ C, 
                  L ~ A, 
                  coords = list(
                    x = c(A = x.exposure, Y = x.outcome, M = x.mediator, C = x.basec, L = x.postc), 
                    y = c(A = y.exposure, Y = y.outcome, M = y.mediator, C = y.basec, L = y.postc)
                  ))
    
    # Generate plot, strip backgrounds, wrap metadata, and style
    p <- ggdag(dag, ...) + 
      theme_dag_blank() + 
      labs(caption = paste0(
        "A (exposure): ", str_wrap(exposure, width = caption.width), "\n",
        "M (mediator): ", str_wrap(paste(mediator, collapse = ", "), width = caption.width), "\n",
        "Y (outcome): ", str_wrap(outcome, width = caption.width), "\n",
        "C (confounders not affected by the exposure): ", str_wrap(paste(basec, collapse = ", "), width = caption.width), "\n",
        "L (confounders affected by the exposure): ", str_wrap(paste(postc, collapse = ", "), width = caption.width)
      )) + 
      theme(plot.caption = element_text(hjust = 0, face = "italic", size = caption.size))
    
    print(p)
  }
  
  # -------------------------------------------------------------------------
  # 4. ENVIRONMENT CLEANUP
  # -------------------------------------------------------------------------
  # To avoid leaving objects behind in the user's workspace, the temporarily
  # defined 'ggproto' object is cleanly removed from the Global Environment.
  rm(ggproto, envir = .GlobalEnv)
}
