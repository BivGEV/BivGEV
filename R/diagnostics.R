

################################################################################################
####################################  diagnostics   ############################################
################################################################################################

diagnostics <- function (x)
{

if (x$BivD == "FGM") {
        cop <- "FGM"
    }
    if (x$BivD == "AMH") {
        cop <- "AMH"
    }
    if (x$BivD == "N") {
        cop <- "Gaussian"
    }
    if (x$BivD == "F") {
        cop <- "Frank"
    }
    if (x$BivD == "C0") {
        cop <- "Clayton"
    }
    if (x$BivD == "C90") {
        cop <- "90 Clayton"
    }
    if (x$BivD == "C180") {
        cop <- "180 Clayton"
    }
    if (x$BivD == "C270") {
        cop <- "270 Clayton"
    }
    if (x$BivD == "J0") {
        cop <- "Joe"
    }
    if (x$BivD == "J90") {
        cop <- "90 Joe"
    }
    if (x$BivD == "J180") {
        cop <- "180 Joe"
    }
    if (x$BivD == "J270") {
        cop <- "270 Joe"
    }
    if (x$BivD == "G0") {
        cop <- "Gumbel"
        lind <- "log( - 1)"
    }
    if (x$BivD == "G90") {
        cop <- "90 Gumbel"
    }
    if (x$BivD == "G180") {
        cop <- "180 Gumbel"
    }
    if (x$BivD == "G270") {
        cop <- "270 Gumbel"
    }
    main.t <- "\nCOPULA:  "
    cat(main.t, cop)
        cat("\nLink function eq.1: GEV")
        cat("\nLink function eq.2: GEV","\n")
        cat("\nTau eq.1:", x$tau.eq1)
        cat("\nTau eq.2:", x$tau.eq2,"\n")

    e.v <- eigen(x$fit$hessian, symmetric = TRUE, only.values = TRUE)$values
    cat("\nLargest absolute gradient value:", max(abs(x$fit$gradient)))
    #if (x$Hes == TRUE)
        mv <- "Observed"
        # else mv <- "Expected"
    if (min(e.v) > 0)
        cat("\n", mv, " information matrix is positive definite\n",
            sep = "")
    else cat("\n", mv, " information matrix is NOT positive definite\n",
        sep = "")
    cat("Eigenvalue range: [", min(e.v), ",", max(e.v), "]\n",
        sep = "")

    cat("\nConverged:", x$fit$converged, "\n")
    cat("\nIterations:", x$fit$iterations, "\n\n")
}


