

############################################################################################
###############################   print.summary.BivGEV     #################################
############################################################################################

print.summary.BivGEV <- function (x, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"),...)
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
    cp <- "  theta = "
    as.p <- x$theta.a
    ct <- "  tau = "
    kt.p <- x$tau.a
    n.sel <- "  n.sel = "
    t1    <- "  tau.eq1 = "
    t2    <- "  tau.eq2 = "
        m1l <- "GEV"
        m2l <- "GEV"

    cat(main.t, cop)
    #    cat("\nMARGIN 1: GEV")
    #    cat("\nMARGIN 2: GEV")

cat("\n\nEQUATION 1")
    #cat("\nLink function:", m1l, "\n")
    cat("\nLink function:", m1l, "\n")
    cat("Formula: ")
    print(x$formula[[1]])
    cat("\n")
    cat("Parametric coefficients:\n")
    printCoefmat(x$tableP1, digits = digits, signif.stars = signif.stars,
        na.print = "NA", ...)
    cat("\n")

cat("\nEQUATION 2")
    #cat("\nLink function for mu.2:", m2l, "\n")
    cat("\nLink function:", m2l, "\n")
    cat("Formula: ")
    print(x$formula[[2]])
    cat("\n")
    cat("Parametric coefficients:\n")
    printCoefmat(x$tableP2, digits = digits, signif.stars = signif.stars,
        na.print = "NA", ...)
    cat("\n")

CIrs <- colMeans(x$CIrs, na.rm = TRUE)
CIkt <- colMeans(x$CIkt, na.rm = TRUE)
decim <- 3


if (x$Model == "BivGEV") {
cat("\nn = ", x$n, cp, format(as.p, digits = decim), "(",
            format(CIrs[1], digits = decim), ",", format(CIrs[2],
                digits = decim), ")", ct, format(kt.p, digits = decim),
            "(", format(CIkt[1], digits = decim), ",", format(CIkt[2],
                digits = decim), ")", "\ntotal edf = ", format(x$t.edf,
                digits = decim), t1, format(x$tau.eq1,digits = decim), t2, format(x$tau.eq2,digits = decim), "\n\n", sep = "")
}


if (x$Model == "SampleSelGEV") {
      cat("\nn = ", x$n, n.sel, x$n.sel, cp, format(as.p, digits = decim), "(",
      format(CIrs[1], digits = decim), ",", format(CIrs[2],
       digits = decim), ")", ct, format(kt.p, digits = decim),
      "(", format(CIkt[1], digits = decim), ",", format(CIkt[2],
      digits = decim), ")", "\ntotal edf = ", format(x$t.edf,
      digits = decim), t1, format(x$tau.eq1,digits = decim), t2, format(x$tau.eq2,digits = decim), "\n\n", sep = "")
}


}



