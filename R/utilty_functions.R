#' @export
print.GFDsurv <- function(x, ...) {
  cat("Call:", "\n")
  print(x$input$formula)

  cat("\n", "P-Value Standard Test", "\n", sep = "")
  print(x$pvalues_stat)
  cat("\n", "P-Value Permutation Test", "\n", sep = "")
  print(x$pvalues_per)
}

#' @export
summary.GFDsurv <- function (x, ...) {
  if ( length(x$rg) == 0 ){
    cat("The chosen weights are linearly independent.", "\n",
        "The test is based on the crossing weight.", "\n","\n")
  }else{
    rg_rep <- paste0("c(", x$rg[[1]][1], ", ", x$rg[[1]][2], ")" )
    if ( length(x$rg) > 1 ){
      for (i in 2:length(x$rg)){
        rg_rep <- paste0( rg_rep, ", c(", x$rg[[i]][1], ", ", x$rg[[i]][2], ")" )
      }
    }
    cat("The chosen weights are", if ( x$indep == FALSE){ " not"},
        " linearly independent.", "\n", "The test is based on ",
        if ( x$cross == TRUE){ "the crossing weight and "},
        length(x$rg), " ","weight", if (length(x$rg) > 1){"s"}, " with exponents ", "\n",
        "     ", "c(r,g)=  ", rg_rep, ".", "\n","\n", sep="")
  }
  print(x)
  }

