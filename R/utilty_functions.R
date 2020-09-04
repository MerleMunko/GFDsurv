#' @export
print.casanova<- function(x, ...) {
  cat("Call:", "\n")
  print(x$input$formula)

  cat("\n", "CASANOVA: Cumulative Aalen survival analyis-of-variance:","\n","\n", sep = "")
  print(x$statistic)
}

#' @export
summary.casanova <- function (x, ...) {
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




#' @export
print.medsanova<- function(x, ...) {
  cat("Call:", "\n")
  print(x$input$formula)

  cat("\n", "medSANOVA: Median survival analyis-of-variance:","\n","\n", sep = "")
  print(x$statistic)
}

#' @export
summary.medsanova <- function (x, ...) {
  print(x)
}



#' @export
print.copsanova <- function(x, ...) {
  cat("Call:", "\n")
  print(x$input$formula)

  cat("\n", "Wildbootstrap multiplier:",x$weights, "\n")

  cat("\n", "copSANOVA: concordance probability SANOVA:","\n","\n", sep = "")
  print(x$statistic)
}

#' @export
summary.copsanova <- function (x, ...) {
  print(x)
}



