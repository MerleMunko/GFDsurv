coeff.check <- function (cross = cross, rg = rg)  
{
  if (is.null(rg)) {
    return(list(rg = rg, cross = cross, indep = TRUE))
  }
  pairs <- rg
  npairs <- length(rg)
  maxpower <- max(sapply(rg, sum))
  if (maxpower == 0) {
    return(list(rg = rg, cross = cross, indep = TRUE))
  }
  coeff <- function(x) {
    r <- x[1]
    g <- x[2]
    k <- 0:g
    c(numeric(r), choose(g, k) * (-1)^k, numeric(max(maxpower - 
                                                       r - g, 0)))
  }
  match.fun(coeff)
  Mcoeff <- sapply(rg, coeff)
  if (cross == TRUE) {
    coeff_cross <- c(1, -2, numeric(max(c(0, maxpower - 1))))
    Mcoeff <- cbind(Mcoeff, coeff_cross)
    npairs <- npairs + 1
    pairs <- c(pairs, "cross")
  }
  change <- FALSE
  rank_coeff <- qr(Mcoeff)$rank
  ndepend_rows <- npairs - rank_coeff
  if (ndepend_rows > 0) {
    change <- TRUE
    for (i in 1:ndepend_rows) {
      index <- sapply(1:npairs, function(x) {
        qr(Mcoeff[, -x])$rank == rank_coeff
      })
      index_remove <- min(which(index))
      Mcoeff <- Mcoeff[, -index_remove]
      pairs <- pairs[-index_remove]
      npairs <- npairs - 1
    }
    if (pairs[npairs] == "cross") {
      return(list(rg = pairs[-npairs], cross = TRUE, indep = FALSE))
    }
    else {
      return(list(rg = pairs, cross = FALSE, indep = FALSE))
    }
  }
  else {
    return(list(rg = rg, cross = cross, indep = TRUE))
  }
}
