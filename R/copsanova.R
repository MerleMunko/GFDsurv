#' CASANOVA: Cumulative Aalen survival analyis-of-variance
#'
#' The function \code{casanova} calculates the Wald-type statistic based on the
#' combination of differently weighted Nelson-Aalen-type integrals. Respective p-values
#' are obtained by a \eqn{\chi^2}-approximation and a permutation approach, respectively.
#' @param formula A model \code{formula} object. The left hand side contains the time variable and the right
#'  hand side contains the factor variables of interest. An interaction term must be
#'  specified.
#' @param event The name of censoring status indicator with values 0=censored and
#' 1=uncensored.
#' The default choice is "event"
#' @param data A data.frame, list or environment containing the variables in formula
#' and the censoring status
#' indicator. Default option is \code{NULL}.
#' @param cross logical. Should the crossing weight w(x) = 1 - 2x be included?
#'  The default is \code{TRUE}.
#' @param rg A list (or \code{NULL}) containing the exponents \code{c(r, g)} of the
#' weights
#'   \eqn{w(x) = x^r (1-x)^g}. Both exponents need to be natural numbers including 0.
#'  Default is \code{list( c(0, 0) )} corresponding to the log-rank weight.
#' @param nperm The number of permutations used for calculating the permuted p-value.
#'   The default option is 1999.
#' @param alpha A number specifying the significance level; the default is 0.05.
#' @param nested.levels.unique A logical specifying whether the levels of the nested
#' factor(s) are labeled uniquely or not.
#'  Default is FALSE, i.e., the levels of the nested factor are the same for each
#'  level of the main factor.
#' @details
#' The \code{casanova} function calculates the Wald-type statistic of weighted
#' Nelson-Aalen type integrals
#' for general factorial survival designs. Crossed as well as hierachically nested designs are
#' implemented. Moreover, the approach allows the combination of
#' different weights into a
#' joint statistic. The user can choose between weights of the following form:
#' w(x) = 1 - 2x (\code{cross = TRUE}) and w(x) = x^r * (1-x)^g for natural numbers
#' r,g (including 0). The function automatically check whether the specified weights
#' fulfill the linear independence assumption and choose a subset of linearly independent
#' weights if the original weights violate the aforemention assumption.
#'
#'   The \code{casanova} function returns the test statistic as well as two
#'   corresponding p-values: the first is based on a \eqn{chi^2} approximation and
#'   the second one is based on a permutation procedure.
#'
#'  @return A \code{casanova} object containing the following components:
#'  \item{pvalues_stat}{The p-values obtained by \eqn{\chi^2}-approximation}
#'  \item{pvalues_per}{The p-values of the permutation approach}
#'  \item{statistics}{The value of the casanova along with degrees of freedom of the
#'  central chi-square distribution and p-value, as well as the p-value of the
#'   permutation procedure.}
#'  \item{rg}{A list containg the exponents of the direction considered in the
#'  statistical analysis
#'  \item{cross}{logical. Was the crossing direction considered in the statistical
#'  analysis}
#'  \item{indep}{logical. Were the directions specified by the user linearly
#'  independent?}
#'  \item{nperm}{The number of permutations used for calculating the permuted p-value.
#' @examples
#' library("survival")
#' data(veteran)
#' out <- casanova(formula ="time ~ trt*celltype",event = "status",
#'  data = veteran)
#'
#' ## Detailed informations:
#' summary(out)
#'
#' @references
#'
#'
#' @export
#'

copsanova <- function(formula, event ="event", data = NULL, BSiter = 1999,
                      weights = "pois",
                     nested.levels.unique = FALSE, tau = NULL){
  input_list <- list(formula = formula,time = time, data = data, BSiter = BSiter,
                     weights = weights)
  #Zeit und in Formel einbinden
  formula2 <-  paste0(formula,"*",event)
  dat <- model.frame(formula2, data)
  #n
  subject <- 1:nrow(dat)
  n_all <- length(subject)

  formula <- as.formula(formula)
  nf <- ncol(dat) - 1 - 1
  nadat <- names(dat)

  names(dat) <- c("time",nadat[2:(1+nf)],"event")

  dat2 <- data.frame(dat, subject = subject)

  nadat2 <- nadat[-c(1,nf+2)]


  fl <- NA
  for (aa in 1:nf) {
    fl[aa] <- nlevels(as.factor(dat[, aa + 1]))
  }
  levels <- list()
  for (jj in 1:nf) {
    levels[[jj]] <- levels(as.factor(dat[, jj + 1]))
  }
  lev_names <- expand.grid(levels)
  if (nf == 1) {
    dat2 <- dat2[order(dat2[, 2]), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(subject),
                     .drop = F)$Measure
    hypo_matrices <- list(diag(fl) - matrix(1/fl, ncol = fl, nrow = fl))
    group <- rep(1:length(n),n)
    dat2$group <- group

    ##########

    event <- dat2[,"event"]
    group <- dat2$group
    diff_groups <- length(unique(group))

    dat2$exit <- dat2$time
    dat2$to <- dat2$event
    data1 <- list()

    print(dat2)
    n <- numeric(0)
    tau_all <- numeric(0)
    for( k in 1:diff_groups){
      dat_tmp <- dat2[dat2$group == k, c("exit","to","group")]
      n <- c(n,length(dat_tmp$to))

      tau_all <- c(tau_all,max(dat_tmp[dat_tmp$to ==1,]$exit))
      # ind_tau <- dat_tmp$exit >= tau
      #
      # dat_tmp$to[ind_tau] <- "1"  # Alles was gr??er oder gleich tau ist, wird als unzensiert angesetzt, damit der Kaplan-Meier-Sch?tzer in tau auf Null f?llt.
      # dat_tmp$exit[ind_tau] <- tau
      data1[[k]] <- dat_tmp
    }

    if(is.null(tau)){
      tau <- min(tau_all)
    }

    for( k in 1:diff_groups){
      dat_tmp <- data1[[k]]
      ind_tau <- dat_tmp$exit >= tau

      dat_tmp$to[ind_tau] <- "1"  # Alles was gr??er oder gleich tau ist, wird als unzensiert angesetzt, damit der Kaplan-Meier-Sch?tzer in tau auf Null f?llt.
      dat_tmp$exit[ind_tau] <- tau
      data1[[k]] <- dat_tmp
    }

    copsanova_erg <- test.data(data1, n, BSiter = BSiter,
                               Gewichte = weights, c.matrix = hypo_matrices)


  }
  else {
    lev_names <- lev_names[do.call(order, lev_names[, 1:nf]),
                           ]
    dat2 <- dat2[do.call(order, dat2[, 2:(nf + 1)]), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    fac_names_original <- fac_names
    perm_names <- t(attr(terms(formula), "factors")[-1, ])
    ###

    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(subject),
                     .drop = F)$Measure
    group <- rep(1:length(n),n)
    dat2$group <- group
    if (length(fac_names) != nf && 2 %in% nr_hypo) {
      stop("A model involving both nested and crossed factors is\n           not impemented!")
    }
    if (length(fac_names) == nf && nf >= 4) {
      stop("Four- and higher way nested designs are\n           not implemented!")
    }
    if (length(fac_names) == nf) {
      TYPE <- "nested"
      if (nested.levels.unique) {
        n <- n[n != 0]
        blev <- list()
        lev_names <- list()
        for (ii in 1:length(levels[[1]])) {
          blev[[ii]] <- levels(as.factor(dat[, 3][dat[,
                                                      2] == levels[[1]][ii]]))
          lev_names[[ii]] <- rep(levels[[1]][ii], length(blev[[ii]]))
        }
        if (nf == 2) {
          lev_names <- as.factor(unlist(lev_names))
          blev <- as.factor(unlist(blev))
          lev_names <- cbind.data.frame(lev_names, blev)
        }
        else {
          lev_names <- lapply(lev_names, rep, length(levels[[3]])/length(levels[[2]]))
          lev_names <- lapply(lev_names, sort)
          lev_names <- as.factor(unlist(lev_names))
          blev <- lapply(blev, rep, length(levels[[3]])/length(levels[[2]]))
          blev <- lapply(blev, sort)
          blev <- as.factor(unlist(blev))
          lev_names <- cbind.data.frame(lev_names, blev,
                                        as.factor(levels[[3]]))
        }
        if (nf == 2) {
          fl[2] <- fl[2]/fl[1]
        }
        else if (nf == 3) {
          fl[3] <- fl[3]/fl[2]
          fl[2] <- fl[2]/fl[1]
        }
      }
      hypo_matrices <- GFD:::HN(fl)
    }
    else {
      TYPE <- "crossed"
      hypo_matrices <- GFD:::HC(fl, perm_names, fac_names)[[1]]
      fac_names <- GFD:::HC(fl, perm_names, fac_names)[[2]]
    }
    if (length(fac_names) != length(hypo_matrices)) {
      stop("Something is wrong: Perhaps a missing interaction term in formula?")
    }
    if (TYPE == "nested" & 0 %in% n & nested.levels.unique ==
        FALSE) {
      stop("The levels of the nested factor are probably labeled uniquely,\n           but nested.levels.unique is not set to TRUE.")
    }
    if (0 %in% n || 1 %in% n) {
      stop("There is at least one factor-level combination\n           with less than 2 observations!")
    }

    ###############################
    event <- dat2[,"event"]
    group <- dat2$group
    diff_groups <- length(unique(group))

    dat2$exit <- dat2$time
    dat2$to <- dat2$event
    data1 <- list()


    n <- numeric(0)
    tau_all <- numeric(0)
      for( k in 1:diff_groups){
        dat_tmp <- dat2[dat2$group == k, c("exit","to","group")]
        n <- c(n,length(dat_tmp$to))

        tau_all <- c(tau_all,max(dat_tmp[dat_tmp$to ==1,]$exit))
        # ind_tau <- dat_tmp$exit >= tau
        #
        # dat_tmp$to[ind_tau] <- "1"  # Alles was gr??er oder gleich tau ist, wird als unzensiert angesetzt, damit der Kaplan-Meier-Sch?tzer in tau auf Null f?llt.
        # dat_tmp$exit[ind_tau] <- tau
         data1[[k]] <- dat_tmp
      }

    if(is.null(tau)){
      tau <- min(tau_all)
    }

        for( k in 1:diff_groups){
      dat_tmp <- data1[[k]]
      ind_tau <- dat_tmp$exit >= tau

      dat_tmp$to[ind_tau] <- "1"  # Alles was gr??er oder gleich tau ist, wird als unzensiert angesetzt, damit der Kaplan-Meier-Sch?tzer in tau auf Null f?llt.
      dat_tmp$exit[ind_tau] <- tau
      data1[[k]] <- dat_tmp
    }

      copsanova_erg <- test.data(data1, n, BSiter = BSiter,
                      Gewichte = weights, c.matrix = hypo_matrices)


  }
  output <- list()
  output$input <- input_list
  output$pvalues <-  copsanova_erg$value
  output$test_statistics <-  copsanova_erg$test_statistics
  output$bsiter <- BSiter
  output$weights <- weights

  output$statistic <- cbind(copsanova_erg$test_statistics,round(copsanova_erg$value,4))
  rownames(output$statistic) <- fac_names
  colnames(output$statistic) <- c("Test statistic","p-value")


  class(output) <- "copsanova"
  return(output)


}
# set.seed(1)
# copsanova("exit ~ sex", "to", data = data, BSiter = 9,
#                        weights = "pois", nested.levels.unique = FALSE)



