#' medSANOVA: Mediann survival analyis-of-variance
#'
#' The function \code{medSANOVA} calculates the p-values of the multiple-direction logrank test based
#' on the \eqn{\chi^2}-approximation and the permutation approach.
#' @param formula A model \code{formula} object. The left hand side contains the time variable and the right
#'  hand side contains the factor variables of interest. An interaction term must be
#'  specified.
#' @param event The name of censoring status indicator with values 0=censored and
#' 1=uncensored.
#' The default choice is "event"
#' @param data A data.frame, list or environment containing the variables in formula
#' and the censoring status
#' indicator. Default option is \code{NULL}.
#' @param nperm The number of permutations used for calculating the permuted p-value.
#'   The default option is 1999.
#' @param alpha A number specifying the significance level; the default is 0.05.
#' @param nested.levels.unique A logical specifying whether the levels of the nested
#' factor(s) are labeled uniquely or not.
#'  Default is FALSE, i.e., the levels of the nested factor are the same for each
#'  level of the main factor.
#'  @param  var_est Kannst hier auch den Namen ändern
#' @details
#' The \code{casanova} function calculates the Wald-type statistic of weighted
#' Nelson-Aalen type integrals
#' for general factorial survival designs. The approach allows the combination of
#' different weights into a
#' joint statistic. The user can choose between weights of the following form:
#' w(x) = 1 - 2x (\code{cross = TRUE}) and w(x) = x^r * (1-x)^g for natural numbers
#' r,g (including 0). The function automatically check whether the specified weights
#' fulfill
#' the linear independence assumption and choose a subset of linearly independent
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
#' @references Ditzhaus, M., Janssen, A. and Pauly, M. (2020). Permutation inference in factorial survival designs with the
#'           CASANOVA. arXiv preprint (arXiv:2004.10818).
#'
