#Sort object of data_gen
#Input:
# values:   matrix - created by data_gen

sort_data <- function(values) values[order(values[, 1]), ]

#create Matrix C, containing the groups
#Input:
# values:   matrix - created by data_gen, possibly sorted by sort_data
#Output:
# Matrix in which a column represents a group and a row an individual
create_c_mat <- function(values) {
  group <- values[, 3]
  c_mat <- matrix(0, nrow = max(group), ncol = nrow(values))

  c_mat[(0:(nrow(values) - 1)) * 3 + group] <- 1
  return(t(c_mat))
}

#Example:
#create_c_mat(sort_data(test_dat))




#Kaplan-Meier estimator and median
#Input:
# values:   matrix created by data_gen and sorted by sort_data
# group:    integer vector containing the group of the observations. Default is
#           the third column of the values, the groups drawn by data_gen

#Output:
# matrix values with additional column containing the KME and the median
# for the corresponding group
KME <- function(values, group = values[, 3]) {
  #Add columns
  values <- cbind(values, 0, 0)
  colnames(values)[c(4, 5)] <- c("KME_VF", "group_med")

  for (i in 1:max(group)) {
    #calculate in group
    values2 <- values[group == i, ]
    n <- nrow(values2)
    vec_temp <- rep(1, n)

    #calculate the KME
    res <- 1 - cumprod(1 - values2[, 2] / (n + 1 - 1:n))

    #find median
    med_ind <- findInterval(0.5-1e-8, res) + 1
    if(med_ind > n) {
      #if median lies not in observations
      med <- Inf
    } else {
      med <- values2[med_ind, 1]
    }

    #add values to matri
    values[group == i, 4] <- res
    values[group == i, 5] <- rep(med, n)
  }
  return(values)
}

#Example
#test_KME <- KME(sort_data(test_dat))


#estimation of standard deviation by formula derived from confidence interval
#two-sided

#Help function to calculate the standard deviation out of the interval. Can
#only be used for one group
#Input:
# values:   - matrix created by dat_gen, sorted by sort_data and with added
#             columns by KME. Contains only one group
# n_all     - integer, total number of observations
# a         - numeric value between 0 and 1, alpha level, default 0.1
#Output:    Standard deviation estimated with the formula derived from interval

int_sd2 <- function(values, a = 0.1, n_all) {
  n <- nrow(values)
  indi <- values[, 1] <= values[, 5]

  Vi <- sum(indi * values[, 2]/(n:1)^2)
  z <- qnorm(1 - a/2)

  li <- max(0, 1/2*(1 - z * sqrt(Vi)))
  #lower quantile
  ui <- min(1, 1/2*(1 + z * sqrt(Vi)))
  #upper quantile

  if(ui > max(values[, 4])) {
    ui <- max(values[, 4])
    a <- 2 * dnorm((2*ui - 1)/sqrt(Vi/n))
    z <- qnorm(1-a/2)
    li <- max(0, 1/2*(1 - z * sqrt(Vi)))

    if(is.infinite(values[1, 5])){
      return(c("sigma" = NA))
    }
  }
  #upper quantile does not work

  fui_ind <- findInterval(ui, values[, 4])
  #index upper quantile, if the upper Interval is not in n, it is the last value
  fli_ind <- findInterval(li, values[, 4])
  #index lower quantile

  sig <- 1/(2*z) * sqrt(n_all) *
    (values[min(fui_ind+1, n), 1] - values[fli_ind+1, 1])

  return(c("sigma" = unname(sig)))
}





#estimation of standard deviation by formula derived from confidence interval
#one-sided

#Help function to calculate the standard deviation out of the interval. Can
#only be used for one group
#Input:
# values:   - matrix created by dat_gen, sorted by sort_data and with added
#             columns by KME. Contains only one group
# n_all     - integer, total number of observations
# a         - numeric value between 0 and 1, alpha level, default 0.1
#Output:    Standard deviation estimated with the formula derived from interval

int_sd3 <- function(values, a = 0.1, n_all) {
  n <- nrow(values)
  indi <- values[, 1] <= values[, 5]

  Vi <- sum(indi * values[, 2]/(n:1)^2)
  z <- qnorm(1 - a/2)

  li <- max(0, 1/2*(1 - z * sqrt(Vi)))
  #lower quantile

  if(is.infinite(values[1, 5])){
    return(c("sigma" = NA))
  }

  fli_ind <- findInterval(li, values[, 4])
  #index lower quantile

  sig <- 1/z * sqrt(n_all) *
    (values[1, 5] - values[fli_ind+1, 1])

  return(c("sigma" = unname(sig)))
}



#wrapper for int_sd several groups in data
#Input:
# values    - matrix created by dat_gen, sorted by sort_data and with added
#             columns by KME. Should contain several groups
# alpha     - numeric value between 0 and 1, alpha level
# group:    - integer vector containing the group of the observations. Default
#             is the third column of the values, the groups drawn by data_gen
#Output:    Standard deviation estimated with interval formula for each group

int_var_groups <- function(values, alpha = 0.1, group = values[, 3],
                           variant) {
  n_all <- nrow(values)
  n_group <- max(group)
  erg <- numeric(n_group)
  med <- numeric(n_group)

  for(i in 1:n_group) {
    values2 <- values[group == i,]

    #calculate for each group

    temp <- do.call(paste0("int_sd",variant),
                    args = list(values = values2, a = alpha, n_all = n_all))
    erg[i] <- temp
    med[i] <- values2[1, 5]

  }

  names(erg) <- 1:n_group
  names(med) <- 1:n_group
  return(list("Median" = med, "Variance" = erg^2))
  #notice that ^2 changes from sd to var...
}

#int_var_groups(test_KME)




#Null hypothesis matrix

#classical null hypothesis
#Input:
# n_group:   - integer, number of groups
#Output:
# Returns the matrix testing the Nullhypothesis for the respective number of
# groups
null_mat_x1 <- function(n_group_a, n_group_b){
  n_group  <- n_group_a * n_group_b
  return(diag(n_group) - matrix(1, n_group, n_group)/n_group)
}
#null_mat_clas(5)



#No main effect of factor a
#Input:
# n_group_a:  - integer, number of levels in factor a
# n_group_b:  - integer, number of levels in factor b
#Output:
# Returns the matrix testing the Nullhypothesis
null_mat_A <- function(n_group_a, n_group_b) {
  pa <- diag(n_group_a) - matrix(1, n_group_a, n_group_a)/n_group_a
  jbb <- matrix(1, n_group_b, n_group_b)/n_group_b
  ha <- pa %x% jbb

  return(ha)
}

#No main effect of factor B
#Input:
# n_group_a:  - integer, number of levels in factor a
# n_group_b:  - integer, number of levels in factor b
#Output:
# Returns the matrix testing the Nullhypothesis
null_mat_B <- function(n_group_a, n_group_b) {
  pb <- diag(n_group_b) - matrix(1, n_group_b, n_group_b)/n_group_b
  jaa <- matrix(1, n_group_a, n_group_a)/n_group_a
  ha <- jaa %x% pb

  return(ha)
}


#null_mat_fac(2, 2)

#No interaction effect
#Input:
# n_group_a:  - integer, number of levels in factor a
# n_group_b:  - integer, number of levels in factor b
#Output:
# Returns the matrix testing the Nullhypothesis
null_mat_AB <- function(n_group_a, n_group_b){
  pa <- diag(n_group_a) - matrix(1, n_group_a, n_group_a)/n_group_a
  pb <- diag(n_group_b) - matrix(1, n_group_b, n_group_b)/n_group_b
  hab <- pa %x% pb

  return(hab)
}

#null_mat_int(2, 2)

#projector
#Input:
# h       - Matrix
#Output: projected matrix

project <- function(h) {
  return(t(h) %*% MASS:::ginv(h %*% t(h)) %*% h)
}


#Teststatistic
#Input:
# n         - integer, total number of observations
# t_mat     - matrix, null hypothesis matrix
# var_out   - Output from one of the variance functions boot_var_groups or
#             int_var_groups, containing median and variance
#
#Output: Teststatistic

test_stat <- function(values, t_mat, var_out) {
  n <- nrow(values)
  med_vec <- var_out$Median
  sig_vec <- var_out$Variance
  if(any(is.infinite(med_vec))) return(NA)
  if(any(is.na(sig_vec))) return(NA)

  return(n * t(t_mat %*% med_vec) %*%
           MASS::ginv(t_mat %*% diag(sig_vec) %*% t(t_mat)) %*%
           (t_mat %*% med_vec))
}


#wrapper to get the teststatistics directly from the sorted data
#Input
# values:   matrix. The data to start with
# group:    integer vector containing the group of the observations. Default is
#           the third column of the values, the groups drawn by data_gen
wrap_sim2 <- function(values, group = values[, 3], t_mat_list, variant){

  C_mat <- function(x){
    t(x) %*% MASS::ginv( x %*% t(x) ) %*% x
  }

  t_mat_list <- lapply(t_mat_list, C_mat)

  values <- sort_data(values)
  values[, 3] <- group
  print(values)
  values_KME <- KME(values, group = values[, 3])

  var_int <- int_var_groups(values_KME, group = values[ ,3], variant = variant)
  erg_int <- list()

  for( i in length(t_mat_list)){
    erg_int[[i]] <- test_stat(values, t_mat_list[[i]], var_out = var_int)
  }

  print(erg_int)
  out <- c( unlist(erg_int))
  return(out)
}


#wrap_sim(values = test_dat,
#         t_mat = null_mat_clas(3))


#Permutations
#Input:
# values    - Matrix, Data to be entered in Simulation
# n_perm    - Integer, Number of permutations
# t_mat     - matrix, test matrix
#Output
# Vector with the quantiles of the test-statistics and number of permutations
#

perm_fun <- function(values, n_perm, t_mat_list = t_mat_list, hyp_list = hyp_list, variant ) {
  values2 <- sort_data(values)
  group_org <- values2[, 3]
  group_new <- replicate(n_perm, sample(group_org))
  test_stat_erg <- apply(group_new, 2,
                         function(x) wrap_sim(values = values2, group = x,
                                              t_mat_list = t_mat_list, hyp_list = hyp_list, variant = variant) )
  q <- list()
  for( hyp in hyp_list){
    q <- unname(quantile(test_stat_erg[paste0("int_", hyp), ], 0.95, na.rm = TRUE))
    q[[hyp]] <- c(q)
    names(q[[hyp]]) <- paste0("q_int")
  }
  return(list(q = q, test_stat_erg = test_stat_erg ) )
}

#perm_fun(test_dat, 1000, null_mat_clas(3))



#wrapper fuer data example
data_test <- function(values, n_perm, t_mat_list, hyp_list, chi_quant, group_s, a, b, variant ) {
  values <- sort_data(values)
  erg_stat <- wrap_sim(values, t_mat_list = t_mat_list, hyp_list = hyp_list, variant = variant)

  out <- list()

  erg_perm <- perm_fun(values, n_perm, t_mat_list = t_mat_list, hyp_list = hyp_list, variant = variant)
  for(hyp in hyp_list){
    q_perm <- erg_perm$test_stat_erg
    t_int_perm <- mean(erg_stat[paste0("int_", hyp)] <= q_perm[paste0("int_", hyp), ], na.rm = TRUE)
    t_int_chi <- 1-pchisq(erg_stat[paste0("int_", hyp)], df = qr(t_mat_list[[hyp]])$rank )

    t_int_perm <- ifelse(is.nan(t_int_perm), NA, t_int_perm)

    out1 <- c("perm" = t_int_perm, "chi" = t_int_chi)
    out[[hyp]] <- out1
  }
  return(out)
}

