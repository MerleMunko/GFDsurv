
data_examp <- function(time, cens, A, B, iter, variant){
  a <- length(unique(A))
  b <- length(unique(B))
  k <- a*b

  group <- rep(0, length(time))
  i <- 0
  for(i_a in unique(A)){
    for(i_b in unique(B)){
      i <- i + 1
      ind <- (A == i_a) & (B == i_b)
      group[ind] <- i
    }
  }
  values <- matrix( c(time, cens, group), ncol = 3, byrow = FALSE)


  ################################################################################
  ###################### Parameter################################################
  ################################################################################
  source("functions_inerle_3.R")



    design = list("A", "B", "A:B")
  hyp_list <- unlist(design)

  design[which(unlist(design) == "A:B")] <- "AB"

  chi_quant <- list()
  t_mat_list <- list()
  for( i_des in 1:length(design)){
    hyp <- hyp_list[i_des]
    t_mat <- do.call(paste0("null_mat_", design[[i_des]]), list(a,b))
    t_mat <- t(t_mat) %*% MASS::ginv( t_mat %*% t(t_mat) ) %*% t_mat
    t_mat_list[[hyp]] <- t_mat
    chi_quant[[hyp]] <- qchisq(p = 0.95, df = qr(t_mat)$rank)
  }

  #settings <- data.frame( Distr = "data_exp", n = paste0(group_s, collapse = "|"), cens = paste0(censp_s, collapse = "|"))

  ################################################################################
  ###################### Hilfsparameter ##########################################
  ################################################################################
 # ngroup <- length(censp_s)



  ################################################################################
  ############# Funktionen #######################################################
  ################################################################################

  #Simulationen
  out <- data_test(values = sort_data(values), n_perm = iter, t_mat_list = t_mat_list, hyp_list = hyp_list, chi_quant = chi_quant, group_s = group_s, a = a, b = b, variant = variant)


  hyp1 <- hyp_list[[1]]
  out_names <- names(out[[hyp1]])
  result <- matrix(1, ncol = length(out[[hyp1]]) , nrow = length(hyp_list))
  colnames(result) <- c(out_names)
  rownames(result) <- unlist(hyp_list)

  i <- 0
  for( hyp in hyp_list){
    i <- i+1
    result[hyp,] <- c(round(100* out[[hyp]],digits=3))
  }

  out <- list(#settings = settings,
              variant = variant, result = result)
  return(out)
}




##################################################################################################
############ Beispiel ##########################################################################
##################################################################################################

library(timereg)

data(csl)

csl2 <- csl[findInterval(unique(csl$id), csl$id), ]
attach(csl2)


## 60 - 69, m
### NUR EFFECT F?R PROP! PASST ABER ZUM BILD
n_perm <- 1 # Permutations- und Bootstrapiterationen
set.seed(1490)
data <- csl2[ (csl2$age >= 0) & (csl2$age < 10) & (csl2$sex == 1),]
data$eventT <- data$eventT +  10^(-4)*runif(length(data[,1]),0,1)
data$prot_groups <- (data$prot.base >= 70)
out_6069_m <- data_examp( data$eventT, data$dc, A= data$treat, B = data$prot_groups, iter = n_perm, variant = 2)
str(data)
out_6069_m
# Ergebnisse f?r variant = 2 (zweiseitiger Intervalsch?tzer)
#$settings
#Distr              n                cens
#1 data_exp 117|103|128|98 28.2|46.6|33.6|53.1
 
#$variant
#[1] 2

#$result
#perm chi.int_A
#A   68.068    65.735
#B    1.802     1.424
#A:B 75.576    74.289


# Ergebnisse f?r variant = 3 (einseitiger Intervalsch?tzer)
#$settings
#Distr              n                cens
#1 data_exp 117|103|128|98 28.2|46.6|33.6|53.1

#$variant
#[1] 3

#$result
#perm chi.int_A
#A   62.362    62.386
#B    1.502     0.673
#A:B 71.972    71.687
