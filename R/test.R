library(condSURV)
data0 <- colonCS

# Achtung bei der codierung f?r die Zensierung. Codierung erfolgt als "1" f?r unzensiert und "cenS" f?r zensiert
data <- data.frame( "exit" = data0$Stime, "to" = ifelse(data0$event == 1, "1", "cens"), "sex" = data0$sex, "treat" = data0$rx)
data1 <- list()
n <- numeric(0)
j <- 1
tau <- 2173

for( sex in unique(data0$sex)){
  for( treat in unique(data0$rx)){
    ind <- (data$sex == sex) & (data$treat == treat)
    n <- c(n,sum(ind))
    data2 <- data[ (data$sex == sex) & (data$treat == treat), ]
    # ind_tau <- data2$exit >= tau
    # data2$to[ind_tau] <- "1"  # Alles was gr??er oder gleich tau ist, wird als unzensiert angesetzt, damit der Kaplan-Meier-Sch?tzer in tau auf Null f?llt.
    # data2$exit[ind_tau] <- tau
    data1[[j]] <- data2
    j <- j + 1
  }
}

for(i in 1:6){
  print(max(data1[[i]][data1[[i]]$to ==1,]$exit))
}
tau

set.seed(1)
system.time( a <-test.data(data1, n, BSiter = 100, alpha = 0.05, Gewichte = "pois", c.matrix = null_mat_A(2,3)))

system.time( a1 <-test.data(data1, n, BSiter = 1000, alpha = 0.05, Gewichte = "pois", c.matrix = null_mat_B(2,3) ))

system.time( a2 <-test.data(data1, n, BSiter = 1000, alpha = 0.05, Gewichte = "pois", c.matrix = null_mat_AB(2,3)))

system.time( a3 <-test.data(data1, n, BSiter = 1000, alpha = 0.05, Gewichte = "pois", c.matrix = null_mat_x1(2,3)))


> a
[1] 0.5984016
> a1
[1] 0.001998002
> a2
[1] 0.002997003
> a3
[1] 0.000999001
