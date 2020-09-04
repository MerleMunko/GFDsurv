#FIXME: - mehrere Hypothesen parallel -ich
#FIXME: output festlegen- besprechung
#FIXME: tau -optionen -besprechung

copsanova <- function(formula, event ="event", data = NULL, BSiter = 1999,
                      weights = "pois",
                     nested.levels.unique = FALSE){
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
    ###############################
    event <- dat2[,"event"]
    group <- dat2$group
    diff_groups <- length(unique(group))

    dat2$exit <- dat2$time
    dat2$to <- dat2$event
    data1 <- list()
    tau <- 2173


    n <- numeric(0)
    for( k in 1:diff_groups){
      dat_tmp <- dat2[dat2$group == k, c("exit","to","group")]
      n <- c(n,length(dat_tmp$to))

      ind_tau <- dat_tmp$exit >= tau

      dat_tmp$to[ind_tau] <- "1"  # Alles was gr??er oder gleich tau ist, wird als unzensiert angesetzt, damit der Kaplan-Meier-Sch?tzer in tau auf Null f?llt.
      dat_tmp$exit[ind_tau] <- tau
      data1[[k]] <- dat_tmp
    }
    tau <- max(data1[[i]][data1[[i]]$to ==1,]$exit)

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
    tau <- 2173


    n <- numeric(0)
      for( k in 1:diff_groups){
        dat_tmp <- dat2[dat2$group == k, c("exit","to","group")]
        n <- c(n,length(dat_tmp$to))

        ind_tau <- dat_tmp$exit >= tau

        dat_tmp$to[ind_tau] <- "1"  # Alles was gr??er oder gleich tau ist, wird als unzensiert angesetzt, damit der Kaplan-Meier-Sch?tzer in tau auf Null f?llt.
        dat_tmp$exit[ind_tau] <- tau
        data1[[k]] <- dat_tmp
      }
    tau <- max(data1[[i]][data1[[i]]$to ==1,]$exit)

      copsanova_erg <- test.data(data1, n, BSiter = BSiter,
                      Gewichte = weights, c.matrix = hypo_matrices)


  }
  output <- list()
  output$input <- input_list
  output$pvalues <-  copsanova_erg$value
  output$bsiter <- BSiter
  output$weights <- weights

  output$statistic <- cbind(copsanova_erg$value,round(copsanova_erg$value,4))
  rownames(output$statistic) <- fac_names
  colnames(output$statistic) <- c("Test statistic","p-value")


  class(output) <- "copsanova"
  return(output)


}
set.seed(1)
datatest1  <- copsanova("exit ~ sex*treat", "to", data = data, BSiter = 9,
                       weights = "pois", nested.levels.unique = FALSE)



copsanova("exit ~ sex*treat", "to", data = data)
data
sort(datatest1[[2]][[6]]$exit)
sort(data1[[1]]$exit)
data[[1]][]

datatset <- data1[[1]]
datatset$exit[-c(2,3,4,5,6)] <- 1222
datatset$exit
1:141 %in% c(2,3,4,5,6)
datatset$to[1:141 %in% c(2,3,4,5,6)] <- "test"
datatset
