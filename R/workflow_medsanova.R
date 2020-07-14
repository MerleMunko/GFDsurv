medSANOVA <- function(formula, event ="event", data = NULL, nperm = 1999, alpha = 0.05,
                      variant= "twosided",
                      nested.levels.unique = FALSE)

#  formel "Var~ FactorA * FactorB"
#  variant c("onesided,twosided")
# variante 2 twosided
# variante 3 onesided
 
  Call:
  [1] "time ~ trt*celltype"

medSANOVA: Cumulative Aalen survival analyis-of-variance:

             Test statistic df p-value p-value perm
trt                6.778400  2  0.0337          0.0
celltype          21.472560  6  0.0015          0.0
trt:celltype       7.439236  6  0.2821          0.5
