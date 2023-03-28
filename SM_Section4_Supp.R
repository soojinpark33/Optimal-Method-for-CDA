# Load package and source file
library(parallel)
source("SM_source.R")

#-----------------------------------------------------------------------------------------------#
# Simulation Study (Section 4 & Supplementary Materials)
#-----------------------------------------------------------------------------------------------#
#' The followings are regression coefficients extracted from the MIDUS data,
#' using (generalized) linear regression, as described in the manuscript:
#' The MIDUS data is not publicly available due to confidentiality concerns,
#' thus the values of regression coefficients are directly provided below.

## For continuous mediator
# intermediate confounder model (lm): X ~ R + C
coef.X <- c(0.4557002, 0.3054280, 0.1454241)
# mediator model (lm): M ~ R + X + C
coef.M <- c(1.25311025, 0.41166946, 0.13201494, -0.07454179)
# outcome model (lm): Y ~ R + X + C + M + R*M
coef.Y <- c(8.7614110, -1.2532416, -0.1691439, -0.4973777, -0.4071892, 0.3876383)

## For binary mediator
# dichotomized mediator model (glm): dichotomized_M ~ R + X + C
coef.M.bin <- c(-0.4042500, 0.5338881, 0.3442321, -0.2116475)
# outcome model (lm): Y ~ R + X + C + dichotomized_M + R*dichotomized_M
coef.Y.bin <- c(8.3331832, -0.9281828, -0.1856586, -0.4859839, -0.2197366, 0.3545898)

### Generate population data for each setting (estimators 2 through 6) -------------------------#
# The following code reproduces Table 2 of the manuscript.
N <- 1e+06  # population size

## For continuous mediator
con.ratio <- c(0.3, 0.5, 1, 2, 3)
con.coefs <- matrix(c(0.244, -1.048, -1.322,  # for ratio = 0.3
                      0.326, -1.048, -1.112,  # for ratio = 0.5
                      0.477, -1.048, -0.900,  # for ratio = 1
                      0.689, -1.049, -0.750,  # for ratio = 2
                      0.852, -1.049, -0.684), # for ratio = 3
                    ncol = 3, byrow = TRUE)
for(i in 1:5){
  # con0.3, con0.5, con1, con2, con3 created.
  assign(paste("con", con.ratio[i], sep = ""),
         popDataGen(N = N, adjustedCoefs = con.coefs[i, ], binary.mediator = FALSE, seed = 10,
                    coef.X = coef.X, coef.M = coef.M, coef.Y = coef.Y))
}

## For binary mediator
bin.ratio <- c(0.1, 0.3, 0.5, 0.7, 0.9)
bin.coefs <- matrix(c(0.552, -0.669, -1.932,  # for ratio = 0.1
                      1.032, -0.644, -1.244,  # for ratio = 0.3
                      1.463, -0.674, -1.063,  # for ratio = 0.5
                      1.906, -0.693, -0.964,  # for ratio = 0.7
                      2.466, -0.709, -0.900), # for ratio = 0.9
                    ncol = 3, byrow = TRUE)
for(i in 1:5){
  # bin0.1, bin0.3, bin0.5, bin0.7, bin0.9 created.
  assign(paste("bin", bin.ratio[i], sep = ""),
         popDataGen(N = N, adjustedCoefs = bin.coefs[i, ], binary.mediator = TRUE, seed = 10,
                    coef.X = coef.X, coef.M = coef.M.bin, coef.Y = coef.Y.bin))
}

# For Estimator 1, the coefficient of R*M interaction is set to 0.
coef.Y[6] <- 0
coef.Y.bin[6] <- 0

### Generate population data for each setting (estimator 1) ------------------------------------#
# The following code reproduces Table 1 of Supplementary Materials.
## For continuous mediator
con.coefs <- matrix(c(0.244, -0.558, -0.931,  # for ratio = 0.3
                      0.325, -0.558, -0.722,  # for ratio = 0.5
                      0.475, -0.558, -0.511,  # for ratio = 1
                      0.686, -0.558, -0.361,  # for ratio = 2
                      0.850, -0.558, -0.295), # for ratio = 3
                    ncol = 3, byrow = TRUE)
for(i in 1:5){
  # con0.3.est1, con0.5.est1, con1.est1, con2.est1, con3.est1 created.
  assign(paste("con", con.ratio[i], ".est1", sep = ""),
         popDataGen(N = N, adjustedCoefs = con.coefs[i, ], binary.mediator = FALSE, seed = 10,
                    coef.X = coef.X, coef.M = coef.M, coef.Y = coef.Y))
}
## For binary mediator
bin.coefs <- matrix(c(0.628, -0.669, -1.762,  # for ratio = 0.1
                      1.189, -0.644, -1.001,  # for ratio = 0.3
                      1.713, -0.674, -0.791,  # for ratio = 0.5
                      2.275, -0.693, -0.673,  # for ratio = 0.7
                      3.184, -0.709, -0.604), # for ratio = 0.9
                    ncol = 3, byrow = TRUE)
for(i in 1:5){
  # bin0.1.est1, bin0.3.est1, bin0.5.est1, bin0.7.est1, bin0.9.est1 created.
  assign(paste("bin", bin.ratio[i], ".est1", sep = ""),
         popDataGen(N = N, adjustedCoefs = bin.coefs[i, ], binary.mediator = TRUE, seed = 10,
                    coef.X = coef.X, coef.M = coef.M.bin, coef.Y = coef.Y.bin))
}

#' The following code reproduces Table 2 of Supplementary Materials.
#' It can be adjusted for any other setting.
pop.dat <- con0.3$data  # population data generated
round(cov(cbind(pop.dat$Y, pop.dat$M, as.numeric(pop.dat$R),
                pop.dat$X, as.numeric(pop.dat$C.grp))), 3)
pop.dat <- bin0.1$data  # population data generated
round(cov(cbind(pop.dat$Y.bin, as.numeric(pop.dat$M.bin), as.numeric(pop.dat$R),
                pop.dat$X, as.numeric(pop.dat$C.grp))), 3)

#' The following code reproduces Table 3 of Supplementary Materials.
#' It can be adjusted for any other setting.
pop.dat <- con0.3.est1$data  # population data generated
round(cov(cbind(pop.dat$Y, pop.dat$M, as.numeric(pop.dat$R),
                pop.dat$X, as.numeric(pop.dat$C.grp))), 3)
pop.dat <- bin0.1.est1$data  # population data generated
round(cov(cbind(pop.dat$Y.bin, as.numeric(pop.dat$M.bin), as.numeric(pop.dat$R),
                pop.dat$X, as.numeric(pop.dat$C.grp))), 3)

### Simulation ---------------------------------------------------------------------------------#
n.sim <- 1000     # no. of simulations
n.boot <- 1000    # no. of bootstrap replicates
n.cores <- 8      # no. of cores to use; 1 for Windows
conf.level <- .95 # confidence level

#' The following code reproduces Figures 1 and 2 of the manuscript (plots),
#' Tables 4 and 5 of Supplementary Materials (numerical results), and also
#' Figure 2 of Supplementary Materials (standardized data).
#' For example, we focus on the case with continuous mediator, ratio r = 0.3, and 
#' sample size n = 1000. Then the following code can be adjusted for any other
#' setting accordingly.
pop.dat <- con0.3$data                # population data generated (estimators 2 through 6)
Red.true <- con0.3$Red.true           # true value of disparity reduction
Rem.true <- con0.3$Rem.true           # true value of disparity remaining
pop.dat.est1 <- con0.3.est1$data      # population data generated (estimator 1)
Red.true.est1 <- con0.3.est1$Red.true # true value of disparity reduction
Rem.true.est1 <- con0.3.est1$Rem.true # true value of disparity remaining
n <- 1000                             # sample size
binary.mediator <- FALSE              # mediator type
std <- FALSE                          # whether or not to standardize data

# print results
system.time(summary0 <- mclapply(1:n.sim, sim.one, pop.dat = pop.dat, Red.true = Red.true, Rem.true = Rem.true,
                     pop.dat.est1 = pop.dat.est1, Red.true.est1 = Red.true.est1, Rem.true.est1 = Rem.true.est1,
                     n = n, n.boot = n.boot, conf.level = conf.level,
                     binary.mediator = binary.mediator, standardize = std, mc.cores = n.cores))
out.full <- simplify2array(summary0)
con.n500.r0.3 <- summarize.results(out.mclapply = out.full, n.sim = n.sim,
                                    Red.true = Red.true, Rem.true = Rem.true,
                                    Red.true.est1 = Red.true.est1, Rem.true.est1 = Rem.true.est1,
                                    binary.mediator = binary.mediator)
cbind(con.n1000.r0.3$disparity.reduction[, -c(3, 4)],
      con.n1000.r0.3$disparity.remaining[, -c(3, 4)])
