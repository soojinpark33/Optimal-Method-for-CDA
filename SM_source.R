# Load package and source file
library(truncnorm)
library(parameters)
library(effectsize)

#-----------------------------------------------------------------------------------------------#
# Simulation Study (Section 4 & Supplementary Materials)
#-----------------------------------------------------------------------------------------------#
# Population data generation
popDataGen <- function(N = 1000000, adjustedCoefs, binary.mediator, seed = NULL,
                       coef.X, coef.M, coef.Y){
  
  require(truncnorm)
  require(parameters)
  require(effectsize)
  
  # Set parameters
  p0 <- 0.5       # target prop. of White men (R = 0)
  lb.C <- 25      # lower bound for C
  ub.C <- 75      # upper bound for C
  avg.C <- 50     # truncated normal mean for C
  sd.C <- 12      # truncated normal sd for C
  shift.C <- 2    # difference in truncated normal means between two groups (R = 0, 1) for C
  mean.X1 <- 0.5  # poisson mean of X for group R = 0
  mean.X2 <- 0.8  # poisson mean of X for group R = 1
  err.fitM <- 1   # for error term of M model
  err.fitX <- 1   # for error term of X models
  err.fitY <- 1   # for error term of Y model
  
  # Create variables for population data
  set.seed(seed)
  R <- as.numeric(rbinom(N, size = 1, prob = p0))
  C <- rep(NA, N)
  C[R == 0] <- truncnorm::rtruncnorm(sum(R == 0), a = lb.C, b = ub.C, mean = avg.C, sd = sd.C)
  C[R == 1] <- truncnorm::rtruncnorm(sum(R == 1), a = lb.C, b = ub.C, mean = avg.C - shift.C, sd = sd.C)
  err.M <- rnorm(N, sd = err.fitM)
  err.X <- rnorm(N, sd = err.fitX)
  err.Y <- rnorm(N, sd = err.fitY)
  C.grp <- cut(C, breaks = c(0, 50, Inf), include.lowest = TRUE, labels = c(1, 2))
  C.grp.ref2 <- relevel(C.grp, ref = 2)	# relevel
  X <- coef.X[1] + coef.X[2]*R + coef.X[3]*ifelse(C.grp == 2, 1, 0) + err.X  
  
  if(binary.mediator){
    
    temp <- coef.M[1] + coef.M[2]*R + coef.M[3]*X + coef.M[4]*ifelse(C.grp == 2, 1, 0) 
    m.prob <- exp(temp)/(1 + exp(temp))
    M.bin <- rbinom(N, 1, m.prob)
    Y.bin <- coef.Y[1] + coef.Y[2]*R + coef.Y[3]*X + coef.Y[4]*ifelse(C.grp == 2, 1, 0) +
      coef.Y[5]*M.bin + coef.Y[6]*R*M.bin + err.Y
    pop.dat <- data.frame(Y.bin, R, X, M.bin, C, C.grp, C.grp.ref2)
    names(pop.dat) <- c("Y.bin", "R", "X", "M.bin", "C", "C.grp", "C.grp.ref2")
    
  } else {
    
    M <- coef.M[1] + coef.M[2]*R + coef.M[3]*X + coef.M[4]*ifelse(C.grp == 2, 1, 0) + err.M
    Y <- coef.Y[1] + coef.Y[2]*R + coef.Y[3]*X + coef.Y[4]*ifelse(C.grp == 2, 1, 0) + coef.Y[5]*M + coef.Y[6]*R*M + err.Y
    pop.dat <- data.frame(Y, R, X, M, C, C.grp, C.grp.ref2)
    names(pop.dat) <- c("Y", "R", "X", "M", "C", "C.grp", "C.grp.ref2")
    
  }
  
  if(binary.mediator){
    
    # Adjust regression coefficients
    coef.M[2] <- adjustedCoefs[1]	# coefficient of R in M model
    coef.Y[2] <- adjustedCoefs[2]	# coefficient of R in Y model
    coef.Y[5] <- adjustedCoefs[3]	# coefficient of M in Y model
    
    # Adjust data accordingly
    temp <- coef.M[1] + coef.M[2]*R + coef.M[3]*X + coef.M[4]*ifelse(C.grp == 2, 1, 0) 
    m.prob <- exp(temp)/(1 + exp(temp))
    M.bin <- rbinom(N, 1, m.prob)
    Y.bin <- coef.Y[1] + coef.Y[2]*R + coef.Y[3]*X + coef.Y[4]*ifelse(C.grp == 2, 1, 0) +
      coef.Y[5]*M.bin + coef.Y[6]*R*M.bin + err.Y
    
    pop.dat$Y.bin <- Y.bin
    pop.dat$M.bin <- M.bin
    
    # Check
    mod1 <- lm(M.bin ~ R + C, data = pop.dat)
    params1 <- parameters::model_parameters(mod1)
    tmp1 <- effectsize::t_to_r(params1$t[- 1], df_error = params1$df_error[- 1])
    RM.effect <- tmp1$r[1]
    mod2 <- lm(Y.bin ~ R + X + M.bin + C + R * M.bin, data = pop.dat)
    params2 <- parameters::model_parameters(mod2)
    tmp2 <- effectsize::t_to_r(params2$t[- 1], df_error = params2$df_error[- 1])
    MY.effect <- tmp2$r[3]
    
    # True values of the disparity reduction and disparity remaining
    set.seed(seed)
    y0 <- mean(pop.dat$Y.bin[pop.dat$R == 0 & pop.dat$C.grp == 2])
    y1 <- mean(pop.dat$Y.bin[pop.dat$R == 1 & pop.dat$C.grp == 2])
    m0 <- pop.dat$M.bin[pop.dat$R == 0 & pop.dat$C.grp == 2]
    m1 <- pop.dat$M.bin[pop.dat$R == 1 & pop.dat$C.grp == 2]
    m01 <- sample(m0, (sum(pop.dat$R == 1 & pop.dat$C.grp == 2)), replace = TRUE)
    ym0 <- mean(coef.Y[1] + coef.Y[2]*1 + coef.Y[3]*pop.dat$X[pop.dat$R == 1 & pop.dat$C.grp == 2] +
                  coef.Y[4]*1 + coef.Y[5]*ifelse(m01 == 1, 1, 0) + coef.Y[6]*1*ifelse(m01 == 1, 1, 0) +
                  err.Y[pop.dat$R == 1 & pop.dat$C.grp == 2])
    
    print.ratio <- (mean(as.numeric(m1)) - mean(as.numeric(m0)))/(coef.Y[5] + coef.Y[6])
    pop.dat$R <- as.factor(pop.dat$R)
    pop.dat$M.bin <- as.factor(pop.dat$M.bin)
    
  } else if(!binary.mediator) {
    
    # Adjust regression coefficients
    coef.M[2] <- adjustedCoefs[1]	# coefficient of R in M model
    coef.Y[2] <- adjustedCoefs[2]	# coefficient of R in Y model
    coef.Y[5] <- adjustedCoefs[3]	# coefficient of M in Y model
    
    # Adjust data accordingly
    M <- coef.M[1] + coef.M[2]*R + coef.M[3]*X + coef.M[4]*ifelse(C.grp == 2, 1, 0) + err.M
    Y <- coef.Y[1] + coef.Y[2]*R + coef.Y[3]*X + coef.Y[4]*ifelse(C.grp == 2, 1, 0) +
      coef.Y[5]*M + coef.Y[6]*R*M + err.Y
    
    pop.dat$M <- M
    pop.dat$Y <- Y
    
    # Check
    mod1 <- lm(M ~ R + C, data = pop.dat)
    params1 <- parameters::model_parameters(mod1)
    tmp1 <- effectsize::t_to_r(params1$t[- 1], df_error = params1$df_error[- 1])
    RM.effect <- tmp1$r[1]
    mod2 <- lm(Y ~ R + X + M + C + R * M, data = pop.dat)
    params2 <- parameters::model_parameters(mod2)
    tmp2 <- effectsize::t_to_r(params2$t[- 1], df_error = params2$df_error[- 1])
    MY.effect <- tmp2$r[3]
    
    # True values of the disparity reduction and disparity remaining
    set.seed(seed)
    y0 <- mean(pop.dat$Y[pop.dat$R == 0 & pop.dat$C.grp == 2])
    y1 <- mean(pop.dat$Y[pop.dat$R == 1 & pop.dat$C.grp == 2])
    m0 <- pop.dat$M[pop.dat$R == 0 & pop.dat$C.grp == 2]
    m1 <- pop.dat$M[pop.dat$R == 1 & pop.dat$C.grp == 2]
    m01 <- sample(m0, (sum(pop.dat$R == 1 & pop.dat$C.grp == 2)), replace = TRUE)
    ym0 <- mean(coef.Y[1] + coef.Y[2]*1 + coef.Y[3]*pop.dat$X[pop.dat$R == 1 & pop.dat$C.grp == 2] +
                  coef.Y[4]*1 + coef.Y[5]*m01 + coef.Y[6]*1*m01 + err.Y[pop.dat$R == 1 & pop.dat$C.grp == 2])
    
    print.ratio <- (mean(m1) - mean(m0))/(coef.Y[5] + coef.Y[6])
    pop.dat$R <- as.factor(pop.dat$R)
    
  }
  
  Red.true <- y1 - ym0
  Rem.true <- ym0 - y0
  
  print.out <- matrix(c(Red.true/(Red.true + Rem.true), print.ratio,
                        RM.effect, MY.effect, adjustedCoefs, Red.true, Rem.true), ncol = 1)
  rownames(print.out) <- c("true effect size (fixed as 0.3)", "true ratio", "R-M effect", "M-Y effect",
                           "adjusted b1", "adjusted c1", "adjusted c3",
                           "disparity reduction", "disparity remaining")
  cat("\n")
  print(round(print.out, 3))
  
  out <- list(data = pop.dat, RM.effect = RM.effect, MY.effect = MY.effect,
              Red.true = Red.true, Rem.true = Rem.true, Ini.true = Red.true + Rem.true,
              ratio.true = print.ratio, true.effectsize = Red.true/(Red.true + Rem.true))
  return(out)
  
}

# Estimator 1. Difference-in-coefficients regression
est1 <- function(samp.data, n.boot, conf.level, Red.true, Rem.true, binary.mediator){
  
  ## Estimator 1. Difference-in-coefficients regression
  
  # Calculate the estimates of the disparity reduction and remaining
  if(binary.mediator){
    mod1 <- lm(Y.bin ~ R + C.grp.ref2, data = samp.data)
    mod2 <- lm(Y.bin ~ R + X + C.grp.ref2, data = samp.data)
    mod3 <- lm(Y.bin ~ R + X + M.bin + C.grp.ref2, data = samp.data)
  } else {
    mod1 <- lm(Y ~ R + C.grp.ref2, data = samp.data)
    mod2 <- lm(Y ~ R + X + C.grp.ref2, data = samp.data)
    mod3 <- lm(Y ~ R + X + M + C.grp.ref2, data = samp.data)
  }
  
  phi <- coef(mod1)
  gamma <- coef(mod2)
  theta <- coef(mod3)
  disp.red <- gamma[2] - theta[2] + (1 - theta[3]/gamma[3]) * (phi[2] - gamma[2])
  disp.rem <- theta[2] + (theta[3]/gamma[3]) * (phi[2] - gamma[2])
  
  # Bootsrapping for coverage
  Red <- Rem <- rep(NA, n.boot)
  
  for(b in 1:n.boot){
    
    set.seed(b)
    boot.ind <- sample(nrow(samp.data), nrow(samp.data), replace = TRUE)
    bt.data <- data.frame(samp.data[boot.ind, ])
    
    if(binary.mediator){
      mod1 <- lm(Y.bin ~ R + C.grp.ref2, data = bt.data)
      mod2 <- lm(Y.bin ~ R + X + C.grp.ref2, data = bt.data)
      mod3 <- lm(Y.bin ~ R + X + M.bin + C.grp.ref2, data = bt.data)
    } else {
      mod1 <- lm(Y ~ R + C.grp.ref2, data = bt.data)
      mod2 <- lm(Y ~ R + X + C.grp.ref2, data = bt.data)
      mod3 <- lm(Y ~ R + X + M + C.grp.ref2, data = bt.data)
    }
    phi <- coef(mod1)
    gamma <- coef(mod2)
    theta <- coef(mod3)
    Red[b] <- gamma[2] - theta[2] + (1 - theta[3]/gamma[3]) * (phi[2] - gamma[2])###theta[2] - previously
    Rem[b] <- theta[2] + (theta[3]/gamma[3]) * (phi[2] - gamma[2])
    
  }
  
  Red.ci <- quantile(Red, prob = c((1 - conf.level)/2, 1/2 + conf.level/2), na.rm = TRUE)
  Rem.ci <- quantile(Rem, prob = c((1 - conf.level)/2, 1/2 + conf.level/2), na.rm = TRUE)
  red.cvg <- Red.true >= Red.ci[1] & Red.true <= Red.ci[2]
  rem.cvg <- Rem.true >= Rem.ci[1] & Rem.true <= Rem.ci[2]
  
  return(c(disp.red, disp.rem, red.cvg, rem.cvg, Red.ci, Rem.ci))
}

# Estimator 2. Product-of-coefficients regression
est2 <- function(samp.data, n.boot, conf.level, Red.true, Rem.true, binary.mediator){
  
  ## Estimator 2. Product-of-coefficients regression
  
  # Calculate the estimates of the disparity reduction and remaining
  if(binary.mediator){
    mod10 <- glm(M.bin ~ C.grp.ref2, data = subset(samp.data, R == 0), family = binomial(logit))
    mod11 <- glm(M.bin ~ C.grp.ref2, data = subset(samp.data, R == 1), family = binomial(logit))
    alpha <- exp(coef(mod11)[1])/(1 + exp(coef(mod11)[1])) - exp(coef(mod10)[1])/(1 + exp(coef(mod10)[1]))
    mod2 <- lm(Y.bin ~ R + X + M.bin + R * M.bin + C.grp.ref2, data = samp.data)
    mod0 <- lm(Y.bin ~ R + C.grp.ref2, data = samp.data)
  } else {
    mod1 <- lm(M ~ R + C.grp.ref2, data = samp.data)
    alpha <- coef(mod1)[2]
    mod2 <- lm(Y ~ R + X + M + R * M + C.grp.ref2, data = samp.data)
    mod0 <- lm(Y ~ R + C.grp.ref2, data = samp.data)
  }
  
  beta <- coef(mod2)
  phi <- coef(mod0)
  disp.red <- alpha * (beta[4] + beta[length(beta)])
  disp.rem <- phi[2] - alpha * (beta[4] + beta[length(beta)])
  
  # Bootsrapping for coverage
  Red <- Rem <- rep(NA, n.boot)
  
  for(b in 1:n.boot){
    
    set.seed(b)
    boot.ind <- sample(nrow(samp.data), nrow(samp.data), replace = TRUE)
    bt.data <- data.frame(samp.data[boot.ind, ])
    
    if(binary.mediator){
      mod10 <- glm(M.bin ~ C.grp.ref2, data = subset(bt.data, R == 0), family = binomial)
      mod11 <- glm(M.bin ~ C.grp.ref2, data = subset(bt.data, R == 1), family = binomial)
      alpha <-exp(coef(mod11)[1])/(1+ exp(coef(mod11)[1])) - exp(coef(mod10)[1])/(1 + exp(coef(mod10)[1]))
      mod2 <- lm(Y.bin ~ R + X + M.bin + R * M.bin + C.grp.ref2, data = bt.data)
      mod0 <- lm(Y.bin ~ R + C.grp.ref2, data = bt.data)
    } else {
      mod1 <- lm(M ~ R + C.grp.ref2, data = bt.data)
      alpha <- coef(mod1)[2]
      mod2 <- lm(Y ~ R + X + M + R * M + C.grp.ref2, data = bt.data)
      mod0 <- lm(Y ~ R + C.grp.ref2, data = bt.data)
    }
    
    beta <- coef(mod2)   
    phi <- coef(mod0)
    Red[b] <- alpha * (beta[4] + beta[length(beta)])
    Rem[b] <- phi[2] - alpha * (beta[4] + beta[length(beta)])
    
  }
  
  Red.ci <- quantile(Red, prob = c((1 - conf.level)/2, 1/2 + conf.level/2), na.rm = TRUE)
  Rem.ci <- quantile(Rem, prob = c((1 - conf.level)/2, 1/2 + conf.level/2), na.rm = TRUE)
  red.cvg <- Red.true >= Red.ci[1] & Red.true <= Red.ci[2]
  rem.cvg <- Rem.true >= Rem.ci[1] & Rem.true <= Rem.ci[2]
  
  return(c(disp.red, disp.rem, red.cvg, rem.cvg, Red.ci, Rem.ci))
}

# Estimator 2. Product-of-coefficients regression (alternative)
est2.v2 <- function(samp.data, n.boot, conf.level, Red.true, Rem.true, binary.mediator){
  
  ## Estimator 2. Product-of-coefficients regression
  
  # Calculate the estimates of the disparity reduction and remaining
  if(binary.mediator){
    mod10 <- glm(M.bin ~ C.grp.ref2, data = subset(samp.data, R == 0), family = binomial(logit))
    mod11 <- glm(M.bin ~ C.grp.ref2, data = subset(samp.data, R == 1), family = binomial(logit))
    alpha1 <- exp(coef(mod11)[1])/(1 + exp(coef(mod11)[1]))- exp(coef(mod10)[1])/(1 + exp(coef(mod10)[1]))
    mod2 <- lm(Y.bin ~ R + X + M.bin + R * M.bin + C.grp.ref2, data = samp.data)
    mod0 <- lm(Y.bin ~ R + C.grp.ref2, data = samp.data)
    
    modX <- lm(X ~ R + C.grp.ref2, data = samp.data)
    gamma1 <- coef(modX)[2]
    alpha0 <- exp(coef(mod10)[1])/(1 + exp(coef(mod10)[1]))
  } else {
    mod1 <- lm(M ~ R + C.grp.ref2, data = samp.data)
    alpha1 <- coef(mod1)[2]
    mod2 <- lm(Y ~ R + X + M + R * M + C.grp.ref2, data = samp.data)
    mod0 <- lm(Y ~ R + C.grp.ref2, data = samp.data)
    
    modX <- lm(X ~ R + C.grp.ref2, data = samp.data)
    gamma1 <- coef(modX)[2]
    alpha0 <- coef(mod1)[1] + coef(mod1)[3] * 2
  }
  beta <- coef(mod2)
  phi <- coef(mod0)
  disp.red <- alpha1 * (beta[4] + beta[length(beta)])
  disp.rem <- beta[2] + beta[3] * gamma1 + beta[length(beta)] * alpha0
  
  ## Bootsrapping for coverage
  Red <- Rem <- rep(NA, n.boot)
  for(b in 1:n.boot){
    
    set.seed(b)
    boot.ind <- sample(nrow(samp.data), nrow(samp.data), replace = TRUE)
    bt.data <- data.frame(samp.data[boot.ind, ])
    
    if(binary.mediator){
      mod10 <- glm(M.bin ~ C.grp.ref2, data = subset(bt.data, R == 0), family = binomial)
      mod11 <- glm(M.bin ~ C.grp.ref2, data = subset(bt.data, R == 1), family = binomial)
      alpha1 <-exp(coef(mod11)[1])/(1 + exp(coef(mod11)[1])) - exp(coef(mod10)[1])/(1 + exp(coef(mod10)[1]))
      mod2 <- lm(Y.bin ~ R + X + M.bin + R * M.bin + C.grp.ref2, data = bt.data)
      mod0 <- lm(Y.bin ~ R + C.grp.ref2, data = bt.data)
      
      modX <- lm(X ~ R + C.grp.ref2, data = bt.data)
      gamma1 <- coef(modX)[2]
      alpha0 <- exp(coef(mod10)[1])/(1+exp(coef(mod10)[1]))
    } else {
      mod1 <- lm(M ~ R + C.grp.ref2, data = bt.data)
      alpha1 <- coef(mod1)[2]
      mod2 <- lm(Y ~ R + X + M + R * M + C.grp.ref2, data = bt.data)
      mod0 <- lm(Y ~ R + C.grp.ref2, data = bt.data)
      
      modX <- lm(X ~ R + C.grp.ref2, data = bt.data)
      gamma1 <- coef(modX)[2]
      alpha0 <- coef(mod1)[1] + coef(mod1)[3] * 2
    }
    beta <- coef(mod2)   
    phi <- coef(mod0)
    Red[b] <- alpha1 * (beta[4] + beta[length(beta)])
    Rem[b] <- beta[2] + beta[3] * gamma1 + beta[length(beta)] * alpha0
    
  }
  Red.ci <- quantile(Red, prob = c((1 - conf.level)/2, 1/2 + conf.level/2), na.rm=TRUE)
  Rem.ci <- quantile(Rem, prob = c((1 - conf.level)/2, 1/2 + conf.level/2), na.rm=TRUE)
  red.cvg <- Red.true >= Red.ci[1] & Red.true <= Red.ci[2]
  rem.cvg <- Rem.true >= Rem.ci[1] & Rem.true <= Rem.ci[2]
  
  return(c(disp.red, disp.rem, red.cvg, rem.cvg, Red.ci, Rem.ci))
  
}

# Estimator 3. RMPW
est3 <- function(samp.data, n.boot, conf.level, Red.true, Rem.true){
  
  ## Estimator 3. RMPW
  
  # Fit regression models for M
  fit.m0 <- glm(M.bin ~ C.grp.ref2, data = subset(samp.data, R == 0), family = binomial(logit))
  fit.m1 <- glm(M.bin ~ X + C.grp.ref2, data = subset(samp.data, R == 1), family = binomial(logit))
  
  # Calculate P(m|r)
  pm1 <- predict(fit.m1, samp.data, type = "response")
  pm0 <- predict(fit.m0, samp.data, type = "response")
  
  # Calculate weights
  WM <- rep(NA, nrow(samp.data))
  ind11 <- samp.data$R == 1 & samp.data$M.bin == 1
  ind10 <- samp.data$R == 1 & samp.data$M.bin == 0
  WM[ind11] <- pm0[ind11]/pm1[ind11]
  WM[ind10] <- (1 - pm0[ind10])/(1 - pm1[ind10])
  WM[samp.data$R == 0] <- 0
  samp.data <- data.frame(samp.data, WM)
  
  # Calculate weighted y
  wmu1xdm <- lm(Y.bin ~ C.grp.ref2, weight = WM, data = subset(samp.data, R == 1))
  wy1 <- lm(Y.bin ~ C.grp.ref2, data = subset(samp.data, R == 1))
  wy0 <- lm(Y.bin ~ C.grp.ref2, data = subset(samp.data, R == 0))
  
  # Calculate the estimates of the disparity reduction and remaining
  disp.red = wy1$coef[1] - wmu1xdm$coef[1]
  disp.rem = wmu1xdm$coef[1] - wy0$coef[1]
  
  # Bootsrapping for coverage
  Red <- Rem <- rep(NA, n.boot)
  
  for(b in 1:n.boot){
    
    set.seed(b)
    boot.ind <- sample(nrow(samp.data), nrow(samp.data), replace = T)
    bt.data <- data.frame(samp.data[boot.ind, ])
    
    # Fit regression models for M
    fit.m0 <- glm(M.bin ~ C.grp.ref2, data = subset(bt.data, R == 0), family = binomial)
    fit.m1 <- glm(M.bin ~ X + C.grp.ref2, data = subset(bt.data, R == 1), family = binomial)
    
    # Calculate P(m|r)
    pm1 <- predict(fit.m1, bt.data, type = "response")
    pm0 <- predict(fit.m0, bt.data, type = "response")
    
    # Calculate weight
    WM <- rep(NA, nrow(bt.data))
    ind11 <- bt.data$R == 1 & bt.data$M.bin == 1
    ind10 <- bt.data$R == 1 & bt.data$M.bin == 0
    WM[ind11] <- pm0[ind11]/pm1[ind11]
    WM[ind10] <- (1 - pm0[ind10])/(1 - pm1[ind10])
    WM[bt.data$R == 0] <- 0
    bt.data <- data.frame(bt.data, WM)
    
    # Calculate weighted y
    wmu1xdm <- lm(Y.bin ~ C.grp.ref2, weight = WM, data = subset(bt.data, R == 1))
    wy1 <- lm(Y.bin ~ C.grp.ref2, data = subset(bt.data, R == 1))
    wy0 <- lm(Y.bin ~ C.grp.ref2, data = subset(bt.data, R == 0))
    
    # Calculate the estimates of the disparity reduction and remaining
    Red[b] <- wy1$coef[1] - wmu1xdm$coef[1]
    Rem[b] <- wmu1xdm$coef[1] - wy0$coef[1]
    
  }
  
  Red.ci <- quantile(Red, prob = c((1 - conf.level)/2, 1/2 + conf.level/2), na.rm = TRUE)
  Rem.ci <- quantile(Rem, prob = c((1 - conf.level)/2, 1/2 + conf.level/2), na.rm = TRUE)
  red.cvg <- Red.true >= Red.ci[1] & Red.true <= Red.ci[2]
  rem.cvg <- Rem.true >= Rem.ci[1] & Rem.true <= Rem.ci[2]
  
  return(c(disp.red, disp.rem, red.cvg, rem.cvg, Red.ci, Rem.ci))
}

# Estimator 4. IORW
est4 <- function(samp.data, n.boot, conf.level, Red.true, Rem.true){
  
  ## Estimator 4. IORW instead of RMPW
  
  # Fit regression models for M
  fit.m.c <- glm(M.bin ~ C.grp.ref2, data = samp.data, family = binomial(logit))
  fit.m.xc <- glm(M.bin ~ X + C.grp.ref2, data = samp.data, family = binomial(logit))
  
  # Fit regression models for R
  fit.r.c <- glm(R ~ C.grp.ref2, data = samp.data, family = binomial(logit))
  fit.r.xc <- glm(R ~ X + C.grp.ref2, data = samp.data, family = binomial(logit))
  fit.r.mc <- glm(R ~ M.bin + C.grp.ref2, data = samp.data, family = binomial(logit))
  fit.r.mxc <- glm(R ~ M.bin + X + C.grp.ref2, data = samp.data, family = binomial(logit))
  
  # Calculate Probabilities
  pr.m.c <- predict(fit.m.c, samp.data, type = "response")
  pr.m.xc <- predict(fit.m.xc, samp.data, type = "response")
  pr.r.c <- predict(fit.r.c, samp.data, type = "response")
  pr.r.xc <- predict(fit.r.xc, samp.data, type = "response")
  pr.r.mc <- predict(fit.r.mc, samp.data, type = "response")
  pr.r.mxc <- predict(fit.r.mxc, samp.data, type = "response")
  
  # Calculate weights
  WM <- rep(NA, nrow(samp.data))
  ind11 <- samp.data$R == 1 & samp.data$M.bin == 1
  ind10 <- samp.data$R == 1 & samp.data$M.bin == 0
  WM[ind11] <- ((1-pr.r.mc[ind11])/pr.r.mxc[ind11]) / ((1 - pr.r.c[ind11])/pr.r.xc[ind11]) * (pr.m.c[ind11]/pr.m.xc[ind11])
  WM[ind10] <- ((1-pr.r.mc[ind10])/pr.r.mxc[ind10]) / ((1 - pr.r.c[ind10])/pr.r.xc[ind10]) * ((1 - pr.m.c[ind10])/(1 - pr.m.xc[ind10]))
  WM[samp.data$R == 0] <- 0
  samp.data <- data.frame(samp.data, WM)
  
  # Calculate weighted y
  wmu1xdm <- lm(Y.bin ~ C.grp.ref2, weight = WM, data = subset(samp.data, R == 1))
  wy1 <- lm(Y.bin ~ C.grp.ref2, data = subset(samp.data, R == 1))
  wy0 <- lm(Y.bin ~ C.grp.ref2, data = subset(samp.data, R == 0))
  
  # Calculate the estimates of the disparity reduction and remaining
  disp.red = wy1$coef[1] - wmu1xdm$coef[1]
  disp.rem = wmu1xdm$coef[1] - wy0$coef[1]
  
  # Bootsrapping for coverage
  Red <- Rem <- rep(NA, n.boot)
  
  for(b in 1:n.boot){
    
    set.seed(b)
    boot.ind <- sample(nrow(samp.data), nrow(samp.data), replace = TRUE)
    bt.data <- data.frame(samp.data[boot.ind, ])
    
    # Fit regression models for M
    fit.m.c <- glm(M.bin ~ C.grp.ref2, data = bt.data, family = binomial)
    fit.m.xc <- glm(M.bin ~ X + C.grp.ref2, data = bt.data, family = binomial)
    
    fit.m.c <- glm(M.bin ~ C.grp.ref2, data = bt.data, family = binomial(logit))
    fit.m.xc <- glm(M.bin ~ X + C.grp.ref2, data = bt.data, family = binomial(logit))
    
    # Fit regression models for R
    fit.r.c <- glm(R ~ C.grp.ref2, data = subset(bt.data), family = binomial(logit))
    fit.r.xc <- glm(R ~ X + C.grp.ref2, data = subset(bt.data), family = binomial(logit))
    fit.r.mc <- glm(R ~ M.bin + C.grp.ref2, data = subset(bt.data), family = binomial(logit))
    fit.r.mxc <- glm(R ~ M.bin + X + C.grp.ref2, data = subset(bt.data), family = binomial(logit))
    
    # Calculate Probabilities
    pr.m.c <- predict(fit.m.c, bt.data, type = "response")
    pr.m.xc <- predict(fit.m.xc, bt.data, type = "response")
    pr.r.c <- predict(fit.r.c, bt.data, type = "response")
    pr.r.xc <- predict(fit.r.xc, bt.data, type = "response")
    pr.r.mc <- predict(fit.r.mc, bt.data, type = "response")
    pr.r.mxc <- predict(fit.r.mxc, bt.data, type = "response")
    
    # Calculate weights
    WM <- rep(NA, nrow(bt.data))
    ind11 <- bt.data$R == 1 & bt.data$M.bin == 1
    ind10 <- bt.data$R == 1 & bt.data$M.bin == 0
    WM[ind11] <- ((1-pr.r.mc[ind11])/pr.r.mxc[ind11]) / ((1 - pr.r.c[ind11])/pr.r.xc[ind11]) * (pr.m.c[ind11]/pr.m.xc[ind11])
    WM[ind10] <- ((1-pr.r.mc[ind10])/pr.r.mxc[ind10]) / ((1 - pr.r.c[ind10])/pr.r.xc[ind10]) * ((1 - pr.m.c[ind10])/(1 - pr.m.xc[ind10]))
    WM[bt.data$R == 0] <- 0
    bt.data <- data.frame(bt.data, WM)
    
    # Calculate weighted y
    wmu1xdm <- lm(Y.bin ~ C.grp.ref2, weight = WM, data = subset(bt.data, R == 1))
    wy1 <- lm(Y.bin ~ C.grp.ref2, data = subset(bt.data, R == 1))
    wy0 <- lm(Y.bin ~ C.grp.ref2, data = subset(bt.data, R == 0))
    
    # Calculate the estimates of the disparity reduction and remaining
    Red[b] <- wy1$coef[1] - wmu1xdm$coef[1]
    Rem[b] <- wmu1xdm$coef[1] - wy0$coef[1]
    
  }
  
  Red.ci <- quantile(Red, prob = c((1 - conf.level)/2, 1/2 + conf.level/2), na.rm = TRUE)
  Rem.ci <- quantile(Rem, prob = c((1 - conf.level)/2, 1/2 + conf.level/2), na.rm = TRUE)
  red.cvg <- Red.true >= Red.ci[1] & Red.true <= Red.ci[2]
  rem.cvg <- Rem.true >= Rem.ci[1] & Rem.true <= Rem.ci[2]
  
  return(c(disp.red, disp.rem, red.cvg, rem.cvg, Red.ci, Rem.ci))
}

# Estimator 5. Single-mediator-imputation
est5 <- function(samp.data, n.boot, conf.level, Red.true, Rem.true, binary.mediator){
  
  ## Estimator 5. Single-mediator-imputation
  
  if(binary.mediator){
    fit.y <- lm(Y.bin ~ R + X + M.bin + R*M.bin + C.grp.ref2, data = samp.data)
    fit.m1 <- glm(M.bin ~ R + C.grp.ref2, data = samp.data, family = binomial(logit))
  } else {
    fit.y <- lm(Y ~ R + X + M + R*M + C.grp.ref2, data = samp.data)
    fit.m1 <- lm(M ~ R + C.grp.ref2, data = samp.data)
  }
  
  # Calculate weighted outcome
  dat <- model.frame(fit.y)
  n <- nrow(dat)
  
  # Predict mediator
  # Create data in which R = 0
  dat_0 <- dat_1 <- dat 
  dat_0[, "R"] <- levels(dat[, "R"])[1]
  
  # Compute predicted values of mediator
  if(binary.mediator){
    wy0 <- lm(Y.bin ~ C.grp.ref2, data = subset(dat, R == 0)) 
    wy1 <- lm(Y.bin ~ C.grp.ref2, data = subset(dat, R == 1)) 
    muM <- predict(fit.m1, newdata = dat_0, type = "response")
  } else {
    wy0 <- lm(Y ~ C.grp.ref2, data = subset(dat, R == 0)) 
    wy1 <- lm(Y ~ C.grp.ref2, data = subset(dat, R == 1)) 
    sigma <- summary(fit.m1)$sigma
    error <- rnorm(n, mean = 0, sd = sigma)
    muM <- predict(fit.m1, type = "response", newdata = dat_0) + error
  }
  
  # Predict outcomes
  if(binary.mediator){
    PredictM <- rbinom(n, size = 1, prob = muM)
    dat_1[, "M.bin"] <- as.factor(PredictM)
    dat$muldm <- predict(fit.y, newdata = dat_1)
  } else {
    dat_1[, "M"] <- muM 
    dat$muldm <- predict(fit.y, newdata = dat_1)
  }
  
  # Compute outcomes after incorporating predicted values of mediator
  a <- lm(muldm ~ C.grp.ref2, data = subset(dat, R == 1))
  wmuldm <- a$coef[1]
  w01 <- mean(as.numeric(wmuldm))
  
  # Calculate the estimates of the disparity reduction and remaining
  disp.red <- wy1$coef[1] - w01
  disp.rem <- w01 - wy0$coef[1]
  
  # Bootsrapping for coverage
  Red <- Rem <- rep(NA, n.boot)
  
  for(b in 1:n.boot){
    
    set.seed(b)
    boot.ind <- sample(nrow(samp.data), nrow(samp.data), replace = TRUE)
    bt.data <- data.frame(samp.data[boot.ind, ])
    
    if(binary.mediator){ 
      fit.y <- lm(Y.bin ~ R + X + M.bin + R*M.bin + C.grp.ref2, data = bt.data)
      fit.m1 <- glm(M.bin ~ R + C.grp.ref2, data = bt.data, family = binomial)
    } else {
      fit.y <- lm(Y ~ R + X + M + R*M + C.grp.ref2, data = bt.data)
      fit.m1 <- lm(M ~ R + C.grp.ref2, data = bt.data)
    }
    
    # Calculate weighted outcome
    dat <- model.frame(fit.y)
    n <- nrow(dat)
    
    # Predict mediator
    # Create data in which R = 0
    dat_0 <- dat_1 <- dat 
    dat_0[, "R"] <- levels(dat[, "R"])[1]
    
    # Compute predicted values of mediator
    if(binary.mediator){
      wy0 <- lm(Y.bin ~ C.grp.ref2, data = subset(dat, R == 0)) 
      wy1 <- lm(Y.bin ~ C.grp.ref2, data = subset(dat, R == 1))
      muM <- predict(fit.m1, newdata = dat_0, type = "response")
    } else {
      wy0 <- lm(Y ~ C.grp.ref2, data = subset(dat, R == 0)) 
      wy1 <- lm(Y ~ C.grp.ref2, data = subset(dat, R == 1)) 
      sigma <- summary(fit.m1)$sigma
      error <- rnorm(n, mean = 0, sd = sigma)
      muM <- predict(fit.m1, type = "response", newdata = dat_0) + error
    }
    
    # Predict outcomes
    if(binary.mediator){
      PredictM <- rbinom(n, size = 1, prob = muM)
      dat_1[, "M.bin"] <- as.factor(PredictM)
      dat$muldm <- predict(fit.y, newdata = dat_1)
    } else {
      dat_1[, "M"] <- muM 
      dat$muldm <- predict(fit.y, newdata = dat_1)
    }
    
    # Compute outcomes after incorporating predicted values of mediator
    a <- lm(muldm ~ C.grp.ref2, data = subset(dat, R == 1))
    wmuldm <- a$coef[1]
    w01 <- mean(as.numeric(wmuldm))
    
    # Calculate the estimates of the disparity reduction and remaining
    Red[b] <- wy1$coef[1] - w01
    Rem[b] <- w01 - wy0$coef[1]
    
  }
  
  Red.ci <- quantile(Red, prob = c((1 - conf.level)/2, 1/2 + conf.level/2), na.rm = TRUE)
  Rem.ci <- quantile(Rem, prob = c((1 - conf.level)/2, 1/2 + conf.level/2), na.rm = TRUE)
  red.cvg <- Red.true >= Red.ci[1] & Red.true <= Red.ci[2]
  rem.cvg <- Rem.true >= Rem.ci[1] & Rem.true <= Rem.ci[2]
  
  if(binary.mediator){
    out <- c(mean(Red), mean(Rem), red.cvg, rem.cvg, Red.ci, Rem.ci)
  } else {
    out <- c(disp.red, disp.rem, red.cvg, rem.cvg, Red.ci, Rem.ci)
  }
  return(out)
  
}

# Estimator 6. Multiple-mediator-imputation
est6 <- function(samp.data, n.boot, conf.level, Red.true, Rem.true, binary.mediator){
  
  ## Estimator 6. Multiple-mediator-imputation
  
  if(binary.mediator){
    fit.y <- lm(Y.bin ~ R + X + M.bin + R*M.bin + C.grp.ref2, data = samp.data)
  } else {
    fit.y <- lm(Y ~ R + X + M + R*M + C.grp.ref2, data = samp.data)
  }
  fit.x <- lm(X ~ R + C.grp.ref2, data = samp.data)
  
  # Calculate weighted outcome
  dat <- model.frame(fit.y)
  n <- nrow(dat)
  
  if(binary.mediator){
    wy0 <- lm(Y.bin ~ C.grp.ref2, data = subset(dat, R == 0)) 
    wy1 <- lm(Y.bin ~ C.grp.ref2, data = subset(dat, R == 1)) 
  } else {
    wy0 <- lm(Y ~ C.grp.ref2, data = subset(dat, R == 0)) 
    wy1 <- lm(Y ~ C.grp.ref2, data = subset(dat, R == 1)) 
  }
  
  # Predict confounder
  # Create data in which R = 2
  dat_1 <- dat 
  dat_1[, "R"] <- levels(dat[, "R"])[2]
  
  sigma <- summary(fit.x)$sigma
  error <- rnorm(n, mean = 0, sd = sigma)
  PredictX <- predict(fit.x, type = "response", newdata = dat_1) + error
  
  # Predict outcomes
  dat_1[, "X"] <- PredictX
  dat$muldm <- predict(fit.y, newdata = dat_1)
  
  # Compute outcomes after incorporating predicted values of mediator
  a <- lm(muldm ~ C.grp.ref2, data = subset(dat, R == 0))
  wmuldm <- a$coef[1]
  w01 <- mean(as.numeric(wmuldm))
  
  # Calculate the estimates of the disparity reduction and remaining
  disp.red <- wy1$coef[1] - w01
  disp.rem <- w01 - wy0$coef[1]
  
  # Bootsrapping for coverage
  Red <- Rem <- rep(NA, n.boot)
  
  for(b in 1:n.boot){
    
    set.seed(b)
    boot.ind <- sample(nrow(samp.data), nrow(samp.data), replace = TRUE)
    bt.data <- data.frame(samp.data[boot.ind, ])
    
    if(binary.mediator){   
      fit.y <- lm(Y.bin ~ R + X + M.bin + R*M.bin + C.grp.ref2, data = bt.data)
    } else {
      fit.y <- lm(Y ~ R + X + M + R*M + C.grp.ref2, data = bt.data)
    }
    fit.x <- lm(X ~ R + C.grp.ref2, data = bt.data)
    
    # Calculate weighted outcome
    dat <- model.frame(fit.y)
    n <- nrow(dat)
    if(binary.mediator){ 
      wy0 <- lm(Y.bin ~ C.grp.ref2, data = subset(dat, R == 0)) 
      wy1 <- lm(Y.bin ~ C.grp.ref2, data = subset(dat, R == 1)) 
    } else {
      wy0 <- lm(Y ~ C.grp.ref2, data = subset(dat, R == 0)) 
      wy1 <- lm(Y ~ C.grp.ref2, data = subset(dat, R == 1)) 
    }
    
    # Predict mediator
    # Create data in which R = 2
    dat_1 <- dat 
    dat_1[, "R"] <- levels(dat[, "R"])[2]
    
    # Compute predicted values of mediator
    sigma <- summary(fit.x)$sigma
    error <- rnorm(n, mean = 0, sd = sigma)
    PredictX <- predict(fit.x, type = "response", newdata = dat_1) + error
    
    # Predict outcomes
    dat_1[, "X"] <- PredictX
    dat$muldm <- predict(fit.y, newdata = dat_1)
    
    # Compute outcomes after incorporating predicted values of mediator
    a <- lm(muldm ~ C.grp.ref2, data = subset(dat, R == 0))
    wmuldm <- a$coef[1]
    w01 <- mean(as.numeric(wmuldm))
    
    # Calculate the estimates of the disparity reduction and remaining
    Red[b] <- wy1$coef[1] - w01
    Rem[b] <- w01 - wy0$coef[1]
    
  }
  
  Red.ci <- quantile(Red, prob = c((1 - conf.level)/2, 1/2 + conf.level/2), na.rm = TRUE)
  Rem.ci <- quantile(Rem, prob = c((1 - conf.level)/2, 1/2 + conf.level/2), na.rm = TRUE)
  red.cvg <- Red.true >= Red.ci[1] & Red.true <= Red.ci[2]
  rem.cvg <- Rem.true >= Rem.ci[1] & Rem.true <= Rem.ci[2]
  
  return(c(disp.red, disp.rem, red.cvg, rem.cvg, Red.ci, Rem.ci))
}

# Estimation functions
sim.one <- function(x, pop.dat, Red.true, Rem.true,
                    pop.dat.est1, Red.true.est1, Rem.true.est1,
                    n = 1000, n.boot = 1000, conf.level = .95, binary.mediator = FALSE,
                    standardize = FALSE){
  
  set.seed(x)
  samp.ind <- sample(nrow(pop.dat), n, replace = FALSE)
  samp.data <- data.frame(pop.dat[samp.ind, ])
  samp.data.est1 <- data.frame(pop.dat.est1[samp.ind, ])
  
  ## standardize
  if(standardize){
    
    if(binary.mediator){
      samp.data$X <- as.numeric(scale(samp.data$X))
      samp.data$Y.bin <- as.numeric(scale(samp.data$Y.bin))
    } else {
      samp.data$X <- as.numeric(scale(samp.data$X))
      samp.data$Y <- as.numeric(scale(samp.data$Y))
      samp.data$M <- as.numeric(scale(samp.data$M))
    }
    
  }
  
  out1 <- tryCatch(est1(samp.data = samp.data.est1, n.boot = n.boot, conf.level = conf.level,
                        Red.true = Red.true.est1, Rem.true = Rem.true.est1, binary.mediator = binary.mediator), 
                   error = function(e) rep(NA, 8), warning = function(w) rep(NA, 8))
  out2 <- tryCatch(est2(samp.data = samp.data, n.boot = n.boot, conf.level = conf.level,
                        Red.true = Red.true, Rem.true = Rem.true, binary.mediator = binary.mediator), 
                   error = function(e) rep(NA, 8), warning = function(w) rep(NA, 8))
  if(binary.mediator){
    out2.alt <- tryCatch(est2.v2(samp.data = samp.data, n.boot = n.boot, conf.level = conf.level,
                          Red.true = Red.true, Rem.true = Rem.true, binary.mediator = binary.mediator), 
                     error = function(e) rep(NA, 8), warning = function(w) rep(NA, 8))
    out3 <- tryCatch(est3(samp.data = samp.data, n.boot = n.boot, conf.level = conf.level,
                          Red.true = Red.true, Rem.true = Rem.true), 
                     error = function(e) rep(NA, 8), warning = function(w) rep(NA, 8))
    out4 <- tryCatch(est4(samp.data = samp.data, n.boot = n.boot, conf.level = conf.level,
                               Red.true = Red.true, Rem.true = Rem.true), 
                     error = function(e) rep(NA, 8), warning = function(w) rep(NA, 8))
  }
  out5 <- tryCatch(est5(samp.data = samp.data, n.boot = n.boot, conf.level = conf.level,
                        Red.true = Red.true, Rem.true = Rem.true, binary.mediator = binary.mediator), 
                   error = function(e) rep(NA, 8), warning = function(w) rep(NA, 8))
  out6 <- tryCatch(est6(samp.data = samp.data, n.boot = n.boot, conf.level = conf.level,
                        Red.true = Red.true, Rem.true = Rem.true, binary.mediator = binary.mediator), 
                   error = function(e) rep(NA, 8), warning = function(w) rep(NA, 8))

  if(binary.mediator){
    out <- c(out1, out2, out3, out4, out5, out6, out2.alt)
  } else {
    out <- c(out1, out2, out5, out6)
  }
  return(out)
  
}

summarize.results <- function(out.mclapply, n.sim = 1000, Red.true, Rem.true,
                              Red.true.est1, Rem.true.est1, binary.mediator){
  
  if(binary.mediator){
    
    no.est <- 7
    
    Red.out <- Rem.out <- matrix(NA, no.est, 6)
    rownames(Red.out) <- rownames(Rem.out) <- c("Est1", "Est2", "Est3", "Est4", "Est5", "Est6", "Est2.alt")
    colnames(Red.out) <- colnames(Rem.out) <-
      c("estimate", paste(conf.level * 100, "% CI Coverage", sep = ""), "Bias", "RMSE", "rBias", "nRMSE")
    
    Red.CI.out <- Rem.CI.out <- matrix(NA, n.sim, no.est * 2)
    colnames(Red.CI.out) <- colnames(Rem.CI.out) <-
      c("Est1.L", "Est1.U", "Est2.L", "Est2.U", "Est3.L", "Est3.U", "Est4.L", "Est4.U", "Est5.L", "Est5.U",
        "Est6.L", "Est6.U", "Est2.alt.L", "Est2.alt.U")
    
  } else {
    
    no.est <- 4
    
    Red.out <- Rem.out <- matrix(NA, no.est, 6)
    rownames(Red.out) <- rownames(Rem.out) <- c("Est1", "Est2", "Est5", "Est6")
    colnames(Red.out) <- colnames(Rem.out) <-
      c("estimate", paste(conf.level * 100, "% CI Coverage", sep = ""), "Bias", "RMSE", "rBias", "nRMSE")
    
    Red.CI.out <- Rem.CI.out <- matrix(NA, n.sim, no.est * 2)
    colnames(Red.CI.out) <- colnames(Rem.CI.out) <-
      c("Est1.L", "Est1.U", "Est2.L", "Est2.U", "Est5.L", "Est5.U", "Est6.L", "Est6.U")
    
  }
  
  Red.res.mat <- Rem.res.mat <- Red.cvg.mat <- Rem.cvg.mat <- NULL
  for(no.sim in 1:n.sim){
    
    Red.res.mat <- cbind(Red.res.mat, out.mclapply[seq(1, by = 8, length.out = no.est), no.sim])
    Rem.res.mat <- cbind(Rem.res.mat, out.mclapply[seq(2, by = 8, length.out = no.est), no.sim])
    Red.cvg.mat <- cbind(Red.cvg.mat, out.mclapply[seq(3, by = 8, length.out = no.est), no.sim])
    Rem.cvg.mat <- cbind(Rem.cvg.mat, out.mclapply[seq(4, by = 8, length.out = no.est), no.sim])
    
    Red.CI.out[no.sim, ] <- out.mclapply[sort(c(seq(5, by = 8, length.out = no.est),
                                                seq(6, by = 8, length.out = no.est))), no.sim]
    Rem.CI.out[no.sim, ] <- out.mclapply[sort(c(seq(7, by = 8, length.out = no.est),
                                                seq(8, by = 8, length.out = no.est))), no.sim]
  }
  
  # Estimate
  Red.out[, 1] <- rowMeans(Red.res.mat, na.rm = TRUE)
  Rem.out[, 1] <- rowMeans(Rem.res.mat, na.rm = TRUE)
  # Coverage
  Red.out[, 2] <- rowMeans(Red.cvg.mat, na.rm = TRUE)
  Rem.out[, 2] <- rowMeans(Rem.cvg.mat, na.rm = TRUE)
  # Bias
  Red.out[- 1, 3] <- Red.true - Red.out[- 1, 1]
  Rem.out[- 1, 3] <- Rem.true - Rem.out[- 1, 1]
  Red.out[1, 3] <- Red.true.est1 - Red.out[1, 1]
  Rem.out[1, 3] <- Rem.true.est1 - Rem.out[1, 1]
  # RMSE
  Red.out[- 1, 4] <- sqrt(rowMeans((Red.res.mat[- 1, ] - Red.true)^2, na.rm = TRUE))
  Rem.out[- 1, 4] <- sqrt(rowMeans((Rem.res.mat[- 1, ] - Rem.true)^2, na.rm = TRUE))
  Red.out[1, 4] <- sqrt(mean((Red.res.mat[1, ] - Red.true.est1)^2, na.rm = TRUE))
  Rem.out[1, 4] <- sqrt(mean((Rem.res.mat[1, ] - Rem.true.est1)^2, na.rm = TRUE))
  # %Bias
  Red.out[- 1, 5] <- Red.out[- 1, 3]/Red.true #* 100
  Rem.out[- 1, 5] <- Rem.out[- 1, 3]/Rem.true #* 100
  Red.out[1, 5] <- Red.out[1, 3]/Red.true.est1 #* 100
  Rem.out[1, 5] <- Rem.out[1, 3]/Rem.true.est1 #* 100
  # %RMSE
  Red.out[- 1, 6] <- Red.out[- 1, 4]/Red.true #* 100
  Rem.out[- 1, 6] <- Rem.out[- 1, 4]/Rem.true #* 100
  Red.out[1, 6] <- Red.out[1, 4]/Red.true.est1 #* 100
  Rem.out[1, 6] <- Rem.out[1, 4]/Rem.true.est1 #* 100
  
  true.value <- matrix(NA, 1, 4)
  colnames(true.value) <- c("disp.red", "disp.rem","disp.red.est1", "disp.rem.est1")
  true.value[1, ] <- c(Red.true, Rem.true, Red.true.est1, Rem.true.est1)
  
  return(list(true.value = true.value,
              disparity.reduction = Red.out, disparity.remaining = Rem.out,
              disparity.reduction.ci = Red.CI.out, disparity.remaining.ci = Rem.CI.out))
  
}
