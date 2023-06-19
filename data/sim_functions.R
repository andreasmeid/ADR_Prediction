# simulation function -----------------------------------------------------
expit <- function(x) { exp(x) / (1 + exp(x)) }
logit <- function(x) { log( x / (1 - x) ) }
range01 <- function(x){(x-min(x))/(max(x)-min(x))}


sim_data <- function(n,          # sample size
                     alpha,      # treatment mean coefficients (vector of length 10)
                     trt_freq,   # frequent treatment (T)
                     pred_freq   # frequent predictors (T)
) {
  
  # sample data 
  # continuous predictors X1:X6
  baseCovars <- mvrnorm(n, numeric(6), diag(6))
  # binary indicators Z1:Z3
  Z1 <- rbinom(n, 1, ifelse(pred_freq,0.5, 0.1)) 
  Z2 <- rbinom(n, 1, ifelse(pred_freq,0.4, 0.05))
  Z3 <- rbinom(n, 1, ifelse(pred_freq,0.6, 0.15))
  
  allCovars <- cbind(baseCovars, Z1, Z2, Z3)
  colnames(allCovars) <- c( paste0("X",1:6), paste0("Z", 1:3) )
  
  # translate to meaningful labels
  baseData <- as.data.frame(baseCovars) 
  baseData$age <- round(qnorm(range01(baseData$V1), mean = 65, sd = 11)) #qnorm(pnorm(X), mean = 72, sd = 8)
  baseData$drugs <- round(qnorm(range01(baseData$V2), mean = 5, sd = 3)) #qnorm(pnorm(X), mean = 72, sd = 8)
  baseData$sex <- ifelse(baseData$V3>median(baseData$V3), 1,0) #qnorm(pnorm(X), mean = 72, sd = 8)
  baseData$prior_hosp <- round(qweibull(range01(baseData$V4), shape = range01(baseData$V4), scale=3)) #qnorm(pnorm(X), mean = 72, sd = 8)
  baseData$score1 <- round(qexp(range01(baseData$V5), rate = range01(baseData$V5)/1.5)) #qnorm(pnorm(X), mean = 72, sd = 8)
  baseData$score2 <- round(qexp(range01(baseData$V6), rate = range01(baseData$V6)/2.5)) #qnorm(pnorm(X), mean = 72, sd = 8)
  baseData$comed1 <- Z1; baseData$comed2 <- Z2; baseData$comed3 <- Z3
  
  
  baseData <- baseData %>% 
    mutate(age = ifelse(age== -Inf, 18, age),
           drugs = ifelse(drugs== -Inf, 0, drugs),
           prior_hosp = ifelse(drugs== -Inf, 0, prior_hosp),
           score1= ifelse(score1>5, 5, score1),
           score2= ifelse(score2<3, 3, score2)) %>% 
    mutate(age = ifelse(is.finite(age)== F, 90, age),
           prior_hosp = ifelse(is.finite(prior_hosp)== F, 15, prior_hosp),
           drugs = ifelse(is.finite(drugs)== F, 15, drugs),
           score1 = ifelse(is.finite(score1)== F, 5, score1),
           score2 = ifelse(is.finite(score2)== F, 7, score2),) %>%
    dplyr::select(age, sex, drugs, prior_hosp, score1, score2, comed1, comed2, comed3)
  
  #PerformanceAnalytics::chart.Correlation(baseData, histogram=TRUE, pch=19)
  
  
  # Treatment model
  pTreat <- expit( cbind(1, allCovars) %*% c(ifelse(trt_freq, -1, round(runif(1, -5,-3))), alpha[-1]) )
  trt <- rbinom(n, 1, pTreat)
  
  ds <- cbind(data.frame(trt, allCovars), baseData)
  
  return( ds )
  
}

#   this function generates one data set, using the parameter values specified in the argument
sim_out <- function(sim_dat,    # simulated data from sim_data function
                    beta       # outcome mean coefficients (vector of length 14)
) {
  
  
  allCovars <- sim_dat %>% dplyr::select(X1:X6, Z1:Z3) %>% as.matrix()
  trt <- as.vector(sim_dat$trt)
  baseData <- sim_dat %>% dplyr::select(-(X1:X6), -(Z1:Z3), -trt)
  
  # outcome model: intercept, treatment, interaction terms
  # no treatment
  Y0mean <- cbind(1, 0, allCovars, 0*allCovars[,"Z1"], 0*allCovars[,"Z2"], 0*allCovars[,"Z3"]) %*% beta
  
  Y0mean <- Y0mean + rnorm(n,0,1)
  p0true=exp(Y0mean)/(1+exp(Y0mean))
  Y0 <-rbinom(n,1,p0true)
  
  # treatment
  Y1mean <- cbind(1, 1, allCovars, 1*allCovars[,"Z1"], 1*allCovars[,"Z2"], 1*allCovars[,"Z3"]) %*% beta
  Y1mean <- Y1mean + rnorm(n,0,1)
  
  p1true=exp(Y1mean)/(1+exp(Y1mean))
  Y1 <-rbinom(n,1,p1true)
  
  # outcome
  Y <- ifelse(trt, Y1, Y0) # observed outcome
  
  trueGrp <- round(Y1mean - Y0mean) # treatment effect group
  
  ds <- cbind(data.frame(trt, Y, allCovars, Y0, Y1, trueGrp), baseData)
  
  return( ds )
}

datagen <- function(n,          # sample size
                    alpha,      # treatment mean coefficients (vector of length 10)
                    beta,        # outcome mean coefficients (vector of length 14)
                    trt_freq, # frequent treatment (T)
                    pred_freq # frequent predictors (T)
) {
  
  # sample data 
  # continuous predictors X1:X6
  baseCovars <- mvrnorm(n, numeric(6), diag(6))
  # binary indicators Z1:Z3
  Z1 <- rbinom(n, 1, ifelse(pred_freq,0.5, 0.1)) 
  Z2 <- rbinom(n, 1, ifelse(pred_freq,0.4, 0.05))
  Z3 <- rbinom(n, 1, ifelse(pred_freq,0.6, 0.15))
  
  allCovars <- cbind(baseCovars, Z1, Z2, Z3)
  colnames(allCovars) <- c( paste0("X",1:6), paste0("Z", 1:3) )
  
  # translate to meaningful labels
  baseData <- as.data.frame(baseCovars) 
  baseData$age <- round(qnorm(range01(baseData$V1), mean = 65, sd = 11)) #qnorm(pnorm(X), mean = 72, sd = 8)
  baseData$drugs <- round(qnorm(range01(baseData$V2), mean = 5, sd = 3)) #qnorm(pnorm(X), mean = 72, sd = 8)
  baseData$sex <- ifelse(baseData$V3>median(baseData$V3), 1,0) #qnorm(pnorm(X), mean = 72, sd = 8)
  baseData$prior_hosp <- round(qweibull(range01(baseData$V4), shape = range01(baseData$V4), scale=3)) #qnorm(pnorm(X), mean = 72, sd = 8)
  baseData$score1 <- round(qexp(range01(baseData$V5), rate = range01(baseData$V5)/1.5)) #qnorm(pnorm(X), mean = 72, sd = 8)
  baseData$score2 <- round(qexp(range01(baseData$V6), rate = range01(baseData$V6)/2.5)) #qnorm(pnorm(X), mean = 72, sd = 8)
  baseData$comed1 <- Z1; baseData$comed2 <- Z2; baseData$comed3 <- Z3
  
  
  baseData <- baseData %>% 
    mutate(age = ifelse(age== -Inf, 18, age),
           drugs = ifelse(drugs== -Inf, 0, drugs),
           prior_hosp = ifelse(drugs== -Inf, 0, prior_hosp),
           score1= ifelse(score1>5, 5, score1),
           score2= ifelse(score2<3, 3, score2)) %>% 
    mutate(age = ifelse(is.finite(age)== F, 90, age),
           prior_hosp = ifelse(is.finite(prior_hosp)== F, 15, prior_hosp),
           drugs = ifelse(is.finite(drugs)== F, 15, drugs),
           score1 = ifelse(is.finite(score1)== F, 5, score1),
           score2 = ifelse(is.finite(score2)== F, 7, score2),) %>%
    dplyr::select(age, sex, drugs, prior_hosp, score1, score2, comed1, comed2, comed3)
  
  #PerformanceAnalytics::chart.Correlation(baseData, histogram=TRUE, pch=19)
  
  
  # Treatment model
  pTreat <- expit( cbind(1, allCovars) %*% c(ifelse(trt_freq, -1, round(runif(1, -5,-3))), alpha[-1]) )
  trt <- rbinom(n, 1, pTreat)
  
  
  
  # outcome model: intercept, treatment, interaction terms
  # no treatment
  Y0mean <- cbind(1, 0, allCovars, 0*allCovars[,"Z1"], 0*allCovars[,"Z2"], 0*allCovars[,"Z3"]) %*% beta
  
  Y0mean <- Y0mean + rnorm(n,0,1)
  p0true=exp(Y0mean)/(1+exp(Y0mean))
  Y0 <-rbinom(n,1,p0true)
  
  # treatment
  Y1mean <- cbind(1, 1, allCovars, 1*allCovars[,"Z1"], 1*allCovars[,"Z2"], 1*allCovars[,"Z3"]) %*% beta
  Y1mean <- Y1mean + rnorm(n,0,1)
  
  p1true=exp(Y1mean)/(1+exp(Y1mean))
  Y1 <-rbinom(n,1,p1true)
  
  # outcome
  Y <- ifelse(trt, Y1, Y0) # observed outcome
  
  trueGrp <- round(Y1mean - Y0mean) # treatment effect group
  
  ds <- cbind(data.frame(trt, Y, allCovars, Y0, Y1, trueGrp), baseData)
  
  return( ds )
}


# evaluation functions ----------------------------------------------------


formel_full <- "Y ~ trt + age + sex + drugs + prior_hosp + score1 + score2 + comed1 + comed2 + comed3 + comed1:trt + comed2:trt + comed3:trt"
preds <- c("trt", "age", "sex", "drugs", "prior_hosp", "score1", "score2", "comed1", "comed2", "comed3", "comed1:trt", "comed2:trt", "comed3:trt")

calc_c <- function(dat, formel="Y ~ trt + age + sex + drugs + prior_hosp + score1 + score2 + comed1 + comed2 + comed3 + comed1:trt + comed2:trt + comed3:trt",
                   boot_n=NULL, boot_from=NULL, old_model=NULL){ # sample size from bootstrapped population drawn from population x
  
  # outcome
  outcome_string <- gsub("\\(|\\)", "", as.formula(formel)[2])# substring(as.formula(formel)[2],1,1); 
  #print(outcome_string)
  
  # split data
  if (is.null(boot_n)) {
    smp_size <- floor(10/12 * nrow(dat))
  } else {
    smp_size <- floor(10/12 * nrow(dat)/4)
  }
  
  # set the seed to make your partition reproducible
  set.seed(saat)
  train_ind <- sample(seq_len(nrow(dat)), size = smp_size)
  
  train <- dat[train_ind, ]
  if (is.null(boot_from)) {
    test <- dat[-train_ind, ]
  } else {
    
    if (is.null(boot_n)) {
      n_new <- nrow(dat) - length(train_ind)
    } else {
      n_new <- boot_n
    }
    bootdat <- dat[-train_ind, ] %>% filter(situation %in% boot_from)
    set.seed(saat)
    boot_ind <- sample(seq_len(nrow(bootdat)), size = n_new, replace = T)
    test <- bootdat[boot_ind, ]
  }
  
  
  # fit model  
  #if (is.null(old_model)) {
    mod <- glm(as.formula(formel), family=binomial(), data=train)
  #} else {
  #  mod <- old_model
  #}
  
  # predictions  
  in_sample_predictions <- predict(mod, type = "response")
  out_of_sample_predictions <- predict(mod, newdata=test, type = "response")
  
  # ROC
  set.seed(saat)
  eval(parse(text=(paste("in_sample_roc <- roc(train$",outcome_string,", in_sample_predictions)",sep=""))))
  set.seed(saat)
  eval(parse(text=(paste("out_of_sample_roc <- roc(test$",outcome_string,", out_of_sample_predictions)",sep=""))))
  
  
  # c-statistics 
  in_sample_c <- as.numeric(in_sample_roc$auc)
  out_of_sample_c <- as.numeric(out_of_sample_roc$auc)
  in_sample_c_sd <- sqrt(var(in_sample_roc))
  out_of_sample_c_sd <- sqrt(var(out_of_sample_roc))
  
  #save(mod, file="data/old_model.RData")

  return(list=c(in_sample_c=in_sample_c, out_of_sample_c=out_of_sample_c,
                in_sample_c_sd=in_sample_c_sd, out_of_sample_c_sd=out_of_sample_c_sd))
  
}



draw_ROC <- function(dat, formel="Y ~ trt + age + sex + drugs + prior_hosp + score1 + score2 + comed1 + comed2 + comed3 + comed1:trt + comed2:trt + comed3:trt",
                             boot_n=NULL, boot_from=NULL, old_model=NULL){ # sample size from bootstrapped population drawn from population x
  
  # outcome
  outcome_string <- gsub("\\(|\\)", "", as.formula(formel)[2])# substring(as.formula(formel)[2],1,1); 
  #print(outcome_string)
  
  # split data
  if (is.null(boot_n)) {
    smp_size <- floor(10/12 * nrow(dat))
  } else {
    smp_size <- floor(10/12 * nrow(dat)/4)
  }
  
  
  # set the seed to make your partition reproducible
  set.seed(saat)
  train_ind <- sample(seq_len(nrow(dat)), size = smp_size)
  
  train <- dat[train_ind, ]
  if (is.null(boot_from)) {
    test <- dat[-train_ind, ]
  } else {
    
    if (is.null(boot_n)) {
      n_new <- nrow(dat) - length(train_ind)
    } else {
      n_new <- boot_n
    }
    bootdat <- dat[-train_ind, ] %>% filter(situation %in% boot_from)
    set.seed(saat)
    boot_ind <- sample(seq_len(nrow(bootdat)), size = n_new, replace = T)
    test <- bootdat[boot_ind, ]
  }
  
  
  # fit model  
  #if (is.null(old_model)) {
    mod <- glm(as.formula(formel), family=binomial(), data=train)
  #} else {
  #  mod <- old_model
  #}
  
  
  # predictions  
  in_sample_predictions <- predict(mod, type = "response")
  out_of_sample_predictions <- predict(mod, newdata=test, type = "response")
  
  # ROC
  set.seed(saat)
  eval(parse(text=(paste("in_sample_roc <- roc(train$",outcome_string,", in_sample_predictions)",sep=""))))
  set.seed(saat)
  eval(parse(text=(paste("out_of_sample_roc <- roc(test$",outcome_string,", out_of_sample_predictions)",sep=""))))
  
  out_of_sample_curve <- ggroc(out_of_sample_roc, alpha = 0.5, colour = "red", linetype = 1, size = 1.5)
  
  return(out_of_sample_curve)
}

  