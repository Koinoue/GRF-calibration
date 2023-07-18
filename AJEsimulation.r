    #############################################################################
    set.seed(123)
    ##Data generation from simulation code by Neal Jawadekar et al. 
    #https://github.com/njawadekar/Simulation-HCF/blob/main/Simulation_HCF.Rmd
    nIndividuals<-10000
    
    # We simulate the individual variables in each dataset
    pid=seq(1, by=1, len=nIndividuals) # this creates a sequential list of "pid" from 1 to nIndividuals   
    A =rbinom(n = nIndividuals, size = 1, prob = 0.5)  # creates a randomly allocated treatment variable, A (each individual has a 50% probability of treatment)
    
    # Here, we used "correlation=low" pattern in their original code
     # devtools::install_github("debruine/faux")
    dat1 <-   rnorm_multi(n = nIndividuals, mu = c(20,23), sd = c(6,9),r = c(0.05), varnames = c("B_cont","C_cont"), empirical = FALSE) 
    dat2 <-   rnorm_multi(n = nIndividuals, mu = c(25,20), sd = c(14,11),r = c(0.15), varnames = c("N","O"), empirical = FALSE)
    
    # Simulate dichotomous variables (B through K)	 
    B = ifelse(dat1$B_cont > median(dat1$B_cont),1,0)
    C = ifelse(dat1$C_cont > median(dat1$C_cont),1,0)
    chisq.test(B, C, correct=FALSE) # This test shows that B and C are associated        
    D = rbinom(n = nIndividuals, size = 1, prob = 0.3)
    E = rbinom(n = nIndividuals, size = 1, prob = 0.7)
    F = rbinom(n = nIndividuals, size = 1, prob = 0.13)
    G = rbinom(n = nIndividuals, size = 1, prob = 0.25)
    H = rbinom(n = nIndividuals, size = 1, prob = 0.30)
    I = rbinom(n = nIndividuals, size = 1, prob = 0.08)
    J = rbinom(n = nIndividuals, size = 1, prob = 0.15)
    K = rbinom(n = nIndividuals, size = 1, prob = 0.25)
    
    # Simulate continuous variables (L through U)
    L = rnorm(n = nIndividuals, mean = 30, sd = 12)
    M = rnorm(n = nIndividuals, mean = 15, sd = 0.1) 
    # Create 2 correlated variables, N and O
    N = dat2$N
    O = dat2$O
    cor(dat2, method = "pearson")[1,2] # This test shows that N and O are correlated
    P = rnorm(n = nIndividuals, mean = 120, sd = 4)
    Q = rnorm(n = nIndividuals, mean = 72, sd = 2.5)
    R = rnorm(n = nIndividuals, mean = 5, sd = 0.5)
    S = rnorm(n = nIndividuals, mean = 22, sd = 2)
    T = rnorm(n = nIndividuals, mean = 50, sd = 3)
    U = rnorm(n = nIndividuals, mean = 100, sd = 10)
    
    # Specify the conditional probabilities of Y within specific strata of treatment/covariate combinations. 
    Yprob <- numeric(nIndividuals) # this creates an empty vector which is used to assign probability of Y for each individual  
    # We are assigning probabilities of Y based on each individual's values of A, B, and N, as well as the sim_cf_rct function's input value of "cates"   
    # Here, we used the situation with high CATE
      Yprob[A == 1 & B == 1 & N >=41]=0.05 
      Yprob[A == 0 & B == 1 & N >= 41]=0.20
      # when B == 1 and N >= 41, CATE of A->Y should be about -0.15
      
      Yprob[A == 1 & B == 0 & N >=41 ]=0.07 
      Yprob[A == 0 & B == 0 & N >= 41]=0.10
      # when B == 0 and N >= 41, CATE of A->Y should be about -0.03
      
      Yprob[A == 1 & B == 0 & N < 41]=0.15
      Yprob[A == 0 & B == 0 & N < 41]=0.10
      # when B == 0 and N < 41, CATE of A->Y should be about 0.05
      
      Yprob[A == 1 & B == 1 & N < 41]=0.15
      Yprob[A == 0 & B == 1 & N < 41]=0.05
      # when B == 1 and N < 41, CATE of A->Y should be about 0.10
    
    # Simulate each individual's outcome as a random draw from the binomial distribution of probabilities assigned above  
    Y=rbinom(n = nIndividuals, size = 1, prob = Yprob) 
    
    # Create a data frame which combines pid, A, Y, and all covariates
    trialdata=data.frame(cbind(pid, A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, Y)) 
    
    #############################################################################
    
    ## Run causal forest model with cross-fitting
    # ref. https://bookdown.org/halflearned/ml-ci-tutorial/hte-i-binary-treatment.html
    # Number of rankings that the predictions will be ranking on 
    # (e.g., 2 for above/below median estimated CATE, 5 for estimated CATE quintiles, etc.)
    library(grf)
    library(sandwich)
    library(lmtest)
    library(ggplot2)
    num.rankings <- 5
    
    # Prepare for data.splitting
    # Assign a fold number to each observation.
    # The argument 'clusters' in the next step will mimick K-fold cross-fitting.
    num.folds <- 10
    
    # Prepare dataset
    X = data.matrix(trialdata[, c(3:22)])  
    Y = trialdata$Y
    W = trialdata$A
    
    #causal forest with honest splitting & 10 folds cross-fitting 
    set.seed(Sys.time())
    set.seed(123) 
    folds <- sample(c(1:num.folds), replace=TRUE, size=nIndividuals)
    forest <- causal_forest(X, Y, W, 
                            honesty = TRUE, 
                            clusters=folds #,
                            #tune.parameters="all"
                            )
    
    predictions <- predict(forest)
    tau.hat <- predictions$predictions
    
    # Rank observations *within each fold* into quintiles according to their CATE predictions.
    ranking <- rep(NA, nIndividuals)
    for (fold in seq(num.folds)) {
      tau.hat.quantiles <- quantile(tau.hat[folds == fold], probs = seq(0, 1, by=1/num.rankings))
      ranking[folds == fold] <- cut(tau.hat[folds == fold], tau.hat.quantiles, include.lowest=TRUE,labels=seq(num.rankings))
    }
    
    # Valid only in randomized settings.
    # Average difference-in-means within each ranking
    
    # Formula y ~ 0 + ranking + ranking:w
    treatment <- "W"
    outcome <- "Y"
    fmla <- paste0(outcome, " ~ 0 + ranking + ranking:", treatment)
    ols.ate <- lm(fmla, data=transform(trialdata, ranking=factor(ranking)))
    ols.ate <- coeftest(ols.ate, vcov=vcovHC(ols.ate, type='HC2'))
    interact <- which(grepl(":", rownames(ols.ate)))
    ols.ate <- data.frame("ols", paste0("Q", seq(num.rankings)), ols.ate[interact, 1:2])
    rownames(ols.ate) <- NULL # just for display
    colnames(ols.ate) <- c("method", "ranking", "estimate", "std.err")
    ols.ate
    
    # Computing AIPW scores.
    tau.hat <- predict(forest)$predictions
    e.hat <- forest$W.hat # P[W=1|X]
    m.hat <- forest$Y.hat # E[Y|X]
    
    # Estimating mu.hat(X, 1) and mu.hat(X, 0) for obs in held-out sample
    # Note: to understand this, read equations 6-8 in this vignette:
    # https://grf-labs.github.io/grf/articles/muhats.html
    mu.hat.0 <- m.hat - e.hat * tau.hat        # E[Y|X,W=0] = E[Y|X] - e(X)*tau(X)
    mu.hat.1 <- m.hat + (1 - e.hat) * tau.hat  # E[Y|X,W=1] = E[Y|X] + (1 - e(X))*tau(X)
    
    # AIPW scores
    aipw.scores <- tau.hat + W / e.hat * (Y -  mu.hat.1) - (1 - W) / (1 - e.hat) * (Y -  mu.hat.0)
    ols <- lm(aipw.scores ~ 0 + factor(ranking))
    forest.ate <- data.frame("aipw", paste0("Q", seq(num.rankings)), coeftest(ols, vcov=vcovHC(ols, "HC2"))[,1:2])
    colnames(forest.ate) <- c("method", "ranking", "estimate", "std.err")
    rownames(forest.ate) <- NULL # just for display
    forest.ate
    
    # Concatenate the two results.
    res <- rbind(forest.ate, ols.ate)
    
    # Plotting the point estimate of average treatment effect 
    # and 95% confidence intervals around it.
    res$estimate<-res$estimate
    ggplot(res) +
      aes(x = ranking, y = estimate, group=method, color=method) + 
      geom_point(position=position_dodge(0.2)) +
      geom_errorbar(aes(ymin=estimate-2*std.err, ymax=estimate+2*std.err), width=.2, position=position_dodge(0.2)) +
      ylab("") + xlab("") +
      ggtitle("Average CATE within each ranking (as defined by predicted CATE)") +
      theme_minimal() +
      theme(legend.position="bottom", legend.title = element_blank())
    #assessing fit: the slope of the calibration line between predicted ARR and observed ARR.
    test_calibration(forest) 
    #hist(data$tau, main="CATE estimates", freq=F)
    #A coefficient of 1 for mean.forest.prediction suggests that the mean forest prediction is correct
    #and a coefficient of 1 for differential.forest.prediction suggests that the forest has captured heterogeneity in the underlying signal.


    #Partial dependence plot
    trialdata$tau.hat = tau.hat
    trialdata %>% 
    ggplot(aes(x = N, y = tau.hat)) +
    geom_point(size = 0.05) +
    theme_bw() +
    ylab("Estimated CATE")

    #Heatmap
    trialdata %>% 
    mutate(N.41 = ifelse(N>=41,1,0),
           N.41 = as.factor(N.41),
           B    = as.factor(B)) %>% 
    group_by(N.41,B) %>% 
    summarise(mean = mean(tau.hat), sd = sd(tau.hat)) %>% 
    mutate(result = paste0("Mean CATE: ",signif(mean,3),"\n","(SD: ",signif(sd,3),")")) %>% 
    ggplot(aes(x = N.41, y = B)) +
    geom_tile(aes(fill = mean)) +
    geom_text(aes(label = result)) +
    scale_fill_gradient2(name = "Mean CATEs") +
    theme_bw() +
    xlab("N ≥ 41")

   #Estimate Rank-Weighted Average Treatment Effect 
    #Please see the original description in the package (https://grf-labs.github.io/grf/reference/rank_average_treatment_effect.html)
    rate <- rank_average_treatment_effect(forest, tau.hat, target = "AUTOC") #Cross-fitting was done by splitting the entire sample into folds and using one fold to make predictions based on the causal forest algorithm trained with the remaining folds.
    #AUTOC and Qini coefficient
    rate
    rate$estimate + data.frame(lower = -1.96 * rate$std.err, upper = 1.96 * rate$std.err, row.names = rate$target)
    # Plot the Targeting Operator Characteristic (TOC) curve.
    plot(rate, xlab="Treated fraction according to CATEs",ylab="Difference in treatment effects", main="Targeting operator characteristics")
    
    #Estimate Rank-Weighted Average Treatment Effect via sample splitting so that treatment prioritization scores 
    #are constructed independently from the evaluation forest training data for valid statistical performance.
    set.seed(123)
    n<-nrow(trialdata)
     train <- sample(1:n, n / 2)
    cf.priority <- causal_forest(X[train, ], Y[train], W[train])
    # Compute a prioritization based on estimated treatment effects.
    priority.cate <-  predict(cf.priority, X[-train, ])$predictions
    # Estimate AUTOC on held out data.
    cf.eval <- causal_forest(X[-train, ], Y[-train], W[-train])
    rate <- rank_average_treatment_effect(cf.eval, priority.cate)
    rate
    # Plot the Targeting Operator Characteristic (TOC) curve.
    plot(rate, xlab="Treated fraction according to CATEs",ylab="Difference in treatment effects", main="Targeting operator characteristics")
