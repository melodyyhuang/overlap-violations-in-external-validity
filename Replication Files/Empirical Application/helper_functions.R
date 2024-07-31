orv<-function(ipw, sigma_w2, C_sigma){
    a = ipw^2/sigma_w2
    return(sqrt(a)/(1+sqrt(a)))
}
var_w<-function(weights, Y){
    return(mean(weights*(Y-mean(Y))^2))
}

cov_w<-function(weights, A, B){
    return(mean(weights*(A-mean(A))*(B-mean(B))))
}
estimate_bias<-function(R2, C_sigma, p, sigma_w2){
    sqrt(R2/(1-R2)*(1+C_sigma*p/(1-p))*p*sigma_w2)
}
estimate_tau_hat<-function(sample_data, population_data, outcome, covariates, omit){
    sample_data$outcome = sample_data[[outcome]]
    population_data$outcome = population_data[[outcome]]

    if(omit %in% c("rural", "urban", "ruralxwealth_low", "ruralxeduc_below9", "ruralxage_under18",
        "urbanxwealth_low", "urbanxeduc_below9", "urbanxage_under18")){
        covariates = covariates[-which(covariates == "urban")]
        covariates = covariates[-which(covariates == "district")]
    }

    if(omit == "female"){
        covariates = covariates[-which(covariates == "female")]
    }


    if(omit %in% c("wealthindex_high", "wealthindex_low")){
        covariates = covariates[-which(covariates == "district")]
    }


    #Check for empty cells:
    excl = sapply(covariates, function(x){
            check_empty = table(population_data$omit, population_data[[x]])
            if(any(check_empty==0)){
                return(TRUE)
            }else{
                return(FALSE)
            }
        })

    covariates = covariates[excl]

    f = as.formula(paste(outcome, "~", paste(covariates, collapse ="+")))

    model_1 = lm(f, data = sample_data[sample_data$assigned==1,])

    model_0 = lm(f, data = sample_data[sample_data$assigned==0,])

    Y1_sample = predict(model_1, newdata = sample_data)
    Y0_sample = predict(model_0, newdata = sample_data)

    Y1_pop = predict(model_1, newdata = population_data[population_data$yop == 0,])
    Y0_pop = predict(model_0, newdata = population_data[population_data$yop == 0,])

    return(list(tau_sample = Y1_sample - Y0_sample,
        tau_pop = Y1_pop - Y0_pop))
}

benchmark<-function(covariate, outcome, df_sample, survey_data, sigma_w2, C_sigma = 1,
    tau_pop = NULL,
    type='weighted'){
    p = mean(survey_data[[covariate]])

    if(type == 'augmented'){
        df_sample$omit_ind =df_sample[[covariate]]
        survey_data$omit_ind =survey_data[[covariate]]
        survey_data_sub = survey_data %>% filter(yop == 0)
        df_sample_sub = df_sample %>% filter(omit_ind == 0)

        tau_hat_omit = estimate_tau_hat(df_sample_sub, survey_data, outcome,
            weighting_covariates, covariate)

        rho = cor(tau_pop - tau_hat_omit$tau_pop, survey_data_sub$omit_ind)
        R2 = rho^2
    }else{
        rho = (cov_w(weights[df_sample$assigned==1],
        df_sample[[outcome]][df_sample$assigned==1],
        df_sample[[covariate]][df_sample$assigned==1])-
        cov_w(weights[df_sample$assigned==0],
        df_sample[[outcome]][df_sample$assigned==0],
        df_sample[[covariate]][df_sample$assigned==0]))/
        (sd(survey_data[[covariate]][survey_data$yop==0])*sqrt(sigma_w2))
        R2 = rho^2#(1+R2^2)
    }


    return(data.frame(R2 = R2, p, rho = rho,
        bias = estimate_bias(R2, 1, p, sigma_w2)))
}


#bootstrap code to obtain standard errors:
run_iter<-function(i, df_sample_orig, survey_data_orig){
    pop_data = survey_data_orig %>% filter(yop == 0)
    sample_data = survey_data_orig %>% filter(yop == 1)
    pop_data = pop_data[sample(1:nrow(pop_data), nrow(pop_data), replace = TRUE),]

    select_index = sample(1:nrow(df_sample_orig), nrow(df_sample_orig), replace = TRUE)
    df_sample = df_sample_orig[select_index,]
    sample_data = sample_data[select_index,]
    survey_data = rbind(pop_data, sample_data)

    model_ps = WeightIt::weightit((1-yop) ~ female + age + urban  + live_together +
    education + wealthindex + hhsize+ district,  method = 'ebal',
    data = survey_data,
    estimand="ATT")

    survey_data$weights = model_ps$weights
    df_sample$weights = model_ps$weights[survey_data$yop==1]

    weights = model_ps$weights[survey_data$yop == 1]

    model_ipw_hrs = estimatr::lm_robust(training_hours_e~assigned,
    data= df_sample, weights=weights,
        cluster = df_sample$groupid)

    sigma_w2_hrs = var_w(weights[df_sample$assigned==1], df_sample$training_hours_e[df_sample$assigned==1])-
               var_w(weights[df_sample$assigned==0], df_sample$training_hours_e[df_sample$assigned==0])

    ipw_hrs = coef(model_ipw_hrs)[2]

    ORV_hrs =  orv(ipw_hrs, sigma_w2_hrs, 1)

    benchmark_covariates_all = c("educ_below9", "educ_above12", "wealthindex_low","wealthindex_high",
         "rural", "urban", "female", "age_under18", "age_18to25", "age_25to35", "age_36to50", "age_50plus",
         "hhsize_high", "hhsize_low", "ruralxwealth_low", "ruralxeduc_below9", "ruralxage_under18", "urbanxwealth_low", 
         "urbanxeduc_below9", "urbanxage_under18")
    df_benchmark_all = data.frame(
        covariate =benchmark_covariates_all,
        lapply(benchmark_covariates_all, benchmark, outcome = "training_hours_e", df_sample, survey_data, sigma_w2_hrs)%>% bind_rows())

    df_benchmark_all$MROB = ipw_hrs/df_benchmark_all$bias

    model_ipw_earnings = estimatr::lm_robust(profits4w_real_p99_e~assigned,
    data= df_sample, weights=weights,
        cluster = df_sample$groupid)

    sigma_w2_earnings = var_w(weights[df_sample$assigned==1], df_sample$profits4w_real_p99_e[df_sample$assigned==1])-
               var_w(weights[df_sample$assigned==0], df_sample$profits4w_real_p99_e[df_sample$assigned==0])

    ipw_earnings = coef(model_ipw_earnings)[2]

    ORV_earnings =  orv(ipw_earnings, sigma_w2_earnings, 1)

    df_benchmark_earnings = data.frame(
        covariate =benchmark_covariates_all,
        lapply(benchmark_covariates_all, benchmark, outcome = "profits4w_real_p99_e",
            df_sample, survey_data, sigma_w2_earnings)%>% bind_rows())


    df_benchmark_earnings$MROB = ipw_earnings/df_benchmark_earnings$bias
    return(list(est_hrs = ipw_hrs,
        ORV_hrs = ORV_hrs, benchmark_hrs = df_benchmark_all,
        est_earnings = ipw_earnings,
        ORV_earnings = ORV_earnings,
        benchmark_earnings = df_benchmark_earnings))

}

run_iter_aug<-function(i, df_sample_orig, survey_data_orig){
    pop_data = survey_data_orig %>% filter(yop == 0)
    sample_data = survey_data_orig %>% filter(yop == 1)
    pop_data = pop_data[sample(1:nrow(pop_data), nrow(pop_data), replace = TRUE),]

    select_index = sample(1:nrow(df_sample_orig), nrow(df_sample_orig), replace = TRUE)
    df_sample = df_sample_orig[select_index,]
    sample_data = sample_data[select_index,]
    survey_data = rbind(pop_data, sample_data)

    model_ps = WeightIt::weightit((1-yop) ~ female + age + urban  + live_together +
    education + wealthindex + hhsize+ district,  method = 'ebal',
    data = survey_data,
    estimand="ATT")

    survey_data$weights = model_ps$weights
    df_sample$weights = model_ps$weights[survey_data$yop==1]

    weights = model_ps$weights[survey_data$yop == 1]

    model_ipw_hrs = estimatr::lm_robust(training_hours_e~assigned,
    data= df_sample, weights=weights,
        cluster = df_sample$groupid)

    ipw_hrs = coef(model_ipw_hrs)[2]

    model_hrs_1 = lm(training_hours_e ~ female + age+ urban + live_together+education+
                                        wealthindex+hhsize+district, data = df_sample[df_sample$assigned==1,])

    model_hrs_0 = lm(training_hours_e ~ female + age+ urban + live_together+education+
                                        wealthindex+hhsize+district, data = df_sample[df_sample$assigned==0,])

    Y1_sample = predict(model_hrs_1, newdata = df_sample)
    Y0_sample = predict(model_hrs_0, newdata = df_sample)

    Y1_pop = predict(model_hrs_1, newdata = survey_data[survey_data$yop == 0,])
    Y0_pop = predict(model_hrs_0, newdata = survey_data[survey_data$yop == 0,])

    aug_ipw_hrs = ipw_hrs - mean(weights*(Y1_sample - Y0_sample))/mean(weights) + mean(Y1_pop - Y0_pop)

    sigma_w2_hrs = var_w(weights[df_sample$assigned==1], (df_sample$training_hours_e - Y1_sample)[df_sample$assigned==1])-
               var_w(weights[df_sample$assigned==0], (df_sample$training_hours_e - Y0_sample)[df_sample$assigned==0])

    ORV_hrs =  orv(aug_ipw_hrs, sigma_w2_hrs, 1)

    df_sample$training_hours_resid = df_sample$training_hours_e - ifelse(df_sample$assigned==1, Y1_sample, Y0_sample)

    benchmark_covariates_all = c("educ_below9", "educ_above12", "wealthindex_low","wealthindex_high",
         "rural", "urban", "female", "age_under18", "age_18to25", "age_25to35", "age_36to50", "age_50plus",
         "hhsize_high", "hhsize_low", "ruralxwealth_low", "ruralxeduc_below9", "ruralxage_under18", "urbanxwealth_low", "urbanxeduc_below9", "urbanxage_under18")
    df_benchmark_all = data.frame(
    covariate =benchmark_covariates_all,
    lapply(benchmark_covariates_all, benchmark,
        outcome = "training_hours_e", df_sample, survey_data, sigma_w2_aug_hrs, tau_pop = Y1_pop-Y0_pop,
        type = 'augmented')%>% bind_rows())

    df_benchmark_all$MROB = aug_ipw_hrs/df_benchmark_all$bias

    model_ipw_earnings = estimatr::lm_robust(profits4w_real_p99_e~assigned,
    data= df_sample, weights=weights,
        cluster = df_sample$groupid)

    ipw_earnings = coef(model_ipw_earnings)[2]

    #Set up aug model:
    model_earnings_1 = lm(profits4w_real_p99_e ~ female + age+ urban + live_together+education+
                                        wealthindex+hhsize, data = df_sample[df_sample$assigned==1,])

    model_earnings_0 = lm(profits4w_real_p99_e ~ female + age+ urban + live_together+education+
                                        wealthindex+hhsize, data = df_sample[df_sample$assigned==0,])

    Y1_earnings_sample = predict(model_earnings_1, newdata = df_sample)
    Y0_earnings_sample = predict(model_earnings_0, newdata = df_sample)

    Y1_earnings_pop = predict(model_earnings_1, newdata = survey_data[survey_data$yop == 0,])
    Y0_earnings_pop = predict(model_earnings_0, newdata = survey_data[survey_data$yop == 0,])

    aug_ipw_earnings = ipw_earnings - mean(weights*(Y1_earnings_sample - Y0_earnings_sample))/mean(weights) +
    mean(Y1_earnings_pop - Y0_earnings_pop)

    df_sample$earnings_resid = df_sample$profits4w_real_p99_e - ifelse(df_sample$assigned==1, Y1_earnings_sample, Y0_earnings_sample)

    sigma_w2_earnings = var_w(weights[df_sample$assigned==1], (df_sample$earnings_resid)[df_sample$assigned==1])-
                   var_w(weights[df_sample$assigned==0], (df_sample$earnings_resid)[df_sample$assigned==0])

    ORV_earnings =  orv(ipw_earnings, sigma_w2_earnings, 1)

    df_benchmark_earnings = data.frame(
        covariate =benchmark_covariates_all,
        lapply(benchmark_covariates_all, benchmark, outcome = "profits4w_real_p99_e",
        df_sample, survey_data, sigma_w2_aug_earnings, tau_pop = Y1_earnings_pop-Y0_earnings_pop,
        type = 'augmented')%>% bind_rows())


    df_benchmark_earnings$MROB = ipw_earnings/df_benchmark_earnings$bias
    return(list(est_hrs = aug_ipw_hrs,
        ORV_hrs = ORV_hrs, benchmark_hrs = df_benchmark_all,
        est_earnings = aug_ipw_earnings,
        ORV_earnings = ORV_earnings,
        benchmark_earnings = df_benchmark_earnings))

}