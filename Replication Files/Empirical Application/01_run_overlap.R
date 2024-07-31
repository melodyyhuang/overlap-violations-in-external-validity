rm(list=ls())
library(tidyverse)
#------------------------------------------------------------------------------------------------
ggMelody<-theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 17, face = "bold"),
                                  axis.text=element_text(size=9), #size = 12
                                  #legend.text=element_text(size=7),
                                  legend.position='bottom', axis.title = element_text(size = 12),
                                  strip.text.x = element_text(size = 12, face='bold'),
                                  strip.text.y = element_text(size = 12, face='bold'),
                                  plot.subtitle=element_text(size = 14, hjust=0.5))
theme_set(ggMelody)

setwd("~/Documents/GitHub/overlap-violations-in-external-validity/Replication Files/Empirical Application")
source('helper_functions.R')
load('Data/uganda_cleaned.Rdata')
source("prep_data_files.R")

bootstrap = FALSE
#------------------------------------------------------------------------------------------
model_ps = WeightIt::weightit((1-yop) ~ female + age + urban  + live_together +
    education + wealthindex + hhsize+ district,  method = 'ebal',
    data = survey_data,
    estimand="ATT")

survey_data$weights = model_ps$weights
df_sample$weights = model_ps$weights[survey_data$yop==1]

weights = model_ps$weights[survey_data$yop == 1]
var_w<-function(weights, Y){
    return(mean(weights*(Y-mean(Y))^2))
}
#------------------------------------------------------------------------------------------
#Hours of vocational training:
#------------------------------------------------------------------------------------------
model_dim_hrs = estimatr::lm_robust(training_hours_e~assigned, data = df_sample)
DiM_hrs = coef(model_dim_hrs)[2]


model_ipw_hrs = estimatr::lm_robust(training_hours_e~assigned,
    data= df_sample, weights=weights,
        cluster = df_sample$groupid)

sigma_w2_hrs = var_w(weights[df_sample$assigned==1], df_sample$training_hours_e[df_sample$assigned==1])-
               var_w(weights[df_sample$assigned==0], df_sample$training_hours_e[df_sample$assigned==0])

ipw_hrs = coef(model_ipw_hrs)[2]

ORV_hrs = orv(ipw_hrs, sigma_w2_hrs, 1)

benchmark_covariates = c("educ_below9", "educ_above12", "wealthindex_low",
     "rural", "female", "age_under18", "hhsize_high")
cleaned_covariate_names = c("Educ. (< 9)", "Educ. (> 12)", "Wealth (Low)",
     "Rural", "Female", "Age (< 18)", "Household Size (Large)")
df_benchmark = data.frame(
    covariate = cleaned_covariate_names,
    lapply(benchmark_covariates, benchmark, outcome = "training_hours_e", df_sample, survey_data, sigma_w2_hrs)%>% bind_rows())

df_benchmark$MROB =ipw_hrs/df_benchmark$bias

print(xtable::xtable(df_benchmark), include.rownames=FALSE)


r2_vals <- seq(0, 0.9, by = 0.025)
p <- seq(0, 0.95, by = 0.025)
data.fit <- expand.grid(r2_vals, p)
names(data.fit) = c("r2", "p")
bias = estimate_bias(data.fit$r2, 1, data.fit$p, sigma_w2_hrs)
df_plot = data.frame(R2 = data.fit$r2, p=data.fit$p,
    bias = bias)
breaks=c(25, 50, seq(0, 1300, by=100), 1500, 2000, 3000)

#FIGURE 2 (a)
df_plot %>% ggplot(aes(x = p, y = R2, z = bias)) +
    geom_contour(col = "slategray", breaks=breaks) + ggtitle("") + #xlab(expression(R[epsilon]^2)) +
    #ylab(expression(rho[epsilon * "," * tau])) +
    metR::geom_text_contour(aes(z = bias), breaks=breaks,
    stroke = 0.2, size=3, color = 'slategray')+ metR::geom_contour_fill(breaks = c(ipw_hrs,
            1000 * ipw_hrs), fill = "skyblue3", alpha = 0.25) +
            geom_contour(breaks = c(ipw_hrs), col = "skyblue3",
                size = 1)+
            geom_hline(yintercept=0, color='slategray')+geom_vline(xintercept=0, color = 'slategray')+
            geom_point(data = df_benchmark, aes(x = p, y = R2,
                color = (covariate == "Educ. (< 9)")))+
                ggrepel::geom_label_repel(data = df_benchmark,
                  aes(x = p, y = (R2), color = covariate == "Educ. (< 9)", label = covariate),
                  fill = "white", nudge_y=0.01)+
                ylab(expression(R[tau*"~"*"V"]^2))+
            theme(legend.position = 'none')+
            scale_color_manual(values = c("black", "red4"))+
            ggtitle("Hours of Vocational Training")+
            annotate('point', x = ORV_hrs, y = ORV_hrs, color = 'navy')+
            annotate('text', x = ORV_hrs-0.01, y = ORV_hrs, label = "ORV= 0.37",
                color = 'navy', vjust=1, hjust=1)
#ggsave("overlap_bias_contour_hrs.pdf", height=5, width=5.5)

#------------------------------------------------------------------------------------------
#Earnings:
model_dim_earnings = estimatr::lm_robust(profits4w_real_p99_e~assigned, data = df_sample)
DiM_earnings = coef(model_dim_earnings)[2]

model_ipw_earnings = estimatr::lm_robust(profits4w_real_p99_e~assigned,
    data= df_sample, weights=weights, cluster = df_sample$groupid)

sigma_w2_earnings = var_w(weights[df_sample$assigned==1], df_sample$profits4w_real_p99_e[df_sample$assigned==1])-
    var_w(weights[df_sample$assigned==0], df_sample$profits4w_real_p99_e[df_sample$assigned==0])

ipw_earnings = coef(model_ipw_earnings)[2]

ORV_earnings = orv(ipw_earnings, sigma_w2_earnings, 1)

df_benchmark_earnings = data.frame(
    covariate = cleaned_covariate_names,
    lapply(benchmark_covariates, benchmark,
        outcome = "profits4w_real_p99_e",
        df_sample, survey_data, sigma_w2_earnings)%>% bind_rows())

df_benchmark_earnings$MROB =ipw_earnings/df_benchmark_earnings$bias
# r2_vals <- seq(0, 0.9, by = 0.025)
# p <- seq(0, 0.95, by = 0.025)
# data.fit <- expand.grid(r2_vals, p)
# names(data.fit) = c("r2", "p")
bias_earnings = estimate_bias(data.fit$r2, 1, data.fit$p, sigma_w2_earnings)
df_plot_earnings = data.frame(R2 = data.fit$r2, p=data.fit$p,
    bias = bias_earnings)

breaks_earnings=c(5, seq(10,150, by=10), 200, 300, 400)
df_benchmark_earnings_plot =df_benchmark_earnings %>% filter(!cleaned_covariate_names %in% c("Educ. (> 12)"))
df_plot_earnings %>% ggplot(aes(x = p, y = R2, z = bias)) +
    geom_contour(col = "slategray", breaks=breaks_earnings) + ggtitle("") + #xlab(expression(R[epsilon]^2)) +
    #ylab(expression(rho[epsilon * "," * tau])) +
    metR::geom_text_contour(aes(z = bias), color = 'slategray', breaks=breaks_earnings,
    stroke = 0.2)+ metR::geom_contour_fill(breaks = c(ipw_earnings,
            1000 * ipw_earnings), fill = "darkolivegreen", alpha = 0.25) +
            geom_contour(breaks = c(ipw_earnings), col = "darkolivegreen",
                size = 1)+
            geom_hline(yintercept=0, color='slategray')+geom_vline(xintercept=0, color='slategray')+
            geom_point(data = df_benchmark_earnings_plot, aes(x = p, p = R2,
                color = (covariate %in% c("Educ. (< 9)", "Age (< 18)", "Female"))))+
                ggrepel::geom_label_repel(data = df_benchmark_earnings_plot,
                  aes(x = p, y = R2,
                    color = (covariate %in% c("Educ. (< 9)", "Female", "Age (< 18)")),
                    label = covariate),
                  fill = "white", nudge_y=0.01)+
                ylab(expression(R[tau*"~"*"V"]^2))+
            scale_color_manual(values = c("black", "red4"))+
            theme(legend.position='none')+
            ggtitle("Cash Earnings")+
            annotate('point', x = ORV_earnings, y = ORV_earnings, color = 'darkslategray')+
            annotate('text', x = ORV_earnings-0.01, y = ORV_earnings, label = "ORV= 0.18",
                color = 'darkslategray', vjust=1, hjust=1)
#ggsave("overlap_bias_contour_profit.pdf", height=5, width=5.5)
#------------------------------------------------------------------------------------------------
if(bootstrap == TRUE){
    options(warn=-1)
    set.seed(331)
    df_bs = parallel::mclapply(1:1000, run_iter, df_sample_orig = df_sample, survey_data_orig = survey_data,
        mc.cores = parallel::detectCores()-1)
    options(warn=0)

    est_hrs_bs = lapply(df_bs, function(x){ return(x$est_hrs)}) %>% as.numeric()
    orv_hrs_bs = lapply(df_bs, function(x){ return(x$ORV_hrs)}) %>% as.numeric()
    est_earnings_bs = lapply(df_bs, function(x){ return(x$est_earnings)}) %>% as.numeric()
    orv_earnings_bs = lapply(df_bs, function(x){ return(x$ORV_earnings)}) %>% as.numeric()

    benchmark_hrs_bs = lapply(df_bs, function(x){ return(x$benchmark_hrs)}) %>% bind_rows()
    benchmark_hrs_bs = benchmark_hrs_bs[-which(benchmark_hrs_bs$R2 == Inf),]
    benchmark_hrs_bs = benchmark_hrs_bs[-which(benchmark_hrs_bs$R2 > 1),]

    df_benchmark_hrs_se = benchmark_hrs_bs %>% group_by(covariate) %>%
    summarize(
        R2 = sd(R2, na.rm=TRUE),
        p = sd(p, na.rm=TRUE),
        bias = sd(bias, na.rm=TRUE)
    )

    pct_confounding = benchmark_hrs_bs %>% group_by(covariate) %>%
    summarize(
        pct= mean(MROB < 1, na.rm=TRUE)*100
    )


    benchmark_earnings_bs = lapply(df_bs, function(x){ return(x$benchmark_earnings)}) %>% bind_rows()
    benchmark_earnings_bs = benchmark_earnings_bs[-which(benchmark_earnings_bs$R2 == Inf),]
    benchmark_earnings_bs = benchmark_earnings_bs[-which(benchmark_earnings_bs$R2 > 1),]
    #benchmark_earnings_bs = benchmark_earnings_bs[-which(benchmark_earnings_bs$R2 < 0),]
    df_benchmark_earnings_se = benchmark_earnings_bs %>% group_by(covariate) %>%
    summarize(
        R2 = sd(R2, na.rm=TRUE),
        p = sd(p, na.rm=TRUE),
        bias = sd(bias, na.rm=TRUE)
    )
    pct_confounding_earnings = benchmark_earnings_bs %>% group_by(covariate) %>%
    summarize(
        pct=mean(MROB < 1, na.rm=TRUE)*100
    )
}

##SUMMARY TABLE:
print(xtable::xtable(
    data.frame(DiM = c(model_dim_hrs$coef[2], model_dim_earnings$coef[2]),
               DiM_se = c(model_dim_hrs$std.error[2], model_dim_earnings$std.error[2]),
               IPW = c(model_ipw_hrs$coef[2], model_ipw_earnings$coef[2]),
               IPW_se = c(model_ipw_hrs$std.error[2], model_ipw_earnings$std.error[2]))), include.rownames=FALSE)


##TABLE 1:
print(xtable::xtable(
    data.frame(DiM = c(model_dim_hrs$coef[2], model_dim_earnings$coef[2]),
               DiM_se = c(model_dim_hrs$std.error[2], model_dim_earnings$std.error[2]),
               IPW = c(model_ipw_hrs$coef[2], model_ipw_earnings$coef[2]),
               IPW_se = c(model_ipw_hrs$std.error[2], model_ipw_earnings$std.error[2]),
               RV = c(orv(ipw_hrs, sigma_w2_hrs, 1), orv(ipw_earnings, sigma_w2_earnings, 1)))), include.rownames=FALSE)


#FULL BENCHMARK TABLE:
benchmark_covariates_all = c("educ_below9", "educ_above12", "wealthindex_low","wealthindex_high",
     "rural", "urban", "female", "age_under18", "age_18to25", "age_25to35", "age_36to50", "age_50plus",
     "hhsize_high", "hhsize_low", "ruralxwealth_low", "ruralxeduc_below9", "ruralxage_under18", "urbanxwealth_low", "urbanxeduc_below9", "urbanxage_under18")
df_benchmark_all = data.frame(
    covariate =benchmark_covariates_all,
    lapply(benchmark_covariates_all, benchmark,
        outcome = "training_hours_e", df_sample, survey_data, sigma_w2_hrs)%>% bind_rows())
df_benchmark_all$MROB = ipw_hrs/df_benchmark_all$bias
df_benchmark_all$pct_confounding = pct_confounding$pct
df_benchmark_all$sign = ifelse(sign(df_benchmark_all$rho)==1, "+", "-")


df_benchmark_all_earnings = data.frame(
    covariate =benchmark_covariates_all,
    lapply(benchmark_covariates_all, benchmark, outcome = "profits4w_real_p99_e",
        df_sample, survey_data, sigma_w2_earnings)%>% bind_rows())
df_benchmark_all_earnings$MROB = ipw_earnings/df_benchmark_all_earnings$bias
df_benchmark_all_earnings$pct_confounding = pct_confounding_earnings$pct
df_benchmark_all_earnings$sign = ifelse(sign(df_benchmark_all_earnings$rho)==1, "+", "-")

##TABLE 2
print(xtable::xtable(cbind(df_benchmark_all %>% dplyr::select(-rho, -pct_confounding),
    df_benchmark_all_earnings %>% dplyr::select(-rho, -covariate, -pct_confounding))), include.rownames=FALSE)

#ADDING IN SE:
df_benchmark_all$R2 = paste(round(df_benchmark_all$R2,2), " (", round(df_benchmark_hrs_se$R2,2), ")", sep="")
df_benchmark_all$p = paste(round(df_benchmark_all$p,2), " (", formatC(df_benchmark_hrs_se$p, format='e',digits=1), ")", sep="")
df_benchmark_all$bias = paste(round(df_benchmark_all$bias,2), " (", round(df_benchmark_hrs_se$bias,2), ")", sep="")
df_benchmark_all_earnings$R2 = paste(round(df_benchmark_all_earnings$R2,2), " (", round(df_benchmark_earnings_se$R2,2), ")", sep="")
df_benchmark_all_earnings$p = paste(round(df_benchmark_all_earnings$p,2), " (", formatC(df_benchmark_earnings_se$p, format = 'e', digits=1), ")", sep="")
df_benchmark_all_earnings$bias = paste(round(df_benchmark_all_earnings$bias,2), " (", round(df_benchmark_earnings_se$bias,2), ")", sep="")

print(xtable::xtable(cbind(df_benchmark_all %>% select(-rho),
    df_benchmark_all_earnings %>% select(-rho, -covariate))), include.rownames=FALSE)


