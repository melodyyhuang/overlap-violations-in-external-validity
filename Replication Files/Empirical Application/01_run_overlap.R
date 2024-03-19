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

setwd('~/Desktop/Overlap/Uganda')

load('Data/uganda_cleaned.Rdata')

drop_index = which(df_sample$ind_found_e1 == 0 | is.na(df_sample$training_hours_e))

survey_data = survey_data[-which(survey_data$partid %in% df_sample$partid[drop_index]),]
df_sample = df_sample[-drop_index,]
df_sample$wealthindex = survey_data$wealthindex[survey_data$yop==1]

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
sigma_hrs = var(df_sample$training_hours_e[df_sample$assigned==1], na.rm=TRUE)-
    var(df_sample$training_hours_e[df_sample$assigned==0], na.rm=TRUE)
ipw_hrs = coef(model_ipw_hrs)[2]
orv<-function(ipw, sigma_w2, C_sigma){
    a = ipw^2/sigma_w2
    return(sqrt(a)/(1+sqrt(a)))
}
ORV_hrs = orv(ipw_hrs, sigma_w2_hrs, 1)

#Code variables for benchmarking
#Indicator 1 = Dropped from Sample
which.max(scale(df_sample$wealthindex) <= qnorm(0.05))

survey_data$educ_below9 = ifelse(survey_data$education <= 9, 1, 0)
survey_data$educ_above12 = ifelse(survey_data$education >=12, 1, 0)
survey_data$wealthindex_low<-ifelse(survey_data$wealthindex < -0.64, 1, 0)
survey_data$wealthindex_high<-ifelse(survey_data$wealthindex > 0.34, 1, 0)
survey_data$rural<-1-survey_data$urban
survey_data$hhsize_high <- ifelse(survey_data$hhsize >= 10, 1, 0)
survey_data$hhsize_low <- ifelse(survey_data$hhsize <= 3, 1, 0)
survey_data$age_under18 <- ifelse(survey_data$age <= 18, 1, 0)
survey_data$age_18to25 <- ifelse(survey_data$age > 18 & survey_data$age <= 25, 1, 0)
survey_data$age_25to35 <- ifelse(survey_data$age > 25 & survey_data$age <= 35, 1, 0)
survey_data$age_36to50 <- ifelse(survey_data$age > 36 & survey_data$age <= 50, 1, 0)
survey_data$age_50plus <- ifelse(survey_data$age >= 50, 1, 0)
survey_data$ruralxwealth_low = survey_data$rural*survey_data$wealthindex_low
survey_data$ruralxeduc_below9 = survey_data$rural*survey_data$educ_below9
survey_data$ruralxage_under18 = survey_data$rural*survey_data$age_under18
survey_data$urbanxwealth_low = survey_data$urban*survey_data$wealthindex_low
survey_data$urbanxeduc_below9 = survey_data$urban*survey_data$educ_below9
survey_data$urbanxage_under18 = survey_data$urban*survey_data$age_under18


df_sample$educ_below9 = ifelse(df_sample$education <= 9, 1, 0)
df_sample$educ_above12 = ifelse(df_sample$education >=12, 1, 0)
df_sample$wealthindex_low = ifelse(df_sample$wealthindex <= -0.64, 1, 0)
df_sample$wealthindex_high = ifelse(df_sample$wealthindex > 0.34, 1, 0)
df_sample$rural = 1-df_sample$urban
df_sample$hhsize_high = ifelse(df_sample$hhsize>=10, 1, 0)
df_sample$hhsize_low = ifelse(df_sample$hhsize<=3, 1, 0)
df_sample$age_under18 <- ifelse(df_sample$age <= 18, 1, 0)
df_sample$age_18to25 <- ifelse(df_sample$age > 18 & df_sample$age <= 25, 1, 0)
df_sample$age_25to35 <- ifelse(df_sample$age > 25 & df_sample$age <= 35, 1, 0)
df_sample$age_36to50 <- ifelse(df_sample$age > 36 & df_sample$age <= 50, 1, 0)
df_sample$age_50plus <- ifelse(df_sample$age >= 50, 1, 0)
df_sample$ruralxwealth_low = df_sample$rural*df_sample$wealthindex_low
df_sample$ruralxeduc_below9 = df_sample$rural*df_sample$educ_below9
df_sample$ruralxage_under18 = df_sample$rural*df_sample$age_under18
df_sample$urbanxwealth_low = df_sample$urban*df_sample$wealthindex_low
df_sample$urbanxeduc_below9 = df_sample$urban*df_sample$educ_below9
df_sample$urbanxage_under18 = df_sample$urban*df_sample$age_under18

cov_w<-function(weights, A, B){
    return(mean(weights*(A-mean(A))*(B-mean(B))))
}
estimate_bias<-function(R2, C_sigma, p, sigma_w2){
    sqrt(R2/(1-R2)*(1+C_sigma*p/(1-p))*p*sigma_w2)
}
benchmark<-function(covariate, outcome, df_sample, survey_data, sigma_w2, C_sigma = 1){
    p = mean(survey_data[[covariate]])

    rho = (cov_w(weights[df_sample$assigned==1],
        df_sample[[outcome]][df_sample$assigned==1],
        df_sample[[covariate]][df_sample$assigned==1])-
    cov_w(weights[df_sample$assigned==0],
        df_sample[[outcome]][df_sample$assigned==0],
        df_sample[[covariate]][df_sample$assigned==0]))/
    (sd(survey_data[[covariate]][survey_data$yop==0])*sqrt(sigma_w2))
    R2 = rho^2#(1+R2^2)


    return(data.frame(R2 = R2, p, rho = rho,
        bias = estimate_bias(R2, 1, p, sigma_w2)))
}
benchmark_covariates = c("educ_below9", "educ_above12", "wealthindex_low",
     "rural", "female", "age_under18", "hhsize_high")
cleaned_covariate_names = c("Educ. (< 9)", "Educ. (> 12)", "Wealth (Low)",
     "Rural", "Female", "Age (< 18)", "Household Size (Large)")
df_benchmark = data.frame(
    covariate = cleaned_covariate_names,
    lapply(benchmark_covariates, benchmark, outcome = "training_hours_e", df_sample, survey_data, sigma_w2_hrs)%>% bind_rows())

# df_benchmark = rbind(
#     data.frame(covariate = "Educ. (< 9)", benchmark('educ_below9', 'training_hours_e', df_sample, survey_data, sigma_w2_hrs)),
#     data.frame(covariate = "Wealth (Low)", benchmark('wealthindex_low', 'training_hours_e',  df_sample, survey_data, sigma_w2_hrs)),
#     data.frame(covariate = 'Rural', benchmark('rural', 'training_hours_e',  df_sample, survey_data, sigma_w2_hrs)),
#     data.frame(covariate = "Household Size (Large)", benchmark('hhsize_high', 'training_hours_e',  df_sample, survey_data, sigma_w2_hrs)),
#     data.frame(covariate="Female", benchmark('female', 'training_hours_e',  df_sample, survey_data, sigma_w2_hrs)))

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
ggsave("overlap_bias_contour_hrs.pdf", height=5, width=5.5)

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
ggsave("overlap_bias_contour_profit.pdf", height=5, width=5.5)
#------------------------------------------------------------------------------------------------
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
    lapply(benchmark_covariates_all, benchmark, outcome = "training_hours_e", df_sample, survey_data, sigma_w2_hrs)%>% bind_rows())
df_benchmark_all$MROB = ipw_hrs/df_benchmark_all$bias
df_benchmark_all$sign = ifelse(sign(df_benchmark_all$rho)==1, "+", "-")

df_benchmark_all_earnings = data.frame(
    covariate =benchmark_covariates_all,
    lapply(benchmark_covariates_all, benchmark, outcome = "profits4w_real_p99_e",
        df_sample, survey_data, sigma_w2_earnings)%>% bind_rows())
df_benchmark_all_earnings$MROB = ipw_earnings/df_benchmark_all_earnings$bias
df_benchmark_all_earnings$sign = ifelse(sign(df_benchmark_all_earnings$rho)==1, "+", "-")

##TABLE 2
print(xtable::xtable(cbind(df_benchmark_all %>% select(-rho),
    df_benchmark_all_earnings %>% select(-rho, -covariate))), include.rownames=FALSE)
