rm(list = ls())
#-------------------------------------------------------------------------------------
#Plotting things:
library(ggplot2)
ggMelody<-theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 17, face = "bold"),
                                  axis.text=element_text(size=9), #size = 12
                                  #legend.text=element_text(size=7),
                                  legend.position='bottom', axis.title = element_text(size = 12),
                                  strip.text.x = element_text(size = 12, face='bold'),
                                  strip.text.y = element_text(size = 12, face='bold'),
                                  plot.subtitle=element_text(size = 14, hjust=0.5))
theme_set(ggMelody)
var_w<-function(weights, Y){
    return(mean(weights*(Y-mean(Y))^2))
}

#N = 10000
library(ggplot2)
library(tidyverse)
library(parallel)
#n=1000
#beta controls for interaction between X and V
generate_data<-function(N=10000, n_sample, is_pop = FALSE, gamma, alpha, beta, rho=0.25, V_binary = TRUE){
    #Selection scores:
    X1 = rnorm(N, 0, 1)
    X2 = rnorm(N, 0, 1)
    prob_X3 = exp(alpha+rho*X1)/(1+exp(alpha+rho*X1))
    X3 = rbinom(N, size = 1, prob = prob_X3)

    #Generate TE:
    u = rnorm(N, 0, 0.25)
    Y0 = 0.5*X1
    tau = 2+X1+X2 + 0.15*(X1^2) - gamma*(1-X3)+(beta*X3*(X1+X2)) + u

    #tau = 2+ X+ gamma*X_overlap
    Y1 = tau+Y0

    if(V_binary){
        V = X3
    }else{
        V = ifelse(tau <= alpha, 1, 0)
    }
    scores =  ifelse(V == 1, 0, plogis(X1+X2))#-0.25*X3
    if(!is_pop) {
      S = rbinom(N, size = 1, prob = scores)
      index = sample(which(S == 1), size=n_sample)
      S[!(1:N %in% index)] = 0
    } else {
      S = rep(0, N)
    }


    return(data.frame(Y0, Y1, X1, X2, X3, V, S, scores))
}


df_pop = generate_data(is_pop=TRUE, gamma = 1, alpha=0.5, beta = 0, rho= 0.25, V_binary = FALSE)
cor(df_pop$V, df_pop$Y1-df_pop$Y0)^2
var((df_pop$Y1-df_pop$Y0)[df_pop$V==0])
mean(df_pop$Y1-df_pop$Y0)
var((df_pop$Y1-df_pop$Y0)[df_pop$V==1])


calculate_orv<-function(ipw, sigma_w2, C_sigma){
    a = ipw^2/sigma_w2
    return(sqrt(a)/(1+sqrt(a)))
}

run_iter<-function(i, n, N, gamma, alpha, beta, rho =0, V_binary = TRUE){
    df_population = generate_data(N, n_sample, is_pop = TRUE, gamma, alpha, beta, rho, V_binary)
    df_sample = generate_data(N*2, n_sample = n, is_pop = FALSE, gamma, alpha, beta, rho, V_binary)

    df_sample = df_sample %>% dplyr::filter(S != 0)
    n_actual = nrow(df_sample)
    Z = rep(0, n_actual)
    Z[sample(1:n_actual, round(n_actual/2))] = 1
    Y = ifelse(Z == 1, df_sample$Y1, df_sample$Y0)

    df_all = rbind(df_population, df_sample)
    model_ps = WeightIt::weightit((1-S)~X1+X2, data = df_all, method='ebal', estimand="ATT")
    # if(length(unique(df_sample$X3)) == 1){
    #     #Implies that there is only one value of X3
    #     model_ps = WeightIt::weightit((1-S)~X1+X2, data = df_all, method='ebal', estimand="ATT")
    # }else{
    #     model_ps = WeightIt::weightit((1-S)~X1+X2+X3, data = df_all, method='ebal', estimand="ATT")
    # }
    weights = model_ps$weights[df_all$S == 1]

    model_ipw = estimatr::lm_robust(Y~Z, weights=weights)
    df_population$tau = df_population$Y1-df_population$Y0

    R2 = cor(df_population$V, df_population$tau)^2
    p = mean(df_population$V)


    df_sample$Z = Z
    df_sample$Y = Y

    bias_decomp = sqrt(R2*var(df_population$tau)*p/(1-p))

    var_tau_bound = var_w(weights[Z==1], Y[Z==1])+var_w(weights[Z==0], Y[Z==0])

    ORV = calculate_orv(coef(model_ipw)[2], var_tau_bound, 1)

    model = lm(tau ~ V, data = df_population)
    C_sigma = var(model$resid[df_population$V==1])/var(model$resid[df_population$V==0])
    bias_decomp_est = sqrt(R2/(1-R2) * (1+C_sigma*p/(1-p))*p*var(df_population$tau[df_population$V==0]))
    bias_decomp_bound = sqrt(R2/(1-R2) * (1+C_sigma*p/(1-p))*p*var_tau_bound)

    true_bias =mean(df_population$tau[df_population$V==0])-mean(df_population$tau)
    ipw_bias = coef(model_ipw)[2] - mean(df_population$tau)


    # #AUGMENTED WEIGHTED ESTIMATOR:
    # model_lm_1 = lm(Y~X1 + X2, data = df_sample[df_sample$Z==1,])
    # model_lm_0 = lm(Y~X1 + X2, data = df_sample[df_sample$Z==0,])

    # df_population_1 = df_population
    # df_population_1$Z = 1
    # df_population_0 = df_population
    # df_population_0$Z = 0


    # Y1_pop = predict(model_lm_1, newdata = df_population_1)
    # Y0_pop = predict(model_lm_0, newdata = df_population_0)
    # tau_pop = Y1_pop-Y0_pop


    # Y1_sample = predict(model_lm_1)
    # Y0_sample = predict(model_lm_0)
    # Yhat_sample = ifelse(df_sample$Z==1, Y1_sample, Y0_sample)

    # model_hat = estimatr::lm_robust(Yhat_sample~Z, weights=weights)
    # aug_ipw = coef(model_ipw)[2] - coef(model_hat)[2] + mean(tau_pop)


    # aug_ipw_bias = aug_ipw - mean(df_population$tau)


    # model_resid = lm(tau ~ V + tau_pop, data = df_population)
    # C_sigma_aug = var(model_resid$resid[df_population$V==1])/var(model_resid$resid[df_population$V==0])
    # R2_aug = cor(df_population$V, df_population$tau - tau_pop)^2
    # p = mean(df_population$V)

    # var_tau_bound_aug = var_w(weights[Z==1], Y[Z==1] - Yhat_sample[Z==1])+var_w(weights[Z==0], Y[Z==0]-Yhat_sample[Z==0])

    # bias_decomp_est_aug = sqrt(R2_aug/(1-R2_aug) * (1+C_sigma_aug*p/(1-p))*p*var(df_population$tau[df_population$V==0]-tau_pop[df_population$V==0]))
    # bias_decomp_bound_aug = sqrt(R2_aug/(1-R2_aug) * (1+C_sigma*p/(1-p))*p*var_tau_bound_aug)


   return(data.frame(gamma, beta, alpha, R2, p, # R2_aug,
    true_bias, ipw = coef(model_ipw)[2], ipw_bias, bias_decomp, bias_decomp_est, C_sigma,
                     bias_decomp_bound, var_tau_bound,
                #aug_ipw_bias, bias_decomp_est_aug, bias_decomp_bound_aug, var_tau_bound_aug,
                     var_tau = var(df_population$tau[df_population$V==0]),
                     ORV))
}

alpha = c(-0.5, -0.25, 0, 0.5)
gamma = c(0.5, 1, 1.5, 2)
beta = 0
set.seed(331)
df_results = c()
for(gamma_val in gamma){
    for(alpha_val in alpha){
        df_results = rbind(
            df_results, mclapply(1:1000, run_iter, n = 1000, N = 10000,
                alpha = alpha_val, gamma = gamma_val, beta = beta,
            mc.cores=detectCores()-1) %>% dplyr::bind_rows())
    }
}

#alpha = c(0.25, 0.5, 0.75, 1)
#gamma = -c(0.5, 1, 2, 5)
beta = 0
set.seed(331)
df_results_cont = c()
for(gamma_val in gamma){
    for(alpha_val in alpha){
        df_results_cont = rbind(
            df_results_cont, mclapply(1:1000, run_iter, n = 1000, N = 10000,
                alpha = alpha_val, gamma = gamma_val, beta = beta, V_binary=FALSE,
            mc.cores=detectCores()-1) %>% dplyr::bind_rows())
    }
}


# df_results %>% ggplot(aes(x = gamma, y = ipw_bias)) + geom_point()
# df_results %>% ggplot(aes(x = as.factor(alpha), y = p, fill = as.factor(gamma))) + geom_boxplot()

# df_results %>% ggplot(aes(x = as.factor(gamma), y = C_sigma, fill = as.factor(alpha))) + geom_boxplot()


# df_results %>% ggplot(aes(x = as.factor(gamma), y = ipw_bias, fill = as.factor(alpha))) + geom_boxplot()
# quartz()
# df_results %>% ggplot(aes(x = as.factor(gamma), y = mean(aug_ipw_bias), fill = as.factor(alpha))) + geom_point()
# df_results %>% ggplot(aes(x = as.factor(gamma), y = mean(ipw_bias), fill = as.factor(alpha))) + geom_point()
# df_results_cont %>% ggplot(aes(x = as.factor(alpha), y = p, fill = as.factor(gamma))) + geom_boxplot()


# df_results_inter %>% ggplot(aes(x = as.factor(gamma), y = R2, fill = as.factor(alpha))) + geom_boxplot()
# df_results_inter %>% ggplot(aes(x = as.factor(alpha), y = p, fill = as.factor(gamma))) + geom_boxplot()
# df_results_inter %>% ggplot(aes(x = as.factor(gamma), y = C_sigma, fill = as.factor(gamma))) + geom_boxplot()


# df_results %>% ggplot(aes(x = abs(true_bias), y = bias_decomp_est, color = as.factor(gamma))) + geom_point()
# df_results %>% ggplot(aes(x = abs(true_bias), y = ORV, color = as.factor(gamma))) + geom_point()
# df_results_cont %>% ggplot(aes(x = abs(true_bias), y = ORV, color = as.factor(gamma))) + geom_point()
# df_results %>% ggplot(aes(x = abs(true_bias), y = mean(ipw), color = as.factor(gamma))) + geom_point()
# df_results_cont %>% ggplot(aes(x = abs(true_bias), y = mean(ORV), color = as.factor(gamma))) + geom_point()


# df_results %>% ggplot(aes(x = as.factor(gamma), y = C_sigma, fill = as.factor(beta))) + geom_boxplot() +
# xlab(expression(alpha)) #+ facet_wrap(~beta)


# df_plot = rbind(data.frame(df_results, type = "Scenario 1 (Binary)", interact = "No Interactions"),
#     data.frame(df_results_cont, type = "Scenario 2 (Continuous)", interact = "No Interactions"),
#     data.frame(df_results_inter, type = "Scenario 3 (Binary)", interact = "w/ Interactions"),
#     data.frame(df_results_inter_cont, type = "Scenario 4 (Continuous)", interact = "w/ Interactions")) %>%
# gather(key = Parameter, value = Quantity, - gamma, -beta, -alpha, -type,-interact) %>%
# filter(Parameter %in% c("R2", "p", "C_sigma"))

#GENERATE FIGURE (a):
df_plot = rbind(data.frame(df_results, type = "Scenario 1 (Binary)", interact = "No Interactions"),
    data.frame(df_results_cont, type = "Scenario 2 (Continuous)", interact = "No Interactions")) %>%
gather(key = Parameter, value = Quantity, - gamma, -beta, -alpha, -type,-interact) %>%
filter(Parameter %in% c("R2", "p", "C_sigma"))


df_plot$Parameter = factor(df_plot$Parameter, levels = c("p", "R2", "C_sigma"),
    labels = c('p', parse(text = latex2exp::TeX("$R^2_{\\tau~V}$")),
    parse(text = latex2exp::TeX("$C_\\sigma$"))))

# df_plot$type = factor(df_plot$type,
#     labels = c(parse(text = latex2exp::TeX("$\\textbf{Scenario\\,1:\\,Binary}$")),
#     parse(text = latex2exp::TeX("$\\textbf{Scenario\\,2:\\,Continuous}$")),
#     parse(text = latex2exp::TeX("$\\textbf{Scenario\\,3:\\,Binary}$")),
#     parse(text = latex2exp::TeX("$\\textbf{Scenario\\,4:\\,Continuous}$"))))

df_plot$type = factor(df_plot$type,
    labels = c(parse(text = latex2exp::TeX("$\\textbf{Scenario\\,1\\,(Binary)}$")),
    parse(text = latex2exp::TeX("$\\textbf{Scenario\\,2\\,(Continuous)}$"))))


# df_plot$interact = factor(df_plot$interact,
#     labels = c(parse(text = latex2exp::TeX("$\\textbf{No\\,Interactions}$")),
#     parse(text = latex2exp::TeX("$\\textbf{Interactions}$"))))
pl1 = c("darkorange2", "#ffc14d", "#3dbbf5", "#88d575")


parameter_plot = df_plot[1:64000,] %>% group_by(gamma = as.factor(abs(gamma)), alpha = as.factor(alpha), Parameter, interact, type) %>%
summarize(
    Quantity = mean(Quantity)
    ) %>%
ggplot(aes(x = gamma, y = Quantity, color = alpha, group = alpha)) +
geom_point(size=2) + facet_grid(Parameter~type, scales='free_y', labeller = label_parsed) + #
geom_line(linewidth=0.9) + xlab(expression(gamma)) +
ylab("Parameter Value")+
scale_color_manual(name = expression(alpha), values = pl1) +
theme(legend.position='right') + ggtitle("(a) Sensitivity Parameters") #+ylim(0,1)

ggsave("parameter_sims.pdf", height=6, width=10.5)


#ORV FIGURE:
p1 = rbind(data.frame(df_results, type = "Scenario 1 (Binary)", interact = "No Interactions"),
    data.frame(df_results_cont, type = "Scenario 2 (Continuous)", interact = "No Interactions")) %>%
ggplot(aes(x = abs(ipw_bias), y = ORV)) +
geom_point(pch = 15, size=0.5, color = 'lightgray', alpha=0.25)


df_plot_orv = rbind(data.frame(df_results, type = "Scenario 1 (Binary)", interact = "No Interactions"),
    data.frame(df_results_cont, type = "Scenario 2 (Continuous)", interact = "No Interactions")) %>%
group_by(gamma, alpha, type) %>%
summarize(
    true_bias = abs(mean(ipw_bias)),
    ORV = mean(ORV)
)

p2 = p1 + geom_point(data = df_plot_orv, aes(x = true_bias, y = ORV, color = as.factor(alpha))) + xlab("True Error (IPW - T-PATE)")+
ggtitle("(b) ORV and Bias")
ggsave("orv_illustration.pdf", height=6, width=6)


library(cowplot)
cowplot::plot_grid(parameter_plot, p2, rel_widths = c(1.5,1))
ggsave("simulation_figure.pdf", height=4.5, width=11)


#------------------------------------------------------------------------------------------------
#APPENDIX RESULTS:


#INTERACTIONS:
# alpha = c(-0.5, -0.25, 0, 0.5)
# gamma = c(0.5, 1, 2, 5)
beta = -0.5
rho = 0.25
set.seed(331)
df_results_inter = c()
for(gamma_val in gamma){
    for(alpha_val in alpha){
        df_results_inter = rbind(
            df_results_inter, mclapply(1:1000, run_iter, n = 1000, N = 10000,
                alpha = alpha_val, gamma = gamma_val, beta = beta, rho = rho, V_binary=TRUE,
            mc.cores=detectCores()-1) %>% dplyr::bind_rows())
    }
}

# alpha = c(0.25, 0.5, 0.75, 1)
# beta = -0.5
set.seed(331)
df_results_inter_cont = c()
for(gamma_val in gamma){
    for(alpha_val in alpha){
        df_results_inter_cont = rbind(
            df_results_inter_cont, mclapply(1:1000, run_iter, n = 1000, N = 10000,
                alpha = alpha_val, gamma = gamma_val, beta = beta, rho=rho, V_binary=FALSE,
            mc.cores=detectCores()-1) %>% dplyr::bind_rows())
    }
}

# #BIAS RECONSTRUCTION PLOT:
# rbind(data.frame(df_results, type = "Scenario 1 (Binary)", interact = "No Interactions"),
#     data.frame(df_results_cont, type = "Scenario 2 (Continuous)", interact = "No Interactions"),
#     data.frame(df_results_inter, type = "Scenario 3 (Binary)", interact = "Interactions"),
#     data.frame(df_results_inter_cont, type = "Scenario 4 (Continuous)", interact = "Interactions")) %>%
# ggplot(aes(x = abs(ipw_bias), y = ORV)) +
# geom_point()#+facet_wrap(~type, scales = 'free')

p1 = rbind(data.frame(df_results, type = "Scenario 1 (Binary)", interact = "No Interactions"),
    data.frame(df_results_cont, type = "Scenario 2 (Continuous)", interact = "No Interactions"),
    data.frame(df_results_inter, type = "Scenario 3 (Binary)", interact = "Interactions"),
    data.frame(df_results_inter_cont, type = "Scenario 4 (Continuous)", interact = "Interactions")) %>%
#gather(key = Parameter, value = Quantity, - gamma, -beta, -alpha, -type,-interact) %>%
#filter(Parameter %in% c("R2", "p", "C_sigma", "true_bias", "ORV")) %>%
ggplot(aes(x = abs(ipw_bias), y = bias_decomp_bound)) +
geom_point(size=0.5, alpha=0.5, color = 'lightgray')+
facet_wrap(~type)+ geom_abline(yintercept=0) + ylim(0, 1.6) + xlim(0, 1.6)

df_summary_plot = rbind(data.frame(df_results, type = "Scenario 1 (Binary)", interact = "No Interactions"),
    data.frame(df_results_cont, type = "Scenario 2 (Continuous)", interact = "No Interactions"),
    data.frame(df_results_inter, type = "Scenario 3 (Binary)", interact = "Interactions"),
    data.frame(df_results_inter_cont, type = "Scenario 4 (Continuous)", interact = "Interactions")) %>%
group_by(gamma, alpha, type) %>%
summarize(
    true_bias = abs(mean(ipw_bias)),
    bias_bound = mean(bias_decomp_bound)
)

p1 + geom_point(data = df_summary_plot, aes(x = true_bias, y = bias_bound))+
xlab("True Error (IPW - T-PATE)") + ylab("Bias Bound")
ggsave("bias_recovery.pdf", height=5, width=6)


rbind(data.frame(df_results, type = "Scenario 1 (Binary)", interact = "No Interactions"),
    data.frame(df_results_cont, type = "Scenario 2 (Continuous)", interact = "No Interactions")) %>%
#gather(key = Parameter, value = Quantity, - gamma, -beta, -alpha, -type,-interact) %>%
#filter(Parameter %in% c("R2", "p", "C_sigma", "true_bias", "ORV")) %>%
ggplot(aes(x = var_tau, y = var_tau_bound)) +
geom_point(size=0.5, alpha=0.5, color = 'lightgray')+
facet_wrap(~type)+ geom_abline(yintercept=0) 

df_bias = rbind(data.frame(df_results, type = "Scenario 1 (Binary)", interact = "No Interactions"),
    data.frame(df_results_cont, type = "Scenario 2 (Continuous)", interact = "No Interactions"),
    data.frame(df_results_inter, type = "Scenario 3 (Binary)", interact = "w/ Interactions"),
    data.frame(df_results_inter_cont, type = "Scenario 4 (Continuous)", interact = "w/ Interactions")) %>%
gather(key = Parameter, value = Quantity, - gamma, -beta, -alpha, -type,-interact) %>%
filter(Parameter %in% c("bias_decomp", "true_bias", "bias_decomp_bound"))

df_bias %>% ggplot(aes(x = true_bias, y = bias_decomp)) + geom_point()
#------------------------------------------------------------------------------------------
#GENERATE FIGURE (a):
df_plot_inter = rbind(data.frame(df_results_inter, type = "Scenario 3"),
    data.frame(df_results_inter_cont, type = "Scenario 4")) %>%
gather(key = Parameter, value = Quantity, - gamma, -beta, -alpha, -type) %>%
filter(Parameter %in% c("R2", "p", "C_sigma"))


#df_plot_inter$Parameter = factor(df_plot_inter$Parameter, levels = c("p", "R2", "C_sigma"))

# df_plot$type = factor(df_plot$type,
#     labels = c(parse(text = latex2exp::TeX("$\\textbf{Scenario\\,1:\\,Binary}$")),
#     parse(text = latex2exp::TeX("$\\textbf{Scenario\\,2:\\,Continuous}$")),
#     parse(text = latex2exp::TeX("$\\textbf{Scenario\\,3:\\,Binary}$")),
#     parse(text = latex2exp::TeX("$\\textbf{Scenario\\,4:\\,Continuous}$"))))

# df_plot$interact = factor(df_plot$interact,
#     labels = c(parse(text = latex2exp::TeX("$\\textbf{No\\,Interactions}$")),
#     parse(text = latex2exp::TeX("$\\textbf{Interactions}$"))))


df_plot_inter$Parameter = factor(df_plot_inter$Parameter, levels = c("p", "R2", "C_sigma"),
    labels = c('p', parse(text = latex2exp::TeX("$R^2_{\\tau~V}$")),
    parse(text = latex2exp::TeX("$C_\\sigma$"))))

# df_plot$type = factor(df_plot$type,
#     labels = c(parse(text = latex2exp::TeX("$\\textbf{Scenario\\,1:\\,Binary}$")),
#     parse(text = latex2exp::TeX("$\\textbf{Scenario\\,2:\\,Continuous}$")),
#     parse(text = latex2exp::TeX("$\\textbf{Scenario\\,3:\\,Binary}$")),
#     parse(text = latex2exp::TeX("$\\textbf{Scenario\\,4:\\,Continuous}$"))))

df_plot_inter$type = factor(df_plot_inter$type,
    labels = c(parse(text = latex2exp::TeX("$\\textbf{Scenario\\,3\\,(Binary)}$")),
    parse(text = latex2exp::TeX("$\\textbf{Scenario\\,4\\,(Continuous)}$"))))

df_plot_inter %>% group_by(gamma = as.factor(abs(gamma)), alpha = as.factor(alpha), Parameter, type) %>%
summarize(
    Quantity = mean(Quantity)
    ) %>%
ggplot(aes(x = gamma, y = Quantity, color = alpha, group = alpha)) +
geom_point(size=2) + facet_grid(Parameter~type, scales='free_y', labeller=label_parsed) + #
geom_line(linewidth=0.9) + xlab(expression(gamma)) +
ylab("Parameter Value")+
scale_color_manual(name = expression(alpha), values = pl1) +
theme(legend.position='right') + ggtitle("(a) Sensitivity Parameters") #+ylim(0,1)
ggsave("parameter_sims_app.pdf", height=6, width=10.5)


#CONDUCT MORE ON INTERACTIONS: 
#Plot C_sigma:
set.seed(331)
df_interactions_binary = c()
beta =  c(-2, -1, 0, 1, 2)
alpha_val = 0
#gamma_val = 1
for(gamma_val in beta){
    for(beta_val in beta){
        df_interactions_binary  = rbind(
        df_interactions_binary , mclapply(1:1000, run_iter, n = 1000, N = 10000,
            alpha = alpha_val, gamma = gamma_val, beta = beta_val, rho =rho, V_binary=TRUE,
        mc.cores=detectCores()-1) %>% dplyr::bind_rows())
    }

}

set.seed(331)
df_interactions_cont = c()
for(gamma_val in beta){
    for(beta_val in beta){
        df_interactions_cont  = rbind(
            df_interactions_cont , mclapply(1:1000, run_iter, n = 1000, N = 10000,
                alpha = alpha_val, gamma = gamma_val, beta = beta_val, rho =rho, V_binary=FALSE,
            mc.cores=detectCores()-1) %>% dplyr::bind_rows())
    }
}

df_plot_inter = rbind(data.frame(df_interactions_binary, type = "Binary"),
    data.frame(df_interactions_cont, type = "Continuous")) %>%
gather(key = Parameter, value = Quantity, - gamma, -beta, -alpha, -type) %>%
filter(Parameter %in% c("R2", "p", "C_sigma"))

    #parse(text = latex2exp::TeX("$C_\\sigma$"))))

df_plot_inter$Parameter = factor(df_plot_inter$Parameter, levels = c("p", "R2", "C_sigma"),
    labels = c('p', parse(text = latex2exp::TeX(r'($R^2_{\tau\sim V}$)')),
    parse(text = latex2exp::TeX("$C_\\sigma$"))))

# df_plot$type = factor(df_plot$type,
#     labels = c(parse(text = latex2exp::TeX("$\\textbf{Scenario\\,1:\\,Binary}$")),
#     parse(text = latex2exp::TeX("$\\textbf{Scenario\\,2:\\,Continuous}$")),
#     parse(text = latex2exp::TeX("$\\textbf{Scenario\\,3:\\,Binary}$")),
#     parse(text = latex2exp::TeX("$\\textbf{Scenario\\,4:\\,Continuous}$"))))

df_plot_inter$type = factor(df_plot_inter$type,
    labels = c(parse(text = latex2exp::TeX("$\\textbf{Scenario\\,3\\,(Binary)}$")),
    parse(text = latex2exp::TeX("$\\textbf{Scenario\\,4\\,(Continuous)}$"))))


pl2 =c("#DF7E66FF", "#E09351FF", "#EDC775FF", "#94B594FF", "#224B5EFF")
df_plot_inter %>% group_by(gamma = as.factor(gamma), alpha = as.factor(alpha), beta = as.factor(beta), Parameter, type) %>%
summarize(
    Quantity = mean(Quantity)
    ) %>%
ggplot(aes(x = beta, y = Quantity, color = gamma, group = gamma)) +
geom_point(size=2) + facet_grid(Parameter~type, scales='free_y', labeller=label_parsed) + #
geom_line(linewidth=0.9) + xlab(expression(varphi)) +
ylab("Parameter Value") +
theme(legend.position='right') +scale_color_manual(name = expression(gamma), values = pl2)
ggsave("parameter_sims_app.pdf", height=6, width=10.5)



df_plot_inter_csigma = df_plot_inter %>% 
group_by(gamma = as.factor(gamma),  beta = as.factor(beta), Parameter, type) %>%
summarize(
    Quantity = mean(Quantity)
    ) %>%
filter(Parameter == "C_sigma") 

df_plot_inter_csigma$Parameter = factor(df_plot_inter_csigma$Parameter, 
    labels = c(parse(text = latex2exp::TeX("$C_{\\sigma}$"))))#,
    #parse(text = latex2exp::TeX("$C_\\sigma$")

df_plot_inter_csigma$type = factor(df_plot_inter_csigma$type,
    labels = c(parse(text = latex2exp::TeX("$\\textbf{Binary}$")),
    parse(text = latex2exp::TeX("$\\textbf{Continuous}$"))))

df_plot_inter_csigma %>% 
ggplot(aes(x = gamma, y = Quantity, color = beta, group = beta)) + geom_hline(yintercept=1, linetype='dotted') + 
geom_point() + geom_line() +
facet_grid(Parameter~type, labeller = label_parsed) + xlab(expression(gamma)) + ylab("Parameter Value") +
scale_color_manual(values = pl2, name = expression(varphi))+
theme(legend.position='right') #+ylim(0,1)
ggsave("csigma_sims.pdf", height=3, width=9)

# "#B75347FF",
df_bias = rbind(data.frame(df_results, Scenario = "(1) Binary"),
    data.frame(df_results_cont, Scenario = "(2) Continuous")) %>%
group_by(gamma = as.factor(gamma), alpha = as.factor(alpha), Scenario) %>%
summarize(
    true_bias = mean(true_bias),
    bias_decomp_bound = mean(bias_decomp_bound),
    bias_decomp = mean(bias_decomp_est)
)
df_bias %>% gather(key = key, value = Bias, -gamma, -alpha, -Scenario, -true_bias) %>%
mutate(label = ifelse(key == "bias_decomp_bound", "(Conservative) Bound", "Oracle Bias Estimate")) %>%
ggplot(aes(x = abs(true_bias), y = Bias, color = label, pch = Scenario)) +
geom_point(data = rbind(data.frame(df_results, Scenario = "(1) Binary"),
    data.frame(df_results_cont, Scenario = "(2) Continuous")),
aes(x = abs(true_bias), y = bias_decomp_bound), color = "skyblue3", pch='.', alpha=0.2)+
geom_abline(slope=1) +
geom_point(size=2) + xlab("(True) Bias") + ylab("Estimated Bias") +
scale_color_manual(values = c(pl2[5], pl2[1]), name ="Est. Type") +
scale_shape_manual(values = c(15, 19)) + theme(legend.position= 'right')
ggsave("bias_recovery.pdf", height=5, width=7)

df_bias %>% gather(key = key, value = Bias, -gamma, -alpha, -Scenario, -true_bias) %>%
mutate(label = ifelse(key == "bias_decomp_bound", "(Conservative) Bound", "Oracle Bias Estimate")) %>%
filter(abs(true_bias) >= Bias)

df_results  %>% ggplot(aes(x = var_tau, y = var_tau_bound)) + geom_point() + xlim(1.9, 5) + ylim(1.9, 5)
df_results_cont %>% ggplot(aes(x = var_tau, y = var_tau_bound)) + geom_point() + xlim(2, 10) + ylim(2,10)