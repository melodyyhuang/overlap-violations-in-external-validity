ggMelody<-theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 17, face = "bold"),
                                  axis.text=element_text(size=9), #size = 12
                                  #legend.text=element_text(size=7),
                                  legend.position='bottom', axis.title = element_text(size = 12),
                                  strip.text.x = element_text(size = 12, face='bold'),
                                  strip.text.y = element_text(size = 12, face='bold'),
                                  plot.subtitle=element_text(size = 14, hjust=0.5))
theme_set(ggMelody)

#Varying subgroup cutoffs:
#Construct
threshold = 9
vary_benchmark<-function(covariate, threshold){
    survey_data$covariate_belowX = ifelse(survey_data[[covariate]] <= threshold, 1, 0)
    df_sample$covariate_belowX = ifelse(df_sample[[covariate]] <= threshold, 1, 0)
    return(
        rbind(
            data.frame(outcome = "Hours of Vocational Training",
                       covariate = covariate,
                       threshold = threshold,
                       benchmark("covariate_belowX", outcome = "training_hours_e",
                        df_sample, survey_data, sigma_w2_hrs)
            ),
            data.frame(outcome = "Cash Earnings",
                       covariate = covariate,
                       threshold = threshold,
                       benchmark("covariate_belowX", outcome = "profits4w_real_p99_e",
                        df_sample, survey_data, sigma_w2_earnings)))
    )
}

df_vary_educ = lapply(2:12, vary_benchmark, covariate = "education") %>% bind_rows()
df_vary_age = lapply(16:27, vary_benchmark, covariate = "age") %>% bind_rows()
df_vary_wealth = lapply(seq(-0.8, -0.1, by = 0.05), vary_benchmark, covariate = "wealthindex") %>% bind_rows()




df_hrs_plot = rbind(df_vary_educ, df_vary_age, df_vary_wealth) %>%
filter(outcome == "Hours of Vocational Training") %>%
mutate(killer_confounder = ifelse(bias > ipw_hrs, "Yes", "No"))
df_hrs_plot$Bias = df_hrs_plot$bias
df_hrs_plot$covariate[df_hrs_plot$covariate == "education"] = "Education"
df_hrs_plot$covariate[df_hrs_plot$covariate == "age"] = "Age"
df_hrs_plot$covariate[df_hrs_plot$covariate == "wealthindex"] = "Wealth Index"
shade_color = 'skyblue3'
alpha_val = 0.2
df_hrs_plot %>%
gather(key = Parameter, value = Value, -c(outcome, covariate, threshold, killer_confounder)) %>%
filter(Parameter != 'rho' & Parameter !="bias")%>%
ggplot(aes(x = threshold, y = Value)) +
#annotate("rect", xmin=, xmax= 45, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "green") +
facet_grid(Parameter~covariate, scales = 'free')+
#Add killer confounder background shading:
geom_rect(data = data.frame(covariate = "Education", Parameter = 'Bias', Value = 1),
    aes(xmin = 3, xmax = 12, ymin = 0, ymax = Inf),
    alpha = alpha_val, fill=shade_color, inherit.aes = FALSE) +
 geom_rect(data = data.frame(covariate = "Age", Parameter = 'Bias', Value = 1),
    aes(xmin = 22, xmax = 27, ymin = 0, ymax = Inf),
    alpha = alpha_val, fill=shade_color, inherit.aes = FALSE) +
 geom_rect(data = data.frame(covariate = "Education", Parameter = 'p', Value = 1),
    aes(xmin = 3, xmax = 12, ymin = 0, ymax = Inf),
    alpha = alpha_val, fill=shade_color, inherit.aes = FALSE) +
 geom_rect(data = data.frame(covariate = "Age", Parameter = 'p', Value = 1),
    aes(xmin = 22, xmax = 27, ymin = 0, ymax = Inf),
    alpha = alpha_val, fill=shade_color, inherit.aes = FALSE) +
 geom_rect(data = data.frame(covariate = "Education", Parameter = 'R2', Value = 1),
    aes(xmin = 3, xmax = 12, ymin = 0, ymax = Inf),
    alpha = alpha_val, fill=shade_color, inherit.aes = FALSE) +
 geom_rect(data = data.frame(covariate = "Age", Parameter = 'R2', Value = 1),
    aes(xmin = 22, xmax = 27, ymin = 0, ymax = Inf),
    alpha = alpha_val, fill=shade_color, inherit.aes = FALSE) +
  #Add original cutoff line:
  geom_vline(data = data.frame(covariate = "Education", Parameter = 'Bias', Value = 1),
    aes(xintercept=9), inherit.aes = FALSE) +
 geom_vline(data = data.frame(covariate = "Education", Parameter = 'p', Value = 1),
    aes(xintercept=9), inherit.aes = FALSE) +
 geom_vline(data = data.frame(covariate = "Education", Parameter = 'R2', Value = 1),
    aes(xintercept=9), inherit.aes = FALSE) +
 geom_vline(data = data.frame(covariate = "Wealth Index", Parameter = 'Bias', Value = 1),
    aes(xintercept=-0.63), inherit.aes = FALSE) +
 geom_vline(data = data.frame(covariate = "Wealth Index", Parameter = 'p', Value = 1),
    aes(xintercept=-0.63), inherit.aes = FALSE) +
 geom_vline(data = data.frame(covariate = "Wealth Index", Parameter = 'R2', Value = 1),
    aes(xintercept=-0.63), inherit.aes = FALSE) +
 geom_vline(data = data.frame(covariate = "Age", Parameter = 'Bias', Value = 1),
    aes(xintercept=18), inherit.aes = FALSE) +
 geom_vline(data = data.frame(covariate = "Age", Parameter = 'p', Value = 1),
    aes(xintercept=18), inherit.aes = FALSE) +
 geom_vline(data = data.frame(covariate = "Age", Parameter = 'R2', Value = 1),
    aes(xintercept=18), inherit.aes = FALSE) +
 geom_point() + geom_line() + xlab("Threshold Cutoff")+ggtitle("Hours of Vocational Training")

ggsave("Figures/vary_cutoff_hrs.pdf", height=6, width=10)


shade_color = 'darkolivegreen'
df_earnings_plot = rbind(df_vary_educ, df_vary_age, df_vary_wealth) %>%
filter(outcome == "Cash Earnings") %>%
mutate(killer_confounder = ifelse(bias > ipw_earnings, "Yes", "No"))
df_earnings_plot$Bias = df_earnings_plot$bias
df_earnings_plot$covariate[df_earnings_plot$covariate == "education"] = "Education"
df_earnings_plot$covariate[df_earnings_plot$covariate == "age"] = "Age"
df_earnings_plot$covariate[df_earnings_plot$covariate == "wealthindex"] = "Wealth Index"
#shade_color = 'red4'
alpha_val = 0.2
df_earnings_plot %>%
gather(key = Parameter, value = Value, -c(outcome, covariate, threshold, killer_confounder)) %>%
filter(Parameter != 'rho' & Parameter !="bias")%>%
ggplot(aes(x = threshold, y = Value)) +
#annotate("rect", xmin=, xmax= 45, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "green") +
facet_grid(Parameter~covariate, scales = 'free')+
#Add killer confounder background shading:
geom_rect(data = data.frame(covariate = "Education", Parameter = 'Bias', Value = 1),
    aes(xmin = 5, xmax = 12, ymin = 0, ymax = Inf),
    alpha = alpha_val, fill=shade_color, inherit.aes = FALSE) +
 geom_rect(data = data.frame(covariate = "Age", Parameter = 'Bias', Value = 1),
    aes(xmin = 16, xmax = 27, ymin = 0, ymax = Inf),
    alpha = alpha_val, fill=shade_color, inherit.aes = FALSE) +
 geom_rect(data = data.frame(covariate = "Education", Parameter = 'p', Value = 1),
    aes(xmin = 5, xmax = 12, ymin = 0, ymax = Inf),
    alpha = alpha_val, fill=shade_color, inherit.aes = FALSE) +
 geom_rect(data = data.frame(covariate = "Age", Parameter = 'p', Value = 1),
    aes(xmin = 16, xmax = 27, ymin = 0, ymax = Inf),
    alpha = alpha_val, fill=shade_color, inherit.aes = FALSE) +
 geom_rect(data = data.frame(covariate = "Education", Parameter = 'R2', Value = 1),
    aes(xmin = 5, xmax = 12, ymin = 0, ymax = Inf),
    alpha = alpha_val, fill=shade_color, inherit.aes = FALSE) +
 geom_rect(data = data.frame(covariate = "Age", Parameter = 'R2', Value = 1),
    aes(xmin = 16, xmax = 27, ymin = 0, ymax = Inf),
    alpha = alpha_val, fill=shade_color, inherit.aes = FALSE) +
 #Add original cutoff line:
 geom_vline(data = data.frame(covariate = "Education", Parameter = 'Bias', Value = 1),
    aes(xintercept=9), inherit.aes = FALSE) +
 geom_vline(data = data.frame(covariate = "Education", Parameter = 'p', Value = 1),
    aes(xintercept=9), inherit.aes = FALSE) +
 geom_vline(data = data.frame(covariate = "Education", Parameter = 'R2', Value = 1),
    aes(xintercept=9), inherit.aes = FALSE) +
 geom_vline(data = data.frame(covariate = "Wealth Index", Parameter = 'Bias', Value = 1),
    aes(xintercept=-0.63), inherit.aes = FALSE) +
 geom_vline(data = data.frame(covariate = "Wealth Index", Parameter = 'p', Value = 1),
    aes(xintercept=-0.63), inherit.aes = FALSE) +
 geom_vline(data = data.frame(covariate = "Wealth Index", Parameter = 'R2', Value = 1),
    aes(xintercept=-0.63), inherit.aes = FALSE) +
 geom_vline(data = data.frame(covariate = "Age", Parameter = 'Bias', Value = 1),
    aes(xintercept=18), inherit.aes = FALSE) +
 geom_vline(data = data.frame(covariate = "Age", Parameter = 'p', Value = 1),
    aes(xintercept=18), inherit.aes = FALSE) +
 geom_vline(data = data.frame(covariate = "Age", Parameter = 'R2', Value = 1),
    aes(xintercept=18), inherit.aes = FALSE) +
 geom_point() + geom_line() + xlab("Threshold Cutoff")+ggtitle("Cash Earnings")

ggsave("Figures/vary_cutoff_earnings.pdf", height=6, width=10)