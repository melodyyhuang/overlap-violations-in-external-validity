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


N = 10000
library(ggplot2)
library(tidyverse)
dgp<-function(p, gamma, noise_scale){
    u0 = rnorm(N, 0, 1)
    u1 = rnorm(N, 0, noise_scale)
    V = rbinom(N, size = 1, prob=p)
    u = ifelse(V==1, u1, u0)
    tau = 2 - gamma*V + u
    C_sigma = var(u1)/var(u0)
    R2 = cor(V, tau)^2
    return(data.frame(tau, V, R2, gamma, p, C_sigma))
}

set.seed(331)
dgp_c0.5 = dgp(0.25, 1.75, 0.5)
dgp_c1 = dgp(0.25, 2, 1)
dgp_c2= dgp(0.25, 2.5, 2)

df_plot = rbind(dgp_c0.5, dgp_c1, dgp_c2) 
df_plot$C_sigma_val = ifelse(df_plot$gamma == 1.75, "Scenario 1: C = 0.25", ifelse(df_plot$gamma==2, "Scenario 2: C = 1", "Scenario 3: C = 4"))
df_plot %>% 
ggplot(aes(x = tau, fill = ifelse(V==1, "Omitted", "Included"), pattern=as.factor(V), pattern_fill = as.factor(V))) + facet_wrap(~C_sigma_val) + 
geom_density_pattern(alpha=0.5, pattern_linetype = 0, pattern_size = 1,
pattern_density = 0.1, pattern_spacing=0.02) + 
scale_fill_manual(values = c( "skyblue3", "slategray"), name = "")+
scale_pattern_manual(values = c("none", "stripe"), guide='none')+
scale_pattern_fill_manual(values = c("skyblue3", "black"), guide='none')+
xlab(expression(tau)) + ylab("Density")
ggsave("c_sigma_example.pdf", height=3.5, width=9)