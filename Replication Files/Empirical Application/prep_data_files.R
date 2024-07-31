
drop_index = which(df_sample$ind_found_e1 == 0 | is.na(df_sample$training_hours_e))

survey_data = survey_data[-which(survey_data$partid %in% df_sample$partid[drop_index]),]
df_sample = df_sample[-drop_index,]
df_sample$wealthindex = survey_data$wealthindex[survey_data$yop==1]
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
