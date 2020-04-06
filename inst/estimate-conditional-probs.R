## estimate conditional probabilities

options(scipen=999)
library(googlesheets4)
library(tidyverse)
library(mgcv)
library(covidSeverity)

age_specific_data <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1pDfQ2SkO--F2WNfjZxx6t9V7K51IZW0bJKV5cwSCZfA/edit#gid=1769840547",
                                               sheet="age risk") %>%
  filter(USE_pipeline==T)

## look at one conditional probability (e.g. proportion symptomatic pSymp_Inf)
raw_params <- age_specific_data %>%
  mutate(N = ifelse(is.na(N),
                    ifelse(!is.na(X) & value>0, round(X/value),
                           ifelse(!is.na(valueR),
                                  round(value*(1-value)*(qnorm(.975)/(valueR-value))^2),
                                  1000)),
                    N),
         X = ifelse(is.na(X),ceiling(N*value),X),
         age_rng = paste0(ageL,"_",ageR))

## reallocate cases to new age groups
## first expand summarized data
expanded_dat <- c()
## do each param separately
params <- unique(raw_params$param)
num_params <- length(params)
for(i in 1:num_params){
  param_dat <- filter(raw_params, param==params[i])
  ## keep the studies separate
  studies <- unique(param_dat$source)
  num_studies <- length(studies)
  for(j in 1:num_studies){
    tmp <- param_dat %>% filter(source==studies[j])
    ## break out each age group
    age_groups <- unique(tmp$age_rng)
    num_ag <- length(age_groups)
    for(k in 1:num_ag){
      tmp_ag <- tmp %>% filter(age_rng==age_groups[k])

      expanded_dat <- bind_rows(expanded_dat,
                                tibble(param=tmp_ag$param[1],
                                       study=j,
                                       age=sample(tmp_ag$ageL[1]:tmp_ag$ageR[1],
                                                  tmp_ag$N[1], replace=T),
                                       x=sample(c(rep(1,tmp_ag$X[1]),
                                                  rep(0,tmp_ag$N[1]-tmp_ag$X[1])),
                                                tmp_ag$N[1], replace=F)))
    }
  }
}

p_symp <- est_age_spec_param(expanded_dat)
ggplot(data=p_symp, aes(x=age_grp)) +
  geom_errorbar(aes(ymin=plogis(logit_mean-1.96*logit_sd), ymax=plogis(logit_mean+1.96*logit_sd)), alpha=0.4) +
  geom_point(aes(y=plogis(logit_mean)))

p_death <- est_age_spec_param(expanded_dat, "pDeath_Symp")
ggplot(data=p_death, aes(x=age_grp)) +
  geom_errorbar(aes(ymin=plogis(logit_mean-1.96*logit_sd), ymax=plogis(logit_mean+1.96*logit_sd)), alpha=0.4) +
  geom_point(aes(y=plogis(logit_mean)))

p_hosp <- est_age_spec_param(expanded_dat, "pHosp_Symp")
ggplot(data=p_hosp, aes(x=age_grp)) +
  geom_errorbar(aes(ymin=plogis(logit_mean-1.96*logit_sd), ymax=plogis(logit_mean+1.96*logit_sd)), alpha=0.4) +
  geom_point(aes(y=plogis(logit_mean)))

p_icu <- est_age_spec_param(expanded_dat, "pICU_Hosp")
ggplot(data=p_icu, aes(x=age_grp)) +
  geom_errorbar(aes(ymin=plogis(logit_mean-1.96*logit_sd), ymax=plogis(logit_mean+1.96*logit_sd)), alpha=0.4) +
  geom_point(aes(y=plogis(logit_mean)))

p_vent <- est_age_spec_param(expanded_dat, "pVent_Symp")
ggplot(data=p_vent, aes(x=age_grp)) +
  geom_errorbar(aes(ymin=plogis(logit_mean-1.96*logit_sd), ymax=plogis(logit_mean+1.96*logit_sd)), alpha=0.4) +
  geom_point(aes(y=plogis(logit_mean)))

## probability symptomatic
p_symp_inf_age
## probability hospitalized given symptomatic
p_hosp_symp_age
## probability death given symptomatic
p_death_symp_age
## probability ICU given hospitalized
p_icu_hosp_age
## probability ventilator given symptomatic
p_vent_symp
