## estimate conditional probabilities

options(scipen=999)
library(googlesheets4)
library(tidyverse)
library(mgcv)
library(covidSeverity)
library(doMC)
registerDoMC(10)
n_sims <- 40
n_preds <- 1e3
age_grps <- c(seq(0,80,by=10),100)

data("US_age_geoid_pct")
data("US_age_geoid_pop")

age_specific_data <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1pDfQ2SkO--F2WNfjZxx6t9V7K51IZW0bJKV5cwSCZfA/edit#gid=1769840547",
                                               sheet="age risk") %>%
  filter(USE_pipeline==T)
1

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
## first expand summarized data, n_sims times
all_params <- foreach(h=1:n_sims, .combine=rbind) %dopar% {
  set.seed(h)
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

  ## make proportion symptomatic given infected
  p_symp <- est_age_spec_param(expanded_dat=expanded_dat,
                               param_to_est="p_symp_inf",
                               age_cats=age_grps,
                               n_preds=n_preds,
                               study_wt="none")
  ggplot(data=p_symp$pred_sum, aes(x=age_grp)) +
    geom_errorbar(aes(ymin=p_symp_inf_lb,
                      ymax=p_symp_inf_ub), alpha=0.4) +
    geom_point(aes(y=p_symp_inf_med)) +
    scale_y_continuous("Probability symptomatic, given infected") +
    scale_x_discrete("Age groups") +
    theme_bw()

  ## CFR
  p_death <- est_age_spec_param(expanded_dat=expanded_dat,
                                param_to_est="p_death_symp",
                                age_cats=age_grps,
                                n_preds=n_preds,
                                study_wt="none")
  ggplot(data=p_death$pred_sum, aes(x=age_grp)) +
    geom_errorbar(aes(ymin=p_death_symp_lb,
                      ymax=p_death_symp_ub), alpha=0.4) +
    geom_point(aes(y=p_death_symp_med)) +
    scale_y_continuous("Case fatality rate") +
    scale_x_discrete("Age groups") +
    theme_bw()

  ## hospitalization rate amongst symptomatic
  p_hosp <- est_age_spec_param(expanded_dat=expanded_dat,
                               param_to_est="p_hosp_symp",
                               age_cats=age_grps,
                               n_preds=n_preds,
                               study_wt="none")
  ggplot(data=p_hosp$pred_sum, aes(x=age_grp)) +
    geom_errorbar(aes(ymin=p_hosp_symp_lb,
                      ymax=p_hosp_symp_ub), alpha=0.4) +
    geom_point(aes(y=p_hosp_symp_med)) +
    scale_y_continuous("Probability hospitalized, given symptomatic") +
    scale_x_discrete("Age groups") +
    theme_bw()

  ## ICU rate amongst hospitalized
  p_icu <- est_age_spec_param(expanded_dat=expanded_dat,
                              param_to_est="p_icu_hosp",
                              age_cats=age_grps,
                              n_preds=n_preds,
                              study_wt="none")
  ggplot(data=p_icu$pred_sum, aes(x=age_grp)) +
    geom_errorbar(aes(ymin=p_icu_hosp_lb,
                      ymax=p_icu_hosp_ub), alpha=0.4) +
    geom_point(aes(y=p_icu_hosp_med)) +
    scale_y_continuous("Probability ICU, given hospitalized") +
    scale_x_discrete("Age groups") +
    theme_bw()

  ## ventilation rate amongst symptomatic
  p_vent <- est_age_spec_param(expanded_dat,
                               param_to_est="p_vent_icu",
                               age_cats=c(0,100),
                               n_preds=n_preds,
                               study_wt="none")
  ggplot(data=p_vent$pred_sum, aes(x=age_grp)) +
    geom_errorbar(aes(ymin=p_vent_icu_lb,
                      ymax=p_vent_icu_ub), alpha=0.4) +
    geom_point(aes(y=p_vent_icu_med)) +
    scale_y_continuous("Probability invasive ventilation, given ICU") +
    scale_x_discrete("Age groups") +
    theme_bw()

  ## get all output parameters
  geoid_params <- est_geoid_params(p_symp, US_age_geoid_pct) %>%
    left_join(est_geoid_params(p_death, US_age_geoid_pct), by="geoid") %>%
    left_join(est_geoid_params(p_hosp, US_age_geoid_pct), by="geoid") %>%
    left_join(est_geoid_params(p_icu, US_age_geoid_pct), by="geoid") %>%
    left_join(est_geoid_params(p_vent,
                               matrix(1, nrow=nrow(US_age_geoid_pct), ncol=1,
                                      dimnames=list(rownames(US_age_geoid_pct),
                                                    "[0,100)"))), by="geoid") %>%
    left_join(est_geoid_rrs(pred_mtx=p_death$pred_mtx,
                            param_to_est="death_symp",
                            geoid_age_mtx=US_age_geoid_pct,
                            geoid_pops=US_age_geoid_pop), by="geoid") %>%
    left_join(est_geoid_rrs(pred_mtx=p_hosp$pred_mtx,
                            param_to_est="hosp_symp",
                            geoid_age_mtx=US_age_geoid_pct,
                            geoid_pops=US_age_geoid_pop), by="geoid") %>%
    mutate(sim=h)

  return(geoid_params)
}

## find the simulation with the median CFR
med_sim <- (all_params %>%
  select(sim, death_symp_overall) %>%
  distinct() %>%
  mutate(dist_from_med=abs(death_symp_overall-median(death_symp_overall))) %>%
  arrange(dist_from_med) %>%
  slice(1))$sim

## median sim
median_sim <- all_params %>%
  filter(sim==med_sim)

## county parameter distributions
median_sim %>%
  summary()

median_sim %>%
  select(-sim) %>%
  write_csv("generated_data/geoid-params.csv")
