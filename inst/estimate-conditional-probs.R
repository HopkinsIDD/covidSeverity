## estimate conditional probabilities

options(scipen=999)
library(googlesheets4)
library(tidyverse)
library(mgcv)
library(covidSeverity)
library(doMC)
registerDoMC(10)
n_sims <- 40

load("data/USpop_geoid_agecat.Rdata")

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
## first expand summarized data, n_sims times
all_params <- foreach(h=1:n_sims, .combine=rbind) %dopar% {
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
  p_symp <- est_age_spec_param(expanded_dat)
  ggplot(data=p_symp, aes(x=age_grp)) +
    geom_errorbar(aes(ymin=plogis(logit_mean-1.96*logit_sd),
                      ymax=plogis(logit_mean+1.96*logit_sd)), alpha=0.4) +
    geom_point(aes(y=plogis(logit_mean)))

  ## CFR
  p_death <- est_age_spec_param(expanded_dat, "pDeath_Symp")
  ggplot(data=p_death, aes(x=age_grp)) +
    geom_errorbar(aes(ymin=plogis(logit_mean-1.96*logit_sd),
                      ymax=plogis(logit_mean+1.96*logit_sd)), alpha=0.4) +
    geom_point(aes(y=plogis(logit_mean)))

  ## hospitalization rate amongst symptomatic
  p_hosp <- est_age_spec_param(expanded_dat, "pHosp_Symp")
  ggplot(data=p_hosp, aes(x=age_grp)) +
    geom_errorbar(aes(ymin=plogis(logit_mean-1.96*logit_sd),
                      ymax=plogis(logit_mean+1.96*logit_sd)), alpha=0.4) +
    geom_point(aes(y=plogis(logit_mean)))

  ## ICU rate amongst hospitalized
  p_icu <- est_age_spec_param(expanded_dat, "pICU_Hosp")
  ggplot(data=p_icu, aes(x=age_grp)) +
    geom_errorbar(aes(ymin=plogis(logit_mean-1.96*logit_sd),
                      ymax=plogis(logit_mean+1.96*logit_sd)), alpha=0.4) +
    geom_point(aes(y=plogis(logit_mean)))

  ## ventilation rate amongst symptomatic
  p_vent <- est_age_spec_param(expanded_dat, "pVent_ICU")
  ggplot(data=p_vent, aes(x=age_grp)) +
    geom_errorbar(aes(ymin=plogis(logit_mean-1.96*logit_sd),
                      ymax=plogis(logit_mean+1.96*logit_sd)), alpha=0.4) +
    geom_point(aes(y=plogis(logit_mean)))

  ## load geoid populations
  geoid_pops <- t(USpop_geoid_agecat$est_age) %>%
    as_tibble() %>%
    setNames(USpop_geoid_agecat$GEOID[,,drop=T]) %>%
    pivot_longer(cols=1:(nrow(USpop_geoid_agecat$GEOID)),
                 names_to="geoid",
                 values_to="age_pop") %>%
    group_by(geoid) %>%
    summarize(pop=sum(age_pop))

  ## geoid symptomatic given
  geoid_symp <- t(USpop_geoid_agecat$p_age) %>%
    as_tibble() %>%
    setNames(USpop_geoid_agecat$GEOID[,,drop=T]) %>%
    bind_cols(p_symp) %>%
    pivot_longer(cols=1:(nrow(USpop_geoid_agecat$GEOID)),
                 names_to="geoid",
                 values_to="wt") %>%
    left_join(geoid_pops) %>%
    group_by(geoid) %>%
    summarize(symp_logit_mean=weighted.mean(logit_mean, wt),
              symp_logit_sd=sqrt(weighted.mean(logit_mean^2+logit_sd^2, wt)-weighted.mean(logit_mean, wt)^2)/sqrt(pop[1]))

  geoid_death <- t(USpop_geoid_agecat$p_age) %>%
    as_tibble() %>%
    setNames(USpop_geoid_agecat$GEOID[,,drop=T]) %>%
    bind_cols(p_death) %>%
    pivot_longer(cols=1:(nrow(USpop_geoid_agecat$GEOID)),
                 names_to="geoid",
                 values_to="wt") %>%
    left_join(geoid_pops) %>%
    group_by(geoid) %>%
    summarize(death_logit_mean=weighted.mean(logit_mean, wt),
              death_logit_sd=sqrt(weighted.mean(logit_mean^2+logit_sd^2, wt)-weighted.mean(logit_mean, wt)^2)/sqrt(mean(pop)))

  geoid_hosp <- t(USpop_geoid_agecat$p_age) %>%
    as_tibble() %>%
    setNames(USpop_geoid_agecat$GEOID[,,drop=T]) %>%
    bind_cols(p_hosp) %>%
    pivot_longer(cols=1:(nrow(USpop_geoid_agecat$GEOID)),
                 names_to="geoid",
                 values_to="wt") %>%
    left_join(geoid_pops) %>%
    group_by(geoid) %>%
    summarize(hosp_logit_mean=weighted.mean(logit_mean, wt),
              hosp_logit_sd=sqrt(weighted.mean(logit_mean^2+logit_sd^2, wt)-weighted.mean(logit_mean, wt)^2)/sqrt(mean(pop)))


  geoid_icu <- t(USpop_geoid_agecat$p_age) %>%
    as_tibble() %>%
    setNames(USpop_geoid_agecat$GEOID[,,drop=T]) %>%
    bind_cols(p_icu) %>%
    pivot_longer(cols=1:(nrow(USpop_geoid_agecat$GEOID)),
                 names_to="geoid",
                 values_to="wt") %>%
    left_join(geoid_pops) %>%
    group_by(geoid) %>%
    summarize(icu_logit_mean=weighted.mean(logit_mean, wt),
              icu_logit_sd=sqrt(weighted.mean(logit_mean^2+logit_sd^2, wt)-weighted.mean(logit_mean, wt)^2)/sqrt(mean(pop)))


  geoid_vent <- t(USpop_geoid_agecat$p_age) %>%
    as_tibble() %>%
    setNames(USpop_geoid_agecat$GEOID[,,drop=T]) %>%
    bind_cols(p_vent) %>%
    pivot_longer(cols=1:(nrow(USpop_geoid_agecat$GEOID)),
                 names_to="geoid",
                 values_to="wt") %>%
    left_join(geoid_pops) %>%
    group_by(geoid) %>%
    summarize(vent_logit_mean=weighted.mean(logit_mean, wt),
              vent_logit_sd=sqrt(weighted.mean(logit_mean^2+logit_sd^2, wt)-weighted.mean(logit_mean, wt)^2)/sqrt(mean(pop)))

  geoid_params <- left_join(geoid_pops, geoid_symp, by="geoid") %>%
    left_join(geoid_death, by="geoid") %>%
    left_join(geoid_hosp, by="geoid") %>%
    left_join(geoid_icu, by="geoid") %>%
    left_join(geoid_vent, by="geoid") %>%
    mutate(p_symp_inf=plogis(symp_logit_mean),
           p_hosp_symp=plogis(hosp_logit_mean),
           p_death_symp=plogis(death_logit_mean),
           p_death_inf=p_death_symp*p_symp_inf,
           p_hosp_inf=p_hosp_symp*p_symp_inf,
           p_icu_hosp=plogis(icu_logit_mean),
           p_vent_icu=plogis(vent_logit_mean),
           sim=h)
  return(geoid_params)
}
all_params %>%
  filter(sim==17) %>%
  select(-sim) %>%
  write_csv("generated_data/geoid-params.csv")
