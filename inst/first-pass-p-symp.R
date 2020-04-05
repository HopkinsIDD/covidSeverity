## set up data for Stan model
##

options(mc.cores=4,
        scipen=999)
library(googlesheets4)
library(tidyverse)
library(rstan)

n_iter <- 1e4
n_warmup <- 1e4-5e3
p_adapt_delta <- 0.99
n_max_treedepth <- 20

age_specific_data <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1pDfQ2SkO--F2WNfjZxx6t9V7K51IZW0bJKV5cwSCZfA/edit#gid=1769840547",
                                               sheet="age risk") %>%
  filter(USE_pipeline==T)

## look at one conditional probability (e.g. proportion symptomatic pSymp_Inf)
raw_params <- age_specific_data %>%
  # filter(param=="pSymp_Inf") %>%
  mutate(N = ifelse(is.na(N),
                    ifelse(!is.na(X), round(X/value),
                           ifelse(!is.na(valueR),
                                  round(value*(1-value)*(qnorm(.975)/(valueR-value))^2),
                                  1000)),
                    N),
         X = ifelse(is.na(X),round(N*value),X))

## reallocate cases to new age groups
## first expand summarized data
expanded_dat <- c()
## do each param separately
params <- unique(raw_params$param)
num_params <- length(params)
for(i in 1:num_params){
  param_dat <- filter(raw_params, param==params[i])
  if(nrow(param_dat)==1){
    next
  }
  ## keep the studies separate
  studies <- unique(param_dat$source)
  num_studies <- length(studies)
  for(j in 1:num_studies){
    tmp <- param_dat %>% filter(source==studies[j])
    ## break out each age group
    age_groups <- unique(tmp$ageL)
    num_ag <- length(age_groups)
    for(k in 1:num_ag){
      tmp_ag <- tmp %>% filter(ageL==age_groups[k])

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
age_params <- expanded_dat %>%
  mutate(age_cat=cut(age,seq(0,100,by=10),right = FALSE),
         age_grp=as.factor(age_cat) %>% as.numeric) %>%
  group_by(param, study, age_cat, age_grp) %>%
  summarize(X=sum(x),
            N=n())

age_stan <- stan_model("R/multilevel_age_prob.stan")

symp_inf_est <- sampling(age_stan,
                         data=list(S=nrow(age_params %>%
                                            filter(param=="pSymp_Inf")),
                                   J=max(age_params$study[age_params$param=="pSymp_Inf"]),
                                   A=max(age_params$age_grp),
                                   agegrp=age_params$age_grp[age_params$param=="pSymp_Inf"],
                                   id=age_params$study[age_params$param=="pSymp_Inf"],
                                   N=age_params$N[age_params$param=="pSymp_Inf"],
                                   Y=age_params$X[age_params$param=="pSymp_Inf"]),
                         pars=c("prob_age"),
                         iter=n_iter,
                         warmup=n_warmup,
                         control=list(adapt_delta=p_adapt_delta,
                                      max_treedepth=n_max_treedepth),
                         save_warmup=F)

death_symp_est <- sampling(age_stan,
                           data=list(S=nrow(age_params %>%
                                              filter(param=="pDeath_Symp")),
                                     J=max(age_params$study[age_params$param=="pDeath_Symp"]),
                                     A=max(age_params$age_grp),
                                     agegrp=age_params$age_grp[age_params$param=="pDeath_Symp"],
                                     id=age_params$study[age_params$param=="pDeath_Symp"],
                                     N=age_params$N[age_params$param=="pDeath_Symp"],
                                     Y=age_params$X[age_params$param=="pDeath_Symp"]),
                           pars=c("prob_age"),
                           iter=n_iter,
                           warmup=n_warmup,
                           control=list(adapt_delta=p_adapt_delta,
                                        max_treedepth=n_max_treedepth),
                           save_warmup=F)

hosp_symp_est <- sampling(age_stan,
                          data=list(S=nrow(age_params %>%
                                             filter(param=="pHosp_Symp")),
                                    J=max(age_params$study[age_params$param=="pHosp_Symp"]),
                                    A=max(age_params$age_grp),
                                    agegrp=age_params$age_grp[age_params$param=="pHosp_Symp"],
                                    id=age_params$study[age_params$param=="pHosp_Symp"],
                                    N=age_params$N[age_params$param=="pHosp_Symp"],
                                    Y=age_params$X[age_params$param=="pHosp_Symp"]),
                          pars=c("prob_age"),
                          iter=n_iter,
                          warmup=n_warmup,
                          control=list(adapt_delta=p_adapt_delta,
                                       max_treedepth=n_max_treedepth),
                          save_warmup=F)

icu_hosp_est <- sampling(age_stan,
                          data=list(S=nrow(age_params %>%
                                             filter(param=="pICU_Hosp")),
                                    J=max(age_params$study[age_params$param=="pICU_Hosp"]),
                                    A=max(age_params$age_grp),
                                    agegrp=age_params$age_grp[age_params$param=="pICU_Hosp"],
                                    id=age_params$study[age_params$param=="pICU_Hosp"],
                                    N=age_params$N[age_params$param=="pICU_Hosp"],
                                    Y=age_params$X[age_params$param=="pICU_Hosp"]),
                          pars=c("prob_age"),
                          iter=n_iter,
                          warmup=n_warmup,
                          control=list(adapt_delta=p_adapt_delta,
                                       max_treedepth=n_max_treedepth),
                          save_warmup=F)

p_symp_inf_age <- extract(symp_inf_est, par="prob_age")[[1]]
apply(p_symp_inf_age, 2, median)

p_death_symp_age <- extract(death_symp_est, par="prob_age")[[1]]
apply(p_death_symp_age, 2, median)

p_hosp_symp_age <- extract(hosp_symp_est, par="prob_age")[[1]]
apply(p_hosp_symp_age, 2, median)

p_icu_hosp_age <- extract(icu_hosp_est, par="prob_age")[[1]]
apply(p_icu_hosp_age, 2, median)

p_vent_symp <- rbinom(n=nrow(p_symp_inf_age),
                      size=raw_params$N[raw_params$param=="pVent_Symp"],
                      prob=raw_params$X[raw_params$param=="pVent_Symp"]/raw_params$N[raw_params$param=="pVent_Symp"])/
  raw_params$N[raw_params$param=="pVent_Symp"]

