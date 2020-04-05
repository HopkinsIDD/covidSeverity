## set up data for Stan model
##

options(mc.cores=4,
        scipen=999)
library(googlesheets4)
library(tidyverse)
library(rstan)

n_iter <- 1e4
n_warmup <- 1e4-5e3
p_adapt_delta <- 0.8
n_max_treedepth <- 10

age_specific_data <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1pDfQ2SkO--F2WNfjZxx6t9V7K51IZW0bJKV5cwSCZfA/edit#gid=1769840547",
                                               sheet="age risk") %>%
  filter(USE_pipeline==T)

## look at one conditional probability (e.g. proportion symptomatic pSymp_Inf)
raw_params <- age_specific_data %>%
  # filter(param=="pSymp_Inf") %>%
  mutate(N = ifelse(is.na(N),
                    ifelse(!is.na(X) & value>0, round(X/value),
                           ifelse(!is.na(valueR),
                                  round(value*(1-value)*(qnorm(.975)/(valueR-value))^2),
                                  1000)),
                    N),
         X = ifelse(is.na(X),round(N*value),X),
         age_cat = paste0(ageL,"_",ageR))

## reallocate cases to new age groups
## first expand summarized data
expanded_dat <- c()
## do each param separately
params <- unique(raw_params$param)
num_params <- length(params)
for(i in 1:num_params){
  param_dat <- filter(raw_params, param==params[i])
  # if(nrow(param_dat)==1){
  #   next
  # }
  ## keep the studies separate
  studies <- unique(param_dat$source)
  num_studies <- length(studies)
  for(j in 1:num_studies){
    tmp <- param_dat %>% filter(source==studies[j])
    ## break out each age group
    age_groups <- unique(tmp$age_cat)
    num_ag <- length(age_groups)
    for(k in 1:num_ag){
      tmp_ag <- tmp %>% filter(age_cat==age_groups[k])

      expanded_dat <- bind_rows(expanded_dat,
                                tibble(param=tmp_ag$param[1],
                                       study=j,
                                       ageL=tmp_ag$ageL,
                                       ageR=tmp_ag$ageR,
                                       x=sample(c(rep(1,tmp_ag$X[1]),
                                                  rep(0,tmp_ag$N[1]-tmp_ag$X[1])),
                                                tmp_ag$N[1], replace=F)))
    }
  }
}
# expanded_dat <- expanded_dat %>%
#   mutate(age_cat=cut(age,seq(0,100,by=10),right = FALSE),
#          age_grp=as.factor(age_cat) %>% as.numeric) %>%
#   group_by(param, study, age_cat, age_grp) %>%
#   summarize(X=sum(x),
#             N=n())


age_stan <- stan_model("R/age_prob_para.stan")

symp_inf_est <- sampling(age_stan,
                         data=list(N=nrow(expanded_dat %>%
                                            filter(param=="pSymp_Inf")),
                                   J=max(expanded_dat$study[expanded_dat$param=="pSymp_Inf"]),
                                   id=expanded_dat$study[expanded_dat$param=="pSymp_Inf"],
                                   id_wt=(expanded_dat %>%
                                            filter(param=="pSymp_Inf") %>%
                                            group_by(study) %>%
                                            summarize(id_wt=n()))$id_wt/
                                     nrow(expanded_dat %>%
                                            filter(param=="pSymp_Inf")),
                                   Y=expanded_dat$x[expanded_dat$param=="pSymp_Inf"],
                                   age_min=expanded_dat$ageL[expanded_dat$param=="pSymp_Inf"],
                                   age_max=expanded_dat$ageR[expanded_dat$param=="pSymp_Inf"]),
                         pars=c("prob_age", "log_lik"),
                         iter=n_iter,
                         warmup=n_warmup,
                         control=list(adapt_delta=p_adapt_delta,
                                      max_treedepth=n_max_treedepth),
                         save_warmup=F)
loo::extract_log_lik(symp_inf_est) %>% loo::loo()

p_symp_inf_age <- extract(symp_inf_est, par="prob_age")[[1]]
apply(p_symp_inf_age, 2, median)
rm(symp_inf_est)

death_dat <- expanded_dat %>%
  filter(param=="pDeath_Symp") %>%
  sample_n(1e4) %>%
  mutate(x=as.logical(x)) %>%
  as.data.frame()

library(brms)
death_stan_code <- make_stancode(x ~ s(ageL, bs="cs") + (1|study), data=death_dat, family=bernoulli(), threshold="flexible")

death_symp_est <- sampling(age_stan,
                           data=list(N=nrow(death_dat),
                                     J=max(death_dat$study),
                                     id=death_dat$study,
                                     id_wt=(death_dat %>%
                                              group_by(study) %>%
                                              summarize(id_wt=n()))$id_wt/
                                       nrow(death_dat),
                                     Y=death_dat$x,
                                     age_min=death_dat$ageL,
                                     age_max=death_dat$ageR),
                           pars=c("prob_age", "log_lik"),
                           iter=n_iter,
                           warmup=n_warmup,
                           control=list(adapt_delta=p_adapt_delta,
                                        max_treedepth=n_max_treedepth),
                           save_warmup=F)
loo::extract_log_lik(death_symp_est) %>% loo::loo()

p_death_symp_age <- extract(death_symp_est, par="prob_age")[[1]]
apply(p_death_symp_age, 2, median)

hosp_symp_est <- sampling(age_stan,
                          data=list(N=nrow(expanded_dat %>%
                                             filter(param=="pHosp_Symp")),
                                    J=max(expanded_dat$study[expanded_dat$param=="pHosp_Symp"]),
                                    id=expanded_dat$study[expanded_dat$param=="pHosp_Symp"],
                                    N=expanded_dat$N[expanded_dat$param=="pHosp_Symp"],
                                    Y=expanded_dat$X[expanded_dat$param=="pHosp_Symp"]),
                          pars=c("prob_age", "log_lik"),
                          iter=n_iter,
                          warmup=n_warmup,
                          control=list(adapt_delta=p_adapt_delta,
                                       max_treedepth=n_max_treedepth),
                          save_warmup=F)

icu_hosp_est <- sampling(age_stan,
                         data=list(S=nrow(expanded_dat %>%
                                            filter(param=="pICU_Hosp")),
                                   J=max(expanded_dat$study[expanded_dat$param=="pICU_Hosp"]),
                                   A=max(expanded_dat$age_grp),
                                   agegrp=expanded_dat$age_grp[expanded_dat$param=="pICU_Hosp"],
                                   id=expanded_dat$study[expanded_dat$param=="pICU_Hosp"],
                                   N=expanded_dat$N[expanded_dat$param=="pICU_Hosp"],
                                   Y=expanded_dat$X[expanded_dat$param=="pICU_Hosp"]),
                         pars=c("prob_age", "log_lik"),
                         iter=n_iter,
                         warmup=n_warmup,
                         control=list(adapt_delta=p_adapt_delta,
                                      max_treedepth=n_max_treedepth),
                         save_warmup=F)

p_hosp_symp_age <- extract(hosp_symp_est, par="prob_age")[[1]]
apply(p_hosp_symp_age, 2, median)

p_icu_hosp_age <- extract(icu_hosp_est, par="prob_age")[[1]]
apply(p_icu_hosp_age, 2, median)

p_vent_symp <- rbinom(n=nrow(p_symp_inf_age),
                      size=raw_params$N[raw_params$param=="pVent_Symp"],
                      prob=raw_params$X[raw_params$param=="pVent_Symp"]/raw_params$N[raw_params$param=="pVent_Symp"])/
  raw_params$N[raw_params$param=="pVent_Symp"]

