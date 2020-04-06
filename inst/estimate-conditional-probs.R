## set up data for Stan model
##

options(mc.cores=4,
        scipen=999)
library(googlesheets4)
library(tidyverse)
library(mgcv)
# library(rstan)
# library(brms)

# n_iter <- 2e3
# n_warmup <- 1e3
# p_adapt_delta <- 0.99
# n_max_treedepth <- 20

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
  # if(nrow(param_dat)==1){
  #   next
  # }
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
                                       # ageL=tmp_ag$ageL,
                                       # ageR=tmp_ag$ageR,
                                       x=sample(c(rep(1,tmp_ag$X[1]),
                                                  rep(0,tmp_ag$N[1]-tmp_ag$X[1])),
                                                tmp_ag$N[1], replace=F)))
    }
  }
}
expanded_dat <- expanded_dat %>%
  mutate(age_cat=cut(age,c(seq(0,80,by=10),100),right = FALSE),
         age_grp=as.factor(age_cat) %>% as.numeric)# %>%
# group_by(param, study, age_cat, age_grp) %>%
# summarize(X=sum(x),
#           N=n())


age_cat_mins <- c(0,10,20,30,40,50,60,70,80)
age_cat_maxs <- c(9,19,29,39,49,59,69,79,99)

# if(filter(expanded_dat, param=="pSymp_Inf") %>% nrow() > 10000){
#   symp_dat <- expanded_dat %>%
#     filter(param=="pSymp_Inf") %>%
#     sample_n(1e4)
# } else{

# }


est_age_spec_param <- function(expanded_dat,
                               param_to_est="pSymp_Inf",
                               age_cats=c(seq(0,80,by=10),100)){
  param_dat <- expanded_dat %>%
    filter(param==param_to_est)
  studies <- unique(param_dat$study)
  param_pred_dat <- c()
  for(i in studies){
    tmp <- filter(param_dat, study==i)
    param_pred_dat <- bind_rows(param_pred_dat,
                                tibble(study=i,
                                       age=min(tmp$age):max(tmp$age),
                                       wt=nrow(tmp)/nrow(param_dat)))
  }

  param_gam <- gam(x~s(age, bs="cs") + s(study, bs="re"),
                   data=param_dat, family=binomial)

  pred_terms <- predict(param_gam, param_pred_dat, type="lpmatrix")

  coef_perms <- rmvn(100, coef(param_gam), param_gam$Vp)

  preds <-  (coef_perms %*% t(pred_terms)) %>%
    as_tibble() %>%
    pivot_longer(cols=everything(),
                 names_to="pred",
                 values_to="est") %>%
    mutate(pred=as.numeric(pred)) %>%
    group_by(pred) %>%
    summarize(logit_mean=mean(est),
              logit_var=var(est))
  # qplot(x=param_pred_dat$age, y=plogis(preds$logit_mean), color=factor(param_pred_dat$study)) + theme(legend.position="none")


  param_preds <- param_pred_dat %>%
    mutate(logit_mean=preds$logit_mean,
           logit_var=preds$logit_var,
           age_grp=cut(age,age_cats,right = FALSE)) %>%
    group_by(age_grp, age) %>%
    summarize(wt_logit_mean=weighted.mean(logit_mean, wt),
              wt_logit_var=weighted.mean(logit_mean^2+logit_var, wt)-weighted.mean(logit_mean, wt)^2) %>%
    group_by(age_grp) %>%
    summarize(logit_mean=mean(wt_logit_mean),
              logit_sd=sqrt(mean(wt_logit_mean^2+wt_logit_var)-mean(wt_logit_mean)^2))
  # ggplot(data=param_preds, aes(x=age_grp)) +
  #   geom_errorbar(aes(ymin=plogis(logit_mean-1.96*logit_sd), ymax=plogis(logit_mean+1.96*logit_sd)), alpha=0.4) +
  #   geom_point(aes(y=plogis(logit_mean)))
  return(param_preds)
}

p_symp <- est_age_spec_param(expanded_dat)
ggplot(data=p_symp, aes(x=age_grp)) +
  geom_errorbar(aes(ymin=plogis(logit_mean-1.96*logit_sd), ymax=plogis(logit_mean+1.96*logit_sd)), alpha=0.4) +
  geom_point(aes(y=plogis(logit_mean)))

## get the code and data from brms to make
# symp_stan_code <- make_stancode(x ~ s(age, bs="cs", k=3) + (1|study),
#                                 data=symp_dat_preds, family=bernoulli())
# symp_stan_data <- make_standata(x ~ s(age, bs="cs", k=3) + (1|study),
#                                 data=symp_dat_preds, family=bernoulli())

## create basis matrix for predictions
# symp_pred_b_matrix <- symp_dat_preds %>%
#   bind_cols(tibble(b1=symp_stan_data$Zs_1_1[,1],
#                    b2=symp_stan_data$Zs_1_1[,2])) %>%
#   slice((nrow(symp_dat)+1):nrow(symp_dat_preds)) %>%
#   select(-study, -x)

# symp_brms <- brm(x ~ s(age, bs="cs", k=3) + (1|study),
#                  data=symp_dat, family=bernoulli(),
#                  pars=c("prob_age"),
#                  # pars=c("prob_age", "log_lik"),
#                  iter=n_iter,
#                  warmup=n_warmup,
#                  control=list(adapt_delta=p_adapt_delta,
#                               max_treedepth=n_max_treedepth),
#                  save_warmup=F)

## set up stan model
# age_stan <- stan_model("R/age_prob_splines_steve.stan")

# symp_inf_est <- sampling(age_stan,
#                          data= list(N=nrow(symp_dat),
#                                    Y=symp_dat$x,
#                                    knots_1=2,
#                                    Bmat=symp_stan_data$Zs_1_1[1:nrow(symp_dat),],
#                                    J=max(symp_dat$study),
#                                    id=symp_dat$study,
#                                    id_wt=(symp_dat %>%
#                                             group_by(study) %>%
#                                             summarize(id_wt=n()))$id_wt/
#                                      nrow(symp_dat),
#                                    preds=nrow(symp_pred_b_matrix),
#                                    Bmat_preds=symp_pred_b_matrix %>%
#                                      select(-age) %>% as.matrix),
#                          iter=n_iter,
#                          warmup=n_warmup,
#                          control=list(adapt_delta=p_adapt_delta,
#                                       max_treedepth=n_max_treedepth),
#                          save_warmup=F)
# loo::extract_log_lik(symp_inf_est) %>% loo::loo()

# p_symp_inf_age <- extract(symp_inf_est, par="age_preds")[[1]]
# apply(p_symp_inf_age, 2, median)
# rm(symp_inf_est)

# if(filter(expanded_dat, param=="pDeath_Symp") %>% nrow() > 1e4){
#   death_dat <- expanded_dat %>%
#     filter(param=="pDeath_Symp") %>%
#     sample_n(1e4)
# } else{

p_death <- est_age_spec_param(expanded_dat, "pDeath_Symp")
ggplot(data=p_death, aes(x=age_grp)) +
  geom_errorbar(aes(ymin=plogis(logit_mean-1.96*logit_sd), ymax=plogis(logit_mean+1.96*logit_sd)), alpha=0.4) +
  geom_point(aes(y=plogis(logit_mean)))

#
# death_stan_data <- make_standata(x ~ s(age, bs="cs", k=3) + (1|study),
#                                  data=death_dat, family=bernoulli())
# death_symp_est <- sampling(age_stan,
#                            data=list(N=nrow(death_dat),
#                                      J=max(death_dat$study),
#                                      id=death_dat$study,
#                                      id_wt=(death_dat %>%
#                                               group_by(study) %>%
#                                               summarize(id_wt=n()))$id_wt/
#                                        nrow(death_dat),
#                                      Y=death_dat$x,
#                                      age_min=death_dat$ageL,
#                                      age_max=death_dat$ageR,
#                                      preds=length(age_cat_maxs),
#                                      pred_min=age_cat_mins,
#                                      pred_max=age_cat_maxs),
#                            pars=c("prob_age"),
#                            # pars=c("prob_age", "log_lik"),
#                            iter=n_iter,
#                            warmup=n_warmup,
#                            control=list(adapt_delta=p_adapt_delta,
#                                         max_treedepth=n_max_treedepth),
#                            save_warmup=F)
# # loo::extract_log_lik(death_symp_est) %>% loo::loo()
#
# p_death_symp_age <- extract(death_symp_est, par="prob_age")[[1]]
# apply(p_death_symp_age, 2, median)
# rm(death_symp_est)

# if(filter(expanded_dat, param=="pHosp_Symp") %>% nrow > 1e4){
#   hosp_dat <- expanded_dat %>%
#     filter(param=="pHosp_Symp") %>%
#     sample_n(1e4)
# } else{

p_hosp <- est_age_spec_param(expanded_dat, "pHosp_Symp")
ggplot(data=p_hosp, aes(x=age_grp)) +
  geom_errorbar(aes(ymin=plogis(logit_mean-1.96*logit_sd), ymax=plogis(logit_mean+1.96*logit_sd)), alpha=0.4) +
  geom_point(aes(y=plogis(logit_mean)))

# hosp_stan_data <- make_standata(x ~ s(age, bs="cs", k=3) + (1|study),
#                                  data=hosp_dat, family=bernoulli())
#
# hosp_symp_est <- sampling(age_stan,
#                           data=list(N=nrow(hosp_dat),
#                                     J=max(hosp_dat$study),
#                                     id=hosp_dat$study,
#                                     id_wt=(hosp_dat %>%
#                                              group_by(study) %>%
#                                              summarize(id_wt=n()))$id_wt/
#                                       nrow(hosp_dat),
#                                     Y=hosp_dat$x,
#                                     age_min=hosp_dat$ageL,
#                                     age_max=hosp_dat$ageR,
#                                     preds=length(age_cat_maxs),
#                                     pred_min=age_cat_mins,
#                                     pred_max=age_cat_maxs),
#                           pars=c("prob_age"),
#                           # pars=c("prob_age", "log_lik"),
#                           iter=n_iter,
#                           warmup=n_warmup,
#                           control=list(adapt_delta=p_adapt_delta,
#                                        max_treedepth=n_max_treedepth),
#                           save_warmup=F)
# # loo::extract_log_lik(hosp_symp_est) %>% loo::loo()
# p_hosp_symp_age <- extract(hosp_symp_est, par="prob_age")[[1]]
# apply(p_hosp_symp_age, 2, median)
# rm(hosp_symp_est)

# if(filter(expanded_dat, param=="pICU_Hosp") %>% nrow > 1e4){
#   icu_dat <- expanded_dat %>%
#     filter(param=="pICU_Hosp") %>%
#     sample_n(1e4)
# } else{

p_icu <- est_age_spec_param(expanded_dat, "pICU_Hosp")
ggplot(data=p_icu, aes(x=age_grp)) +
  geom_errorbar(aes(ymin=plogis(logit_mean-1.96*logit_sd), ymax=plogis(logit_mean+1.96*logit_sd)), alpha=0.4) +
  geom_point(aes(y=plogis(logit_mean)))

#
# icu_stan_data <- make_standata(x ~ s(age, bs="cs", k=3) + (1|study),
#                                 data=icu_dat, family=bernoulli())
#
# icu_hosp_est <- sampling(age_stan,
#                          data=list(N=nrow(icu_dat),
#                                    J=max(icu_dat$study),
#                                    id=icu_dat$study,
#                                    id_wt=(icu_dat %>%
#                                             group_by(study) %>%
#                                             summarize(id_wt=n()))$id_wt/
#                                      nrow(icu_dat),
#                                    Y=icu_dat$x,
#                                    age_min=icu_dat$ageL,
#                                    age_max=icu_dat$ageR,
#                                    preds=length(age_cat_maxs),
#                                    pred_min=age_cat_mins,
#                                    pred_max=age_cat_maxs),
#                          pars=c("prob_age"),
#                          # pars=c("prob_age", "log_lik"),
#                          iter=n_iter,
#                          warmup=n_warmup,
#                          control=list(adapt_delta=p_adapt_delta,
#                                       max_treedepth=n_max_treedepth),
#                          save_warmup=F)
# # loo::extract_log_lik(icu_hosp_est) %>% loo::loo()
# p_icu_hosp_age <- extract(icu_hosp_est, par="prob_age")[[1]]
# apply(p_icu_hosp_age, 2, median)
# rm(icu_hosp_est)

p_vent_symp <- rbinom(n=nrow(p_symp_inf_age),
                      size=raw_params$N[raw_params$param=="pVent_Symp"],
                      prob=raw_params$X[raw_params$param=="pVent_Symp"]/
                        raw_params$N[raw_params$param=="pVent_Symp"])/
  raw_params$N[raw_params$param=="pVent_Symp"]

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
