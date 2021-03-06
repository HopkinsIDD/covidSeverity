## estimate conditional probabilities

options(scipen=999)
library(covidSeverity)
library(doMC)
registerDoMC(10)
n_sims <- 100
n_preds <- 1e3
age_grps <- c(seq(0,80,by=10),100)

data("US_age_geoid_pct")
data("US_age_geoid_pop")
## if have access to private files, use this
# raw_age_estimates <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1pDfQ2SkO--F2WNfjZxx6t9V7K51IZW0bJKV5cwSCZfA/edit#gid=1769840547",
#                                                sheet="age risk")
# 1
## otherwise use this
data("raw_age_estimates")

## look at one conditional probability (e.g. proportion symptomatic pSymp_Inf)
raw_params <- raw_age_estimates %>%
  filter(USE_pipeline==T) %>%
  mutate(p = ifelse(!is.na(p), p,
                    ifelse(!is.na(X) & !is.na(N), X/N,
                           ifelse(!is.na(pR) & !is.na(pL),
                                  (pL+pR)/2, NA))),
         N = ifelse(!is.na(N), N,
                    ifelse(!is.na(X) & p>0, round(X/p),
                           ifelse(!is.na(pR),
                                  round(p*(1-p)*(qnorm(.975)/(pR-p))^2),
                                  1000))),
         X = ifelse(is.na(X),ceiling(N*p),X))

geoid_symp_inf <- foreach(h=1:n_sims, .combine=rbind) %dopar% {
  set.seed(h)
  ## make proportion symptomatic given infected
  tmp <- est_age_param(age_cat_data=raw_params,
                          param_to_est="p_symp_inf",
                          age_cats=age_grps,
                          n_preds=n_preds,
                          study_wt="none")
  est_geoid_rrs(tmp$pred_mtx,
                param_to_est="p_symp_inf",
                geoid_age_mtx=US_age_geoid_pct,
                geoid_pops=US_age_geoid_pop) %>%
    mutate(sim=h)
}
symp_sim_med <- (geoid_symp_inf %>%
                   select(sim, p_symp_inf_overall) %>%
                   distinct() %>%
                   mutate(dist_from_med=abs(p_symp_inf_overall-median(p_symp_inf_overall))) %>%
                   arrange(dist_from_med) %>%
                   slice(1))$sim

## CFR
geoid_death_symp <- foreach(h=1:n_sims, .combine=rbind) %dopar% {
  set.seed(h)
  ## make proportion symptomatic given infected
  tmp <- est_age_param(age_cat_data=raw_params,
                       param_to_est="p_death_symp",
                       age_cats=age_grps,
                       n_preds=n_preds,
                       study_wt="none")
  est_geoid_rrs(tmp$pred_mtx,
                param_to_est="p_death_symp",
                geoid_age_mtx=US_age_geoid_pct,
                geoid_pops=US_age_geoid_pop) %>%
    mutate(sim=h)
}
death_sim_med <- (geoid_death_symp %>%
                    select(sim, p_death_symp_overall) %>%
                    distinct() %>%
                    mutate(dist_from_med=abs(p_death_symp_overall-median(p_death_symp_overall))) %>%
                    arrange(dist_from_med) %>%
                    slice(1))$sim

## hospitalization rate amongst symptomatic
geoid_hosp_symp <- foreach(h=1:n_sims, .combine=rbind) %dopar% {
  set.seed(h)
  ## make proportion symptomatic given infected
  tmp <- est_age_param(age_cat_data=raw_params,
                       param_to_est="p_hosp_symp",
                       age_cats=age_grps,
                       n_preds=n_preds,
                       study_wt="none")
  est_geoid_rrs(tmp$pred_mtx,
                param_to_est="p_hosp_symp",
                geoid_age_mtx=US_age_geoid_pct,
                geoid_pops=US_age_geoid_pop) %>%
    mutate(sim=h)
}
hosp_sim_med <- (geoid_hosp_symp %>%
                    select(sim, p_hosp_symp_overall) %>%
                    distinct() %>%
                    mutate(dist_from_med=abs(p_hosp_symp_overall-median(p_hosp_symp_overall))) %>%
                    arrange(dist_from_med) %>%
                    slice(1))$sim


## ICU rate amongst hospitalized
geoid_icu_hosp <- foreach(h=1:n_sims, .combine=rbind) %dopar% {
  set.seed(h)
  ## make proportion symptomatic given infected
  tmp <- est_age_param(age_cat_data=raw_params,
                       param_to_est="p_icu_hosp",
                       age_cats=age_grps,
                       n_preds=n_preds,
                       study_wt="none")
  est_geoid_rrs(tmp$pred_mtx,
                param_to_est="p_icu_hosp",
                geoid_age_mtx=US_age_geoid_pct,
                geoid_pops=US_age_geoid_pop) %>%
    mutate(sim=h)
}
icu_sim_med <- (geoid_icu_hosp %>%
                   select(sim, p_icu_hosp_overall) %>%
                   distinct() %>%
                   mutate(dist_from_med=abs(p_icu_hosp_overall-median(p_icu_hosp_overall))) %>%
                   arrange(dist_from_med) %>%
                   slice(1))$sim

## ventilation rate amongst symptomatic
geoid_vent_icu <- foreach(h=1:n_sims, .combine=rbind) %dopar% {
  set.seed(h)
  ## make proportion symptomatic given infected
  tmp <- est_age_param(age_cat_data=raw_params,
                       param_to_est="p_vent_icu",
                       age_cats=age_grps,
                       n_preds=n_preds,
                       study_wt="none")
  est_geoid_rrs(tmp$pred_mtx,
                param_to_est="p_vent_icu",
                geoid_age_mtx=US_age_geoid_pct,
                geoid_pops=US_age_geoid_pop) %>%
    mutate(sim=h)
}
vent_sim_med <- (geoid_vent_icu %>%
                   select(sim, p_vent_icu_overall) %>%
                   distinct() %>%
                   mutate(dist_from_med=abs(p_vent_icu_overall-median(p_vent_icu_overall))) %>%
                   arrange(dist_from_med) %>%
                   slice(1))$sim

## get all the median sims
## make proportion symptomatic given infected
set.seed(symp_sim_med)
p_symp_inf <- est_age_param(age_cat_data=raw_params,
                        param_to_est="p_symp_inf",
                        age_cats=age_grps,
                        n_preds=n_preds,
                        study_wt="none")

## make proportion death giv0000en symptomatic
set.seed(death_sim_med)
p_death_symp <- est_age_param(age_cat_data=raw_params,
                        param_to_est="p_death_symp",
                        age_cats=age_grps,
                        n_preds=n_preds,
                        study_wt="none")

## make proportion hospitalized given symptomatic
set.seed(hosp_sim_med)
p_hosp_symp <- est_age_param(age_cat_data=raw_params,
                        param_to_est="p_hosp_symp",
                        age_cats=age_grps,
                        n_preds=n_preds,
                        study_wt="none")

## make proportion ICU given hospitalized
set.seed(icu_sim_med)
p_icu_hosp <- est_age_param(age_cat_data=raw_params,
                       param_to_est="p_icu_hosp",
                       age_cats=age_grps,
                       n_preds=n_preds,
                       study_wt="none")

## make proportion ventilated given ICU
set.seed(vent_sim_med)
p_vent_icu <- est_age_param(age_cat_data=raw_params,
                       param_to_est="p_vent_icu",
                       age_cats=age_grps,
                       n_preds=n_preds,
                       study_wt="none")


## get all output parameters
US_geoid_params <- est_geoid_params(p_symp_inf$pred_mtx,
                                    param_to_est=p_symp_inf$param_to_est,
                                    US_age_geoid_pct) %>%
  left_join(est_geoid_params(p_death_symp$pred_mtx * p_symp_inf$pred_mtx,
                             param_to_est="p_death_inf",
                             US_age_geoid_pct), by="geoid") %>%
  left_join(est_geoid_params(p_hosp_symp$pred_mtx * p_symp_inf$pred_mtx,
                             param_to_est="p_hosp_inf",
                             US_age_geoid_pct), by="geoid") %>%
  left_join(est_geoid_params(p_icu_hosp$pred_mtx * p_hosp_symp$pred_mtx *
                               p_symp_inf$pred_mtx,
                             param_to_est="p_icu_inf",
                             US_age_geoid_pct), by="geoid") %>%
  left_join(est_geoid_params(p_vent_icu$pred_mtx * p_icu_hosp$pred_mtx *
                               p_hosp_symp$pred_mtx * p_symp_inf$pred_mtx,
                             param_to_est="p_vent_inf",
                             US_age_geoid_pct), by="geoid") %>%
  mutate(p_death_symp=p_death_inf/p_symp_inf,
         p_hosp_symp=p_hosp_inf/p_symp_inf,
         p_mild_symp=1-p_hosp_symp,
         p_icu_hosp=p_icu_inf/p_hosp_inf,
         p_vent_icu=p_vent_inf/p_icu_inf) %>%
  left_join(est_geoid_rrs(pred_mtx=p_death_symp$pred_mtx * p_symp_inf$pred_mtx,
                          param_to_est="death_inf",
                          geoid_age_mtx=US_age_geoid_pct,
                          geoid_pops=US_age_geoid_pop), by="geoid") %>%
  left_join(est_geoid_rrs(pred_mtx=p_hosp_symp$pred_mtx * p_symp_inf$pred_mtx,
                          param_to_est="hosp_inf",
                          geoid_age_mtx=US_age_geoid_pct,
                          geoid_pops=US_age_geoid_pop), by="geoid") %>%
  left_join(est_geoid_rrs(pred_mtx=p_symp_inf$pred_mtx,
                          param_to_est="symp_inf",
                          geoid_age_mtx=US_age_geoid_pct,
                          geoid_pops=US_age_geoid_pop), by="geoid") %>%
  left_join(est_geoid_rrs(pred_mtx=p_icu_hosp$pred_mtx * p_hosp_symp$pred_mtx *
                            p_symp_inf$pred_mtx,
                          param_to_est="icu_inf",
                          geoid_age_mtx=US_age_geoid_pct,
                          geoid_pops=US_age_geoid_pop), by="geoid") %>%
  left_join(est_geoid_rrs(pred_mtx=p_vent_icu$pred_mtx * p_icu_hosp$pred_mtx *
                            p_hosp_symp$pred_mtx * p_symp_inf$pred_mtx,
                          param_to_est="vent_inf",
                          geoid_age_mtx=US_age_geoid_pct,
                          geoid_pops=US_age_geoid_pop), by="geoid") %>%
  mutate(p_ihiv=p_vent_inf,
         p_ihi=p_icu_inf-p_vent_inf,
         p_ih=p_hosp_inf-p_icu_inf,
         p_i=1-p_hosp_inf)

## county parameter distributions
US_geoid_params %>%
  summary()

## save a csv file
US_geoid_params %>%
  write_csv("generated_data/geoid-params.csv")

## save as a data file that can be called from package
usethis::use_data(US_geoid_params, overwrite=T)

## rename some variables for use with outputs.py
US_output_geoid_params <- US_geoid_params %>% 
  left_join(est_geoid_rrs(pred_mtx=p_vent_icu$pred_mtx,
                          param_to_est="vent_icu",
                          geoid_age_mtx=US_age_geoid_pct,
                          geoid_pops=US_age_geoid_pop), by="geoid") %>%
  left_join(est_geoid_rrs(pred_mtx=p_icu_hosp$pred_mtx,
                          param_to_est="icu_hosp",
                          geoid_age_mtx=US_age_geoid_pct,
                          geoid_pops=US_age_geoid_pop), by="geoid") %>%
                          dplyr::rename(`PincidH|incidence` = p_hosp_inf,
                                        `PincidICU|incidH` = p_icu_hosp,
                                        `PincidVent|incidICU` = p_vent_icu,
                                        `PincidD|incidence` = p_death_inf,
                                        `RincidH|incidence` = rr_hosp_inf,
                                        `RincidICU|incidH` = rr_icu_hosp,
                                        `RincidVent|incidICU` = rr_vent_icu,
                                        `RincidD|incidence` = rr_death_inf) %>%
                          dplyr::select(geoid, starts_with("P", ignore.case = F), starts_with("R", ignore.case = F))
## save a csv file
US_output_geoid_params %>%
  write_csv("generated_data/geoid-params-output.csv")

## save as a data file that can be called from package
usethis::use_data(US_output_geoid_params, overwrite=T)

## Make a long version of the output params
US_output_geoid_params_long <- US_output_geoid_params %>%
                               pivot_longer(-geoid, names_to = "varname", values_to = "value") %>%
                               mutate(quantity=substr(varname, 1, 1),
                                      varname=substr(varname, 2, 20)) %>%
                               separate(varname, into = c("outcome", "source"))
## save a csv file
US_output_geoid_params_long %>%
  write_csv("generated_data/geoid-params-output-long.csv")

## save as a data file that can be called from package
usethis::use_data(US_output_geoid_params_long, overwrite=T)

## Symptomatic plot
ggplot(data=p_symp_inf$pred_sum, aes(x=age_grp)) +
  geom_errorbar(aes(ymin=p_symp_inf_lb,
                    ymax=p_symp_inf_ub), alpha=0.4) +
  geom_point(aes(y=p_symp_inf_med)) +
  scale_y_continuous("Probability symptomatic, given infected") +
  scale_x_discrete("Age groups") +
  theme_bw()

## CFR
ggplot(data=p_death_symp$pred_sum, aes(x=age_grp)) +
  geom_errorbar(aes(ymin=p_death_symp_lb,
                    ymax=p_death_symp_ub), alpha=0.4) +
  geom_point(aes(y=p_death_symp_med)) +
  scale_y_continuous("Case fatality rate") +
  scale_x_discrete("Age groups") +
  theme_bw()

## IFR
p_death_inf <- list(pred_mtx=p_death_symp$pred_mtx * p_symp_inf$pred_mtx,
                   param_to_est="p_death_inf")
p_death_inf$pred_sum <- p_death_inf$pred_mtx %>%
  as_tibble() %>%
  mutate(age_grp=cut(age_grps[1:(length(age_grps)-1)], age_grps, right=F)) %>%
  pivot_longer(cols=starts_with("V"),
               names_to="pred",
               values_to="est") %>%
  group_by(age_grp) %>%
  summarize(est_mean=mean(est),
            est_med=median(est),
            est_lb=quantile(est, probs=.025),
            est_ub=quantile(est, probs=.975))

## hospitalization rate amongst symptomatic
ggplot(data=p_hosp_symp$pred_sum, aes(x=age_grp)) +
  geom_errorbar(aes(ymin=p_hosp_symp_lb,
                    ymax=p_hosp_symp_ub), alpha=0.4) +
  geom_point(aes(y=p_hosp_symp_med)) +
  scale_y_continuous("Probability hospitalized, given symptomatic") +
  scale_x_discrete("Age groups") +
  theme_bw()

p_mild_symp <- list(pred_mtx=1-p_hosp_symp$pred_mtx,
                    param_to_est="p_mild_symp")
p_mild_symp$pred_sum <- p_mild_symp$pred_mtx %>%
  as_tibble() %>%
  mutate(age_grp=cut(age_grps[1:(length(age_grps)-1)], age_grps, right=F)) %>%
  pivot_longer(cols=starts_with("V"),
               names_to="pred",
               values_to="est") %>%
  group_by(age_grp) %>%
  summarize(est_mean=mean(est),
            est_med=median(est),
            est_lb=quantile(est, probs=.025),
            est_ub=quantile(est, probs=.975))

## hospitalization rate amongst infections
p_hosp_inf <- list(pred_mtx=p_hosp_symp$pred_mtx * p_symp_inf$pred_mtx,
                   param_to_est="p_hosp_inf")
p_hosp_inf$pred_sum <- p_hosp_inf$pred_mtx %>%
  as_tibble() %>%
  mutate(age_grp=cut(age_grps[1:(length(age_grps)-1)], age_grps, right=F)) %>%
  pivot_longer(cols=starts_with("V"),
               names_to="pred",
               values_to="est") %>%
  group_by(age_grp) %>%
  summarize(est_mean=mean(est),
            est_med=median(est),
            est_lb=quantile(est, probs=.025),
            est_ub=quantile(est, probs=.975))
ggplot(data=p_hosp_inf$pred_sum, aes(x=age_grp)) +
  geom_errorbar(aes(ymin=est_lb,
                    ymax=est_ub), alpha=0.4) +
  geom_point(aes(y=est_med)) +
  scale_y_continuous("Probability hospitalized, given infection") +
  scale_x_discrete("Age groups") +
  theme_bw()

## ICU rate amongst hospitalized
ggplot(data=p_icu_hosp$pred_sum, aes(x=age_grp)) +
  geom_errorbar(aes(ymin=p_icu_hosp_lb,
                    ymax=p_icu_hosp_ub), alpha=0.4) +
  geom_point(aes(y=p_icu_hosp_med)) +
  scale_y_continuous("Probability ICU, given hospitalized") +
  scale_x_discrete("Age groups") +
  theme_bw()

## ventilation rate amongst symptomatic
ggplot(data=p_vent_icu$pred_sum, aes(x=age_grp)) +
  geom_errorbar(aes(ymin=p_vent_icu_lb,
                    ymax=p_vent_icu_ub), alpha=0.4) +
  geom_point(aes(y=p_vent_icu_med)) +
  scale_y_continuous("Probability invasive ventilation, given ICU") +
  scale_x_discrete("Age groups") +
  theme_bw()

## save as list
saveRDS(list(p_symp_inf,
             p_death_symp,
             p_death_inf,
             p_hosp_symp,
             p_hosp_inf,
             p_icu_hosp,
             p_vent_icu),
        file="generated_data/param-age-dist.rds")
