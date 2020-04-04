## set up data for Stan model
##

options(mc.cores=4,
        scipen=999)
library(googlesheets4)
library(tidyverse)
library(rstan)

n_iter <- 1e3
n_warmup <- 1e3-5e2
p_adapt_delta <- 0.8
n_max_treedepth <- 10

age_specific_data <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1pDfQ2SkO--F2WNfjZxx6t9V7K51IZW0bJKV5cwSCZfA/edit#gid=1769840547",
                                               sheet="age risk")

## look at one conditional probability (e.g. proportion symptomatic pSymp_Inf)
raw_symp_inf <- age_specific_data %>%
  filter(USE_pipeline==T, param=="pSymp_Inf") %>%
  mutate(N = ifelse(is.na(N),
                    round(value*(1-value)*(qnorm(.975)/(valueR-value))^2),
                    N),
         X = ifelse(is.na(X),round(N*value),X))

## reallocate cases to new age groups
## first expand summarized data
expanded_symp_inf <- c()
## keep the studies separate
studies <- unique(raw_symp_inf$source)
num_studies <- length(studies)
for(i in 1:num_studies){
  tmp <- raw_symp_inf %>% filter(source==studies[i])
  ## break out each age group
  age_groups <- unique(tmp$ageL)
  num_ag <- length(age_groups)
  for(j in 1:num_ag){
    tmp_ag <- tmp %>% filter(ageL==age_groups[j])

    expanded_symp_inf <- bind_rows(expanded_symp_inf,
                                   tibble(study=i,
                                          param=tmp$param[1],
                                          age=sample(tmp$ageL[j]:tmp$ageR[j],
                                                     tmp$N[j], replace=T),
                                          x=sample(c(rep(1,tmp$X[j]),
                                                     rep(0,tmp$N[j]-tmp$X[j])),
                                                   tmp$N[j], replace=F)))
  }
}

age_symp_inf <- expanded_symp_inf %>%
  mutate(age_cat=cut(age,seq(0,100,by=10),right = FALSE),
         age_grp=as.factor(age_cat) %>% as.numeric) #%>%
  # group_by(study, age_cat) %>%
  # summarize(X = sum(x),
  #           N = n())

age_stan <- stan_model("R/multilevel_age_prob.stan")

symp_inf_est <- sampling(age_stan,
                         data=list(N=nrow(age_symp_inf),
                                   J=max(age_symp_inf$study),
                                   A=max(age_symp_inf$age_grp),
                                   agegrp=age_symp_inf$age_grp,
                                   id=age_symp_inf$study,
                                   Y=age_symp_inf$x
                         ),
                         iter=n_iter,
                         warmup=n_warmup,
                         control=list(adapt_delta=p_adapt_delta,
                                      max_treedepth=n_max_treedepth),
                         save_warmup=F)

p_symp_inf_age <- extract(symp_inf_est, par="prob_age")[[1]]
apply(p_symp_inf_age, 2, mean)
