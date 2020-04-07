
##' Function to standardize a given probability by age for each
##' GEOID - assuming a single vector of age-specific probabilities
##'
##' @param age_county_pop age specific population by county
##' @param p_vec age specific probability of event
##' @param var_name variable name of standardized pop
##'
##' @return p, standardized probability for each GEOID
##'
##' @import tidyverse
##'
##' @export

get_p_standardized <- function(age_county_pop, p_vec, var_name){

    p_age <- age_county_pop %>%
             select(GEOID, cat_l, page) %>%
             pivot_wider(names_from = cat_l, values_from = page)
    GEOID <- p_age[,1]
    p_age <- as.matrix(p_age[,-1])

    p_tmp <- sweep(p_age, 2, p_vec, "*")
    p_stand <- rowSums(p_tmp)

    return(p_severe_)
}


#### OLD ####


##' estimates probability of being a severe case by age from Shenzen data
##' results from this STAN model saved for easy use later on
# get_severe_age_shenzhen <- function( ){
#
    # sym_dat <- data.frame(
    #     age_cat=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70+"),
    #     mild=c(7,3,13,22,15,21,17,4),
    #     moderate=c(13,9,21,64,40,46,49,12),
    #     severe=c(0,0,0,1,5,7,20,2)
    # )
    # fev_dat <- data.frame(
    #     age_cat=c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70+"),
    #     no=c(6,3,3,13,6,10,16,4),
    #     yes=c(14,9,31,74,54,64,70,14)
    # )
#
#
    # dat <- full_join(sym_dat, fev_dat)
    # dat <- dat %>% as.data.frame() %>% mutate(tot = yes+no)
    # dat <- dat %>% mutate(not_severe = tot-severe)
    # dat <- dat %>% mutate(p_severe = severe/tot)
    #
    # # Use stan
    #
    # library(rstanarm)
    #
    # t_prior <- student_t(df = 7, location = 0, scale = 2.5, autoscale = FALSE)
    # fit1 <- stan_glm(cbind(severe, tot-severe) ~ age_cat, data = dat,
    #                  family = binomial(link = "logit"),
    #                  prior = t_prior, prior_intercept = t_prior,
    #                  cores = 4, seed = 12345)
    #
    # PPD <- posterior_predict(fit1)
    # prob <- PPD
    # for(i in 1:nrow(PPD)){
    #     prob[i,] <- PPD[i,] / dat$tot
    # }
#
#     write_csv(prob %>% as.data.frame(), "data/severe_age_prob.csv")
#     return(prob)
#
# }
#
#

##'
##' Get population distribution and aggregate it to 10 year age groups
##' - this is set up to use population estimates from the World Populaiton Prospects data
##'
##' @param country country of interest
##'
# get_age_pop <- function(country){
#
#     require(stringi)
#     #require(globaltoolbox)
#
#     pop_data <- read_csv("data/WPP2019_POP.csv")
#     pop_data <- pop_data %>%
#         mutate(country_clean = stringi::stri_trans_general(location, "Latin-ASCII")) %>%
#         filter(tolower(country_clean) == tolower(country)) %>% filter(year==max(year))
#
#     # print for a double check
#     print(pop_data$location)
#     pop_data <- pop_data[,-(1:4)] %>% select(-country_clean)
#     dat <- as.numeric(pop_data)
#     names(dat) <- colnames(pop_data)
#     return(dat)
# }





##'
##' Get population distribution and aggregate it to 10 year age groups
##'  - this is set up to use population estimates from the World Populaiton Prospects data
##'
##' @param country country of interest
##'
# get_p_severe <- function(country="China"){
#
#     # Load prob(severe | age) from shenzhen
#     prob <- read_csv("data/severe_age_prob.csv")
#
#
#     #  population by age
#     nage_ <- get_age_pop(country) * 1000
#     nage_[8] <- sum(nage_[8:11])
#     nage_ <- nage_[1:8]
#     pr_age10_ <- nage_ / sum(nage_)
#
#     p_severe_tmp <- prob
#     for(i in 1:nrow(prob)){
#         p_severe_tmp[i,] <- prob[i,] * pr_age10_
#     }
#     p_severe_ <- rowSums(p_severe_tmp)
#
#     fit_ <- fitdistrplus::fitdist(p_severe_, "gamma", "mle")
#
#
#     p_severe_ <- list(ests = p_severe_,
#                       mean=mean(p_severe_),
#                       ll=quantile(p_severe_, .025),
#                       ul=quantile(p_severe_, .975),
#                       q25=quantile(p_severe_, .25),
#                       q75=quantile(p_severe_, .75),
#                       shape = coef(fit_)["shape"],
#                       rate = coef(fit_)["rate"])
#
#     return(p_severe_)
# }



##' estimates proportion of severe cases adjusted for population structure
##'
##' @param pr_age10 proportion of population in 10 year age bins
##'
# get_p_severe_pop <- function(pr_age10){
#
#     # Load prob(severe | age) from shenzhen
#     prob <- read_csv("data/severe_age_prob.csv")
#
#     #  sum all proportion of age old than 70
#     pr_age10[8] <- sum(pr_age10[8:length(pr_age10)])
#
#     p_severe_tmp <- prob
#     for(i in 1:nrow(prob)){
#         p_severe_tmp[i,] <- prob[i,] * pr_age10
#     }
#     p_severe_ <- rowSums(p_severe_tmp)
#
#     fit_ <- fitdistrplus::fitdist(p_severe_, "gamma", "mle")
#
#
#     p_severe_ <- list(ests = p_severe_,
#                       mean=mean(p_severe_),
#                       ll=quantile(p_severe_, .025),
#                       ul=quantile(p_severe_, .975),
#                       q25=quantile(p_severe_, .25),
#                       q75=quantile(p_severe_, .75),
#                       shape = coef(fit_)["shape"],
#                       rate = coef(fit_)["rate"])
#
#     return(p_severe_)
# }


#' Estimate age-specific parameter
#'
#' @param expanded_dat data.frame with following columns:
#' \itemize{
#'   \item \code{param} name of the parameter to estimate
#'   \item \code{study} a number representing the study from which the observation came
#'   \item \code{x} outcome of the outcome (for now a binary 0 or 1)
#'   }
#' @param param_to_est character, one of the parameter names in expanded_dat$param
#' @param age_cats numeric vector, the cutoff values for each age group
#'
#' @return dataframe with the following columns
#' \itemize{
#'   \item \code{age_grp} range of ages (derived from age_cats)
#'   \item \code{logit_mean} the mean of a logit-normal distribution for param_to_est
#'   \item \code{logit_sd} the sd of a logit-normal distribution for param_to_est
#' }
#' @export
#'
#' @examples
est_age_spec_param <- function(expanded_dat,
                               param_to_est="p_symp_inf",
                               age_cats=c(seq(0,80,by=10),100),
                               n_preds,
                               study_wt="none" #c("n", "root_n", "equal"))
                               ) {
  if(!(param_to_est %in% expanded_dat$param)){
    stop("param_to_est must exist within expanded_dat$param")
  }
  param_dat <- expanded_dat %>%
    filter(param==param_to_est)

  studies <- unique(param_dat$study)
  num_studies <- length(studies)
  if(study_wt=="none"){
    param_pred_dat <- tibble(study=as.factor(1),
                             age=0:99,
                             wt=1)
  } else{
    param_pred_dat <- c()
    for(i in 1:num_studies){
      tmp <- filter(param_dat, study==studies[i])
      param_pred_dat <- bind_rows(param_pred_dat,
                                  tibble(study=studies[i],
                                         age=min(tmp$age):max(tmp$age),
                                         wt=ifelse(study_wt=="n",
                                                   nrow(tmp)/nrow(param_dat),
                                                   ifelse(study_wt=="root_n",
                                                          sqrt(nrow(tmp)),
                                                          1/num_studies))))
    }
  }
  if(num_studies > 1){
    param_gam <- gam(x~ s(age, bs="cs") + s(study, bs="re"),
                     data=param_dat %>% mutate(study=as.factor(study)),
                     family=binomial())

    pred_terms <- predict(param_gam,
                          param_pred_dat %>% mutate(study=as.factor(study)),
                          type="lpmatrix")
  } else{
    param_gam <- gam(x~s(age, bs="cs"),
                     data=param_dat, family=binomial())

    pred_terms <- predict(param_gam, param_pred_dat, type="lpmatrix")
  }

  coef_perms <- rmvn(n_preds, coef(param_gam), param_gam$Vp)

  if(study_wt=="none"){
    preds <-  t(coef_perms[,1:(ncol(coef_perms)-num_studies)] %*%
                  t(pred_terms[,1:(ncol(coef_perms)-num_studies)])) %>%
      as_tibble() %>%
      bind_cols(param_pred_dat) %>%
      pivot_longer(cols=starts_with("V"),
                   names_to="pred",
                   values_to="est") %>%
      group_by(age, pred) %>%
      summarize(wt_est=weighted.mean(plogis(est), wt)) %>%
      mutate(age_grp=cut(age, age_cats, right=F)) %>%
      group_by(age_grp,pred) %>%
      summarize(est=mean(wt_est)) %>%
      ungroup() %>%
      pivot_wider(names_from=pred,
                  values_from=est) %>%
      select(-age_grp) %>%
      as.matrix()
  } else{
    preds <-  t(coef_perms %*% t(pred_terms)) %>%
      as_tibble() %>%
      bind_cols(param_pred_dat) %>%
      pivot_longer(cols=starts_with("V"),
                   names_to="pred",
                   values_to="est") %>%
      group_by(age, pred) %>%
      summarize(wt_est=weighted.mean(plogis(est), wt)) %>%
      mutate(age_grp=cut(age, age_cats, right=F)) %>%
      group_by(age_grp,pred) %>%
      summarize(est=mean(wt_est)) %>%
      ungroup() %>%
      pivot_wider(names_from=pred,
                  values_from=est) %>%
      select(-age_grp) %>%
      as.matrix()
  }

  pred_summary <- preds %>%
    as_tibble() %>%
    mutate(age_grp=cut(age_cats[1:(length(age_cats)-1)], age_cats, right=F)) %>%
    pivot_longer(cols=starts_with("V"),
                 names_to="pred",
                 values_to="est") %>%
    group_by(age_grp) %>%
    summarize(est_med=median(est),
              est_lb=quantile(est, probs=.025),
              est_ub=quantile(est, probs=.975))
  colnames(pred_summary) <- gsub("est", param_to_est, colnames(pred_summary))
  # preds <-  (coef_perms %*% t(pred_terms)) %>%
  #   as_tibble() %>%
  #   pivot_longer(cols=everything(),
  #                names_to="pred",
  #                values_to="est") %>%
  #   mutate(pred=as.numeric(pred)) %>%
  #   group_by(pred) %>%
  #   summarize(logit_mean=mean(est),
  #             logit_sd=sd(est))
  # qplot(x=param_pred_dat$age, y=plogis(preds$logit_mean), color=factor(param_pred_dat$study)) + theme(legend.position="none")

  # param_preds <- param_pred_dat %>%
  #   mutate(logit_mean=preds$logit_mean,
  #          logit_sd=preds$logit_sd,
  #          age_grp=cut(age,age_cats,right = FALSE)) %>%
  #   group_by(age_grp, age) %>%
  #   summarize(wt_logit_mean=weighted.mean(logit_mean, wt),
  #             wt_logit_sd=sqrt(weighted.mean(logit_mean^2+logit_sd^2, wt)-weighted.mean(logit_mean, wt)^2)) %>%
  #   group_by(age_grp) %>%
  #   summarize(logit_mean=mean(wt_logit_mean),
  #             logit_sd=sqrt(mean(wt_logit_mean^2+wt_logit_sd^2)-mean(wt_logit_mean)^2))
  # ggplot(data=param_preds, aes(x=age_grp)) +
  #   geom_errorbar(aes(ymin=plogis(logit_mean-1.96*logit_sd), ymax=plogis(logit_mean+1.96*logit_sd)), alpha=0.4) +
  #   geom_point(aes(y=plogis(logit_mean)))
  return(list(pred_mtx=preds, pred_sum=pred_summary, param_to_est=param_to_est))
}

#' Find parameter values for GEOIDs
#'
#' @param age_params list, output of est_age_spec_param()
#' @param geoid_age_mtx matrix of population by age and GEOID (e.g. data(US_age_geoid_pct))
#'
#' @return dataframe with GEOIDs with the median probability for the given parameter
#' @export
#'
#' @examples
est_geoid_params <- function(age_params,
                             geoid_age_mtx=US_age_geoid_pct){
  geoid_preds <- geoid_age_mtx %*% age_params$pred_mtx %>%
    as_tibble(rownames="geoid") %>%
    pivot_longer(cols=-geoid,
                 names_to="pred",
                 values_to="est") %>%
    group_by(geoid) %>%
    summarize(p_est=median(est))
  colnames(geoid_preds)[2] <- age_params$param_to_est
  return(geoid_preds)
}

#' Estimate relative rates for GEOIDs
#'
#' @param pred_mtx matrix of predicted parameter values by age category and simulation from est_age_spec_param()$pred_mtx
#' @param param_to_est character string to name output column
#' @param geoid_age_mtx matrix of proportion of population by age category and GEOID (e.g. data(US_age_geoid_pct))
#' @param geoid_pops matrix of total population by age category and GEOID (e.g. data(US_age_geoid_pop))
#'
#' @return dataframe with GEOIDs and relative rates
#' @export
#'
#' @examples
est_geoid_rrs <- function(pred_mtx,
                          param_to_est,
                          geoid_age_mtx=US_age_geoid_pct,
                          geoid_pops=US_age_geoid_pop){
  ## calc nationwide rate per sim
  nationwide_rate <- geoid_pops %*% pred_mtx %>%
    as_tibble(rownames="geoid") %>%
    pivot_longer(cols=-geoid,
                 names_to="pred",
                 values_to="est") %>%
    group_by(pred) %>%
    summarize(total_est=sum(est)) %>%
    mutate(total_pop=sum(US_age_geoid_pop),
           overall_rate=total_est/total_pop)

  ## get geoid relative rates per sim
  geoid_rrs <- geoid_age_mtx %*% pred_mtx %>%
    as_tibble(rownames="geoid") %>%
    pivot_longer(cols=-geoid,
                 names_to="pred",
                 values_to="est") %>%
    left_join(nationwide_rate, by="pred") %>%
    group_by(geoid) %>%
    summarize(rr_est=median(est/overall_rate),
              overall_est=median(overall_rate))
  if(!missing(param_to_est))
    colnames(geoid_rrs)[2:3] <- c(paste0("rr_", param_to_est),
                                  paste0(param_to_est, "_overall"))
  return(geoid_rrs)
}
