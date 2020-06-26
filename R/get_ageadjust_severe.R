
#' Age-specific parameter estimation
#'
#' Fits a cubic spline across multiple studies to estimate the parameter of interest by specified age categories.
#'
#' @param age_cat_data
#' @param param_to_est character, one of the parameter names in age_cat_data$param
#' @param age_cats numeric vector, the cutoff values for each age group
#' @param n_preds numeric, number of predictions for the gam model to make
#' @param study_wt character string, how to weight the random effects for the studies after the model is fit, there are multiple options:
#' \itemize{
#'   \item "none" (default) fits the random effects, but doesn't use them in the predictions of age-specific conditional probabilities
#'   \item "n" weights each random effect by the number of observations in each study
#'   \item "root_n" weights each random effect by the square root of the number of observations
#'   \item "equal" weights each random effect equally
#' }
#' @return list with 3 items
#' \itemize{
#'   \item \code{pred_mtx} matrix, rows are age categories, columns are simulations, values are estimates of the conditional probabilities
#'   \item \code{pred_sum} tbl, with the median, lower and upper bounds of a 95\% CI for each age category
#'   \item \code{param_to_est} character_string, same as in Arguments
#' }
#' @import mgcv
#' @export
#'
#' @examples
#' ## load param data
#' data("raw_age_estimates")
#' ## filter to only those for pipeline
#' raw_params <- filter(raw_age_estimates, USE_pipeline==T)
#'
#' ## run function
#' p_vent <- est_age_param(age_cat_data=raw_params,
#'                         param_to_est="p_vent_icu",
#'                         n_preds=1000)
#'
#' ## view summary
#' p_vent$pred_sum
est_age_param <- function(age_cat_data,
                          param_to_est="p_symp_inf",
                          age_cats=c(seq(0,80,by=10),100),
                          n_preds,
                          study_wt="none" #c("n", "root_n", "equal"))
) {
  if(!(param_to_est %in% age_cat_data$param)){
    stop("param_to_est must exist within expanded_dat$param")
  }
  ## expand categorical age data into individual observations
  ## with randomly assigned ages
  param_dat <- expand_age_data(age_cat_data=age_cat_data,
                               param_to_est=param_to_est)

  ## create a prediction data frame that has the age range for each study
  ## and a weight to apply for each study
  studies <- unique(param_dat$study)
  num_studies <- length(studies)
  ## if study_wt is "none" then the only important aspect of the prediction
  ## set is covering the entire age range
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
    summarize(est_mean=mean(est),
              est_med=median(est),
              est_lb=quantile(est, probs=.025),
              est_ub=quantile(est, probs=.975))
  #colnames(pred_summary) <- gsub("est", param_to_est, colnames(pred_summary))
  return(list(param_to_est=param_to_est, pred_mtx=preds, pred_sum=pred_summary, param_to_est=param_to_est))
}

#' Find parameter values for GEOIDs
#'
#' @param age_params list, output of est_age_param()
#' @param geoid_age_mtx matrix of population by age and GEOID (e.g. data(US_age_geoid_pct))
#'
#' @return tbl, median probabilities for the given parameter by GEOID
#' @export
#'
#' @examples
est_geoid_params <- function(pred_mtx,
                             param_to_est,
                             geoid_age_mtx=US_age_geoid_pct){
  geoid_preds <- geoid_age_mtx %*% pred_mtx %>%
    as_tibble(rownames="geoid") %>%
    pivot_longer(cols=-geoid,
                 names_to="pred",
                 values_to="est") %>%
    group_by(geoid) %>%
    summarize(p_est=median(est))
  colnames(geoid_preds)[2] <- param_to_est
  return(geoid_preds)
}



#' Estimate relative rates for GEOIDs
#'
#' @param pred_mtx matrix of predicted parameter values by age category and simulation from est_age_param()$pred_mtx
#' @param param_to_est character string to name output column
#' @param geoid_age_mtx matrix of proportion of population by age category and GEOID (e.g. data(US_age_geoid_pct))
#' @param geoid_pops matrix of total population by age category and GEOID (e.g. data(US_age_geoid_pop))
#'
#' @return tbl, with relative rates by GEOID and the overall estimate across all GEOIDs
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
    mutate(total_pop=sum(geoid_pops),
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



#' Expanded age category data
#'
#' draw ages from categories for better model fitting
#'
#' @param age_cat_data data.frame
#' @param param_to_est character string of parameter of interest
#'
#' @return data.frame with following columns:
#' \itemize{
#'   \item \code{param} name of the parameter to estimate
#'   \item \code{study} integer, representing the study from which the observation came
#'   \item \code{age} numeric, age of the observed individual
#'   \item \code{x} numeric, binary 0 or 1 for whether outcome was observed
#'   }
#' @export
#'
#' @examples
#' ## load param data
#' data("raw_age_estimates")
#' ## filter to only those for pipeline
#' raw_params <- filter(raw_age_estimates, USE_pipeline==T)
#' ## expand data for "p_symp_inf"
#' expand_age_data(raw_params, "p_symp_inf")
expand_age_data <- function(age_cat_data,
                            param_to_est="p_symp_inf"){
  ## filter to only parameter of interest
  param_dat <- filter(age_cat_data, param==param_to_est) %>%
    mutate(age_rng = paste0(ageL,"_",ageR))
  expanded_dat <- c()
  ## break out each study
  studies <- unique(param_dat$source)
  num_studies <- length(studies)
  for(i in 1:num_studies){
    tmp <- param_dat %>% filter(source==studies[i])
    ## break out each age group
    age_groups <- unique(tmp$age_rng)
    num_ag <- length(age_groups)
    for(j in 1:num_ag){
      tmp_ag <- tmp %>% filter(age_rng==age_groups[j])
      expanded_dat <- bind_rows(expanded_dat,
                                tibble(param=param_to_est,
                                       study=i,
                                       age=sample(tmp_ag$ageL:tmp_ag$ageR,
                                                  tmp_ag$N, replace=T),
                                       x=sample(c(rep(1,tmp_ag$X),
                                                  rep(0,tmp_ag$N-tmp_ag$X)),
                                                tmp_ag$N, replace=F)))
    }
  }
  return(expanded_dat)
}
