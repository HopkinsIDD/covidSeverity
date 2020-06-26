




#' Multiply Pr(age|var) * Pr(age) for each geoid
#' instead of inputting p_hosp to normalize to 1, we just normalize to 1
#'
#' @param p_var 
#' @param p_age_long 
#' @param params_age 
#'
#' @return
#' @export
#' 
#' @import dplyr
#'
#' @examples
get_p_age_given_outcome <- function(p_var = "p_hosp_inf", p_age_long, params_age, age_groups){
    
    
    p_age_ <- p_age_long %>% 
        dplyr::select(-pop) %>% 
        dplyr::rename(p_age = prop) %>%
        dplyr::left_join(params_age[[p_var]] %>% 
                             dplyr::mutate(age = age_groups) %>% 
                             dplyr::select(-age_grp), by=c("age")) %>%
        dplyr::mutate(est_mean = est_mean * p_age,
                      est_med = est_med * p_age,
                      est_ub = est_ub * p_age,
                      est_lb = est_lb * p_age) %>%
        dplyr::group_by(geoid) %>%
        # Normalizing the outputs
        dplyr::mutate(est_mean = est_mean / sum(est_mean),
                      est_med = est_med / sum(est_med),
                      est_ub = est_ub / sum(est_ub),
                      est_lb = est_lb / sum(est_lb)) %>%
        dplyr::ungroup()
    
    return(p_age_)
}




# # Load age adjustment data used in the CSP model
# geoid_params <- readr::read_csv(file.path("output", "geoid-params.csv"))
# params_age_dist <- readRDS(file=file.path("output", "param-age-dist.rds"))
# 
# param_name <- params_age <- list()
# for (i in 1:length(params_age_dist)){
#     param_name[[i]] <- params_age_dist[[i]]$param_to_est
#     params_age[[i]] <- params_age_dist[[i]]$pred_sum
# }
# param_name <- unlist(param_name)
# names(params_age) <- param_name
# 



#' Title
#'
#' @param p_age_long 
#' @param params_age 
#'
#' @return
#' @export
#' 
#' @import dplyr
#'
#' @examples
get_prob_ages_all <- function(p_age_long, params_age){

    prob_data <- p_age_long %>% dplyr::select(geoid, age, p_inf = prop) %>%
        dplyr::left_join(
            get_p_age_given_outcome(p_var = "p_symp_inf", p_age_long, params_age) %>% 
                dplyr::select(geoid, age, p_symp_inf=est_med)) %>%
        dplyr::left_join(
            get_p_age_given_outcome(p_var = "p_mild_inf", p_age_long, params_age) %>% 
                dplyr::select(geoid, age, p_mild_inf=est_med)) %>%
        dplyr::left_join(
            get_p_age_given_outcome(p_var = "p_hosp_inf", p_age_long, params_age) %>% 
                dplyr::select(geoid, age, p_hosp_inf=est_med)) %>%
        dplyr::left_join(
            get_p_age_given_outcome(p_var = "p_icu_inf", p_age_long, params_age) %>% 
                dplyr::select(geoid, age, p_icu_inf=est_med)) %>%
        dplyr::left_join(
            get_p_age_given_outcome(p_var = "p_vent_inf", p_age_long, params_age) %>% 
                dplyr::select(geoid, age, p_vent_inf=est_med)) %>%
        dplyr::left_join(
            get_p_age_given_outcome(p_var = "p_death_inf", p_age_long, params_age) %>% 
                dplyr::select(geoid, age, p_death_inf=est_med)) 

    # map probabilities to outputs
    prob_data <- prob_data %>%
        dplyr::select(geoid, age,
                      incidSymp = p_symp_inf,
                      incidI = p_inf,
                      incidMild = p_mild_inf,
                      incidH = p_hosp_inf,
                      incidICU = p_icu_inf,
                      incidVent = p_vent_inf,
                      incidD = p_death_inf)
    
    return(prob_data)
    
}










#' Title
#'
#' @param ... 
#'
#' @return
#' @export
#' 
#' @import dplyr
#' @import tibble
#'
#' @examples
rmultinom_t <- function(...){
    rmultinom(...) %>% 
        t %>% 
        tibble::as_tibble(.name_repair = "minimal")
}



#' Title
#'
#' @param data 
#' @param outcome_var 
#' @param prob_data 
#'
#' @return
#' @export
#' 
#' @import dplyr
#' @importFrom purrr map
#'
#' @examples
prob_fun <- function(data,
                     outcome_var = "incidH",
                     prob_data, 
                     incl_outcome=TRUE){
    
    probs <- prob_data[[outcome_var]]
    names(probs) <- prob_data$age
    
    if (incl_outcome){
        purrr::map(data[[outcome_var]], rmultinom_t, n=1, prob = probs) %>%
            dplyr::bind_rows() %>% setNames(paste0(outcome_var, "_", prob_data$age))
    }else{
        purrr::map(data[[outcome_var]], rmultinom_t, n=1, prob = probs) %>%
            dplyr::bind_rows()
    }
}



#' #' Title
#' #'
#' #' @param data 
#' #' @param outcome_var 
#' #' @param prob_data 
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @import dplyr
#' #' @importFrom purrr map
#' #' 
#' #' @examples
#' prob_fun_strat <- function(data,
#'                            outcome_var = "incidH",
#'                            prob_data){
#'     
#'     purrr:::map(data[[outcome_var]], rmultinom_t, n=1, prob = prob_data[[outcome_var]]) %>%
#'         dplyr::bind_rows() %>% setNames(paste0(outcome_var, "_", prob_data$age))
#'     
#' }







#' Age Stratify an Outcome Variable
#'
#' @param outcome_var 
#' @param outcome_data 
#' @param prob_data_list 
#' @param add 
#'
#' @return
#' @export
#' 
#' @import dplyr
#' @importFrom purrr map2
#'
#' @examples
age_strat_outcome_wide <- function(outcome_var="incidI", outcome_data, prob_data_list, add=TRUE){
    
    output <- outcome_data %>% split(.$geoid) %>%
        purrr::map2(.x=., .y=prob_data_list, ~prob_fun(data=.x, prob_data = .y, outcome_var=outcome_var, incl_outcome = TRUE)) %>%
        dplyr::bind_rows()
    
    if (add){
        return(dplyr::bind_cols(outcome_data, output))
    } else {
        return(dplyr::bind_cols(outcome_data %>% dplyr::select(geoid, date), output))
    }
}


#' Age Stratify an Outcome Variable, long form
#'
#' @param outcome_var 
#' @param outcome_data 
#' @param prob_data_list 
#' @param add 
#'
#' @return
#' @export
#' 
#' @import dplyr
#' @importFrom purrr map2
#' @importFrom tidyr pivot_longer
#'
#' @examples
age_strat_outcome_long <- function(outcome_var="incidI", outcome_data, prob_data_list){
    
    output <- outcome_data %>% split(.$geoid) %>%
        purrr::map2(.x=., .y=prob_data_list, ~prob_fun(data=.x, prob_data = .y, outcome_var=outcome_var, incl_outcome = FALSE)) %>%
        dplyr::bind_rows() %>% 
        dplyr::bind_cols(outcome_data %>% dplyr::select(geoid, date)) %>% 
        tidyr::pivot_longer(cols = -c(geoid, date), names_to = "age", values_to = outcome_var)
    
    return(output)
}










#' Title
#'
#' @param outcomes 
#' @param outcome_data 
#' @param prob_data 
#'
#' @return
#' @export
#' 
#' @import dplyr
#'
#' @examples
get_age_strat_outcomes <- function(outcomes=c("incidI","incidMild","incidH","incidICU","incidVent","incidD"), 
                                   outcome_data, 
                                   prob_data,
                                   wide=FALSE){
    
    prob_data_list <- prob_data %>% split(.$geoid)
    
    if (wide){
        for (i in seq_along(outcomes)){
            outcome_data <- age_strat_outcome_wide(outcome_var=outcomes[i], outcome_data, prob_data_list, add=TRUE)
        }
        return(outcome_data)
        
    } else {
        res <- age_strat_outcome_long(outcome_var=outcomes[1], outcome_data, prob_data_list)
        for (i in 2:length(outcomes)){
            res1 <- age_strat_outcome_long(outcome_var=outcomes[i], outcome_data, prob_data_list)
            res <- dplyr::bind_cols(res, res1 %>% dplyr::select(outcomes[[i]]))
        }
        res <- res %>% dplyr::bind_rows(outcome_data %>% dplyr::mutate(age="all"))
        return(res)
    }
}






