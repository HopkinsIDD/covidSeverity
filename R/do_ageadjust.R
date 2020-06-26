

#' convert age matrix to long data
#'
#' Fits a cubic spline across multiple studies to estimate the parameter of interest by specified age categories.
#'
#' @param age_matrix wide format matrix of age data, with row.name for geoid, and only columns for age groups. Column names need to be in the format of 10-19, 10_19, 10:19, or 10to19.
#' @return long data.frame of age data 
#' 
#' @import dplyr
#' @importFrom stringr str_pad
#' 
#' @export
#'
convert_age_matrix_long <- function(age_matrix){
    
    age_separator <- ifelse(any(grepl("_",colnames(age_matrix))), "_",
                            ifelse(any(grepl("-",colnames(age_matrix))), "-",
                                   ifelse(any(grep(":",colnames(age_matrix))), ":",
                                          ifelse(any(grep("to",colnames(age_matrix))), "to", NA))))
    
    if (is.na(age_separator)) return("ERROR: Age categories need to be separated by '-','_',':', or 'to'.")
    
    # pop by age group
    age_data_long <- age_matrix %>% tibble::as_tibble() %>% 
        dplyr::mutate(geoid = row.names(age_matrix)) %>%
        tidyr::pivot_longer(cols=c(-geoid), names_to = "age_range", values_to = "pop")
    age_data_long <- age_data_long %>% 
        tidyr::separate(age_range, c("age_l", "age_r"), sep=age_separator, remove=FALSE) %>%
        dplyr::group_by(geoid) %>%
        dplyr::mutate(prop = round(pop / sum(pop),4)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(age_cat = paste0("age_", 
                                       stringr::str_pad(age_l, width = 2, pad = "0"), 
                                       stringr::str_pad(age_r, width = 2, pad = "0")),
                      age_l = as.integer(age_l),
                      age_r = as.integer(age_r)) %>%
        dplyr::mutate(age_mid = round((age_l + age_r)/2, 1))
    
    return(age_data_long)
}
    






#' Take age data for non-uniform age groups and make it into 10-year age groups
#'
#' @param age_data long format age data with columns: [geoid, age_l, age_r, pop]
#' @param age_range_pred ages to predict (default: 0:89)
#' @param pop_zeros population to give age groups with 0 individuals. breaks because of log if left at 0.
#' 
#' @return 10-year age data
#' 
#' @import dplyr
#' @importFrom tibble as_tibble
#' 
#' @export
#'
make_10yr_age_data <- function(age_data, age_range_pred=0:89, pop_zeros=0.0001){    
    
    # population per year and mid age of group
    age_data <- age_data %>% 
        dplyr::mutate(pop_per_year = pop / (age_r - age_l),
                      age_mid = round((age_r + age_l) / 2, 1))
    
    # Smooth to 1-year age groups, by geoid then predict to 10-yr age groups
    fit_1yr_ages <- function(age_mid, pop_per_year, age_range_pred=age_range_pred){
        age_smth <- smooth.spline(age_mid, log(pop_per_year), df=length(age_mid)-1)
        tibble::as_tibble(predict(age_smth, age_range_pred)) %>%
            `colnames<-`(c("age","pop")) %>% 
            dplyr::mutate(pop = round(exp(pop)))
    }

    suppressWarnings( age_data_smoothed <- age_data %>% 
        dplyr::select(geoid, age_mid, pop_per_year) %>%
        tidyr::nest(-geoid) %>%
        dplyr::mutate(smth = purrr::map(data, ~fit_1yr_ages(.$age_mid, .$pop_per_year, age_range_pred))) %>%
        dplyr::select(-data) %>% tidyr::unnest(smth) %>% 
        tibble::as_tibble() %>% 
        dplyr::group_by(geoid) %>%
        dplyr::mutate(prop = round(pop / sum(pop),4)) %>%
        dplyr::ungroup() %>% 
        dplyr::mutate(age10 = floor(age/10)*10) )
    
    # Get age by 10yr groups
    age_data10 <- age_data_smoothed %>% 
        dplyr::group_by(geoid, age10) %>% 
        dplyr::summarise(prop = sum(prop),
                         pop = sum(pop)) %>%
        dplyr::mutate(age = paste0(age10, "_", age10+9)) %>%
        tibble::as_tibble()

    return(age_data10)
}





#' Setup age data for the parameter estimation
#'
#' @param age_data long format age data in 10 year age groups, with columns: [geoid, age, pop, prob]
#' 
#' @return list with 3 items
#' \itemize{
#'   \item \code{age_matrix} matrix, rows are age categories, columns are simulations, values are estimates of the conditional probabilities
#'   \item \code{age_matrix_pct} tbl, with the median, lower and upper bounds of a 95\% CI for each age category
#'   \item \code{age_data_long} character_string, same as in Arguments
#' }
#' 
#' @import dplyr
#' @importFrom tibble as_tibble
#' @importFrom tidyselect everything
#' @importFrom stringr str_pad
#' 
#' @export
#'
setup_age_data <- function(age_data){    
    
    
    # Make sure it matches the age groups of the probabilities
    age_data <- age_data %>%
        dplyr::group_by(geoid) %>%
        tidyr::separate(age, into=c("age_l", "age_r"), remove=FALSE) %>%
        dplyr::mutate(age_l = ifelse(as.integer(age_l)>=80, 80, age_l),
                      age_r = ifelse(as.integer(age_r)>=81, 100, age_r)) %>%
        dplyr::group_by(geoid, age_l, age_r) %>% 
        dplyr::summarise(pop = sum(pop),
                         prop = sum(prop)) %>%
        dplyr::mutate(age = paste(age_l, age_r, sep="_")) %>%
        tidyr::as_tibble() %>%
        dplyr::select(-age_l, -age_r)
    
    
    #convert long data to matrix
    age_matrix <- age_data %>% tidyr::pivot_wider(id_cols=geoid, 
                                    names_from = age, 
                                    values_from = pop)
    geoids_ <- age_matrix$geoid
    age_matrix <- age_matrix %>% 
        dplyr::select(-geoid) %>% as.matrix()
    row.names(age_matrix) <- geoids_
    
    #convert long data to matrix
    age_matrix_pct <- age_data %>% tidyr::pivot_wider(id_cols=geoid, 
                                                  names_from = age, 
                                                  values_from = prop)
    geoids_ <- age_matrix_pct$geoid
    age_matrix_pct <- age_matrix_pct %>% 
        dplyr::select(-geoid) %>% as.matrix()
    row.names(age_matrix_pct) <- geoids_
    
    return(list(age_data_long=age_data, age_matrix=age_matrix, age_matrix_pct=age_matrix_pct))
}
    






#' Do age-specific parameter estimation
#'
#' runs the full age estimation process. requires clean, setup age data setup using setup_age_data
#'
#' @param param_to_est parameter of interest
#' @param age_matrix wide format of age population
#' @param age_matrix_pct wide format of age population proportions
#' @param raw_params raw probabilities
#' @param n_sims numeric vector, the cutoff values for each age group
#' @param n_preds numeric, number of predictions for the gam model to make
#' @param age_grps age groups included in the estimation
#' @param study_wt study weight
#'
#' @return list
#' \itemize{
#'   \item \code{geoid_est}
#'   \item \code{sim_med} median simulation, to use for seeding
#' }
#' 
#' @import dplyr doParallel
#' @importFrom tibble as_tibble
#' @importFrom readr write_csv
#' @importFrom googlesheets4 read_sheet
#' @importFrom tidyselect starts_with
#' 
#' @export
#' 
do_severity_sampling <- function(
    param_to_est="p_symp_inf",
    age_matrix,
    age_matrix_pct,
    raw_params,
    n_sims=40,
    n_preds=1000,
    age_grps=c(seq(0,80,by=10),100),
    study_wt="none"
    ){
    
    #symptomatic given infection
    geoid_est <- foreach::foreach(h=1:n_sims, .combine=rbind, 
                                       .export = c("est_age_param", "est_geoid_rrs"),
                                       .packages = c("dplyr","tidyr")) %dopar% {
                                           ## make proportion symptomatic given infected
                                           tmp <- est_age_param(age_cat_data=raw_params,
                                                                param_to_est=param_to_est,
                                                                age_cats=age_grps,
                                                                n_preds=n_preds,
                                                                study_wt=study_wt)
                                           est_geoid_rrs(tmp$pred_mtx,
                                                         param_to_est=param_to_est,
                                                         geoid_age_mtx=age_matrix_pct,
                                                         geoid_pops=age_matrix) %>%
                                               dplyr::mutate(sim=h)
                                       }
    
    sim_med <- (geoid_est %>% 
                         dplyr::select(sim, prob_var = paste0(param_to_est, "_overall")) %>%
                         dplyr::distinct() %>%
                         dplyr::mutate(dist_from_med = abs( prob_var - median(prob_var))) %>%
                         dplyr::arrange(dist_from_med) %>%
                         dplyr::slice(1))$sim
    
    return(list(geoid_est=geoid_est, 
                sim_med=sim_med))
}

    


#' Do age-specific parameter estimation
#'
#' runs the full age estimation process. requires clean, setup age data setup using setup_age_data
#'
#' @param age_data long form age data with 10-year age groups and with columns: [geoid, age, pop, prop]
#' @param cores integer, number of cores to use for parallel processing
#' @param n_sims numeric vector, the cutoff values for each age group
#' @param n_preds numeric, number of predictions for the gam model to make
#' @param age_grps age groups included in the estimation
#' @param googlesheet_access whether the user has access to the raw data on google docs; default is FALSE
#' @param output_dir where the user wants the output data to be saved
#' @return list
#' 
#' @import dplyr doParallel foreach
#' @importFrom tibble as_tibble
#' @importFrom readr write_csv
#' @importFrom googlesheets4 read_sheet
#' @importFrom tidyselect starts_with
#' 
#' @export
#' 
get_ageadjustments <- function( age_data,
                                cores=4, 
                                n_sims=40,
                                n_preds=1000,
                                age_grps=c(seq(0,80,by=10),100),
                                googlesheet_access=FALSE,
                                output_dir="output",
                                pop_name = NULL) {
    
    # make sure there is an output directory
    dir.create(output_dir, recursive = TRUE, showWarnings=FALSE)
    print(paste0("Output saved in ",output_dir))
    
    # get cleaned age data
    age_data_clean <- setup_age_data(age_data = age_data)
    

    # Get data
    if (googlesheet_access){
        ## if have access to private files, use this
        raw_age_estimates <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1pDfQ2SkO--F2WNfjZxx6t9V7K51IZW0bJKV5cwSCZfA/edit#gid=1769840547",
                                                       sheet="age risk") #%>% dplyr::filter(publicly_available==T)
    } else {
        data("raw_age_estimates") 
        raw_age_estimates <- raw_age_estimates #%>% dplyr::filter(publicly_available==T)
    }
    
    # Get raw parameters
    raw_params <- raw_age_estimates %>%
        dplyr::filter(USE_pipeline==T) %>%
        dplyr::mutate(p = ifelse(!is.na(p), p,
                          ifelse(!is.na(X) & !is.na(N), X/N,
                                 ifelse(!is.na(pR) & !is.na(pL),
                                        (pL+pR)/2, NA))),
               N = ifelse(!is.na(N), N,
                          ifelse(!is.na(X) & p>0, round(X/p),
                                 ifelse(!is.na(pR),
                                        round(p*(1-p)*(qnorm(.975)/(pR-p))^2),
                                        1000))),
               X = ifelse(is.na(X),ceiling(N*p),X))
    
    
    # Set up parallel
    doParallel::registerDoParallel(cores)
    
    
    #symptomatic given infection
    p_symp_inf <- do_severity_sampling(param_to_est="p_symp_inf", 
                                       raw_params = raw_params,
                                       n_sims=n_sims,
                                       n_preds=n_preds, 
                                       age_matrix = age_data_clean$age_matrix,
                                       age_matrix_pct = age_data_clean$age_matrix_pct,
                                       age_grps=age_grps)
    
    # Death given symptomatic
    p_death_symp <- do_severity_sampling(param_to_est="p_death_symp", 
                                       raw_params = raw_params,
                                       n_sims=n_sims,
                                       n_preds=n_preds, 
                                       age_matrix = age_data_clean$age_matrix,
                                       age_matrix_pct = age_data_clean$age_matrix_pct,
                                       age_grps=age_grps)
    
    ## hospitalization rate amongst symptomatic
    p_hosp_symp <- do_severity_sampling(param_to_est="p_hosp_symp", 
                                       raw_params = raw_params,
                                       n_sims=n_sims,
                                       n_preds=n_preds, 
                                       age_matrix = age_data_clean$age_matrix,
                                       age_matrix_pct = age_data_clean$age_matrix_pct,
                                       age_grps=age_grps)
    
    ## ICU rate amongst hospitalized
    p_icu_hosp <- do_severity_sampling(param_to_est="p_icu_hosp", 
                                        raw_params = raw_params,
                                        n_sims=n_sims,
                                        n_preds=n_preds, 
                                        age_matrix = age_data_clean$age_matrix,
                                        age_matrix_pct = age_data_clean$age_matrix_pct,
                                        age_grps=age_grps)
    
    ## ventilation rate amongst symptomatic
    p_vent_icu <- do_severity_sampling(param_to_est="p_vent_icu", 
                                       raw_params = raw_params,
                                       n_sims=n_sims,
                                       n_preds=n_preds, 
                                       age_matrix = age_data_clean$age_matrix,
                                       age_matrix_pct = age_data_clean$age_matrix_pct,
                                       age_grps=age_grps)
    
    doParallel::stopImplicitCluster()
    
    ## Get Medians...................................................................
    ## get all the median sims
    
    ## make proportion symptomatic given infected
    set.seed(p_symp_inf$sim_med)
    p_symp_inf <- est_age_param(age_cat_data=raw_params,
                                param_to_est="p_symp_inf",
                                age_cats=age_grps,
                                n_preds=n_preds,
                                study_wt="none")
    
    ## make proportion death given symptomatic
    set.seed(p_death_symp$sim_med)
    p_death_symp <- est_age_param(age_cat_data=raw_params,
                                  param_to_est="p_death_symp",
                                  age_cats=age_grps,
                                  n_preds=n_preds,
                                  study_wt="none")
    
    ## make proportion hospitalized given symptomatic
    set.seed(p_hosp_symp$sim_med)
    p_hosp_symp <- est_age_param(age_cat_data=raw_params,
                                 param_to_est="p_hosp_symp",
                                 age_cats=age_grps,
                                 n_preds=n_preds,
                                 study_wt="none")
    
    ## make proportion symptomatic given infected
    set.seed(p_icu_hosp$sim_med)
    p_icu_hosp <- est_age_param(age_cat_data=raw_params,
                                param_to_est="p_icu_hosp",
                                age_cats=age_grps,
                                n_preds=n_preds,
                                study_wt="none")
    
    ## make proportion symptomatic given infected
    set.seed(p_vent_icu$sim_med)
    p_vent_icu <- est_age_param(age_cat_data=raw_params,
                                param_to_est="p_vent_icu",
                                age_cats=age_grps,
                                n_preds=n_preds,
                                study_wt="none")
    
    
    # Summarization function
    pred_sum_fn <- function(data, age_grps){
        tibble::as_tibble(data) %>%
            dplyr::mutate(age_grp=cut(age_grps[1:(length(age_grps)-1)], age_grps, right=F)) %>%
            tidyr::pivot_longer(cols=tidyselect::starts_with("V"),
                                names_to="pred",
                                values_to="est") %>%
            dplyr::group_by(age_grp) %>%
            dplyr::summarize(est_mean=mean(est),
                             est_med=median(est),
                             est_lb=quantile(est, probs=.025),
                             est_ub=quantile(est, probs=.975))
    }
    
    
    
    ## hospitalization rate amongst infections
    p_hosp_inf <- list(param_to_est="p_hosp_inf", pred_mtx = p_hosp_symp$pred_mtx * p_symp_inf$pred_mtx)
    p_hosp_inf$pred_sum <- p_hosp_inf$pred_mtx %>% pred_sum_fn(age_grps = age_grps)
    
    ## IFR
    p_death_inf <- list(param_to_est="p_death_inf", pred_mtx = p_death_symp$pred_mtx * p_symp_inf$pred_mtx)
    p_death_inf$pred_sum <- p_death_inf$pred_mtx %>% pred_sum_fn(age_grps = age_grps)
    
    ## mild rate amongst symptomatic
    p_mild_symp <- list(param_to_est="p_mild_symp", pred_mtx = 1-p_hosp_symp$pred_mtx)
    p_mild_symp$pred_sum <- p_mild_symp$pred_mtx %>% pred_sum_fn(age_grps = age_grps)
    
    ## mild illneess rate amongst infected
    p_mild_inf <- list(param_to_est="p_mild_inf", pred_mtx = p_mild_symp$pred_mtx * p_symp_inf$pred_mtx)
    p_mild_inf$pred_sum <- p_mild_inf$pred_mtx %>% pred_sum_fn(age_grps = age_grps)
    
    ## ventilation rate amongst infected
    p_vent_inf <- list(param_to_est="p_vent_inf", 
                       pred_mtx = p_vent_icu$pred_mtx * p_icu_hosp$pred_mtx * p_hosp_symp$pred_mtx * p_symp_inf$pred_mtx)
    p_vent_inf$pred_sum <- p_vent_inf$pred_mtx %>% pred_sum_fn(age_grps = age_grps)
    
    ## icu rate amongst infected
    p_icu_inf <- list(param_to_est="p_icu_inf", 
                       pred_mtx = p_icu_hosp$pred_mtx * p_hosp_symp$pred_mtx * p_symp_inf$pred_mtx)
    p_icu_inf$pred_sum <- p_icu_inf$pred_mtx %>% pred_sum_fn(age_grps = age_grps)



    ## get all output parameters
    geoid_params <- est_geoid_params(p_symp_inf$pred_mtx,
                                     param_to_est=p_symp_inf$param_to_est,
                                     age_data_clean$age_matrix_pct) %>%
        dplyr::left_join(est_geoid_params(p_death_inf$pred_mtx,
                                   param_to_est="p_death_inf",
                                   age_data_clean$age_matrix_pct), by="geoid") %>%
        dplyr::left_join(est_geoid_params(p_hosp_inf$pred_mtx,
                                   param_to_est="p_hosp_inf",
                                   age_data_clean$age_matrix_pct), by="geoid") %>%
        dplyr::left_join(est_geoid_params(p_mild_inf$pred_mtx,
                                   param_to_est="p_mild_inf",
                                   age_data_clean$age_matrix_pct), by="geoid") %>%
        dplyr::left_join(est_geoid_params(p_icu_inf$pred_mtx,
                                   param_to_est="p_icu_inf",
                                   age_data_clean$age_matrix_pct), by="geoid") %>%
        dplyr::left_join(est_geoid_params(p_vent_inf$pred_mtx,
                                   param_to_est="p_vent_inf",
                                   age_data_clean$age_matrix_pct), by="geoid") %>%
        dplyr::mutate(p_death_symp=p_death_inf/p_symp_inf,
               p_hosp_symp=p_hosp_inf/p_symp_inf,
               p_mild_symp=1-p_hosp_symp,
               p_icu_hosp=p_icu_inf/p_hosp_inf,
               p_vent_icu=p_vent_inf/p_icu_inf) %>%
        dplyr::left_join(est_geoid_rrs(pred_mtx=p_death_symp$pred_mtx * p_symp_inf$pred_mtx,
                                param_to_est="death_inf",
                                geoid_age_mtx=age_data_clean$age_matrix_pct,
                                geoid_pops=age_data_clean$age_matrix), by="geoid") %>%
        dplyr::left_join(est_geoid_rrs(pred_mtx=p_hosp_symp$pred_mtx * p_symp_inf$pred_mtx,
                                param_to_est="hosp_inf",
                                geoid_age_mtx=age_data_clean$age_matrix_pct,
                                geoid_pops=age_data_clean$age_matrix), by="geoid") %>%
        dplyr::left_join(est_geoid_rrs(pred_mtx=p_symp_inf$pred_mtx,
                                param_to_est="symp_inf",
                                geoid_age_mtx=age_data_clean$age_matrix_pct,
                                geoid_pops=age_data_clean$age_matrix), by="geoid") %>%
        dplyr::left_join(est_geoid_rrs(pred_mtx=p_icu_hosp$pred_mtx * p_hosp_symp$pred_mtx *
                                p_symp_inf$pred_mtx,
                                param_to_est="icu_inf",
                                geoid_age_mtx=age_data_clean$age_matrix_pct,
                                geoid_pops=age_data_clean$age_matrix), by="geoid") %>%
        dplyr::left_join(est_geoid_rrs(pred_mtx=p_vent_icu$pred_mtx * p_icu_hosp$pred_mtx *
                                p_hosp_symp$pred_mtx * p_symp_inf$pred_mtx,
                                param_to_est="vent_inf",
                                geoid_age_mtx=age_data_clean$age_matrix_pct,
                                geoid_pops=age_data_clean$age_matrix), by="geoid")

    ## Summarize and Save...........................................................   
    
    ## save a csv file
    geoid_params %>%
        readr::write_csv(file.path(output_dir, 
                                   ifelse(!is.null(pop_name), paste0(pop_name, "_geoid-params.csv"),
                                          "geoid-params.csv")))
    ## save as list
    saveRDS(list(p_symp_inf = p_symp_inf,
                 p_mild_inf = p_mild_inf,
                 p_hosp_inf = p_hosp_inf,
                 p_icu_inf = p_icu_inf,
                 p_vent_inf = p_vent_inf,
                 p_death_inf = p_death_inf,
                 p_mild_symp = p_mild_symp,
                 p_hosp_symp = p_hosp_symp,
                 p_death_symp = p_death_symp,
                 p_icu_hosp = p_icu_hosp,
                 p_vent_icu = p_vent_icu),
            file=file.path(output_dir, 
                           ifelse(!is.null(pop_name), paste0(pop_name, "_param-age-dist.rds"),
                                  "param-age-dist.rds")))
            
    return(geoid_params)
}
        









##' Get population distribution from WPP and aggregate it to 10 year age groups
##' - this is set up to use population estimates from the World Populaiton Prospects data
#'
#' @param country country of interest
#' @param ISO3 Optional; ISO3 code for country of interest.
#'
#' @import dplyr
#' 
#' @return 10-year age group populations for country of interest, long form
#' 
#' 
#' @export
#'
get_country_age_pop <- function(countryname, ISO3=NULL){
    
    data("wpp2019", package="covidSeverity")
    
    if (is.null(ISO3)){
        dat <- wpp2019 %>% dplyr::filter(country %in% countryname)
    } else {
        dat <- wpp2019 %>% dplyr::filter(ISO %in% ISO3)
    }
    
    dat <- dat %>% 
        dplyr::mutate(age_l = as.integer(age_l),
                      age_r = as.integer(age_r))
    
    datover80 <- dat %>% filter(age_l>=80) %>% 
        dplyr::group_by(country, ISO) %>%
        dplyr::summarise(pop = sum(pop),
                        prop = sum(prop),
                        age_l = 80,
                        age_r = 100,
                        age = "80_100") %>% tibble::as_tibble()
    dat <- dat %>% 
        dplyr::filter(age_l<80) %>%
        full_join(datover80, by = c("country", "age", "age_l", "age_r", "pop", "ISO", "prop")) %>% 
        tibble::as_tibble()
        
    return(dat)
}


##' Get population distribution from WPP and aggregate it to 10 year age groups
##' - this is set up to use population estimates from the World Populaiton Prospects data
#'
#' @param country country of interest
#' @param ISO3 Optional; ISO3 code for country of interest.
#'
#' @import dplyr
#' 
#' @return 10-year age group populations for country of interest, long form
#' 
#' 
#' @export
#'
get_allcountries_age_pop <- function(){
    
    data("wpp2019", package="covidSeverity")
    dat <- wpp2019 %>% dplyr::filter(!is.na(ISO))

    dat <- dat %>% 
        dplyr::mutate(age_l = as.integer(age_l),
                      age_r = as.integer(age_r))
    
    datover80 <- dat %>% filter(age_l>=80) %>% 
        dplyr::group_by(country, ISO) %>%
        dplyr::summarise(pop = sum(pop),
                         prop = sum(prop),
                         age_l = 80,
                         age_r = 100,
                         age = "80_100") %>% tibble::as_tibble()
    dat <- dat %>% 
        dplyr::filter(age_l<80) %>%
        full_join(datover80, by = c("country", "age", "age_l", "age_r", "pop", "ISO", "prop")) %>% 
        tibble::as_tibble()
    
    return(dat)
}





    