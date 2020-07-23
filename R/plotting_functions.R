## Plotting functions



# library(sf)
# library(tidyverse)
# library(viridis)
# 
# # First generate a geodata.csv like used in the CSP
# geodata <- readr::read_csv("../raw_data/BGD/BGD_age_pop_10yr.csv")
# geodata <- geodata %>% 
#     dplyr::group_by(geoid, adm2) %>%
#     dplyr::summarize(pop = round(sum(pop)),
#                      admin0 = "BGD") %>%
#     dplyr::select(geoid, admin2=adm2, admin0, pop)
# readr::write_csv(geodata, "../raw_data/BGD/geodata.csv")
# 
# 


#' #' map_severity
#' #' 
#' #' Map the outcome probabilities by geoid.
#' #'
#' #' @param severity_var variable to plot 
#' #' @param severity_var_label label for legend title
#' #' @param geoid_params path to geoid_params.csv file
#' #' @param geodata path to geodata.csv
#' #' @param shp path to shapefile
#' #' @param shp_name_var name of variable for names of geo-locations
#' #' @param label_locs Locations to label (i.e., districts or counties)
#' #'
#' #' @return
#' #' @export
#' #' 
#' #' @import tidyverse readr sf
#' #' @importFrom viridis scale_fill_viridis
#' #'
#' #' @examples
#' map_severity <- function(severity_var = "p_death_inf", 
#'                          severity_var_label = "IFR",
#'                          geoid_params = "../raw_data/BGD/geoid-params.csv",
#'                          geodata = "../raw_data/BGD/geodata.csv",
#'                          shp = "../raw_data/BGD/bgd_admbnda_adm2_bbs_20180410/bgd_admbnda_adm2_bbs_20180410.shp",
#'                          shp_name_var = "ADM2_EN",
#'                          label_locs = NA   ){
#'     suppressMessages({
#'         geoid_params_ <- readr::read_csv(geoid_params)
#'         mapdat <- read_sf(shp) %>%
#'             dplyr::mutate(NAME = get(shp_name_var))
#'         geodat <- readr::read_csv(geodata)
#'     })
#'     
#'     if ("pop" %in% colnames(mapdat)){
#'         mapdat <- mapdat %>% dplyr::select(-pop)
#'     }
#'     if ("geoid" %in% colnames(mapdat)){
#'         mapdat <- mapdat %>% dplyr::rename(geoid_v0 = geoid)
#'     }
#'     
#'     # Combine data into shapefile
#'     mapdat <- mapdat %>% 
#'         dplyr::left_join(geodat, by=c("NAME"="admin2")) %>%
#'         dplyr::left_join(geoid_params_, by=c("geoid"))
#'     
#'     
#'     # Plot it
#'     plt <- ggplot(mapdat) +
#'         geom_sf(aes(fill = get(severity_var))) +
#'         theme_void() +
#'         viridis::scale_fill_viridis() +
#'         labs(fill=severity_var_label)
#'     if (!is.na(label_locs)){
#'         if (tolower(label_locs)=="all"){
#'             plt <- plt +
#'                 geom_text(mapdat, aes(X,Y, label=NAME), size=5)
#'         } else if (!is.na(label_locs)){
#'             label_dat <-  mapdat %>% filter(NAME %in% label_locs)
#'             plt <- plt + geom_sf_text(data=label_dat, aes(label=NAME), size=3, hjust=0.5)
#'         }
#'     }
#'     
#'     return(plt)
#' }


# Example
# map_severity(severity_var = "p_death_inf", 
#              severity_var_label = "IFR",
#              geoid_params = "../raw_data/BGD/geoid-params.csv",
#              geodata = "../raw_data/BGD/geodata.csv",
#              shp = "../raw_data/BGD/bgd_admbnda_adm2_bbs_20180410/bgd_admbnda_adm2_bbs_20180410.shp",
#              shp_name_var = "ADM2_EN",
#              label_locs = NA)








library(sf)
library(tidyverse)
library(viridis)



#' Title
#'
#' @param severity_var 
#' @param severity_var_label 
#' @param geoid_params 
#' @param geodata 
#' @param shp 
#' @param shp_name_var 
#' @param label_locs 
#'
#' @return
#' @export
#'
#' @examples
map_severity <- function(severity_var = "p_death_inf", 
                         severity_var_label = "IFR",
                         geoid_params = "../raw_data/BGD/geoid-params.csv",
                         geodata = "../raw_data/BGD/geodata.csv",
                         shp = "../raw_data/BGD/bgd_admbnda_adm2_bbs_20180410/bgd_admbnda_adm2_bbs_20180410.shp",
                         shp_name_var = "ADM2_EN",
                         label_locs = NA   ){
    suppressMessages({
        geoid_params_ <- readr::read_csv(geoid_params)
        mapdat <- read_sf(shp) %>%
            dplyr::mutate(NAME = get(shp_name_var))
        geodat <- readr::read_csv(geodata)
    })
    
    if ("pop" %in% colnames(mapdat)){
        mapdat <- mapdat %>% dplyr::select(-pop)
    }
    if ("geoid" %in% colnames(mapdat)){
        mapdat <- mapdat %>% dplyr::rename(geoid_v0 = geoid)
    }
    
    # Combine data into shapefile
    mapdat <- mapdat %>% 
        dplyr::left_join(geodat, by=c("NAME"="admin2")) %>%
        dplyr::left_join(geoid_params_, by=c("geoid"))
    
    
    # Plot it
    plt <- ggplot(mapdat) +
        geom_sf(aes(fill = get(severity_var))) +
        theme_void() +
        viridis::scale_fill_viridis() +
        labs(fill=severity_var_label)
    if (!is.na(label_locs)){
        if (tolower(label_locs)=="all"){
            plt <- plt +
                geom_text(mapdat, aes(X,Y, label=NAME), size=5)
        } else if (!is.na(label_locs)){
            label_dat <-  mapdat %>% filter(NAME %in% label_locs)
            plt <- plt + geom_sf_text(data=label_dat, aes(label=NAME), size=3, hjust=0.5)
        }
    }
    
    return(plt)
}






#' Title
#'
#' @param severity_var 
#' @param severity_var_label 
#' @param .region 
#' @param .country 
#' @param geoid_params 
#' @param geodata 
#' @param shp_adm2 
#' @param shp_adm0 
#' @param shp2_name_var 
#' @param legend_pos 
#' @param label_locs 
#'
#' @return
#' @export
#'
#' @examples
map_severity_gadm <- function(severity_var = "p_death_inf", 
                              severity_var_label = "IFR",
                              .region = "Africa",
                              .country = NULL,
                              region_var = "region",
                              geoid_params = "../raw_data/global_severity/GLOBAL_geoid_params.csv",
                              geodata = "../raw_data/global_severity/GLOBAL_geodata.csv",
                              shp_adm2 = "raw_data/gadm36_levels_shp/gadm36_2.shp",
                              shp_adm0 = "raw_data/gadm36_levels_shp/gadm36_0.shp",
                              shp2_name_var = "GID_2",
                              legend_pos = c(.9,.9),
                              leg_breaks = c(0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.007, 0.01, 0.02, 0.03),
                              log_trans = TRUE,
                              leg_range = c(0,1),
                              color_pal = "plasma",
                              label_locs = NA ){
    
    suppressMessages({
        if (is.character(geoid_params)){
            geoid_params <- readr::read_csv(geoid_params)
        }
        if (is.character(geodata)){
            geodata <- readr::read_csv(geodata)
        }
        if (is.character(shp_adm2)){
            shp_adm2 <- sf::read_sf(shp_adm2)
        }
        if (is.character(shp_adm0)){
            shp_adm0 <- sf::read_sf(shp_adm0)
        }
    })
    
    shp_adm2 <- shp_adm2 %>% dplyr::mutate(geoid = get(shp2_name_var))
    
    
    
    if ("pop" %in% colnames(shp_adm2)){
        shp_adm2 <- shp_adm2 %>% dplyr::select(-pop)
    }
    
    # Add regions
    shp_adm0 <- shp_adm0 %>% 
        dplyr::mutate(region = globaltoolboxlite::get_region(GID_0),
                      who_region = globaltoolboxlite::get_whoregion(GID_0),
                      subregion = globaltoolboxlite::get_subregion(GID_0))
    
    # Combine data into shapefile
    shp_adm2 <- shp_adm2 %>% 
        dplyr::left_join(geodata, by=c("geoid")) %>%
        dplyr::left_join(geoid_params, by=c("geoid"))
    
    # Filter to desired locations
    if (!is.null(.country)){
        shp_adm2 <- shp_adm2 %>% dplyr::filter(GID_0 %in% .country)
        shp_adm0 <- shp_adm0 %>% dplyr::filter(GID_0 %in% .country)
    }
    if (!is.null(.region)){
        shp_adm2 <- shp_adm2 %>% 
            dplyr::mutate(reg = get(region_var)) %>%
            dplyr::filter(reg %in% .region)
        shp_adm0 <- shp_adm0 %>% 
            dplyr::mutate(reg = get(region_var)) %>%
            dplyr::filter(reg %in% .region)
    }
    
    
    
    # Plot it
    plt <- ggplot() +
        geom_sf(data=shp_adm0, fill="grey", size = .75, color = "black") +
        geom_sf(data=shp_adm2, aes(fill = get(severity_var)), size = .25, color = "black") +
        geom_sf(data=shp_adm0, fill=NA, size = .75, color = "black") +
        theme_void() +
        labs(fill=severity_var_label) + 
        theme(legend.position = legend_pos) 
    
    plt
    
    if (log_trans){
        plt <- plt + scale_fill_viridis_c(breaks = leg_breaks, 
                                          labels = leg_breaks, 
                                          trans = scales::pseudo_log_trans(), 
                                          na.value="grey", 
                                          option = color_pal)
    } else {
        plt <- plt + viridis::scale_fill_viridis(na.value="grey",   
                                                 limits = leg_range,
                                                 option = color_pal)
        #begin=leg_range[1], end=leg_range[2])
    }
    
    plt
    
    # Add labels if needed
    
    if (!is.na(label_locs)){
        if (tolower(label_locs)=="all"){
            plt <- plt +
                geom_text(shp_adm2, aes(X,Y, label=NAME), size=5)
        } else if (!is.na(label_locs)){
            label_dat <-  shp_adm2 %>% filter(NAME %in% label_locs)
            plt <- plt + geom_sf_text(data=label_dat, aes(label=NAME), size=3, hjust=0.5)
        }
    }
    
    return(plt)
}










#' Title
#'
#' @param severity_var 
#' @param severity_var_label 
#' @param .region 
#' @param .country 
#' @param geoid_params 
#' @param geodata 
#' @param shp_adm2 
#' @param shp_adm0 
#' @param shp2_name_var 
#' @param legend_pos 
#' @param label_locs 
#'
#' @return
#' @export
#'
#' @examples
merge_severity_gadm_data <- function(
                              geoid_params = "raw_data/global_severity/GLOBAL_geoid_params.csv",
                              geodata = "raw_data/global_severity/GLOBAL_geodata.csv",
                              shp_adm2 = "raw_data/gadm36_levels_shp/gadm36_2.shp",
                              shp_adm0 = "raw_data/gadm36_levels_shp/gadm36_0.shp",
                              shp2_name_var = "GID_2"){
    
    
    suppressMessages({
        if (is.character(geoid_params)){
            geoid_params <- readr::read_csv(geoid_params)
        }
        if (is.character(geodata)){
            geodata <- readr::read_csv(geodata)
        }
        if (is.character(shp_adm2)){
            shp_adm2 <- sf::read_sf(shp_adm2)
        }
        if (is.character(shp_adm0)){
            shp_adm0 <- sf::read_sf(shp_adm0)
        }
    })
    
    shp_adm2 <- shp_adm2 %>% dplyr::mutate(geoid = get(shp2_name_var)) 
    
    
    if ("pop" %in% colnames(shp_adm2)){
        shp_adm2 <- shp_adm2 %>% dplyr::select(-pop)
    }
    
    # Add regions
    shp_adm0 <- shp_adm0 %>% 
        dplyr::mutate(region = globaltoolboxlite::get_region(GID_0),
                      who_region = globaltoolboxlite::get_whoregion(GID_0),
                      subregion = globaltoolboxlite::get_subregion(GID_0))
    
    # Combine data into shapefile
    shp_adm2 <- shp_adm2 %>% 
        dplyr::left_join(geodata, by=c("geoid")) %>%
        dplyr::left_join(geoid_params, by=c("geoid"))

 return(list(shp_adm2=shp_adm2, shp_adm0=shp_adm0))   
    
}








# Fast version with all data already merged so we can experiment (merging happens in function above)


#' Title
#'
#' @param severity_var 
#' @param severity_var_label 
#' @param .region 
#' @param .country 
#' @param geoid_params 
#' @param geodata 
#' @param shp_adm2 
#' @param shp_adm0 
#' @param shp2_name_var 
#' @param legend_pos 
#' @param label_locs 
#'
#' @return
#' @export
#'
#' @examples
map_severity_gadm_fast <- function(severity_var = "p_death_inf", 
                              severity_var_label = "IFR",
                              .region = "Africa",
                              .country = NULL,
                              region_var = "region",
                              merged_data = mapdata,
                              shp2_name_var = "GID_2",
                              legend_pos = c(.9,.9),
                              leg_breaks = c(0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.007, 0.01, 0.02, 0.03),
                              log_trans = TRUE,
                              leg_range = c(0,1),
                              color_pal = "plasma",
                              border_size = c(.25, .05),
                              label_locs = NA ){
    
    
    shp_adm0 <- merged_data$shp_adm0
    shp_adm2 <- merged_data$shp_adm2
        
    # Filter to desired locations
    if (!is.null(.country)){
        shp_adm2 <- shp_adm2 %>% dplyr::filter(GID_0 %in% .country)
        shp_adm0 <- shp_adm0 %>% dplyr::filter(GID_0 %in% .country)
    }
    
    if (!is.null(.region)){
        shp_adm2 <- shp_adm2 %>% 
            dplyr::mutate(reg = get(region_var)) %>%
            dplyr::filter(reg %in% .region)
        shp_adm0 <- shp_adm0 %>% 
            dplyr::mutate(reg = get(region_var)) %>%
            dplyr::filter(reg %in% .region)
    }
    
    
    # Plot it
    plt <- ggplot() +
        geom_sf(data=shp_adm0, fill="grey", size = border_size[1], color = "black") +
        geom_sf(data=shp_adm2, aes(fill = get(severity_var)), size = border_size[2], color = "black") +
        geom_sf(data=shp_adm0, fill=NA, size = border_size[1], color = "black") +
        theme_void() +
        labs(fill=severity_var_label) + 
        theme(legend.position = legend_pos) 
    
    plt
    
    if (log_trans){
        plt <- plt + scale_fill_viridis_c(breaks = leg_breaks, 
                                          labels = leg_breaks, 
                                          trans = scales::pseudo_log_trans(), 
                                          na.value="grey", 
                                          option = color_pal)
    } else {
        plt <- plt + viridis::scale_fill_viridis(na.value="grey",   
                                                 limits = leg_range,
                                                 option = color_pal)
        #begin=leg_range[1], end=leg_range[2])
    }
    
    plt
    
    # Add labels if needed
    
    if (!is.na(label_locs)){
        if (tolower(label_locs)=="all"){
            plt <- plt +
                geom_text(shp_adm2, aes(X,Y, label=NAME), size=5)
        } else if (!is.na(label_locs)){
            label_dat <-  shp_adm2 %>% filter(NAME %in% label_locs)
            plt <- plt + geom_sf_text(data=label_dat, aes(label=NAME), size=3, hjust=0.5)
        }
    }
    
    return(plt)
}











#' Title
#'
#' @param severity_var 
#' @param severity_var_label 
#' @param .region 
#' @param .country 
#' @param geoid_params 
#' @param geodata 
#' @param shp_adm2 
#' @param shp_adm0 
#' @param shp2_name_var 
#' @param legend_pos 
#' @param label_locs 
#'
#' @return
#' @export
#'
#' @examples
map_ageprop_gadm <- function(severity_var = "p_death_inf", 
                              severity_var_label = "IFR",
                              .region = "Africa",
                              .country = NULL,
                              region_var = "region",
                              age_pop_data = "raw_data/global_severity/GLOBAL_age_pop_10yr.csv",
                              geodata = "raw_data/global_severity/GLOBAL_geodata.csv",
                              shp_adm2 = "raw_data/gadm36_levels_shp/gadm36_2.shp",
                              shp_adm0 = "raw_data/gadm36_levels_shp/gadm36_0.shp",
                              shp2_name_var = "GID_2",
                              legend_pos = c(.9,.9),
                              leg_breaks = c(0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.007, 0.01, 0.02, 0.03),
                              log_trans = TRUE,
                              leg_range = c(0,1),
                              color_pal = "plasma",
                              label_locs = NA ){
    
    # Load the data
    suppressMessages({
        if (is.character(age_pop_data)){
            geodata <- readr::read_csv(age_pop_data)
        }
        if (is.character(geodata)){
            geodata <- readr::read_csv(geodata)
        }
        mapdat <- sf::read_sf(shp_adm2) %>%
            dplyr::mutate(geoid = get(shp2_name_var))
        mapdat0 <- sf::read_sf(shp_adm0)
    })
    
    
    if ("pop" %in% colnames(mapdat)){
        mapdat <- mapdat %>% dplyr::select(-pop)
    }
    
    # Add regions
    mapdat0 <- mapdat0 %>% 
        dplyr::mutate(region = globaltoolboxlite::get_region(GID_0),
                      who_region = globaltoolboxlite::get_whoregion(GID_0),
                      subregion = globaltoolboxlite::get_subregion(GID_0))
    
    # Combine data into shapefile
    mapdat <- mapdat %>% 
        dplyr::left_join(geodata, by=c("geoid")) %>%
        dplyr::left_join(geoid_params, by=c("geoid"))
    
    # Filter to desired locations
    if (!is.null(.country)){
        mapdat <- mapdat %>% dplyr::filter(GID_0 %in% .country)
        mapdat0 <- mapdat0 %>% dplyr::filter(GID_0 %in% .country)
    }
    if (!is.null(.region)){
        mapdat <- mapdat %>% 
            dplyr::mutate(reg = get(region_var)) %>%
            dplyr::filter(reg %in% .region)
        mapdat0 <- mapdat0 %>% 
            dplyr::mutate(reg = get(region_var)) %>%
            dplyr::filter(reg %in% .region)
    }
    
    
    
    # Plot it
    plt <- ggplot() +
        geom_sf(data=mapdat0, fill="grey", size = 1.1, color = "black") +
        geom_sf(data=mapdat, aes(fill = get(severity_var)), size = .25, color = "black") +
        geom_sf(data=mapdat0, fill=NA, size = 1.1, color = "black") +
        theme_void() +
        labs(fill=severity_var_label) + 
        theme(legend.position = legend_pos) 
    
    plt
    
    if (log_trans){
        plt <- plt + scale_fill_viridis_c(breaks = leg_breaks, 
                                          labels = leg_breaks, 
                                          trans = scales::pseudo_log_trans(), 
                                          na.value="grey", 
                                          option = color_pal)
    } else {
        plt <- plt + viridis::scale_fill_viridis(na.value="grey",   
                                                 limits = leg_range,
                                                 option = color_pal)
        #begin=leg_range[1], end=leg_range[2])
    }
    
    plt
    
    # Add labels if needed
    
    if (!is.na(label_locs)){
        if (tolower(label_locs)=="all"){
            plt <- plt +
                geom_text(mapdat, aes(X,Y, label=NAME), size=5)
        } else if (!is.na(label_locs)){
            label_dat <-  mapdat %>% filter(NAME %in% label_locs)
            plt <- plt + geom_sf_text(data=label_dat, aes(label=NAME), size=3, hjust=0.5)
        }
    }
    
    return(plt)
}

