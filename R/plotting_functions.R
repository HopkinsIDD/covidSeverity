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


#' map_severity
#' 
#' Map the outcome probabilities by geoid.
#'
#' @param severity_var variable to plot 
#' @param severity_var_label label for legend title
#' @param geoid_params path to geoid_params.csv file
#' @param geodata path to geodata.csv
#' @param shp path to shapefile
#' @param shp_name_var name of variable for names of geo-locations
#' @param label_locs Locations to label (i.e., districts or counties)
#'
#' @return
#' @export
#' 
#' @import tidyverse readr sf
#' @importFrom viridis scale_fill_viridis
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


# Example
# map_severity(severity_var = "p_death_inf", 
#              severity_var_label = "IFR",
#              geoid_params = "../raw_data/BGD/geoid-params.csv",
#              geodata = "../raw_data/BGD/geodata.csv",
#              shp = "../raw_data/BGD/bgd_admbnda_adm2_bbs_20180410/bgd_admbnda_adm2_bbs_20180410.shp",
#              shp_name_var = "ADM2_EN",
#              label_locs = NA)