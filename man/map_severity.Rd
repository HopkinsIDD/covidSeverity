% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{map_severity}
\alias{map_severity}
\title{map_severity}
\usage{
map_severity(
  severity_var = "p_death_inf",
  severity_var_label = "IFR",
  geoid_params = "../raw_data/BGD/geoid-params.csv",
  geodata = "../raw_data/BGD/geodata.csv",
 
    shp = "../raw_data/BGD/bgd_admbnda_adm2_bbs_20180410/bgd_admbnda_adm2_bbs_20180410.shp",
  shp_name_var = "ADM2_EN",
  label_locs = NA
)
}
\arguments{
\item{severity_var}{variable to plot}

\item{severity_var_label}{label for legend title}

\item{geoid_params}{path to geoid_params.csv file}

\item{geodata}{path to geodata.csv}

\item{shp}{path to shapefile}

\item{shp_name_var}{name of variable for names of geo-locations}

\item{label_locs}{Locations to label (i.e., districts or counties)}
}
\value{

}
\description{
Map the outcome probabilities by geoid.
}
