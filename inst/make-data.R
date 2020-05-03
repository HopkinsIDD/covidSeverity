## make data
library(covidSeverity)


## Adding US territories
## This is currently a rather bootleg operation

## mdb file has the column names for the txt file, though nothing really lines up
# AS_idx <- Hmisc::mdb.get(file="raw_data/DP_AS_Census2010_Access2003_v7.mdb")
## American Samoa
AS <- read_csv("raw_data/as0000120101.txt", col_names=F) %>%
  slice(1) %>%
  select(a_0_4=X7, a_5_9=X8, a_10_14=X9, a_15_19=X10, a_20_24=X11, a_25_29=X12,
         a_30_34=X13, a_35_39=X14, a_40_44=X15, a_45_49=X16, a_50_54=X17,
         a_55_59=X18, a_60_64=X19, a_65_69=X20, a_70_74=X21, a_75_79=X22,
         a_80_84=X23, a_85_99=X24) %>%
  pivot_longer(cols=everything(),
               names_to="age_cat",
               values_to="estimate") %>%
  separate(age_cat, c("a", "age_l", "age_r"), "_") %>%
  mutate(GEOID="60000",
         NAME="American Samoa",
         age_l=as.numeric(age_l),
         age_r=as.numeric(age_r)) %>%
  format_age_county_pop()

## Guam
GU <- read_csv("raw_data/gu0000120101.txt", col_names=F) %>%
  slice(1) %>%
  select(a_0_4=X7, a_5_9=X8, a_10_14=X9, a_15_19=X10, a_20_24=X11, a_25_29=X12,
         a_30_34=X13, a_35_39=X14, a_40_44=X15, a_45_49=X16, a_50_54=X17,
         a_55_59=X18, a_60_64=X19, a_65_69=X20, a_70_74=X21, a_75_79=X22,
         a_80_84=X23, a_85_99=X24) %>%
  pivot_longer(cols=everything(),
               names_to="age_cat",
               values_to="estimate") %>%
  separate(age_cat, c("a", "age_l", "age_r"), "_") %>%
  mutate(GEOID="66000",
         NAME="Guam",
         age_l=as.numeric(age_l),
         age_r=as.numeric(age_r)) %>%
  format_age_county_pop()

## Virgin Islands
VI <- read_csv("raw_data/vi0000120101.txt", col_names=F) %>%
  slice(1) %>%
  select(a_0_4=X7, a_5_9=X8, a_10_14=X9, a_15_19=X10, a_20_24=X11, a_25_29=X12,
         a_30_34=X13, a_35_39=X14, a_40_44=X15, a_45_49=X16, a_50_54=X17,
         a_55_59=X18, a_60_64=X19, a_65_69=X20, a_70_74=X21, a_75_79=X22,
         a_80_84=X23, a_85_99=X24) %>%
  pivot_longer(cols=everything(),
               names_to="age_cat",
               values_to="estimate") %>%
  separate(age_cat, c("a", "age_l", "age_r"), "_") %>%
  mutate(GEOID="78000",
         NAME="Virgin Islands",
         age_l=as.numeric(age_l),
         age_r=as.numeric(age_r)) %>%
  format_age_county_pop()

## Virgin Islands
MP <- read_csv("raw_data/mp0000120101.txt", col_names=F) %>%
  slice(1) %>%
  select(a_0_4=X7, a_5_9=X8, a_10_14=X9, a_15_19=X10, a_20_24=X11, a_25_29=X12,
         a_30_34=X13, a_35_39=X14, a_40_44=X15, a_45_49=X16, a_50_54=X17,
         a_55_59=X18, a_60_64=X19, a_65_69=X20, a_70_74=X21, a_75_79=X22,
         a_80_84=X23, a_85_99=X24) %>%
  pivot_longer(cols=everything(),
               names_to="age_cat",
               values_to="estimate") %>%
  separate(age_cat, c("a", "age_l", "age_r"), "_") %>%
  mutate(GEOID="69000",
         NAME="Commonwealth of the Northern Mariana Islands",
         age_l=as.numeric(age_l),
         age_r=as.numeric(age_r)) %>%
  format_age_county_pop()

## PR aggregated
PR <- load_age_county_pop(census_api_key = "c235e1b5620232fab506af060c5f8580604d89c1",
                          modeled_states = "PR",
                          census_year = 2010) %>%
  group_by(age_l, age_r) %>%
  summarize(estimate=sum(estimate),
            GEOID="72000",
            NAME="Puerto Rico") %>%
  format_age_county_pop()

## make geoid data
USpop_geoid_agecat <- load_age_county_pop(census_api_key = "c235e1b5620232fab506af060c5f8580604d89c1",
                                          modeled_states = c(state.abb, "PR", "DC"),
                                          census_year = 2010) %>%
  format_age_county_pop() %>%
  bind_rows(AS, GU, VI, MP, PR) %>%
  make_age_matrix()

rownames(USpop_geoid_agecat$p_age) <- USpop_geoid_agecat$GEOID[,1, drop=T]

US_age_geoid_pct <- USpop_geoid_agecat$p_age

usethis::use_data(US_age_geoid_pct, overwrite=T)

rownames(USpop_geoid_agecat$est_age) <- USpop_geoid_agecat$GEOID[,1, drop=T]

US_age_geoid_pop <- USpop_geoid_agecat$est_age

usethis::use_data(US_age_geoid_pop, overwrite=T)

## pull data for age-specific parameters
raw_age_estimates <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1pDfQ2SkO--F2WNfjZxx6t9V7K51IZW0bJKV5cwSCZfA/edit#gid=1769840547",
                                               sheet="age risk") %>%
  filter(publicly_available==T)
1
usethis::use_data(raw_age_estimates, overwrite=T)
