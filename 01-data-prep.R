#  Description: 
#     Prep MESA data for exposure analysis with air pollution and MRI outcomes 
#
#  Author(s): 
#     Cam Reimer 
#
#  Last Edited: 
#     2024-10-03
#
#  Additional Notes: 
#
# ------------------------------------------------------------------------------

library(dplyr)
library(labelled)

# set working directory to where you've stored RDS objects. 
# this is highly recommended and much faster than reading sas7bdat files.
setwd("C:/Users/camer/OneDrive/Documents/SPH/MESA_exposome/data/Robjects")


### 1. Covariates --------------------------------------------------------------

# Exam 1: race/ethnicity, gender/sex, education,
mesa1 <- readRDS("MESAE1FinalLabel20240307.rds") |> 
  select(idno, race1c, gender1, educ1)

# Exam 5: age, household income, marital status, mesa site, (hypertension, hyperlipidemia, diabetes)
mesa5 <- readRDS("MESAe5_FinalLabel_20240213.rds") |> 
  select(idno, age5c, marital5, income5, site5c)
#,htn5c) #hypertension

phys <- readRDS("MESAe5_PhysAct_20130327.rds") |> 
  select(idno, pamvcm5c)

covs <- full_join(mesa1, mesa5, by = "idno") |> 
  full_join(phys, by = "idno")


### 2. Air Pollution Exposures -------------------------------------------------

air <- read.csv("../AirPollution/dr1055_wide.csv") 

# Get average exposures 
air$PM25_estimate <- air[, grepl("^PM25_estimate", names(air))] |> apply(MARGIN = 1, FUN = mean, na.rm = T)
air$NO2_estimate <- air[, grepl("^NO2_estimate", names(air))] |> apply(MARGIN = 1, FUN = mean, na.rm = T)
air$NOx_estimate <- air[, grepl("^NOx_estimate", names(air))] |> apply(MARGIN = 1, FUN = mean, na.rm = T)
air$O3_estimate <- air[, grepl("^O3_estimate", names(air))] |> apply(MARGIN = 1, FUN = mean, na.rm = T)
air$LAC_estimate <- air[, grepl("^LAC_estimate", names(air))] |> apply(MARGIN = 1, FUN = mean, na.rm = T)


### 3. MRI Outcomes ------------------------------------------------------------

# Brain MRI Region of Interest Volume 
vol <- readRDS("MESAe6_1stBMRIROIVol_20240513.rds") |>
  filter(is.nan(qc_code)) |>
  select(idno, agebrainmri6c,icv, gm) |>
  rename(gm_vol = gm)

# Fractional Anisotropy
fa <- readRDS("MESAe6_1stBMRIFAMUSE_20230801.rds") |> 
  filter(qc_flag == 0) |>
  select(idno, wm) |>
  rename(wm_fa = wm)

# White Matter Hyperintensity Volume
wmh <- readRDS("MESAe6_1stBMRIWMHVMUSE_20230801.rds") |> 
  filter(qc_flag == 0) |>
  select(idno, wm) |>
  rename(wm_wmh = wm)

# Merge all outcomes (n=1251 complete cases)
outcomes <- vol |>
  full_join(fa, by = "idno") |>
  full_join(wmh, by = "idno") |> 
  mutate(wm_wmh_adj = log(wm_wmh / icv))


### 4. Merge -------------------------------------------------------------------

out <- outcomes |> 
  left_join(covs, by = 'idno') |>
  left_join(air, by = 'idno')


### 5. Missingness + Exclusions ------------------------------------------------

missing <- colSums(is.na(out[names(covs)])) |> as.data.frame()
names(missing) <- "n"
missing <- missing %>% mutate(percent = round((n / nrow(out)) * 100, 1))


out <- out |> tidyr::drop_na(all_of(names(outcomes))) # n = 1251


### 6. Assign factor labels ----------------------------------------------------


out$gender1 <- factor(out$gender1, levels = c(0, 1), 
                      labels = c("Female", "Male"))
out$race1c <- factor(out$race1c, levels = c(1, 2, 3, 4),
                     labels = c("White", "Chinese", "Black", "Hispanic/Latino"))
out$site5c <- factor(out$site5c, levels = c(3, 4, 5, 6, 7, 8), 
                     labels = c("WFU", "COL", "JHU", "MN", "NWU", "UCLA"))
out$income5 <- factor(out$income5, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
                      labels = c("< $25,000", "< $25,000", "< $25,000", "< $25,000", "< $25,000", "< $25,000",
                                 "$25,000-$50,000", "$25,000-$50,000", "$25,000-$50,000", "$25,000-$50,000",
                                 "$50,000-$75,000",
                                 "$75,000 +", "$75,000 +", "$75,000 +", "$75,000 +"))
out$marital5 <- factor(out$marital5, levels = c(1, 2, 3, 4, 5, 6), 
                       labels = c("Married/Living with Partner", "Widowed", "Divorced", "Separated", "Never Married", "Prefer not to answer"))
out$educ1 <- factor(out$educ1, levels = c(0, 1, 2, 3, 4, 5, 6, 7, 8), 
                    labels = c("No Schooling", "Grades 1-8", "Grades 9-11", "High School/GED", "Some college but no degree", 
                               "Technical School Certificate", "Associate Degree", "Bachelor's Degree", "Graduate or Professional School"), ordered = TRUE)

var_label(out$agebrainmri6c) <- "Age at Exam 6 (MRI)"
var_label(out$gender1) <- "Gender/Sex"
var_label(out$race1c) <- "Race/Ethnicity"
var_label(out$marital5) <- "Marital Status"
var_label(out$income5) <- "Household Income"
var_label(out$educ1) <- "Highest Level of Education Completed"
var_label(out$site5c) <- "Site of Exam 5"
var_label(out$pamvcm5c) <- "Total Moderate and Vigorous Activity MET score"
var_label(out$gm_vol) <- "Gray Matter Volume"
var_label(out$wm_wmh) <- "White Matter Hyperintensity Volume"
var_label(out$wm_wmh_adj) <- "White Matter Hyperintensity Volume - Adjusted"
var_label(out$wm_fa) <- "White Matter Fractional Anisotropy"
var_label(out$PM25_estimate) <- "PM 2.5"
var_label(out$NO2_estimate) <- "NO2"
var_label(out$NOx_estimate) <- "NOx"
var_label(out$LAC_estimate) <- "LAC"
var_label(out$O3_estimate) <- "O3"


### 7. Hot deck Imputation -----------------------------------------------------

imp_vars <- c("educ1", "marital5", "income5")
imputed <- VIM::hotdeck(out,
                        variable = imp_vars,
                        ord_var = "agebrainmri6c",
                        domain_var = c("gender1" ,"race1c"),
                        imp_var = FALSE) 

missing <- colSums(is.na(imputed[names(imputed)])) |> as.data.frame()
names(missing) <- "n"
missing <- missing |> mutate(percent = round((n / nrow(imputed)) * 100, 1))

saveRDS(imputed, "C:/Users/camer/OneDrive/Documents/SPH/MESA_exposome/air_pollution/data/mesa_air_pollution_dat.rds")





