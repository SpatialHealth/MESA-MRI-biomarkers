#  Description: 
#     Generate descriptive tables for MESA participants and exposures
#
#  Author(s): 
#     Cam Reimer 
#
#  Last Edited: 
#     2024-10-21
#
#  Additional Notes: 
#
# ------------------------------------------------------------------------------

library(tableone)
library(tidyverse)
library(dplyr)
library(labelled)
library(corrplot)

data <- readRDS("C:/Users/camer/OneDrive/Documents/SPH/MESA_exposome/air_pollution/data/mesa_air_pollution_dat.rds") 


### 1. Table 1: participant demographics ---------------------------------------

myVars <- c("agebrainmri6c", "gender1", "race1c", "marital5", "income5", "educ1", "site5c", 
            "gm_vol", "wm_wmh", "wm_wmh_adj", "wm_fa")
catVars <- c("gender1", "race1c", "marital5", "income5", "educ1", "site5c")

# no strata for now
t1 <- CreateTableOne(vars = myVars, data = data, factorVars = catVars, addOverall = FALSE)
print(t1, formatOptions = list(big.mark = ","), varLabels = TRUE, dropEqual = TRUE)

#Export
t1_out <- print(t1, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, varLabels = TRUE, dropEqual = TRUE)
write.csv(t1_out, file="C:/Users/camer/OneDrive/Documents/SPH/MESA_exposome/air_pollution/output/descriptive/table1.csv")


### 2. Table 2: exposures stratified by race/ethnicity -------------------------

exVars <- c("PM25_estimate", "NO2_estimate", "NOx_estimate", "LAC_estimate", "O3_estimate")
t2 <- CreateTableOne(vars = exVars, data = data, strata = "race1c", addOverall = TRUE)
print(t2, formatOptions = list(big.mark = ","), varLabels = TRUE, dropEqual = TRUE)

t2_out <- print(t2, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, varLabels = TRUE, dropEqual = TRUE)
write.csv(t2_out, file="C:/Users/camer/OneDrive/Documents/SPH/MESA_exposome/air_pollution/output/descriptive/table2.csv")


### 3. Table 3: exposure correlations ------------------------------------------

exposures <- data |> select(all_of(exVars))

names(exposures) <- c("PM 2.5", "NO2", "NOx", "LAC", "O3")

m <- round(cor(exposures, method = "spearman", use = "na.or.complete"), 2)
write.csv(m, file="C:/Users/camer/OneDrive/Documents/SPH/MESA_exposome/air_pollution/output/descriptive/exposure_correlations.csv")

png("C:/Users/camer/OneDrive/Documents/SPH/MESA_exposome/air_pollution/output/descriptive/exposure_correlations.png", 
    width = 9, height = 9, units = 'in', res = 600, pointsize = 14, bg = NA)
corrplot(m, method="color", 
         type="upper", 
         #order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         outline = F,
         tl.cex = 1.3, 
         number.cex = 1.2,
         cl.cex = 1.2,
         oma = c(0,0,0,0),
         bg = NA
)
dev.off()
