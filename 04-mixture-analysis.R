#  Description: 
#     Assess mixed effects of one group of exposures with qgcomp and WQS regression.
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

library(ggplot2)
library(qgcomp)
library(dplyr)
library(patchwork)
library(gWQS)

source("C:/Users/camer/OneDrive/Documents/SPH/MESA_exposome/scripts/get-weights-plot.R")


### 1. Prep data for analysis --------------------------------------------------

data_labels <- read.csv("C:/Users/camer/OneDrive/Documents/SPH/MESA_exposome/data_labels.csv") 
data <- readRDS("C:/Users/camer/OneDrive/Documents/SPH/MESA_exposome/air_pollution/data/mesa_air_pollution_dat.rds")

Xnm <- c("PM25_estimate", "NO2_estimate", "NOx_estimate", "LAC_estimate", "O3_estimate")

covars <- c("agebrainmri6c", "gender1", "race1c", "marital5", "income5", "educ1")

df <- dplyr::select(data, all_of(c("idno", Xnm, covars, "icv", "gm_vol", "wm_wmh", "wm_wmh_adj", "wm_fa"))) |> 
  na.omit() |> 
  mutate(gm_vol = gm_vol * 0.001)
df[Xnm] <- apply(df[Xnm], 2, FUN = scale, center = T, scale = T) |> as.data.frame()


### 2. Qgcomp ------------------------------------------------------------------

# Not controlling for other exposures
gm_vol_fit <- qgcomp.glm.noboot(as.formula(paste0("gm_vol", " ~ ", paste0(c(Xnm, covars, "icv"), collapse = " + "))),
                                 expnms = Xnm,
                                 df, family = gaussian(), q = 4, bayes = TRUE)



### 3. gWQS --------------------------------------------------------------------

gm_vol_wqs <- gwqs(as.formula(paste0("gm_vol", " ~ ", paste0(c("wqs", covars, "icv"), collapse = " + "))), 
                    mix_name = Xnm, data = df, 
                    q = 10, validation = 0.6, b = 5, b1_pos = FALSE, rh = 5,
                    family = "gaussian", seed = 2016)


### 4. Combined Plots ----------------------------------------------------------

p1 <- plot.weights(gm_vol_fit, data_labels[,1:2], group = "") + 
  ggtitle(label = "Qgcomp weights",
          subtitle = paste0("PSI: ", round(gm_vol_fit$psi, 2), "  (-PSI: ", round(gm_vol_fit$neg.psi, 2), ", +PSI: ", round(gm_vol_fit$pos.psi, 2), ")"))

wqs_results <- summary(gm_vol_wqs)$coefficients[2, ] |> round(2)
p2 <- gwqs_barplot(gm_vol_wqs) + 
  ggtitle(label = "WQS weights (estimate)",
          subtitle = paste0("WQS coefficient: ", wqs_results[1], ", 95% CI (", wqs_results[3], ", ", wqs_results[4], ")"))

p3 <- gwqs_boxplot(gm_vol_wqs) + 
  ggtitle(label = "WQS weights",
          subtitle = paste0("WQS coefficient: ", wqs_results[1], ", 95% CI (", wqs_results[3], ", ", wqs_results[4], ")"))

p1 + p2 + p3 + plot_layout(ncol = 3) + plot_annotation(title = "Air Pollution: associations with Grey Matter Volume (ml)")

ggsave("C:/Users/camer/OneDrive/Documents/SPH/MESA_exposome/air_pollution/output/mixture_analysis/gm_vol.png", 
       dpi = 300, 
       height = 6, 
       width = 16, 
       units = "in")


### 5. gWQS Positive and Negative ----------------------------------------------

wqs_fit <- gm_vol_wqs <- gwqs(gm_vol ~ pwqs + nwqs, 
                              mix_name = Xnm, data = df, 
                              q = 10, validation = 0.6, b = 5, b1_pos = FALSE, rh = 5,
                              family = "gaussian", seed = 2016)

results_pos <- summary(wqs_fit)$coefficients[2, ] |> round(2)
results_neg <- summary(wqs_fit)$coefficients[3, ] |> round(2)

out <- gwqs_barplot(wqs_fit) 

out2 <- gwqs_boxplot(wqs_fit)

out + out2 + plot_layout(ncol = 2, widths = c(4,5)) + plot_annotation(title = "Air Pollution: associations with Grey Matter Volume (ml)",
                                                                      subtitle = paste0("POSITIVE: [", results_pos[1], ", 95% CI (", results_pos[3], ", ", results_pos[4], ")]   ",
                                                                                        "NEGATIVE: [", results_neg[1], ", 95% CI (", results_neg[3], ", ", results_neg[4], ")]"))
ggsave("C:/Users/camer/OneDrive/Documents/SPH/MESA_exposome/air_pollution/output/mixture_analysis/gm_vol_wqs_pos_neg.png", 
       dpi = 300, 
       height = 6, 
       width = 12, 
       units = "in")
