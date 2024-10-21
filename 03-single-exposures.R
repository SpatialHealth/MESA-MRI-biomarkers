#  Description: 
#     Assess associations with single exposures
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
library(kableExtra)

### 1. Prep data for analysis --------------------------------------------------

data_labels <- read.csv("C:/Users/camer/OneDrive/Documents/SPH/MESA_exposome/data_labels.csv") 
data <- readRDS("C:/Users/camer/OneDrive/Documents/SPH/MESA_exposome/air_pollution/data/mesa_air_pollution_dat.rds")

Xnm <- c("PM25_estimate", "NO2_estimate", "NOx_estimate", "LAC_estimate", "O3_estimate")

covars <- c("agebrainmri6c", "gender1", "race1c", "marital5", "income5", "educ1")

df <- dplyr::select(data, all_of(c("idno", Xnm, covars, "icv", "gm_vol", "wm_wmh", "wm_wmh_adj", "wm_fa"))) |> 
  na.omit() |> 
  mutate(gm_vol = gm_vol * 0.001)
df[Xnm] <- apply(df[Xnm], 2, FUN = scale, center = T, scale = T) |> as.data.frame()



### 2. Linear regression handler function --------------------------------------

linreg <- function(df, outcome, exvar, covars, sig_figs, null_aic) {
  xname <- paste0(c(exvar, covars), collapse = " + ")
  
  formula1 <- as.formula(paste(outcome,"~", xname))
  
  lm.fit <- do.call("lm", list(data = quote(df), formula1))
  
  main <- as.data.frame(cbind(coef(lm.fit), confint(lm.fit)))[2,] |> round(sig_figs)
  main <- cbind(main, 
                round(summary(lm.fit)$coefficients[2, 4], sig_figs),
                # paste0(main$V1," (", main$`2.5 %`, ", ", main$`97.5 %`, ")"),
                round(AIC(lm.fit) - null_aic, 2),
                round(summary(lm.fit)$adj.r.squared, 3))
  
  
  # colnames(main) = c("beta", "lowerCI", "upperCI", "pvalue", "full", "\u0394AIC", "Adj. R\U00B2")
  colnames(main) = c("beta", "lowerCI", "upperCI", "p-value", "\u0394AIC", "Adj. R\U00B2")
  
  return(main) 
}


### 3. Linear regression -------------------------------------------------------

gm_vol <- wm_wmh_adj <-  wm_fa <- data.frame()

gm_vol_null_aic <- AIC(lm(gm_vol ~ agebrainmri6c + gender1 + race1c + marital5 + income5 + educ1 + icv, df))
wm_wmh_adj_null_aic <- AIC(lm(wm_wmh_adj ~ agebrainmri6c + gender1 + race1c + marital5 + income5 + educ1, df))
wm_fa_null_aic <- AIC(lm(wm_fa ~ agebrainmri6c + gender1 + race1c + marital5 + income5 + educ1, df))

for (exvar in Xnm) {
  gm_vol <- linreg(df = df, 
                   outcome = "gm_vol", 
                   exvar = exvar, 
                   covars = c(covars, "icv"),
                   sig_figs = 2,
                   null_aic = gm_vol_null_aic) %>%
    cbind(exvar, .) %>%
    rbind(gm_vol, .)
  
  wm_wmh_adj <- linreg(df = df, 
                       outcome = "wm_wmh_adj", 
                       exvar = exvar, 
                       covars = covars,
                       sig_figs = 4,
                       null_aic = wm_wmh_adj_null_aic) %>%
    cbind(exvar, .) %>%
    rbind(wm_wmh_adj, .)
  
  wm_fa <- linreg(df = df, 
                  outcome = "wm_fa", 
                  exvar = exvar, 
                  covars = covars,
                  sig_figs = 5,
                  null_aic = wm_fa_null_aic) %>%
    cbind(exvar, .) %>%
    rbind(wm_fa, .)
  
}

gm_vol <- left_join(gm_vol, data_labels, by = join_by("exvar" == "var"))
wm_wmh_adj <- left_join(wm_wmh_adj, data_labels, by = join_by("exvar" == "var"))
wm_fa <- left_join(wm_fa, data_labels, by = join_by("exvar" == "var"))

write.csv(gm_vol, "C:/Users/camer/OneDrive/Documents/SPH/MESA_exposome/air_pollution/output/single_exposures/linreg_gm_vol.csv")
write.csv(wm_wmh_adj, "C:/Users/camer/OneDrive/Documents/SPH/MESA_exposome/air_pollution/output/single_exposures/linreg_wm_wmh_adj.csv.csv")
write.csv(wm_fa, "C:/Users/camer/OneDrive/Documents/SPH/MESA_exposome/air_pollution/output/single_exposures/linreg_wm_fa.csv")


### 4. Format Tables -----------------------------------------------------------

for (out_name in c("gm_vol", "wm_wmh_adj", "wm_fa")) {
  
  out_var <- get(out_name)
  out <- out_var[, c("label", "beta", "p-value", "\u0394AIC", "Adj. R\U00B2")]
  colnames(out)[1] <- c("Exposures")
  alpha <- 0.05
  bolded <- which(out$`p-value` < alpha)
  
  
  ### FIX UNITS ######
  
  # k <- kbl(out, caption = paste0(out_name, ' (ml)')) %>%
  k <- kbl(out, caption = out_name) %>%
    kable_paper("striped", full_width = F) %>%
    kable_styling() %>%
    row_spec(bolded,bold=T,hline_after = T) 
  
  print(k)
  
  save_kable(k,
             paste0("C:/Users/camer/OneDrive/Documents/SPH/MESA_exposome/air_pollution/output/single_exposures/results_", out_name, ".jpeg"),
             bs_theme = "flatly")
  
}



