##### ISD Medits tutorial #######################################################
# Purpose: build a tidy workflow to filter species, standardize counts/biomass,
#          fit a truncated-Pareto (isdbayes + brms) model with an interaction
#          between fishing effort and temperature, visualize conditional effects,
#          and (optionally) summarize the interaction with a simple LM + 95% CIs.
################################################################################

# ---- 0) Libraries -------------------------------------------------------------
# Core wrangling/plotting
library(dplyr)
library(tidyr)
library(readxl)
library(tidyverse)  # if preferred; includes dplyr/ggplot2 etc.
library(ggplot2)
library(ggthemes)
library(scales)

# Spatial & utilities (keep if needed elsewhere)
library(gridExtra)
library(automap)
library(sf)
library(terra)
library(sp)
library(raster)
library(spdep)

# Modeling & post-processing
library(interactions)  # for interact_plot on the LM surface (viz helper)
library(isdbayes)      # truncated Pareto for size spectra
library(brms)          # Bayesian modeling
library(tidybayes)     # posterior tidy helpers
library(broom)         # tidy() with conf.int for LM
library(metR)          # geom_contour_fill
library(viridis)       # scale_fill_viridis_c

# ---- 1) Data objects ----------------------------------------------------------
# Assumes these two data.frames already exist in your workspace:
#   TA_GSA16_treated_new       = haul-level ancillary info (gear, area, swept area, coords…)
#   Medits_TC_GSA16_1994_2023  = length-class tallies per haul/species (TC table)
TA <- TA_GSA16_treated_new
TC <- Medits_TC_GSA16_1994_2023

# ---- 2) Build Haul ID + swept area on TA -------------------------------------
# NOTE: SWEPT_AREA in km^2 as in your formula: ((WING_OPENING*0.1) * DISTANCE)/1e6
TA <- TA %>%
  mutate(
    SWEPT_AREA = ((WING_OPENING * 0.1) * DISTANCE) / 1e6,
    HAUL_ID    = paste0("N", HAUL_NUMBER, "_Y", YEAR, "_Z", AREA)
  )
summary(TA$SWEPT_AREA)

# Quick sanity check (optional)
# plot(TA$lon, TA$lat)

# ---- 3) Build Haul ID on TC + species name -----------------------------------
TC <- TC %>%
  mutate(
    HAUL_ID  = paste0("N", HAUL_NUMBER, "_Y", YEAR, "_Z", AREA),
    Fullname = paste0(GENUS, "_", SPECIES)
  )
length(unique(TC$Fullname))
# print(unique(TC$Fullname))  # uncomment to inspect all species codes

# ---- 4) Merge TC + TA at haul level ------------------------------------------
TAC <- TC %>% left_join(TA, by = "HAUL_ID")

# ---- 5) Choose target species + filter years ---------------------------------
# You were repeatedly overwriting TAC_filter for multiple species.
# Use a single, explicit choice here. Change 'species_code' as needed.
# Examples:
#   "MERL_MER" (fish), "MULL_BAR", "MULL_SUR"
#   "PAPE_LON", "ARIS_FOL", "NEPR_NOR" (crustacea)
#   "LOLI_VUL", "ILLE_COI", "OCTO_VUL" (cephalopods)
species_code <- "MULL_BAR"

# Years to exclude (as in your original code)
years_exclude <- c(2013, 2014, 2017, 2020, 2021)

TAC_filter <- TAC %>%
  filter(Fullname == species_code) %>%
  filter(!(YEAR.y %in% years_exclude)) %>%  # YEAR.y comes from join; keep if present
  filter(YEAR.x > 2000)                      # or use YEAR if only one year field exists

# ---- 6) Standardization to counts per swept area ------------------------------
# Subsampling factor: fraction -> whole sample
TAC_filter <- TAC_filter %>%
  mutate(
    subcamp = WEIGHT_OF_THE_FRACTION / WEIGHT_OF_THE_SAMPLE_MEASURED,
    # Standardized number per km^2 (or per SWEPT_AREA unit)
    num_st  = (NUMBER_OF_INDIVIDUALS_IN_THE_LENGTH_CLASS_AND_MATURITY_STAGE * subcamp) /
      SWEPT_AREA
  )
summary(TAC_filter$subcamp)
summary(TAC_filter$num_st)

# ---- 7) Length–weight (example constants; set SPECIES-SPECIFIC 'a' and 'b') ---
# For "red shrimp" you had: a = 0.0025, b = 2.48.
# Keep/adjust per species_code if needed.
a <- 0.0025
b <- 2.48

# LENGTH_CLASS appears to be in mm; convert to cm for L^b if needed.
TAC_filter <- TAC_filter %>%
  mutate(
    ww_g   = a * ( (LENGTH_CLASS / 10) ^ b ),  # weight in grams
    num_st2 = floor(num_st)                    # integerize for uncount()
  ) %>%
  filter(!is.na(num_st2) & num_st2 > 0)

# Expand rows so each individual is one row (for isdbayes counts interface)
TAC_repeated <- TAC_filter %>%
  uncount(num_st2)

# ---- 8) Global xmin & site-level xmax for truncated Pareto --------------------
# If 'Site' is not present in your data, either:
#   - Use HAUL_ID as site, or
#   - Create a site definition (e.g., paste(AREA, YEAR, sep = "_"))
site_col <- if ("Site" %in% names(TAC_filter)) "Site" else "HAUL_ID"

# Compute xmax by site; xmin is global minimum
xmax_by_site <- TAC_filter %>%
  group_by(.data[[site_col]]) %>%
  summarise(xmax = max(ww_g, na.rm = TRUE), .groups = "drop")

xmin_global <- min(TAC_filter$ww_g, na.rm = TRUE)

# Attach xmin/xmax to each row (for vreal())
TAC_minmax <- TAC_filter %>%
  left_join(xmax_by_site, by = setNames("Site", site_col) %>% replace("Site", site_col)) %>%
  mutate(xmin = xmin_global)

# ---- 9) Center covariates (center-only; no scaling by sd) ---------------------
# IMPORTANT: center-only means: centered_var = original - mean(original)
# To back-transform later: original = centered_var + mean(original)
TAC_minmax <- TAC_minmax %>%
  mutate(
    # Responses/covariates used in model
    fish_c   = scale(Fishing,               center = TRUE, scale = FALSE),
    thetao_c = scale(thetao_mean,           center = TRUE, scale = FALSE),
    depth_c  = scale(log(HAULING_DEPTH),    center = TRUE, scale = FALSE),
    temp_c   = scale(BOTTOM_TEMPERATURE_BEGINNING, center = TRUE, scale = FALSE)
  )

# Store the original means for back-transforms in plots
mean_fish_orig  <- mean(TAC_minmax$Fishing,     na.rm = TRUE)
mean_theta_orig <- mean(TAC_minmax$thetao_mean, na.rm = TRUE)

# ---- 10) Fit truncated Pareto model (isdbayes + brms) -------------------------
# NOTE: You must define `stanvars` as per isdbayes docs BEFORE running brm().
#       (e.g., the helper that injects the required Stan functions/data for paretocounts()).
#       If you already have `stanvars` in your environment, this will just work.
#       Family must be isdbayes::paretocounts() and the response uses vreal(count, xmin, xmax).
#
# Threads/cores: adjust to your machine. `threads=12, cores=32` is aggressive.
# Initial values randomized to avoid problematic inits.

fit_nor <- brm(
  formula = ww_g | vreal(num_st2, xmin, xmax) ~ fish_c * thetao_c,
  data    = TAC_minmax,
  stanvars = stanvars,                  # <-- ensure this exists (see NOTE above)
  family  = paretocounts(),
  prior   = c(
    prior(normal(-2, 0.2), class = "Intercept"),
    prior(normal( 0, 0.1), class = "b")
  ),
  iter    = 4000,
  chains  = 2,
  threads = 12,
  cores   = 32,
  init    = "random",
  init_r  = 0.1
)

# Inspect model + data
print(fit_nor)
summary(TAC_minmax)

# ---- 11) Conditional effects on fish_c : thetao_c -----------------------------
# Extract the conditional surface for the interaction (on the model scale).
int_plot <- conditional_effects(fit_nor, effects = "fish_c:thetao_c")

# Turn into a data.frame and create un-centered covariates to label axes nicely
int_plot_data <- as_tibble(int_plot$`fish_c:thetao_c`) %>%
  mutate(
    # Back-transform because we used center-only:
    thetao_unscaled = thetao_c + mean_theta_orig,
    fish_unscaled   = fish_c   + mean_fish_orig
  )

# Optional: remove extreme predictions (your original used estimate__ > -3)
int_plot_data <- int_plot_data %>%
  filter(estimate__ > -3)

# ---- 12) Continuous interaction plot (contour of lambda over covariates) ------
ggplot(int_plot_data, aes(x = thetao_unscaled, y = fish_unscaled, z = estimate__)) +
  metR::geom_contour_fill() +                                       # filled contours
  geom_contour(aes(z = estimate__), color = "black",
               linetype = "dotted", linewidth = 0.3) +              # contour lines
  scale_fill_viridis_c(direction = -1) +
  theme_classic() +
  labs(
    x    = "Temperature (θᵒ; un-centered)",
    y    = "Fishing effort (un-centered)",
    fill = "λ (estimate__)"
  )

# ---- 13) Summarize the surface with a simple LM + 95% CI ----------------------
# NOTE: This LM is *only* a descriptive summary of the conditional-effects surface
#       (not a replacement for Bayesian inference). It’s useful for clean CIs and
#       helper plots like interact_plot().
mod_final <- lm(estimate__ ~ fish_c * thetao_c, data = int_plot_data)
summary(mod_final)

# 95% CI for coefficients (high precision printing to avoid identical lwr/upr)
options(pillar.sigfig = 12)
broom::tidy(mod_final, conf.int = TRUE, conf.level = 0.95)

# ---- 14) Categorical interaction lines (LM helper viz) ------------------------
# Choose moderator values (quantiles are usually more informative than 1/2/3)
modx_vals <- quantile(int_plot_data$thetao_c, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)

interact_plot(
  mod_final,
  pred         = fish_c,
  modx         = thetao_c,
  modx.values  = modx_vals,
  plot.points  = FALSE,
  point.size   = 2,
  interval     = FALSE
) +
  theme_few() +
  scale_color_viridis(discrete = TRUE) +
  labs(
    x = "Fishing effort (centered)",
    y = "λ (estimate__)",
    color = "θᵒ (centered)"
  )

# ---- 15) OPTIONAL: Confidence bands for the LM predictions --------------------
# This is optional and is separate from the brms model; useful for quick ribbons.
newdat <- expand.grid(
  fish_c   = seq(min(int_plot_data$fish_c, na.rm = TRUE),
                 max(int_plot_data$fish_c, na.rm = TRUE), length.out = 200),
  thetao_c = as.numeric(modx_vals)
)
pred_ci <- cbind(newdat,
                 predict(mod_final, newdata = newdat,
                         interval = "confidence", level = 0.95))

# Example ribbon plot (centered scales)
ggplot(pred_ci,
       aes(x = fish_c, y = fit, ymin = lwr, ymax = upr,
           color = factor(thetao_c), fill = factor(thetao_c))) +
  geom_ribbon(alpha = 0.15, linewidth = 0) +
  geom_line() +
  theme_few() +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "Fishing effort (centered)",
       y = "λ (estimate__ from LM)",
       color = "θᵒ (centered)", fill = "θᵒ (centered)")

# ---- 16) Troubleshooting notes ------------------------------------------------
# * If you see "object 'fish_c' not found": ensure you built int_plot_data from
#   conditional_effects(fit_nor, "fish_c:thetao_c") and are using THAT data.frame.
# * If you see "object 'int_plot_data_SUL' not found": that object never existed;
#   this script standardizes on `int_plot_data`.
# * If broom::tidy shows identical CI bounds, increase printed precision:
#     options(pillar.sigfig = 12)
# * If 'stanvars' is missing, define it per isdbayes doc for paretocounts().
# * Back-transforming centered covariates:
#     original = centered + mean(original)
#   Do not multiply by sd because we used center-only (scale = FALSE).
