#### DEPTH ANALYSES####

##### ISD Medits tutorial 

#library#
library(dplyr)
library(tidyr)
library(readxl)
library(tidyverse)
library(gridExtra)
library(scales)
library(automap)
library(sf)
library(terra)
library(sp)
library(raster)
library(spdep)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(interactions)
library(isdbayes)
library(brms)
library(tidybayes)

# Organization of data #

#data
TA <- TA_GSA16_treated_new
TC <- Medits_TC_GSA16_1994_2023

#ID code
#TA
TA$SWEPT_AREA <- ((TA$WING_OPENING*0.1)*TA$DISTANCE)/1000000
summary(TA$SWEPT_AREA)
TA$HAUL_ID <- paste0("N",TA$HAUL_NUMBER,"_Y",TA$YEAR,"_Z",TA$AREA)

plot(TA$lon,TA$lat)

#TC
TC$HAUL_ID <- paste0("N",TC$HAUL_NUMBER,"_Y",TC$YEAR,"_Z",TC$AREA)
single_haul_TC <- TC %>% dplyr::filter(HAUL_ID == "N10_Y2015_Z10")
single_haul_TC

#Species ID
TC$Fullname <- paste0(TC$GENUS,"_",TC$SPECIES)
length(unique(TC$Fullname))
print(unique(c$Fullname))

#Unique
TAC <- TC %>% left_join(TA, by = "HAUL_ID")

TAC_filter <- TAC %>%
  filter(YEAR.x == "2000")

TAC_filter <- TAC %>%
  filter(YEAR.x == "2010")

TAC_filter <- TAC %>%
  filter(YEAR.x == "2023")

TAC_filter <- TAC %>%
  filter(Fullname == "MERL_MER") %>%
  filter(., YEAR.y!= '2020') %>%
  filter(., YEAR.y!= '2014') %>%
  filter(., YEAR.y!= '2013') %>%
  filter(., YEAR.y!= '2017') %>%
  filter(., YEAR.y!= '2021') %>%
  filter(., YEAR.x > 2000)

TAC_filter <- TAC %>%
  filter(Fullname == "PAPE_LON") %>%
  filter(., YEAR.y!= '2020') %>%
  filter(., YEAR.y!= '2014') %>%
  filter(., YEAR.y!= '2013') %>%
  filter(., YEAR.y!= '2017') %>%
  filter(., YEAR.y!= '2021') %>%
  filter(., YEAR.x > 2000)


TAC_filter <- TAC %>%
  filter(Fullname == "LOLI_VUL")

TAC_filter <- TAC %>%
  filter(Fullname == "NEPR_NOR")


TAC_minmax <- case_when(
  HAULING_DEPTH < 400 ~ "Low",
  HAULING_DEPTH >= 400 & bmi < 600 ~ "Medium",
  HAULING_DEPTH >= 600 ~ "High"
)


#Standartization
TAC_filter$subcamp <-  TAC_filter$WEIGHT_OF_THE_FRACTION/TAC_filter$WEIGHT_OF_THE_SAMPLE_MEASURED
summary(TAC_filter$subcamp)

TAC_filter$num_st<-(TAC_filter$NUMBER_OF_INDIVIDUALS_IN_THE_LENGTH_CLASS_AND_MATURITY_STAGE*TAC_filter$subcamp)/TAC_filter$SWEPT_AREA
summary(TAC_filter$num_st)

#1) Lenght to Weight Merluzzo

a= 0.00501

b= 3.12

#2)Lenght to Weight Trilia

a=0.00871

b=3.09

#3)Lenght to Weight Red Shrimp

a= 0.0025

b= 2.48

TAC_filter$ww_g<-a*(TAC_filter$LENGTH_CLASS/10)^b

TAC_filter$num_st2<-floor(TAC_filter$num_st)

TAC_filter <- TAC_filter %>%
  filter(., num_st2 != 'NA')

TAC_repeated <- TAC_filter %>%
  uncount(num_st2)

# use the posterior lambdas to simulate body sizes at each site/sample/etc. Use the global xmin and site-level xmax values to do this instead of
# the culled values.

TAC_minmax = TAC_filter %>% 
  group_by(Site) %>%
  # ungroup %>% 
  mutate(xmax = max(ww_g)) %>%  # site xmax
  ungroup %>%
  mutate(xmin = min(ww_g))  

#fit model

#1) fit model

TAC_minmax = TAC_minmax %>% 
  mutate(ww_c= scale(ww_g, center = T, scale = F),
         thetao_c = scale(thetao_mean, center = T, scale = F),
         fish_c = scale(Fishing, center = T, scale = F),
         depth_c = scale(log(HAULING_DEPTH), center = T, scale = F))

#### FIT MODEL ####


fit_TAC = brm(ww_g| vreal(num_st2, xmin, xmax) ~ fish_c*depth_c + (1 | Site) + (1| HAUL_ID) + (1| YEAR.x),
              data = TAC_minmax,
              stanvars = stanvars,    # required for truncated Pareto via isdbayes package
              family = paretocounts(),# required for truncated Pareto via isdbayes package
              prior = c(prior(normal(-2, 0.2), class = "Intercept"),
                        prior(normal(0, 0.1), class = "b"),
                        prior(exponential(7), class = "sd")),
              iter = 300,
              chains = 1)


fit_TAC$data
summary(fit_TAC)

save.image("Savina_merluzzoISD.RData")
ggsave("grafico_thetao_fish.tiff",last_plot(),device = "tiff",dpi = "retina")

#### - PLOTING FIGURES - ###
##Figure 1##
thetao_mean = mean(unique(TAC_minmax$thetao_mean))
sd_thetao = sd(unique(TAC_minmax$thetao_mean))
mean_fish = mean(unique(TAC_minmax$Fishing))
sd_fish = sd(unique(TAC_minmax$Fishing))

fish_c = quantile(fit_TAC$data$fish_c, probs = c(0.25, 0.5, 0.75)) %>% 
  as_tibble() %>% rename(fish_c = value)

thetao_c = quantile(fit_TAC$data$thetao_c, probs = c(0.25, 0.5, 0.75)) %>% 
  as_tibble() %>% rename(thetao_c = value)

temp_c = quantile(fit_TAC$data$temp_c, probs = c(0.25, 0.5, 0.75)) %>% 
  as_tibble() %>% rename(temp_c = value)

int_plot = conditional_effects(fit_TAC, effects = "fish_c", conditions = thetao_c)

int_plot

int_plot = conditional_effects(fit_TAC, effects = "fish_c:thetao_c")

int_plot

int_plot_fish = plot(int_plot, plot = FALSE)[[1]]

# Customize the plot with ggplot2 themes
int_plot_fish + labs(title = "Fishing effort effect on Merluzzo merluzzo", x = "Fishing effort", y = "λ", color = fish_c) +
  theme_few()

int_plot_data = int_plot$`fish_c` %>% as_tibble() %>% 
  mutate(mat = (thetao_c*sd_thetao) + thetao_mean,
         fish = (mean_fish*sd_fish) + mean_fish)

data$depth_c <- as.numeric(data$depth_c)

#%>% 
mutate(quantile_fish = case_when(fish == min(fish) ~ "Low Fishing",
                                 fish == max(fish) ~ "High Fishing",
                                 TRUE ~ "Median Fishing"))
