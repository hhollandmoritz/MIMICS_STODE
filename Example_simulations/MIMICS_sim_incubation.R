library(dplyr)
library(ggplot2)
library(DT)
library(purrr)
library(Metrics)
library(here)
## Set working drive
setwd(here())

#----------------------------------------
# Load necessary model components
#----------------------------------------

# Load reverse Michaelis-Menton MIMICS function
source("functions/RXEQ.R")

# Load function to calculate model input variables
source("functions/calc_Tpars.R")

# Set MIMICS parameters
source("Parameters/MIMICS_parameters_MSBio_Incubation.R")  #--> Default "Sandbox" - Wieder et al. 2015

# Bring in MIMICS incubation simulation function
source("functions/MIMICS_sim_incubation.R")

#-----------------------------------------

# Bring in MSBio forcing data
MSB_data <- read.csv("Example_simulations/Data/MSBio_MIM_forcings.csv")
head(MSB_data)
LTER_SITE_DATA <- read.csv("Example_simulations/Data/LTER_SITE_1.csv")

# read incubation data 
incubation_simulation_summary <- read.table("Example_simulations/Data/incubation_summary_simulation.tsv",
                                            sep = "\t", header = T)


# Lidet data from wieder et al. 2015; units are as follows:
# MAT = degC
# MAP = mm
# SOC = kg C/m^2
# Sand% = percent
# clay = percent
# ANPP = gC m^-2y^-1
# Lignin = % 
# N = soil N conctent; units? I think %
# Litter C:N = unitless; ratio of 
# Lig_N = ratio of lignin percent to N percent


# For my data; is Permafrost C from fernando in mgC/g soil? or mgC/mg soil?

ArcticLTER <- LTER_SITE_DATA %>%
  filter(Site == "ARC")

# Litter input rate is calculated as:
#EST_LIT <- (ANPP / (365*24)) * 1e3 / 1e4
#ANPP is assumed to be in gC per year
# So I added 0.07 mg C per g dws per week
# Gonna use a bulk density estimate of 1.3 g/cm^3 based on this paper: https://www-sciencedirect-com.unh.idm.oclc.org/science/article/pii/S0034425720303771#f0010
# ml in 10cm x 1m^2 plot = 10*100*100 = 100000 cm^3 = 100000 mL
# = 1.3*100000 = 130000 g soil in m^2
# so 0.07 mgC * 130000 g soil = 9100 mgC/m^2 soil per week
# = 9100*52 = 473200 mg C/m^2 Soil per year
# = 473.2 g C/m^2 per year
# 



# final numbers
incubation_mimics_end <- incubation_simulation_summary %>%
  group_by(Carbon.Source) %>%
  select(Carbon.Source, Microbial.Community, mgC_g_dws_per_week, LitCNRatio, gwc, PermafrostCN, PermafrostC, PermafrostN, mean, se, Date) %>%
  filter(Date == max(Date)) %>%
  ungroup() %>%
  rename(ug_Carbon_per_g_soil_mean = `mean`,
         ug_Carbon_per_g_soil_se = `se`) %>%
  mutate(ID = row_number(),
         ANPP =  mgC_g_dws_per_week*130000*52/1000, # gC per m^2 per year; assuming bulk density of 1.3 still have to figure out the scaling
         SITE = ArcticLTER$Site,
         TSOI = ArcticLTER$TSOI, # Site soil temperature
         CLAY = ArcticLTER$CLAY, # may want to modify with NEON data?
         LIG = ifelse(Carbon.Source == "Lignin/Cellulose",
                      23.75098091,
                      0), #; % lignin; not sure if 0 is really the right option for the other treatments; since there is probalby some lignin in the soil
         N = ifelse(Carbon.Source == "Lignin/Cellulose", 3.309898167,
                    ifelse(Carbon.Source == "Exudate Cocktail", 2.001428571, 
                           0)),#??? will use litter N % but not sure if this is correct
         CN = ifelse(Carbon.Source == "Lignin/Cellulose", 36.78762706,
                     ifelse(Carbon.Source == "Exudate Cocktail", 26.57458958, 
                            0)), # taken from c-addition ratios; note that the LTER CN is nearly identical to Lignin/cellulose cocktail
         LIG_N = ifelse(Carbon.Source == "Lignin/Cellulose", 7.175743696,
                        ifelse(Carbon.Source == "Exudate Cocktail", 0, 
                               0)),
         MAT = -7, # from table in Wieder et al. 
         TINC = 4, # incubation temperature
         soilC_perc = ArcticLTER$SOC,
         GWC = gwc,
  ) %>%
  mutate(DAY = as.Date(Date) - as.Date("2024-10-31")) %>%
  mutate(CO2C_prop_mean = ug_Carbon_per_g_soil_mean,# unsure if I'm off by factors of 10, need to chekc units
         CO2C_prop_se = ug_Carbon_per_g_soil_se)

incubation_mimics_end %>%
  filter(Microbial.Community == "Permafrost Alone") %>%
  select(Carbon.Source, ID:GWC)

# Over time
incubation_mimics_dat <- incubation_simulation_summary %>%  group_by(Carbon.Source) %>%
  select(Carbon.Source, Microbial.Community, mgC_g_dws_per_week, LitCNRatio, gwc, PermafrostCN, PermafrostC, PermafrostN, Date, mean, se) %>% 
  ungroup() %>%
  rename(ug_Carbon_per_g_soil_mean = `mean`,
         ug_Carbon_per_g_soil_se = `se`) %>%
  left_join(incubation_mimics_end %>% select(Carbon.Source, Microbial.Community, ID),
            by = c("Carbon.Source", "Microbial.Community")) %>%
  mutate(ANPP =  mgC_g_dws_per_week*130000*52/1000, # gC per m^2 per year; assuming bulk density of 1.3 still have to figure out the scaling
         SITE = ArcticLTER$Site,
         TSOI = ArcticLTER$TSOI, # Site soil temperature
         CLAY = ArcticLTER$CLAY, # may want to modify with NEON data?
         LIG = ifelse(Carbon.Source == "Lignin/Cellulose",
                      23.75098091,
                      0), #; % lignin; not sure if 0 is really the right option for the other treatments; since there is probalby some lignin in the soil
         N = ifelse(Carbon.Source == "Lignin/Cellulose", 3.309898167,
                     ifelse(Carbon.Source == "Exudate Cocktail", 2.001428571, 
                            0)),#??? will use litter N % but not sure if this is correct
         CN = ifelse(Carbon.Source == "Lignin/Cellulose", 36.78762706,
                     ifelse(Carbon.Source == "Exudate Cocktail", 26.57458958, 
                            0)), # taken from c-addition ratios; note that the LTER CN is nearly identical to Lignin/cellulose cocktail
         LIG_N = ifelse(Carbon.Source == "Lignin/Cellulose", 7.175743696,
                             ifelse(Carbon.Source == "Exudate Cocktail", 0, 
                                    0)),
         MAT = -7, # from table in Wieder et al. 
         TINC = 4, # incubation temperature
         soilC_perc = ArcticLTER$SOC,
         GWC = gwc,
         ) %>%
  mutate(DAY = as.Date(Date) - as.Date("2024-10-31")) %>%
  mutate(CO2C_prop_mean = ug_Carbon_per_g_soil_mean, # unsure if I'm off by factors of 10, need to chekc units
         CO2C_prop_se = ug_Carbon_per_g_soil_se)
head(incubation_mimics_dat)



forcing_df <- MSB_data

forcing_df <- incubation_mimics_dat


# write out simulation end to incubation data
write.table(x = incubation_mimics_end, "Example_simulations/Data/incubation_end.tsv",
                                            sep = "\t", row.names = F)

# forcing data for east toolik
# need: incubation ID, SITE (TOOL), TSOI (soil temp?), clay content;
# lig N, CN, LIG_N, MAT, TINC <- Incubationt temperature?; Treatmetn; 
# MAT

#---------------------------------------------------------
# Run single incubation
#---------------------------------------------------------
#MIMinc_out <- MIMICS_INCUBATION(forcing_df[1,], days=105, step="hourly", output_type = 1)


#---------------------------------------------------------
#>> Run all rows in forcing dataset
#---------------------------------------------------------
forcing_df <- MSB_data
forcing_df <- incubation_mimics_end
MIMrun <- forcing_df %>% split(1:nrow(forcing_df)) %>% 
            map(~MIMICS_INCUBATION(df=., days=max(forcing_df$DAY), step="daily", output_type = 2)) %>% 
            bind_rows() 

MC_MIMICS <- MIMrun %>% left_join(incubation_mimics_end %>% select(-SITE), by="ID")


MIMrun_1 <- forcing_df %>% split(1:nrow(forcing_df)) %>% 
  map(~MIMICS_INCUBATION(df=., days=max(forcing_df$DAY), step="daily", output_type = 1)) %>% 
  bind_rows() 

MC_MIMICS_1 <- MIMrun_1 %>% left_join(incubation_mimics_dat %>% 
                                        select(-SITE, -DAY, -CO2C_prop_mean, -CO2C_prop_se), by=c("ID")) %>%
  left_join(incubation_mimics_dat %>% select(ID, DAY, CO2C_prop_mean, CO2C_prop_se) %>%
              mutate(DAY = as.numeric(DAY)), by = c("ID", "DAY"))


##########################################################
# Explore MIMICS incubation output
##########################################################

#############################
#Plot settings
# Set ggplot theme
theme_set(theme_bw())

# Palette
toolik <- c("lingonberry" ="#B72F40",
            "cloudberry" = "#FFA07A", 
            "sphagnum fimbriatum"="#C6EC91", 
            "prudhoe bay" ="#35756B", 
            "koyukuk river" = "#14B49B",
            "arctic blue" ="#94E1FF", 
            "crowberry" = "#4B0082", 
            "O horizon" = "#3E2723",
            "arctic sunrise" = "#FFBF29",
            "clay" = "#A6611A",
            "permafrost"="#A9A9A9")

Inocula_colors <- c("#94E1FF", "#FFBF29","#C6EC91")
Inocula_levels <- c("Permafrost Alone", "Active Layer + Permafrost", "Rhizosphere + Permafrost")

carbon_colors <- c("#3E2723","#A6611A","#4B0082")
carbon_levels <- c("No Carbon Added", "Lignin/Cellulose", "Exudate Cocktail")
carbon_names <- c("No Carbon Added", "Lignin/Cellulose", "Exudate Cocktail")
library(tidyverse)
#############################
Side_by_side.plot <- MC_MIMICS_1 %>%
  select(Carbon.Source, Microbial.Community, DAY, contains("prop")) %>%
  mutate(CO2C_prop_mean = CO2C_prop_mean,
         CO2C_prop_se = CO2C_prop_se) %>%
  pivot_longer(all_of(c("CO2_prop_totC", "CO2C_prop_mean")), names_to = "CO2Source", values_to = "CO2_prop") %>% 
  mutate(CO2C_prop_se = ifelse(CO2Source == "CO2C_prop_mean", CO2C_prop_se, NA)) %>%
  mutate(CO2Source = ifelse(CO2Source == "CO2_prop_totC", "MIMICS - CO2_prop_totC", "Incubation - Cum_ug_CO2_mean *10^-6")) %>%
  mutate(Microbial.Community = ifelse(grepl("MIMICS", CO2Source), "MIMICS Simulated\nCommunity", Microbial.Community)) %>%
  distinct() %>%
  mutate(CO2_prop = 100*CO2_prop,
         CO2C_prop_se = 100*CO2C_prop_se) %>%
 ggplot(aes(x = DAY, y = CO2_prop)) +
  geom_point(aes(color = Carbon.Source, shape = Microbial.Community), size = 3, alpha = 0.5) +
  geom_errorbar(aes(ymin = CO2_prop - CO2C_prop_se, 
                    ymax = CO2_prop + CO2C_prop_se, 
                    color = Carbon.Source), width = 0.3) +
  scale_color_manual(name = "Carbon Source", values = carbon_colors, breaks = carbon_levels) +
  scale_shape_manual(name = "Microbial Community", values = c(15, 16,17, 21), breaks = c(Inocula_levels, "MIMICS Simulated\nCommunity")) +
  #facet_grid(Carbon.Source ~ .) +
  ylab("Percent of total carbon respired") +
  guides(col = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom", legend.title.position = "top", legend.box.margin = margin(r = 0.2, l = 0.2))

Side_by_side.plot

Side_by_side.plot_agu <- Side_by_side.plot +
  guides(col = guide_none(), shape = guide_legend(nrow = 2)) +
  theme(aspect.ratio = 1,
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1)))
ggsave(Side_by_side.plot_agu, path = "~/Downloads", file = "side_by_side_mimicsInc.png", width = 5, height = 5)
ggsave(Side_by_side.plot, path = "~/Downloads", file = "All_ontop_side_by_side_mimicsInc.png", width = 6, height = 8)
  
MC_MIMICS_1 %>% select(Carbon.Source, DAY, MICr, MICK) %>%
  distinct() %>%
ggplot(aes(x = DAY, y = MICr/MICK)) +
  geom_point(aes(color = Carbon.Source, shape = Carbon.Source), size = 3, alpha = 0.5) +
  scale_color_manual(values = carbon_colors, breaks = carbon_levels) +
  scale_shape_manual(values = c(3, 16, 21))

ggplot(MC_MIMICS_1, aes(x = DAY, y = CO2_MICK)) +
  geom_point(aes(color = Carbon.Source, shape = Carbon.Source), size = 3, alpha = 0.5) +
  scale_color_manual(values = carbon_colors, breaks = carbon_levels)+
  scale_shape_manual(values = c(3, 16, 21))

ggplot(MC_MIMICS_1, aes(x = DAY, y = CO2_MICr)) +
  geom_point(aes(color = Carbon.Source, shape = Carbon.Source), size = 3, alpha = 0.5) +
  scale_color_manual(values = carbon_colors, breaks = carbon_levels)+
  scale_shape_manual(values = c(3, 16, 21))

##########################################################
# Plot lab vs MIMICS incubation cumulative respiration
##########################################################

plot_df <- MC_MIMICS #%>% filter(moisture.trt != 20)

MC_CO2Cp_cor <- round(cor(plot_df$CO2_prop_totC, plot_df$CO2C_prop), 4)
MC_CO2Cp_rmse <- round(rmse(plot_df$CO2_prop_totC, plot_df$CO2C_prop), 4)
MC_CO2Cp_fit <- lm(plot_df$CO2_prop_totC ~ plot_df$CO2C_prop)

# Extract coefficient values
intercept <- coef(MC_CO2Cp_fit)[1]
slope <- coef(MC_CO2Cp_fit)[2]

# Print the equation
line_eqn = paste0("y = ", round(slope, 6), "x + ", round(intercept, 4))

# Plot the relationship
plt <- ggplot(plot_df,
              aes(y=CO2_prop_totC, x=CO2C_prop, color=factor(moisture.trt))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth=1, alpha=0.6) +
  geom_point(size=3) +
  #geom_smooth(method = "lm", se = FALSE, alpha=0.8) +
  geom_line(stat="smooth",method = "lm", formula = y ~ x,
            linewidth = 1,
            #linetype ="dashed",
            alpha = 0.3) +
  ylab("MIMICS (CO2C / total C)") +
  xlab("Laboratory (CO2C / total C)") +
  labs(#title = "Calibrated kinetic terms & fW curve\nCumulative respiration after 105 day incubation
  #      ",#\nUsing site MAT and moisture controls on decomposition rate", #Default MIMICS parameters",
        subtitle = paste0("r^2 = ", MC_CO2Cp_cor, "   ", "RMSE = ", MC_CO2Cp_rmse),
        caption = "Dashed line = 1:1",
        color = "Moisture") +
  #ylim(0, 0.13) + xlim(0, 0.13) +
  theme_minimal()
plt


plot_df %>%
  ggplot(aes(x = DAY, y = CO2C_prop)) +
  geom_point(aes(color = SITE))
