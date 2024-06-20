########################################################################################################################################################
# This file will produce the plots of the climate trajectories shown in the paper
# 
#
# Author: Julian Rogger
# Last update: 30.05.24
########################################################################################################################################################

# Libraries
library(ggplot2)
library(raster)
library(scales)
library(ggpubr)
library(RColorBrewer)
library(mgcv)
library(cowplot)
library(viridis)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(chronosphere)
library(khroma)

# Import cluster analysis results ------------------------------------------------------------------------------------------------
out_df_0 <- readRDS("../output/exp_A_out/summary_measures_out.rds")
out_df_56 <- readRDS("../output/exp_B_out/summary_measures_out.rds")
out_df_200 <- readRDS("../output/exp_C_out/summary_measures_out.rds")
out_df_250 <- readRDS("../output/exp_D_out/summary_measures_out.rds") 

t_out_0 <- readRDS("../output/exp_A_out/appended_summary_table_out.rds") # temporal out
t_out_56 <- readRDS("../output/exp_B_out/appended_summary_table_out.rds") # temporal out
t_out_200 <- readRDS("../output/exp_C_out/appended_summary_table_out.rds") # temporal out
t_out_250 <- readRDS("../output/exp_D_out/appended_summary_table_out.rds") # temporal out

# Get experimental type assigned in 'degassing_experiment_analysis.R'
t_out_0$exp_type = NA
t_out_56$exp_type = NA
t_out_250$exp_type = NA
t_out_200$exp_type = NA
for(i in unique(out_df_0$exp_ID)){
  t_out_0$exp_type[t_out_0$exp_ID == i] <- out_df_0$exp_type[out_df_0$exp_ID == i]
}
for(i in unique(out_df_56$exp_ID)){
  t_out_56$exp_type[t_out_56$exp_ID == i] <- out_df_56$exp_type[out_df_56$exp_ID == i]
}
for(i in unique(out_df_200$exp_ID)){
  t_out_200$exp_type[t_out_200$exp_ID == i] <- out_df_200$exp_type[out_df_200$exp_ID == i]
}
for(i in unique(out_df_250$exp_ID)){
  t_out_250$exp_type[t_out_250$exp_ID == i] <- out_df_250$exp_type[out_df_250$exp_ID == i]
}

t_out_0$year <- 0
t_out_56$year <- 56
t_out_200$year <- 200
t_out_250$year <- 250

out_df_0$year <- 0
out_df_56$year <- 56
out_df_200$year <- 200
out_df_250$year <- 250


# Import cluster results of alternative boundary conditions ----------------------------------------
# Triassic-Jurassic: 6-fold weathering enhancement
out_df_200_S1 <- readRDS("../output/exp_S1_out/summary_measures_out.rds")
t_out_200_S1 <- readRDS("../output/exp_S1_out/appended_summary_table_out.rds")# temporal out
t_out_200_S1$exp_type = NA
for(i in unique(out_df_200_S1$exp_ID)){
  t_out_200_S1$exp_type[t_out_200_S1$exp_ID == i] <- out_df_200_S1$exp_type[out_df_200_S1$exp_ID == i]
}
t_out_200_S1$year <- 200
out_df_200_S1$year <- 200

# Permian-Triassic: 6-fold weathering enhancement 
out_df_250_S2 <- readRDS("../output/exp_S2_out/summary_measures_out.rds")
t_out_250_S2 <- readRDS("../output/exp_S2_out/appended_summary_table_out.rds")# temporal out
t_out_250_S2$exp_type = NA
for(i in unique(out_df_250_S2$exp_ID)){
  t_out_250_S2$exp_type[t_out_250_S2$exp_ID == i] <- out_df_250_S2$exp_type[out_df_250_S2$exp_ID == i]
}
t_out_250_S2$year <- 250
out_df_250_S2$year <- 250


# Permian-Triassic: 800 ppm starting CO2 & 50'000 Gt degassing
out_df_250_S3 <- readRDS("../output/exp_S3_out/summary_measures_out.rds")
t_out_250_S3 <- readRDS("../output/exp_S3_out/appended_summary_table_out.rds")# temporal out
t_out_250_S3$exp_type = NA
for(i in unique(out_df_250_S3$exp_ID)){
  t_out_250_S3$exp_type[t_out_250_S3$exp_ID == i] <- out_df_250_S3$exp_type[out_df_250_S3$exp_ID == i]
}
t_out_250_S3$year <- 250
out_df_250_S3$year <- 250


# Paleocene-Eocene: 500 ppm starting CO2, 8'000 Gt C degassing
out_df_56_S4 <- readRDS("../output/exp_S4_out/summary_measures_out.rds")
t_out_56_S4 <- readRDS("../output/exp_S4_out/appended_summary_table_out.rds")# temporal out
t_out_56_S4$exp_type = NA
for(i in unique(out_df_56_S4$exp_ID)){
  t_out_56_S4$exp_type[t_out_56_S4$exp_ID == i] <- out_df_56_S4$exp_type[out_df_56_S4$exp_ID == i]
}
t_out_56_S4$year <- 56
out_df_56_S4$year <- 56



########################################################################################################################################################

# 1) Climate trajectories Permian-Triassic  ----------------------------------------------------------------------
index <- t_out_250$exp_ID[t_out_250$exp_type =="new_steady" & t_out_250$time == 8.0 & t_out_250$locb < 1e+12]
t_out_250$exp_type[t_out_250$exp_ID %in% index] <- "die_out"

GAST_dens_250 <- ggplot(subset(t_out_250, preplant_par == 0.25), aes(time-1, GAST, group = exp_ID)) + 
  geom_line(alpha=0.025, linewidth = 0.2, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  geom_point(data = subset(t_out_250, exp_ID == 577), aes(time-1, GAST), shape = 2, col = "lightblue", size = 0.5) +
  geom_point(data = subset(t_out_250, exp_ID == 12122), aes(time-1, GAST), shape = 4, col = "pink", size = 0.5) +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#972e1a","#39259E","#BFA037","black"), labels = c("Die out", "New steady\n state", "Transient", "Return")) +
  xlab("Years relative to LIP onset (Ma)") +
  ylab("GAST (°C)") +
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "none")
ggsave("climate_trajectories_PT_250.png", plot = GAST_dens_250, device = "png", path = "./", scale = 1, width = 12, height = 10, units = "cm", dpi = 400)


# 2) Max warming ~ biological parameters -----------------------------------------------------------------------------------------
out_df_250_split_disp <- split(out_df_250, f=list(out_df_250$dispersal_par, out_df_250$preplant_par))
disp_summary <- lapply(out_df_250_split_disp, FUN = function(list_entry){
  df <- data.frame(
    "dispersal_par" = mean(list_entry$dispersal_par),
    "preplant_par" = mean(list_entry$preplant_par),
    "med_GAST_max" = median(list_entry$GAST_max),
    "up_GAST_max" = quantile(list_entry$GAST_max, 1),
    "low_GAST_max" = quantile(list_entry$GAST_max, 0),
    "med_CO2_ppm_max" = median(list_entry$CO2_ppm_max ),
    "up_CO2_ppm_max" = quantile(list_entry$CO2_ppm_max, 1),
    "low_CO2_ppm_max" = quantile(list_entry$CO2_ppm_max, 0),
    "med_pre_degass_average_temp" = median(list_entry$pre_degass_average_temp),
    "low_pre_degass_average_temp" = quantile(list_entry$pre_degass_average_temp, 1),
    "up_pre_degass_average_temp" = quantile(list_entry$pre_degass_average_temp, 0)
  )
  return(df)
})
disp_summary <- do.call(rbind.data.frame, disp_summary)

out_df_250_split_therm <- split(out_df_250, f=list(out_df_250$thermal_adapt_par, out_df_250$preplant_par))
therm_summary <- lapply(out_df_250_split_therm, FUN = function(list_entry){
  df <- data.frame(
    "thermal_adapt_par" = mean(list_entry$thermal_adapt_par),
    "preplant_par" = mean(list_entry$preplant_par),
    "med_GAST_max" = median(list_entry$GAST_max),
    "up_GAST_max" = quantile(list_entry$GAST_max, 1),
    "low_GAST_max" = quantile(list_entry$GAST_max, 0),
    "med_CO2_ppm_max" = median(list_entry$CO2_ppm_max ),
    "up_CO2_ppm_max" = quantile(list_entry$CO2_ppm_max, 1),
    "low_CO2_ppm_max" = quantile(list_entry$CO2_ppm_max, 0),
    "med_pre_degass_average_temp" = median(list_entry$pre_degass_average_temp),
    "low_pre_degass_average_temp" = quantile(list_entry$pre_degass_average_temp, 1),
    "up_pre_degass_average_temp" = quantile(list_entry$pre_degass_average_temp, 0)
  )
  return(df)
})
therm_summary <- do.call(rbind.data.frame, therm_summary)
out_df_250_split_therm <- split(out_df_250, f=list(out_df_250$aridity_adapt_par, out_df_250$preplant_par))
arid_summary <- lapply(out_df_250_split_therm, FUN = function(list_entry){
  df <- data.frame(
    "aridity_adapt_par" = mean(list_entry$aridity_adapt_par),
    "preplant_par" = mean(list_entry$preplant_par),
    "med_GAST_max" = median(list_entry$GAST_max),
    "up_GAST_max" = quantile(list_entry$GAST_max, 1),
    "low_GAST_max" = quantile(list_entry$GAST_max, 0),
    "med_CO2_ppm_max" = median(list_entry$CO2_ppm_max ),
    "up_CO2_ppm_max" = quantile(list_entry$CO2_ppm_max, 1),
    "low_CO2_ppm_max" = quantile(list_entry$CO2_ppm_max, 0),
    "med_pre_degass_average_temp" = median(list_entry$pre_degass_average_temp),
    "low_pre_degass_average_temp" = quantile(list_entry$pre_degass_average_temp, 1),
    "up_pre_degass_average_temp" = quantile(list_entry$pre_degass_average_temp, 0)
  )
  return(df)
})
arid_summary <- do.call(rbind.data.frame, arid_summary)


# Temperature increase
disp_pl <- ggplot(data = disp_summary) + 
  geom_ribbon(aes(x = (dispersal_par/1e+5)*1e+6, ymin = low_GAST_max-med_pre_degass_average_temp, ymax = up_GAST_max-med_pre_degass_average_temp), fill = "grey", alpha = 0.6) + # Convert dispersal units to km Myr-1
  theme_classic()  + 
  xlab(expression("Scale of dispersal kernel (km"~Ma^-1*")")) + 
  ylab(expression(Delta*"GAST"~"(°C)")) +
  theme(text = element_text(size = 6)) + 
  theme(panel.grid.minor = element_blank()) + 
  scale_x_log10(labels = label_number(drop0trailing = TRUE)) + 
  annotation_logticks(sides="b", size = 0.2, short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm"))
  
therm_pl <- ggplot(data = therm_summary) + 
  geom_ribbon(aes(x = thermal_adapt_par*1e+6, ymin = low_GAST_max-med_pre_degass_average_temp, ymax = up_GAST_max-med_pre_degass_average_temp), fill = "grey", alpha = 0.6) +
  theme_classic()   + 
  xlab(expression("Temperature niche adaptation rate (°C"~Ma^-1*")")) + 
  ylab(expression(Delta*"GAST"~"(°C)")) +
  scale_x_log10(labels = label_number(drop0trailing = TRUE)) + 
  annotation_logticks(sides="b", size = 0.2, short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
  theme(text = element_text(size = 6))  + 
  theme(panel.grid.minor = element_blank())


arid_pl <- ggplot(data = arid_summary) + 
  geom_ribbon(aes(x = (aridity_adapt_par/1e+5)*1e+6, ymin = low_GAST_max-med_pre_degass_average_temp, ymax = up_GAST_max-med_pre_degass_average_temp), fill = "grey", alpha = 0.6) + # Convert dispersal units to - Myr-1
  theme_classic()  + 
  xlab(expression("Aridity niche adaptation rate (BAI units"~Ma^-1*")")) + 
  ylab(expression(Delta*"GAST"~"(°C)")) +
  scale_x_log10(labels = label_number(drop0trailing = TRUE)) + 
  annotation_logticks(sides="b", size = 0.2, short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
  theme(text = element_text(size = 6))  + 
  theme(panel.grid.minor = element_blank()) 

param_sens <- ggarrange(disp_pl, therm_pl + ylab(""), arid_pl + ylab(""), labels=c("A"), ncol=3, nrow=1, common.legend = TRUE, legend="bottom", font.label = list(size=10), hjust = -0.5)
ggsave("GAST_vs_parameters_250Ma.png", plot = param_sens, device = "png", path = "./", scale = 1, width = 18, height = 6, units = "cm", dpi = 400)

# 3) Percentage of scenarios per parameter value ------------------------------------------------------------------------------------------------------------
index <- t_out_250$exp_ID[(t_out_250$exp_type =="new_steady" | t_out_250$exp_type == "die_out") & t_out_250$time == 8.0 & t_out_250$locb < 1e+12]
t_out_250$exp_type[t_out_250$exp_ID %in% index] <- "die_out"
out_df_250$exp_type[out_df_250$exp_ID %in% index] <- "die_out"


out_df_250_split <- split(out_df_250, f=out_df_250$dispersal_par)
disp_num_exp <- lapply(out_df_250_split, FUN = function(entry){
  num_die_out <- length(entry[entry$exp_type == "die_out", 1])
  num_return <- length(entry[entry$exp_type == "return", 1])
  num_new_st <- length(entry[entry$exp_type == "new_steady", 1])
  num_otw <- length(entry[entry$exp_type == "otw", 1])
  sum_exp <- sum(num_die_out, num_return, num_new_st, num_otw)
  df <- data.frame("dispersal_par" = mean(entry$dispersal_par), "num_die_out" = num_die_out, "num_return" = num_return, "num_new_st"=num_new_st, "num_otw" = num_otw, "sum_exp" = sum_exp)
  return(df)
})
disp_num_exp_df <- do.call(rbind.data.frame, disp_num_exp)

out_df_250_split_therm <- split(out_df_250, f=out_df_250$thermal_adapt_par)
therm_num_exp <- lapply(out_df_250_split_therm, FUN = function(entry){
  num_die_out <- length(entry[entry$exp_type == "die_out", 1])
  num_return <- length(entry[entry$exp_type == "return", 1])
  num_new_st <- length(entry[entry$exp_type == "new_steady", 1])
  num_otw <- length(entry[entry$exp_type == "otw", 1])
  sum_exp <- sum(num_die_out, num_return, num_new_st, num_otw)
  df <- data.frame("thermal_adapt_par" = mean(entry$thermal_adapt_par), "num_die_out" = num_die_out, "num_return" = num_return, "num_new_st"=num_new_st, "num_otw" = num_otw, "sum_exp" = sum_exp)
  return(df)
})
therm_num_exp_df <- do.call(rbind.data.frame, therm_num_exp)

out_df_250_split_arid <- split(out_df_250, f=out_df_250$aridity_adapt_par)
arid_num_exp <- lapply(out_df_250_split_arid, FUN = function(entry){
  num_die_out <- length(entry[entry$exp_type == "die_out", 1])
  num_return <- length(entry[entry$exp_type == "return", 1])
  num_new_st <- length(entry[entry$exp_type == "new_steady", 1])
  num_otw <- length(entry[entry$exp_type == "otw", 1])
  sum_exp <- sum(num_die_out, num_return, num_new_st, num_otw)
  df <- data.frame("aridity_adapt_par" = mean(entry$aridity_adapt_par), "num_die_out" = num_die_out, "num_return" = num_return, "num_new_st"=num_new_st, "num_otw" = num_otw, "sum_exp" = sum_exp)
  return(df)
})
arid_num_exp_df <- do.call(rbind.data.frame, arid_num_exp)


disp_per <- ggplot(data = disp_num_exp_df) + 
  geom_line(aes(x = (dispersal_par/1e+5)*1e+6, y = num_die_out/(sum_exp)*100, col = "Die\nout"), size=2, alpha = 0.7) + 
  geom_line(aes(x = (dispersal_par/1e+5)*1e+6, y = num_new_st/(sum_exp)*100, col = "New steady\nstate"), size=2, alpha = 0.7) + 
  geom_line(aes(x = (dispersal_par/1e+5)*1e+6, y = num_otw/(sum_exp)*100, col = "Transient"), size=2, alpha = 0.7) + 
  geom_line(aes(x = (dispersal_par/1e+5)*1e+6, y = num_return/(sum_exp)*100, col = "Return to\ninitial conditions"), size=2, alpha = 0.7) + 
  scale_color_manual(values=c("#972e1a","#074578","#BFA037","black"), 
                     name="Climate trajectory", 
                     limits = c("Die\nout", "New steady\nstate", "Transient", "Return to\ninitial conditions")) + 
  theme_classic() + 
  xlab(expression("Scale of dispersal kernel (km"~Ma^-1*")")) + 
  ylab("Outcome frequency (%)") + 
  scale_x_log10(labels = label_number(drop0trailing = TRUE))+
  annotation_logticks(sides="b", size = 0.2, short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
  theme(text = element_text(size = 6), 
        legend.position = "bottom", 
        legend.justification = 0.25) +
  ylim(c(0, 100))

therm_per <- ggplot(data = therm_num_exp_df) + 
  geom_line(aes(x = thermal_adapt_par*(1e+6)  , y = num_die_out/(sum_exp)*100, col = "Die\nout"), size=2, alpha = 0.7) + 
  geom_line(aes(x = thermal_adapt_par*(1e+6) , y = num_new_st/(sum_exp)*100, col = "New steady\nstate"), size=2, alpha = 0.7) + 
  geom_line(aes(x = thermal_adapt_par*(1e+6), y = num_otw/(sum_exp)*100, col = "Transient"), size=2, alpha = 0.7) + 
  geom_line(aes(x = thermal_adapt_par*(1e+6) , y = num_return/(sum_exp)*100, col = "Return to\ninitial conditions"), size=2, alpha = 0.7) + 
  scale_color_manual(values=c("#972e1a","#074578","#BFA037","black"), 
                     name="Climate trajectory", 
                     limits = c("Die\nout", "New steady\nstate", "Transient", "Return to\ninitial conditions")) + 
  theme_classic() + 
  xlab(expression("Temperature niche adaptation rate (°C"~Ma^-1*")")) + 
  ylab("") + 
  scale_x_log10(labels = label_number(drop0trailing = TRUE))+
  annotation_logticks(sides="b", size = 0.2, short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
  theme(text = element_text(size = 6), 
        legend.position = "bottom", 
        legend.justification = 0.25) +
  ylim(c(0, 100))

arid_per <- ggplot(data = arid_num_exp_df) + 
  geom_line(aes(x = (aridity_adapt_par/1e+5)*1e+6 , y = num_die_out/(sum_exp)*100, col = "Die\nout"), size=2, alpha = 0.7) + 
  geom_line(aes(x = (aridity_adapt_par/1e+5)*1e+6  , y = num_new_st/(sum_exp)*100, col = "New steady\nstate"), size=2, alpha = 0.7) + 
  geom_line(aes(x = (aridity_adapt_par/1e+5)*1e+6 , y = num_otw/(sum_exp)*100, col = "Transient"), size=2, alpha = 0.7) + 
  geom_line(aes(x = (aridity_adapt_par/1e+5)*1e+6  , y = num_return/(sum_exp)*100, col = "Return to\ninitial conditions"), size=2, alpha = 0.7) + 
  scale_color_manual(values=c("#972e1a","#074578","#BFA037","black"), 
                     name="Climate trajectory", 
                     limits = c("Die\nout", "New steady\nstate", "Transient", "Return to\ninitial conditions")) + 
  theme_classic() + 
  xlab(expression("Aridity niche adaptation rate (BAI units"~Ma^-1*")")) + 
  ylab("") + 
  scale_x_log10(labels = label_number(drop0trailing = TRUE))+
  annotation_logticks(sides="b", size = 0.2, short = unit(0.1, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
  theme(text = element_text(size = 6), 
        legend.position = "bottom", 
        legend.justification = 0.25) +
  ylim(c(0, 60))

outcome_frequency <- ggarrange(disp_per, therm_per, arid_per, labels=c("D"),hjust = -0.3, vjust = 1.5, ncol=3, nrow=1, common.legend = TRUE, legend = "top",  font.label = list(size=10))
ggsave("Trajectories_percentage_250Ma.png", plot = outcome_frequency, device = "png", path = "./",scale = 1, width = 18, height = 6, units = "cm", dpi = 400, bg = "white")


# 4) Spatial productivity plots for illustration -----------------------------------------------------------------------------------
# Take plots only from the plausible range, not the entire explored parameter space

t_out_250$plausibility <- "plausible"
index <- t_out_250$exp_ID[t_out_250$exp_type == "new_steady" & t_out_250$time == 8.0 & t_out_250$NPP < 1e+12]
t_out_250$plausibility[t_out_250$exp_ID %in% index] <- "not_plausible"
t_out_250_plausible <- t_out_250[t_out_250$plausibility == "plausible" & t_out_250$dispersal_par >= 600 & t_out_250$aridity_adapt_par > 0 & t_out_250$thermal_adapt_par > 0, ]

# Check out fast recovery example (called recovery plot)
out_min <- readRDS("../output/exp_D_out/experiment_run_12122/model_run_output_detail.rds")
out_final_rec <- rasterFromXYZ(out_min[[31]]) # Timestep 31: two million years after LIP onset
plot(out_final_rec$productivity_rate)

# Check out slow recovery example (called disrupted plot)
out_max <- readRDS("../output/exp_D_out/experiment_run_577/model_run_output_detail.rds")
out_final_dis <- rasterFromXYZ(out_max[[31]]) 
plot(out_final_dis$productivity_rate)

robinson <- CRS("+proj=robin +over")
bb <- sf::st_union(sf::st_make_grid(
  st_bbox(c(xmin = -180,
            xmax = 180,
            ymax = 90,
            ymin = -90), crs = st_crs(4326)),
  n = 100))
bb_robinson <- st_transform(bb, as.character(robinson))

# Recovery plot
prod_rec <- out_final_rec$productivity_rate
crs(prod_rec) <- "+proj=longlat +datum=WGS84 +no_defs"

prod_rec <- raster::projectRaster(prod_rec, crs = robinson)
prod_rec_df <- as.data.frame(prod_rec, xy =T)

plot_recovery <- ggplot() +
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = "#add8e6", size = 0.2) + 
  geom_raster(data = prod_rec_df, aes(x = x, y = y, fill = productivity_rate)) + 
  scale_fill_gradient2(expression("NPP ("*g~C~m^-2*")  "), low = "#E5D0B1", mid = "#729E39", high = "#084113", na.value = "transparent", limits=c(0, 1050), midpoint = 500, breaks = c(0, 500, 1000), guide = guide_colorbar(barheight = 0.8, barwidth = 4)) +
  theme_void() + 
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 6)) + 
  theme(plot.margin=unit(c(-0.5, 0, 0, 0), "cm"))  + 
  theme(plot.margin = margin(b = 0.5, unit = "cm"))


# Disrupted plot
prod_dis <- out_final_dis$productivity_rate
crs(prod_dis) <- "+proj=longlat +datum=WGS84 +no_defs"

prod_dis <- raster::projectRaster(prod_dis, crs = robinson)
prod_dis_df <- as.data.frame(prod_dis, xy =T)

plot_disrupted <- ggplot() +
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = "#add8e6", size = 0.2) + 
  geom_raster(data = prod_dis_df, aes(x = x, y = y, fill = productivity_rate)) + 
  scale_fill_gradient2(expression("NPP ("*g~C~m^-2*")  "), low = "#E5D0B1", mid = "#729E39", high = "#084113", na.value = "transparent", limits=c(0, 1050), midpoint = 500, breaks = c(0, 500, 1000), guide = guide_colorbar(barheight = 0.8, barwidth = 4)) +
  theme_void() + 
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 6)) + 
  theme(plot.margin = margin(b = 0.5, unit = "cm"))

productivity_comp <- ggarrange(plot_disrupted, plot_recovery, ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")
middle_arranged <- ggarrange(GAST_dens_250, productivity_comp, ncol = 2, widths = c(2, 1), labels = c("B", "C"), font.label = list(size=10))
comp_fig1 <- ggarrange(param_sens, middle_arranged, outcome_frequency + theme(legend.text = element_text(size = 12)), ncol = 1, nrow = 3, heights = c(0.55, 1.25, 0.75))
ggsave("Figure_1_250Ma.png", plot = comp_fig1, device = "png", path = "./",scale = 1, width = 18, height = 21, units = "cm", dpi = 400, bg = "white")


# 5) Comparison of trajectories between continental configurations ----------------------------------------------------------------------------------------------------------------------------------------------
# Calculate delta GAST 
Tinit_250 <- mean(t_out_250$GAST[t_out_250$time == 0])
Tinit_200 <- mean(t_out_200$GAST[t_out_200$time == 0]) 
Tinit_56 <- mean(t_out_56$GAST[t_out_56$time == 0]) 
Tinit_0 <- mean(t_out_0$GAST[t_out_0$time == 0]) 
t_out_250$DGAST <- t_out_250$GAST - Tinit_250
t_out_200$DGAST <- t_out_200$GAST - Tinit_200
t_out_56$DGAST <- t_out_56$GAST - Tinit_56
t_out_0$DGAST <- t_out_0$GAST - Tinit_0

t_out_250$plausibility <- "plausible"
index <- t_out_250$exp_ID[(t_out_250$exp_type == "new_steady" | t_out_250$exp_type == "die_out") & t_out_250$time == 8.0 & t_out_250$NPP < 1e+12]
t_out_250$plausibility[t_out_250$exp_ID %in% index] <- "not_plausible"
t_out_250_plausible <- t_out_250[t_out_250$plausibility == "plausible" & t_out_250$dispersal_par >= 600 & t_out_250$aridity_adapt_par > 0 & t_out_250$aridity_adapt_par > 0, ]

DGAST_dens_250 <- ggplot(subset(t_out_250_plausible, preplant_par == 0.25), aes(time-1, DGAST, group = exp_ID)) + 
  geom_line(alpha=0.05, linewidth = 0.4, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("") +
  ylab(expression(Delta*"GAST"~"(°C)")) +
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top") + 
  ylim(c(-0.25, 10))

t_out_200$plausibility <- "plausible"
index <- t_out_200$exp_ID[(t_out_200$exp_type == "new_steady" | t_out_200$exp_type == "die_out") & t_out_200$time == 8.0 & t_out_200$NPP < 1e+12]
t_out_200$plausibility[t_out_200$exp_ID %in% index] <- "not_plausible"
t_out_200_plausible <- t_out_200[t_out_200$plausibility == "plausible" & t_out_200$dispersal_par >= 600 & t_out_200$aridity_adapt_par > 0 & t_out_200$thermal_adapt_par > 0, ]

DGAST_dens_200 <- ggplot(subset(t_out_200_plausible, preplant_par == 0.25), aes(time-1, DGAST, group = exp_ID)) + 
  geom_line(alpha=0.05, linewidth = 0.4, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("") +
  ylab(expression(Delta*"GAST"~"(°C)")) +
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top") + 
  ylim(c(-0.25,10))


t_out_56$plausibility <- "plausible"
index <- t_out_56$exp_ID[(t_out_56$exp_type == "new_steady" | t_out_56$exp_type == "die_out") & t_out_56$time == 8.0 & t_out_56$NPP < 1e+12]
t_out_56$plausibility[t_out_56$exp_ID %in% index] <- "not_plausible"
t_out_56_plausible <- t_out_56[t_out_56$plausibility == "plausible"  & t_out_56$dispersal_par >= 600 & t_out_56$thermal_adapt_par > 0 & t_out_56$aridity_adapt_par > 0, ]

DGAST_dens_56 <- ggplot(subset(t_out_56_plausible, preplant_par == 0.25), aes(time-1, DGAST, group = exp_ID)) + 
  geom_line(alpha=0.05, linewidth = 0.4, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Return")) +
  xlab("") +
  ylab(expression(Delta*"GAST"~"(°C)")) +
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top") + 
  ylim(c(-0.3, 10))


t_out_0$plausibility <- "plausible"
index <- t_out_0$exp_ID[(t_out_0$exp_type == "new_steady" | t_out_0$exp_type == "die_out") & t_out_0$time == 8.0 & t_out_0$NPP < 1e+12]
t_out_0$plausibility[t_out_0$exp_ID %in% index] <- "not_plausible"
t_out_0_plausible <- t_out_0[t_out_0$plausibility == "plausible" & t_out_0$dispersal_par >= 600 & t_out_0$aridity_adapt_par > 0 & t_out_0$thermal_adapt_par > 0, ]

DGAST_dens_0 <- ggplot(subset(t_out_0_plausible, preplant_par == 0.25), aes(time-1, DGAST, group = exp_ID)) + 
  geom_line(alpha=0.05, linewidth = 0.4, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("Years relative to LIP onset (Ma)") +
  ylab(expression(Delta*"GAST"~"(°C)")) +
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top")

plausible_trajectories <- ggarrange(DGAST_dens_250, DGAST_dens_200, DGAST_dens_56, DGAST_dens_0, 
                                    labels = c("A", "B", "C", "D"), ncol = 1, nrow = 4, font.label = list(size=10), 
                                    common.legend = TRUE, legend = "bottom")

ggsave("plausible_trajectories_combined.png", plot = plausible_trajectories, device = "png", path = "./",
       width = 8, height = 22.7, units = "cm", dpi = 400, bg = "white")


# Combination with temperature plots in inkscape - separate analysis and plotting directory for proxies in Matlab


# 6) All carbon fluxes for 250 Ma -------------------------------------------------------------
index <- t_out_250$exp_ID[(t_out_250$exp_type =="new_steady" | t_out_250$exp_type == "die_out") & t_out_250$time == 8.0 & t_out_250$locb < 1e+12]
t_out_250$exp_type[t_out_250$exp_ID %in% index] <- "die_out"

CO2_dens <- ggplot(subset(t_out_250, preplant_par == 0.25), aes(time-1, CO2, group = exp_ID)) + 
  geom_line(alpha=0.025, linewidth = 0.2, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#972e1a","#39259E","#BFA037","black"), labels = c("Die out", "New steady\n state", "Transient", "Return")) +
  xlab("") +
  ylab(expression(CO[2]~"(ppm)"))  +
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "bottom")


silw_dens <- ggplot(subset(t_out_250, preplant_par == 0.25), aes(time-1, silw, group = exp_ID)) + 
  geom_line(alpha=0.025, linewidth = 0.2, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#972e1a","#39259E","#BFA037","black"), labels = c("Die out", "New steady\n state", "Transient", "Return")) +
  xlab("") +
  ylab(expression(F[silw]~"("*mol~C~yr^-1*")"))+
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "bottom")

sfw_dens <- ggplot(subset(t_out_250, preplant_par == 0.25), aes(time-1, sfw, group = exp_ID)) + 
  geom_line(alpha=0.025, linewidth = 0.2, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#972e1a","#39259E","#BFA037","black"), labels = c("Die out", "New steady\n state", "Transient", "Return")) +
  xlab("") +
  ylab(expression(F[sfw]~"("*mol~C~yr^-1*")"))+
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "bottom")


locb_dens <- ggplot(subset(t_out_250, preplant_par == 0.25), aes(time-1, locb, group = exp_ID)) + 
  geom_line(alpha=0.025, linewidth = 0.2, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#972e1a","#39259E","#BFA037","black"), labels = c("Die out", "New steady\n state", "Transient", "Return")) +
  xlab("Years relative to LIP onset (Ma)") +
  ylab(expression(F[locb]~"("*mol~C~yr^-1*")")) +
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "bottom")

mocb_dens <- ggplot(subset(t_out_250, preplant_par == 0.25), aes(time-1, mocb, group = exp_ID)) + 
  geom_line(alpha=0.025, linewidth = 0.2, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#972e1a","#39259E","#BFA037","black"), labels = c("Die out", "New steady\n state", "Transient", "Return")) +
  xlab("Years relative to LIP onset (Ma)") +
  ylab(expression(F[mocb]~"("*mol~C~yr^-1*")"))+
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "bottom")

degass_dens <- ggplot(subset(t_out_250, preplant_par == 0.25), aes(time-1, degass, group = exp_ID)) + 
  geom_line(alpha=0.025, linewidth = 0.2, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#972e1a","#39259E","#BFA037","black"), labels = c("Die out", "New steady\n state", "Transient", "Return")) +
  xlab("Model years (Myr)") +
  ylab(expression(F[degassing]~"("*mol~C~yr^-1*")"))+
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "bottom")

all_fluxes <- ggarrange(CO2_dens, silw_dens, locb_dens, mocb_dens, degass_dens, nrow = 3, ncol = 2, labels = c("A", "B", "C", "D", "E"), font.label = list(size=10), common.legend = TRUE, legend = "bottom")
ggsave("all_C_fluxes_PT_250Ma.png", plot = all_fluxes, device = "png", path = "./", scale = 1, width = 16, height = 18, units = "cm", dpi = 400, bg = "white")

# 7) Carbonate proxy approximation  -----------------------------------------------------------------------------
# Proxy calculation 
t_out_250$d13carb <- -6 + ((t_out_250$locb + t_out_250$mocb)/(t_out_250$locb + t_out_250$mocb + t_out_250$silw + t_out_250$sfw + t_out_250$carbw)) * 27
t_out_200$d13carb <- -6 + ((t_out_200$locb + t_out_200$mocb)/(t_out_200$locb + t_out_200$mocb + t_out_200$silw + t_out_200$sfw + t_out_200$carbw)) * 27
t_out_56$d13carb <- -6 + ((t_out_56$locb + t_out_56$mocb)/(t_out_56$locb + t_out_56$mocb + t_out_56$silw + t_out_56$sfw + t_out_56$carbw)) * 27
Carbonate_init_250 <- mean(t_out_250$d13carb[t_out_250$time == 0])
Carbonate_init_200 <- mean(t_out_200$d13carb[t_out_200$time == 0]) 
Carbonate_init_56 <- mean(t_out_56$d13carb[t_out_56$time == 0]) 
t_out_250$delta_d13carb <- t_out_250$d13carb - Carbonate_init_250
t_out_200$delta_d13carb <- t_out_200$d13carb - Carbonate_init_200
t_out_56$delta_d13carb <- t_out_56$d13carb - Carbonate_init_56


# Shift expID numbering for plotting
t_out_200$exp_ID <- t_out_200$exp_ID + 2*12122
t_out_56$exp_ID <- t_out_56$exp_ID + 4*12122


# Categorization 
t_out_complete <- rbind(t_out_250, t_out_200, t_out_56)
t_out_complete$plausibility <- "plausible"
index <- t_out_complete$exp_ID[t_out_complete$exp_type == "new_steady" & t_out_complete$time == 6.0 & t_out_complete$NPP < 1e+12]
t_out_complete$plausibility[t_out_complete$exp_ID %in% index] <- "not_plausible"
t_out_complete_plausible <- t_out_complete[t_out_complete$plausibility == "plausible" & t_out_complete$thermal_adapt_par > 0 & t_out_complete$aridity_adapt_par > 0 & t_out_complete$dispersal_par >= 600, ]

d13Ccarb_dens <- ggplot(subset(t_out_complete_plausible, year %in% c(250, 56)) , aes(time-1, delta_d13carb, group = exp_ID)) + 
  geom_line(alpha=0.03, size = 0.3, aes(color = as.factor(year)), key_glyph = "rect") + 
  theme_classic() +
  scale_color_manual(values=c("#FD9A52","#812B92")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  xlab("Years relative to LIP onset (Ma)") +
  ylab(expression(delta^13*C ~ "(\u2030)")) +
  theme(text = element_text(size = 8, family = "sans"), 
        panel.grid.minor = element_blank()) + 
  theme(legend.position = "none")
ggsave("model_d13C_approximation.png", plot = d13Ccarb_dens, device = "png", path = "./", scale = 1, width = 8, height = 5, units = "cm", dpi = 400, bg = "white")



# 8) Alternative degassing boundary conditions --------------------------------------------------------------------

t_out_200_S1$plausibility <- "plausible"
index <- t_out_200_S1$exp_ID[(t_out_200_S1$exp_type == "new_steady" | t_out_200_S1$exp_type == "die_out") & t_out_200_S1$time == 8.0 & t_out_200_S1$NPP < 1e+12]
t_out_200_S1$plausibility[t_out_200_S1$exp_ID %in% index] <- "not_plausible"
t_out_200_S1_plausible <- t_out_200_S1[t_out_200_S1$plausibility == "plausible" & t_out_200_S1$dispersal_par >= 600 & t_out_200_S1$aridity_adapt_par > 0 & t_out_200_S1$aridity_adapt_par > 0, ]

GAST_dens_200_S1 <- ggplot(subset(t_out_200_S1_plausible), aes(time, GAST, group = exp_ID)) + 
  geom_line(alpha=0.03, linewidth = 0.3, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("Years relative to LIP onset (Ma)") +
  ylab(expression("GAST (°C)")) +
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top")


t_out_250_S2$plausibility <- "plausible"
index <- t_out_250_S2$exp_ID[(t_out_250_S2$exp_type == "new_steady" | t_out_250_S2$exp_type == "die_out") & t_out_250_S2$time == 8.0 & t_out_250_S2$NPP < 1e+12]
t_out_250_S2$plausibility[t_out_250_S2$exp_ID %in% index] <- "not_plausible"
t_out_250_S2_plausible <- t_out_250_S2[t_out_250_S2$plausibility == "plausible" & t_out_250_S2$dispersal_par >= 600 & t_out_250_S2$aridity_adapt_par > 0 & t_out_250_S2$aridity_adapt_par > 0, ]

GAST_dens_250_S2 <- ggplot(subset(t_out_250_S2_plausible), aes(time, GAST, group = exp_ID)) + 
  geom_line(alpha=0.03, linewidth = 0.3, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("Years relative to LIP onset (Ma)") +
  ylab(expression("GAST (°C)")) +
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top")



t_out_250_S3$plausibility <- "plausible"
index <- t_out_250_S3$exp_ID[(t_out_250_S3$exp_type == "new_steady" | t_out_250_S3$exp_type == "die_out") & t_out_250_S3$time == 8.0 & t_out_250_S3$NPP < 1e+12]
t_out_250_S3$plausibility[t_out_250_S3$exp_ID %in% index] <- "not_plausible"
t_out_250_S3_plausible <- t_out_250_S3[t_out_250_S3$plausibility == "plausible" & t_out_250_S3$dispersal_par >= 600 & t_out_250_S3$aridity_adapt_par > 0 & t_out_250_S3$aridity_adapt_par > 0, ]

GAST_dens_250_S3 <- ggplot(subset(t_out_250_S3_plausible), aes(time, GAST, group = exp_ID)) + 
  geom_line(alpha=0.03, linewidth = 0.3, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("Years relative to LIP onset (Ma)") +
  ylab(expression("GAST (°C)")) +
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top")


t_out_56_S4$plausibility <- "plausible"
index <- t_out_56_S4$exp_ID[(t_out_56_S4$exp_type == "new_steady" | t_out_56_S4$exp_type == "die_out") & t_out_56_S4$time == 8.0 & t_out_56_S4$NPP < 1e+12]
t_out_56_S4$plausibility[t_out_56_S4$exp_ID %in% index] <- "not_plausible"
t_out_56_S4_plausible <- t_out_56_S4[t_out_56_S4$plausibility == "plausible" & t_out_56_S4$dispersal_par >= 600 & t_out_56_S4$aridity_adapt_par > 0 & t_out_56_S4$aridity_adapt_par > 0, ]

GAST_dens_56_S4 <- ggplot(subset(t_out_56_S4_plausible), aes(time, GAST, group = exp_ID)) + 
  geom_line(alpha=0.03, linewidth = 0.3, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("Years relative to LIP onset (Ma)") +
  ylab(expression("GAST (°C)")) +
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top")

supplementary_bc <- ggarrange(GAST_dens_250_S2, GAST_dens_200_S1, GAST_dens_56_S4, GAST_dens_250_S3, nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"), font.label = list(size=10), common.legend = TRUE, legend = "bottom")
ggsave("supp_fig_alternative_BC.png", plot = supplementary_bc, device = "png", path = "./", scale = 1, width = 16, height = 14, units = "cm", dpi = 400, bg = "white")


# 9) Modelled present-day productivity  ------------------------------------------
out_present <- readRDS("../output/exp_A_out/experiment_run_12122/model_run_output_detail.rds")
present_data <- rasterFromXYZ(out_present[[1]])
productivity_present <- present_data$productivity_rate

robinson <- CRS("+proj=robin +over")
bb <- sf::st_union(sf::st_make_grid(
  st_bbox(c(xmin = -180,
            xmax = 180,
            ymax = 90,
            ymin = -90), crs = st_crs(4326)),
  n = 100))
bb_robinson <- st_transform(bb, as.character(robinson))

crs(productivity_present) <- "+proj=longlat +datum=WGS84 +no_defs"
productivity_present <- raster::projectRaster(productivity_present, crs = robinson)
productivity_present_df <- as.data.frame(productivity_present, xy =T)

plot_present_model <- ggplot() +
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = "#add8e6", size = 0.2) + 
  geom_raster(data = productivity_present_df, aes(x = x, y = y, fill = productivity_rate)) + 
  scale_fill_gradient2(expression("NPP ("*g~C~m^-2*")  "), low = "#E5D0B1", mid = "#729E39", high = "#084113", na.value = "transparent", limits=c(0, 1300), midpoint = 650, breaks = c(0, 500, 1000), guide = guide_colorbar(barheight = 0.8, barwidth = 4)) +
  theme_void() + 
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 6)) + 
  theme(plot.margin=unit(c(-0.5, 0, 0, 0), "cm"))  + 
  theme(plot.margin = margin(b = 0.5, unit = "cm"))
plot_present_model

ggsave("NPP_present.png", plot = NPP_comp, device = "png", path = "./", scale = 1, width = 8, height = 7, units = "cm", dpi = 300, bg = "white")



# 10) Comparison silicate weathering and primary productivity for all paleogeographies ------------------------------------------
present <- readRDS("../output/exp_A_out/experiment_run_12122/model_run_output_detail.rds")
present <- rasterFromXYZ(present[[1]])

PE <- readRDS("../output/exp_B_out/experiment_run_12122/model_run_output_detail.rds")
PE <- rasterFromXYZ(PE[[1]])

PT <- readRDS("../output/exp_C_out/experiment_run_12122/model_run_output_detail.rds")
PT <- rasterFromXYZ(PT[[1]])

TJ <- readRDS("../output/exp_D_out/experiment_run_12122/model_run_output_detail.rds")
TJ <- rasterFromXYZ(TJ[[1]])


par(mfrow = c(2, 2))
plot(present$productivity_rate, zlim=c(0, 1200))
plot(PE$productivity_rate, zlim=c(0, 1200))
plot(PT$productivity_rate, zlim=c(0, 1200))
plot(TJ$productivity_rate, zlim=c(0, 1200))


silw_conversion <- 1e+13 / 339296153 # calibraftion factor from degassing_main.R file
plot(present$silw * silw_conversion, zlim=c(0, 1))
plot(PE$silw * silw_conversion, zlim=c(0, 1))
plot(PT$silw * silw_conversion, zlim=c(0, 1))
plot(TJ$silw * silw_conversion, zlim=c(0, 1))


robinson <- CRS("+proj=robin +over")
bb <- sf::st_union(sf::st_make_grid(
  st_bbox(c(xmin = -180,
            xmax = 180,
            ymax = 90,
            ymin = -90), crs = st_crs(4326)),
  n = 100))
bb_robinson <- st_transform(bb, as.character(robinson))

crs(present) <- "+proj=longlat +datum=WGS84 +no_defs"
crs(PE) <- "+proj=longlat +datum=WGS84 +no_defs"
crs(PT) <- "+proj=longlat +datum=WGS84 +no_defs"
crs(TJ) <- "+proj=longlat +datum=WGS84 +no_defs"


present <- raster::projectRaster(present, crs = robinson)
present_df <- as.data.frame(present, xy =T)
PE <- raster::projectRaster(PE, crs = robinson)
PE_df <- as.data.frame(PE, xy =T)
PT <- raster::projectRaster(PT, crs = robinson)
PT_df <- as.data.frame(PT, xy =T)
TJ <- raster::projectRaster(TJ, crs = robinson)
TJ_df <- as.data.frame(TJ, xy =T)


productivity_present <- ggplot() +
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = "#add8e6", size = 0.2) + 
  geom_raster(data = present_df, aes(x = x, y = y, fill = productivity_rate)) + 
  scale_fill_gradient2(expression("NPP ("*g~C~m^-2*")  "), low = "#E5D0B1", mid = "#729E39", high = "#084113", na.value = "transparent", limits=c(0, 1200), midpoint = 600, breaks = c(0, 500, 1000), guide = guide_colorbar(barheight = 0.8, barwidth = 4)) +
  theme_void() + 
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 6)) + 
  theme(plot.margin=unit(c(-0.5, 0, 0, 0), "cm"))  + 
  theme(plot.margin = margin(b = 0.5, unit = "cm"))

productivity_PE <- ggplot() +
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = "#add8e6", size = 0.2) + 
  geom_raster(data = PE_df, aes(x = x, y = y, fill = productivity_rate)) + 
  scale_fill_gradient2(expression("NPP ("*g~C~m^-2*")  "), low = "#E5D0B1", mid = "#729E39", high = "#084113", na.value = "transparent", limits=c(0, 1200), midpoint = 600, breaks = c(0, 500, 1000), guide = guide_colorbar(barheight = 0.8, barwidth = 4)) +
  theme_void() + 
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 6)) + 
  theme(plot.margin=unit(c(-0.5, 0, 0, 0), "cm"))  + 
  theme(plot.margin = margin(b = 0.5, unit = "cm"))

productivity_TJ <- ggplot() +
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = "#add8e6", size = 0.2) + 
  geom_raster(data = TJ_df, aes(x = x, y = y, fill = productivity_rate)) + 
  scale_fill_gradient2(expression("NPP ("*g~C~m^-2*")  "), low = "#E5D0B1", mid = "#729E39", high = "#084113", na.value = "transparent", limits=c(0, 1200), midpoint = 600, breaks = c(0, 500, 1000), guide = guide_colorbar(barheight = 0.8, barwidth = 4)) +
  theme_void() + 
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 6)) + 
  theme(plot.margin=unit(c(-0.5, 0, 0, 0), "cm"))  + 
  theme(plot.margin = margin(b = 0.5, unit = "cm"))

productivity_PT <- ggplot() +
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = "#add8e6", size = 0.2) + 
  geom_raster(data = PT_df, aes(x = x, y = y, fill = productivity_rate)) + 
  scale_fill_gradient2(expression("NPP ("*g~C~m^-2*")  "), low = "#E5D0B1", mid = "#729E39", high = "#084113", na.value = "transparent", limits=c(0, 1200), midpoint = 600, breaks = c(0, 500, 1000), guide = guide_colorbar(barheight = 0.8, barwidth = 4)) +
  theme_void() + 
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 6)) + 
  theme(plot.margin=unit(c(-0.5, 0, 0, 0), "cm"))  + 
  theme(plot.margin = margin(b = 0.5, unit = "cm"))

left <- ggarrange(productivity_present, productivity_PE, productivity_TJ, productivity_PT, ncol = 1, nrow = 4, 
                  labels = c("A", "C", "E", "G"), font.label = list(size = 10), common.legend = T, legend = "bottom")


silw_present <- ggplot() +
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = "#add8e6", size = 0.2) + 
  geom_raster(data = present_df, aes(x = x, y = y, fill = silw * silw_conversion * 12)) + 
  # scale_fill_viridis(expression("silw ("*g~C~m^-2*")  "), option = "B",oob = squish, na.value = "transparent", guide = guide_colorbar(barheight = 0.8, barwidth = 4)) +
  binned_scale(expression("silw ("*g~C~m^-2*")  "),
               aesthetics = "fill",
               scale_name = "stepsn", 
               palette = function(x) c("#000004FF", "#330A5FFF", "#781C6DFF", "#BB3754FF", "#ED6925FF", "#FCB519FF", "#FCFFA4FF"),
               breaks = c(0.01, 0.05, 0.1, 0.5, 1, 5, 10),
               limits = c(0, 10),
               show.limits = TRUE, 
               guide = "colorsteps"
  ) +  
  theme_void() + 
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 6)) + 
  theme(plot.margin=unit(c(-0.5, 0, 0, 0), "cm"))  + 
  theme(plot.margin = margin(b = 0.5, unit = "cm"))


silw_PE <- ggplot() +
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = "#add8e6", size = 0.2) + 
  geom_raster(data = PE_df, aes(x = x, y = y, fill = silw * silw_conversion * 12)) + 
  # scale_fill_viridis(expression("silw ("*g~C~m^-2*")  "), option = "B",oob = squish, na.value = "transparent", guide = guide_colorbar(barheight = 0.8, barwidth = 4)) +
  binned_scale(expression("silw ("*g~C~m^-2*")  "),
               aesthetics = "fill",
               scale_name = "stepsn", 
               palette = function(x) c("#000004FF", "#330A5FFF", "#781C6DFF", "#BB3754FF", "#ED6925FF", "#FCB519FF", "#FCFFA4FF"),
               breaks = c(0.01, 0.05, 0.1, 0.5, 1, 5, 10),
               limits = c(0, 10),
               show.limits = TRUE, 
               guide = "colorsteps"
  ) + 
  theme_void() + 
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 6)) + 
  theme(plot.margin=unit(c(-0.5, 0, 0, 0), "cm"))  + 
  theme(plot.margin = margin(b = 0.5, unit = "cm"))

silw_TJ <- ggplot() +
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = "#add8e6", size = 0.2) + 
  geom_raster(data = TJ_df, aes(x = x, y = y, fill = silw * silw_conversion * 12)) + 
  # scale_fill_viridis(expression("silw ("*g~C~m^-2*")  "), option = "B",oob = squish, na.value = "transparent", guide = guide_colorbar(barheight = 0.8, barwidth = 4)) +
  binned_scale(expression("silw ("*g~C~m^-2*")  "),
               aesthetics = "fill",
               scale_name = "stepsn", 
               palette = function(x) c("#000004FF", "#330A5FFF", "#781C6DFF", "#BB3754FF", "#ED6925FF", "#FCB519FF", "#FCFFA4FF"),
               breaks = c(0.01, 0.05, 0.1, 0.5, 1, 5, 10),
               limits = c(0, 10),
               show.limits = TRUE, 
               guide = "colorsteps"
  ) + 
  theme_void() + 
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 6)) + 
  theme(plot.margin=unit(c(-0.5, 0, 0, 0), "cm"))  + 
  theme(plot.margin = margin(b = 0.5, unit = "cm"))

silw_PT <- ggplot() +
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = "#add8e6", size = 0.2) + 
  geom_raster(data = PT_df, aes(x = x, y = y, fill = silw * silw_conversion * 12)) + 
  # scale_fill_viridis(expression("silw ("*g~C~m^-2*")  "), option = "B",oob = squish, na.value = "transparent", guide = guide_colorbar(barheight = 0.8, barwidth = 4)) +
  binned_scale(expression("silw ("*g~C~m^-2*")  "),
               aesthetics = "fill",
               scale_name = "stepsn", 
               palette = function(x) c("#000004FF", "#330A5FFF", "#781C6DFF", "#BB3754FF", "#ED6925FF", "#FCB519FF", "#FCFFA4FF"),
               breaks = c(0.01, 0.05, 0.1, 0.5, 1, 5, 10),
               limits = c(0, 10),
               show.limits = TRUE, 
               guide = "colorsteps"
  ) + 
  theme_void() + 
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 6)) + 
  theme(plot.margin=unit(c(-0.5, 0, 0, 0), "cm"))  + 
  theme(plot.margin = margin(b = 0.5, unit = "cm"))

right <- ggarrange(silw_present, silw_PE, silw_TJ, silw_PT, ncol = 1, nrow = 4, 
                   labels = c("B", "D", "F", "H"), font.label = list(size = 10), common.legend = T, legend = "bottom")


arranged <- ggarrange(left, right, ncol = 2)
ggsave("spatial_silw_locb.png", plot = arranged, device = "png", path = "./", scale = 1, width = 16, height = 18.5, units = "cm", dpi = 400, bg = "white")


# 11) temporal silw, locb and dA for all relevant scenarios  ------------------------------------------

t_out_250$plausibility <- "plausible"
index <- t_out_250$exp_ID[(t_out_250$exp_type == "new_steady" | t_out_250$exp_type == "die_out") & t_out_250$time == 8.0 & t_out_250$NPP < 1e+12]
t_out_250$plausibility[t_out_250$exp_ID %in% index] <- "not_plausible"
t_out_250_plausible <- t_out_250[t_out_250$plausibility == "plausible" & t_out_250$dispersal_par >= 600 & t_out_250$aridity_adapt_par > 0 & t_out_250$aridity_adapt_par > 0, ]

silw_dens_250 <- ggplot(subset(t_out_250_plausible, preplant_par == 0.25), aes(time-1, silw, group = exp_ID)) + 
  geom_line(alpha=0.05, linewidth = 0.4, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("Time relative to LIP onset (Ma)") +
  ylab(expression(F[silw]~"("*mol~C~yr^-1*")"))+
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top")

locb_dens_250 <- ggplot(subset(t_out_250_plausible, preplant_par == 0.25), aes(time-1, locb, group = exp_ID)) + 
  geom_line(alpha=0.05, linewidth = 0.4, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("Time relative to LIP onset (Ma)") +
  ylab(expression(F[locb]~"("*mol~C~yr^-1*")"))+
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top")

dA_dens_250 <- ggplot(subset(t_out_250_plausible, preplant_par == 0.25), aes(time-1, degass - silw - locb - mocb, group = exp_ID)) + 
  geom_line(alpha=0.05, linewidth = 0.4, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("Time relative to LIP onset (Ma)") +
  ylab(expression(dAO~"("*mol~C~yr^-1*")"))+
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top")




t_out_200$plausibility <- "plausible"
index <- t_out_200$exp_ID[(t_out_200$exp_type == "new_steady" | t_out_200$exp_type == "die_out") & t_out_200$time == 8.0 & t_out_200$NPP < 1e+12]
t_out_200$plausibility[t_out_200$exp_ID %in% index] <- "not_plausible"
t_out_200_plausible <- t_out_200[t_out_200$plausibility == "plausible" & t_out_200$dispersal_par >= 600 & t_out_200$aridity_adapt_par > 0 & t_out_200$aridity_adapt_par > 0, ]

silw_dens_200 <- ggplot(subset(t_out_200_plausible, preplant_par == 0.25), aes(time-1, silw, group = exp_ID)) + 
  geom_line(alpha=0.05, linewidth = 0.4, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("") +
  ylab(expression(F[silw]~"("*mol~C~yr^-1*")"))+
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top")

locb_dens_200 <- ggplot(subset(t_out_200_plausible, preplant_par == 0.25), aes(time-1, locb, group = exp_ID)) + 
  geom_line(alpha=0.05, linewidth = 0.4, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("") +
  ylab(expression(F[locb]~"("*mol~C~yr^-1*")"))+
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top")

dA_dens_200 <- ggplot(subset(t_out_200_plausible, preplant_par == 0.25), aes(time-1, degass - silw - locb - mocb, group = exp_ID)) + 
  geom_line(alpha=0.05, linewidth = 0.4, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("Time relative to LIP onset (Ma)") +
  ylab(expression(dAO~"("*mol~C~yr^-1*")"))+
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top")



t_out_56$plausibility <- "plausible"
index <- t_out_56$exp_ID[(t_out_56$exp_type == "new_steady" | t_out_56$exp_type == "die_out") & t_out_56$time == 8.0 & t_out_56$NPP < 1e+12]
t_out_56$plausibility[t_out_56$exp_ID %in% index] <- "not_plausible"
t_out_56_plausible <- t_out_56[t_out_56$plausibility == "plausible" & t_out_56$dispersal_par >= 600 & t_out_56$aridity_adapt_par > 0 & t_out_56$aridity_adapt_par > 0, ]

silw_dens_56 <- ggplot(subset(t_out_56_plausible, preplant_par == 0.25), aes(time-1, silw, group = exp_ID)) + 
  geom_line(alpha=0.05, linewidth = 0.4, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("") +
  ylab(expression(F[silw]~"("*mol~C~yr^-1*")"))+
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top")

locb_dens_56 <- ggplot(subset(t_out_56_plausible, preplant_par == 0.25), aes(time-1, locb, group = exp_ID)) + 
  geom_line(alpha=0.05, linewidth = 0.4, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("") +
  ylab(expression(F[locb]~"("*mol~C~yr^-1*")"))+
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top")

dA_dens_56 <- ggplot(subset(t_out_56_plausible, preplant_par == 0.25), aes(time-1, degass - silw - locb - mocb, group = exp_ID)) + 
  geom_line(alpha=0.05, linewidth = 0.4, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("Time relative to LIP onset (Ma)") +
  ylab(expression(dAO~"("*mol~C~yr^-1*")"))+
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top")



t_out_0$plausibility <- "plausible"
index <- t_out_0$exp_ID[(t_out_0$exp_type == "new_steady" | t_out_0$exp_type == "die_out") & t_out_0$time == 8.0 & t_out_0$NPP < 1e+12]
t_out_0$plausibility[t_out_0$exp_ID %in% index] <- "not_plausible"
t_out_0_plausible <- t_out_0[t_out_0$plausibility == "plausible" & t_out_0$dispersal_par >= 600 & t_out_0$aridity_adapt_par > 0 & t_out_0$aridity_adapt_par > 0, ]

silw_dens_0 <- ggplot(subset(t_out_0_plausible, preplant_par == 0.25), aes(time-1, silw, group = exp_ID)) + 
  geom_line(alpha=0.05, linewidth = 0.4, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("") +
  ylab(expression(F[silw]~"("*mol~C~yr^-1*")"))+
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top")

locb_dens_0 <- ggplot(subset(t_out_0_plausible, preplant_par == 0.25), aes(time-1, locb, group = exp_ID)) + 
  geom_line(alpha=0.05, linewidth = 0.4, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("") +
  ylab(expression(F[locb]~"("*mol~C~yr^-1*")"))+
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top")

dA_dens_0 <- ggplot(subset(t_out_0_plausible, preplant_par == 0.25), aes(time-1, degass - silw - locb - mocb, group = exp_ID)) + 
  geom_line(alpha=0.05, linewidth = 0.4, aes(color = as.factor(exp_type)), key_glyph = "rect") + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  scale_color_manual("Climate trajectory", values=c("#074578","#BFA037","black"), labels = c("New steady\nstate", "Transient", "Return")) +
  xlab("Time relative to LIP onset (Ma)") +
  ylab(expression(dAO~"("*mol~C~yr^-1*")"))+
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "top")


silw_locb_plausible <- ggarrange(silw_dens_0, locb_dens_0,
                                 silw_dens_56, locb_dens_56,
                                 silw_dens_200, locb_dens_200, 
                                 silw_dens_250, locb_dens_250, 
                                 ncol = 2, nrow = 4, common.legend = T, legend = "bottom",
                                 labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), font.label = list(size = 10))

ggsave("silw_locb_figerprints.png", plot = silw_locb_plausible, device = "png", path = "./", scale = 1, width = 18, height = 18.5, units = "cm", dpi = 400, bg = "white")



