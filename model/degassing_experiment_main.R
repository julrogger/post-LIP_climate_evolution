########################################################################################################################################################
# Main function file for degassing model experiments
# 
# Author: Julian Rogger
# Last update: 30.05.24
########################################################################################################################################################

# Vegetation parameters and degassing boundary conditions are defined in frontend file: 
#
# thermal_adapt_par: rate at which floras can adapt temperature niche, in °C per year
# aridity_adapt_par rate at which floras can adapt aridity niche, in Budyko aridity units per time step
# dispersal_par: rate at which floras can disperse in space, in km per time step 
# PREPLANT, defines maximum weathering enhancement factor, either (1/4) or (1/6) for 4-fold and 6-fold weathering enhancement respectively
#
# CO2_ppm_start: This is the starting CO2 concentration for the model experiment
# t: this is the geologic time period and paleogeography considered, current options are: -250e+6 for the Permian-Triassic, -200e+6 for the Triassic Jurassic, -55e+6 for the Paleocene-Eocene, 0 for the present day configuration
# degass_mass: this is the amount of C degassing during the LIP Period, this should be given in grams
# degass_duration: this is the number of time steps over which degassing occurs
# seafloor_weathering_feedback: should seafloor weathering be included as an additional carbon cycle feedback? Parametrization is associated with uncertainty and can only be considered as a global flux (in contrast to all other fluxes that are sptially explicit)

degassing_experiment_main_function <- function(thermal_adapt_par, aridity_adapt_par, dispersal_par, PREPLANT, CO2_ppm_start, t, degass_mass, degass_duration, seafloor_weathering_feedback){

# Create empty directories for model output
dir.create("./biology/model_output", showWarnings = F)
dir.create("./biology/model_landscape", showWarnings = F)
dir.create("./biology/model_output/biology_configuration_file", showWarnings = F)
dir.create("./biology/model_output/biology_configuration_file/species", showWarnings = F)
dir.create("./biology/model_output/biology_configuration_file/landscapes", showWarnings = F)
dir.create("./biology/model_output/biology_configuration_file/world_trait_map", showWarnings = F)


# Locate topographic maps - should be located in main folder
topo_maps_prefix <- '../paleogeography/'

# Locate climate maps - should be located in main folder
climate_prefix <- "../climate_data/"

# Initial conditions 
CO2_ppm <- CO2_ppm_start # Starting CO2 concentration defined in frontend
# Atmospheric CO2 follows relation from Kump & Arthur 1999, CO2 = CO2,0 * (A_t / A_0)^2 (all references in paper)
# Where A is the carbon content of the atmosphere-ocean carbon pool
# References for A_0 at a CO2 concentration of 280 ppm from Lenton et al. 2018
A0 <- sqrt(CO2_ppm_start/280) * (3.193 * 10^18) # Assess atmosphere-ocean carbon content at starting CO2
A <- A0 # Initialise A_t
RCO2 <- 1
t_geol <- t/1e+6 # time in Ma

# Climate LUT information 
CO2 <- c(100, seq(250, 1000, 50), seq(1125, 3000, 125), seq(3500, 8000, 500)) # these are the CO2 concentrations for which climate model runs are available that can be used for climate interpolation

# Initialize time stepping 
model_period <- 8 # model years in Ma
total_time_steps <- model_period*10 # time steps per Myr 
DT <- (model_period * 1e+6)/(model_period * 10) # 100'000 years timesteps
timestep <- 0 


# Initialize output structure
out <- data.frame("time" = rep(NA, total_time_steps), # time step 
                 "CO2" = rep(NA, total_time_steps), # Atmospheric CO2 in ppm 
                 "locb" = rep(NA, total_time_steps), # land derived organic carbon burial
                 "mocb" = rep(NA, total_time_steps), # marine derived organic carbon burial
                 "silw" = rep(NA, total_time_steps), # silicate weathering carbon burial
                 "carbw" = rep(NA, total_time_steps), # carbonate weathering carbon burial
                 "sfw" = rep(NA, total_time_steps), # seafloor weathering (optional)
                 "degass" = rep(NA, total_time_steps), # degassing carbon flux
                 "GAST" = rep(NA, total_time_steps)) # Global average surface temperature
output_collection <- list() # for collection of all spatial rasters


for (timestep in c(0:(total_time_steps))){
print(paste("Continental configuration of", t, "Ma", sep=" "))

########################################################################################################
# 1.) Climate look up 
#######################################################################################################

  # Find lower and upper frames of time and CO2 
if(CO2_ppm >= 8000){
  lower_CO2 <- 8000
  upper_CO2 <- 8000
  CO2_contribution_lower <- 1
  CO2_contribution_upper <- 0
} else {
  lower_CO2 <- CO2[which(CO2_ppm - CO2 == min(CO2_ppm - CO2[CO2_ppm - CO2 > 0]))]
  upper_CO2 <- CO2[which(CO2_ppm - CO2 == max(CO2_ppm - CO2[CO2_ppm - CO2 <= 0]))]
  CO2_contribution_lower <- 1 - (CO2_ppm - lower_CO2)/(upper_CO2 - lower_CO2)
  CO2_contribution_upper <- 1 - (upper_CO2 - CO2_ppm)/(upper_CO2 - lower_CO2)
}


###### TOPOGRAPHY
template <- raster(nrows=48, ncols=96, crs="+proj=longlat +datum=WGS84 +no_defs", resolution = c(3.75, 3.75)) # raster template for interpolation (resolution of climate model output)

# read and resample topography to model resolution
file_list <- list.files(topo_maps_prefix)
index <- grep(paste('_', t_geol/-1, 'Ma.nc', sep=""), file_list)
topo <- raster(paste(topo_maps_prefix,file_list[index], sep=""))
topo <- resample(topo, template, method="bilinear")
topo_land <- topo ; topo_land[topo_land < 0] = NA


###### TEMPERATURE
# Interpolate climate from climate stack
Tair_raster_lower_co2 <- rotate(mean(brick(paste(climate_prefix,"run_",t_geol/(-1),"_Ma/ymonmean_lonlat_run_",lower_CO2,"_ppm.nc", sep=""), var="tas")))
Tair_raster_upper_co2 <- rotate(mean(brick(paste(climate_prefix,"run_",t_geol/(-1),"_Ma/ymonmean_lonlat_run_",upper_CO2,"_ppm.nc", sep=""), var="tas")))
Tair_raster_current_co2 <- (CO2_contribution_lower * Tair_raster_lower_co2 + CO2_contribution_upper * Tair_raster_upper_co2)-273.15 #[K] -> [°C]

GAST <- weighted.mean(Tair_raster_current_co2, area(Tair_raster_current_co2), na.rm=TRUE)
print(paste("At current",CO2_ppm, "ppm conc, GAST is",round(GAST, digits=2), sep=" "))

###### RUNOFF
precip_raster_lower_co2 <- rotate(mean(brick(paste(climate_prefix,"run_",t_geol/(-1),"_Ma/ymonmean_lonlat_run_",lower_CO2,"_ppm.nc", sep=""), var="prl"))) +
  rotate(mean(brick(paste(climate_prefix,"run_",t_geol/(-1),"_Ma/ymonmean_lonlat_run_",lower_CO2,"_ppm.nc", sep=""), var="prc")))
evap_raster_lower_co2 <- rotate(mean(brick(paste(climate_prefix,"run_",t_geol/(-1),"_Ma/ymonmean_lonlat_run_",lower_CO2,"_ppm.nc", sep=""), var="evap")))
runoff_raster_lower_co2 <- precip_raster_lower_co2 + evap_raster_lower_co2

precip_raster_upper_co2 <- rotate(mean(brick(paste(climate_prefix,"run_",t_geol/(-1),"_Ma/ymonmean_lonlat_run_",upper_CO2,"_ppm.nc", sep=""), var="prl"))) +
  rotate(mean(brick(paste(climate_prefix,"run_",t_geol/(-1),"_Ma/ymonmean_lonlat_run_",upper_CO2,"_ppm.nc", sep=""), var="prc")))
evap_raster_upper_co2 <- rotate(mean(brick(paste(climate_prefix,"run_",t_geol/(-1),"_Ma/ymonmean_lonlat_run_",upper_CO2,"_ppm.nc", sep=""), var="evap")))
runoff_raster_upper_co2 <- precip_raster_upper_co2 + evap_raster_upper_co2

runoff_raster_current_co2 <- (CO2_contribution_lower * runoff_raster_lower_co2 + CO2_contribution_upper * runoff_raster_upper_co2)*(365*24*60*60*1000) # [m/s2] -> [mm/year]

# crop to continents
runoff_raster_current_co2[topo < 0] = NA
# negative net water flux means net evaporation, no ruoff
runoff_raster_current_co2[runoff_raster_current_co2 <= 0] = 0


###### budyko aridity = Rn / (lambda * precip) -> Three raster needed: 1) Temperature 2) Net radiation 3) Precipitation 
# 1) Temperature 
Tair_raster_current_co2_K <- Tair_raster_current_co2 + 273.15 # [K]

# lambda is latent heat of evaporation in MJ m-3 and is a function of temperature
latent_heat <- (3.146 - 0.002361*Tair_raster_current_co2_K) * 10^3

# 2) Radiation components
rss_raster_lower_co2 <- rotate(mean(brick(paste(climate_prefix,"run_",t_geol/(-1),"_Ma/ymonmean_lonlat_run_",lower_CO2,"_ppm.nc", sep=""), var="rss"))) # Net shorwave radiation
rls_raster_lower_co2 <- rotate(mean(brick(paste(climate_prefix,"run_",t_geol/(-1),"_Ma/ymonmean_lonlat_run_",lower_CO2,"_ppm.nc", sep=""), var="rls"))) # Net longwave radiation
rnet_lower_co2 <- rss_raster_lower_co2 + rls_raster_lower_co2 # Net radiation

rss_raster_upper_co2 <- rotate(mean(brick(paste(climate_prefix,"run_",t_geol/(-1),"_Ma/ymonmean_lonlat_run_",upper_CO2,"_ppm.nc", sep=""), var="rss")))
rls_raster_upper_co2 <- rotate(mean(brick(paste(climate_prefix,"run_",t_geol/(-1),"_Ma/ymonmean_lonlat_run_",upper_CO2,"_ppm.nc", sep=""), var="rls")))
rnet_upper_co2 <- rss_raster_upper_co2 + rls_raster_upper_co2

rnet_raster_current_co2 <- CO2_contribution_lower * rnet_lower_co2 + CO2_contribution_upper * rnet_upper_co2
rss_raster_current_co2 <- CO2_contribution_lower * rss_raster_lower_co2 + CO2_contribution_upper * rss_raster_upper_co2

# crop to continents
rnet_raster_current_co2[topo < 0] <- NA
rss_raster_current_co2[topo < 0] <- NA

# Unit conversion 
rss_raster_current_co2 <- rss_raster_current_co2*(365*24*60*60*1e-6) # [J s-1 m-2] -> [MJ year-1 m-2]
rnet_raster_current_co2 <- rnet_raster_current_co2*(365*24*60*60*1e-6) # [J s-1 m-2] -> [MJ year-1 m-2]

# 3) Total precipitation [m s-1]
precip_raster_current_co2 <- CO2_contribution_lower * precip_raster_lower_co2 + CO2_contribution_upper * precip_raster_upper_co2
# crop to continents
precip_raster_current_co2[topo < 0] <- NA

# Unit conversion
precip_raster_current_co2 <- precip_raster_current_co2 * 365*24*60*60 # [m s-1] -> [m year-1]

# Budyko aridity index
budyko_aridity <- rnet_raster_current_co2/(latent_heat*precip_raster_current_co2)
budyko_aridity[budyko_aridity <= 0] <- 0
budyko_aridity[budyko_aridity >= 6] <- 6 # declare everything above 6 as super arid, following Takeshima et al. 2020



########################################################################################################
# 2.) Biology
########################################################################################################

# Initialise gen3sis model 
if (timestep == 0){
  # Create landscape for gen3sis - repeat initial topography and climate for number of time steps - will be updated for each following time step (but gen3sis needs a complete landscape to start)
  # This transfers climate data to the gen3sis model
  landscapes_list <- list()
  for (i in c(0:(total_time_steps))){
    landscapes_list$Tair <- c(landscapes_list$Tair, Tair_raster_current_co2)
    landscapes_list$topo <- c(landscapes_list$topo, topo_land)
    landscapes_list$budyko_aridity <- c(landscapes_list$budyko_aridity, budyko_aridity)
    landscapes_list$radiation <- c(landscapes_list$radiation, rss_raster_current_co2)
  }

  # Standard cost function to travel grid cells - used in dispersal function of gen3sis 
  cost_function_dispersal <- function(source, habitable_src, dest, habitable_dest) {
    if (!all(habitable_src, habitable_dest)) {
      return(4/1000)
    } else {
      return(1/1000)
    }
  }

  create_input_landscape(landscapes = landscapes_list,
                          timesteps = as.character(0:(total_time_steps)),
                          cost_function = cost_function_dispersal,
                          directions = 8,
                          output_directory = "biology/model_landscape",
                          overwrite = T,
                          crs = "+proj=longlat +datum=WGS84 +no_defs",
                          calculate_full_distance_matrices = T,
                          verbose = T)

  print("starting gen3sis part")
  
  # Transfer important parameters to config file
  parameter_list <- list()
  parameter_list[["total_time_steps"]] <- total_time_steps
  parameter_list[["thermal_adapt_par"]] <- thermal_adapt_par*DT # As temperature adaptation rates were initialized in °C per year, multiplication with time step results in °C per timestep (not needed for other parameters as initialized in rate / time step)
  parameter_list[["aridity_adapt_par"]] <- aridity_adapt_par
  parameter_list[["dispersal_par"]] <- dispersal_par
  parameter_list[["timestep"]] <- timestep
  saveRDS(parameter_list, "./biology/parameter_list.rds")
  
  
  # This runs the gen3sis vegetation module
  vegetation_simulation <- run_simulation(config = "biology/biology_configuration_file.R",
                           landscape = "biology/model_landscape",
                           verbose=T, #  progress printed
                           output_directory="biology/model_output")
  
  unlink("./biology/parameter_list.rds")
  
} else if (timestep > 0){
  
  # Read and modify landscape according to current atmospheric CO2
  
  landscapes <- readRDS("./biology/model_landscape/landscapes.rds")

  Tair_raster_current_co2_cropped <- Tair_raster_current_co2
  Tair_raster_current_co2_cropped[is.na(topo_land)] <-  NA
  Tair_column <- as.data.frame(Tair_raster_current_co2_cropped, xy=TRUE)[, -c(1:2)]
  landscapes$Tair[, as.character(timestep)] <- Tair_column
  
  budyko_aridity_column <- as.data.frame(budyko_aridity, xy=TRUE)[, -c(1:2)]
  landscapes$budyko_aridity[, as.character(timestep)] <- budyko_aridity_column
  
  radiation_column <- as.data.frame(rss_raster_current_co2, xy=TRUE)[, -c(1:2)]
  landscapes$radiation[, as.character(timestep)] <- radiation_column
  
  # Save updated landscape in gen3sis landscape folder
  saveRDS(landscapes, file="./biology/model_landscape/landscapes.rds")

  print("starting gen3sis part")
  # Transfer important parameters to config file
  parameter_list <- list()
  parameter_list[["total_time_steps"]] <- total_time_steps
  parameter_list[["thermal_adapt_par"]] <- thermal_adapt_par*DT
  parameter_list[["aridity_adapt_par"]] <- aridity_adapt_par
  parameter_list[["dispersal_par"]] <- dispersal_par
  parameter_list[["timestep"]] <- timestep
  saveRDS(parameter_list, "./biology/parameter_list.rds")
  vegetation_simulation <- run_simulation(config = "biology/biology_configuration_file.R",
                                         landscape = "biology/model_landscape",
                                         verbose=T, #  progress printed
                                         output_directory="biology/model_output")
  unlink("./biology/parameter_list.rds")
}



########################################################################################################
# 3.) Carbon cycle model - calculate all relevant fluxes and atmosphere-ocean carbon mass balance
########################################################################################################

# Template
template$area <- area(template) * 1e+6 # in m2
template_df <- as.data.frame(template, xy=T)

# Import all results from gen3sis vegetation model
landscape <- readRDS(paste("./biology/model_output/biology_configuration_file/landscapes/landscape_t_",timestep,".rds", sep=""))
landscape_df <- as.data.frame(landscape$coordinates, xy=TRUE)
landscape_df$Tair <- landscape$environment[, "Tair"]
landscape_df$budyko_aridity <- landscape$environment[, "budyko_aridity"]
landscape_df$radiation <- landscape$environment[, "radiation"]
landscape_df$topo <- landscape$environment[, "topo"]
output <- left_join(template_df, landscape_df, by=c("x", "y"))

# Import vegetation state 
vegetation_traits <- readRDS(paste("./biology/model_output/biology_configuration_file/world_trait_map/map_t_",timestep,".rds", sep=""))
vegetation_traits_df <- as.data.frame(vegetation_traits, xy=TRUE)
output <- left_join(output, vegetation_traits_df, by=c("x", "y"))


######################################################################
################ Productivity calculations 
#####

NPP_init_adapt <- 0.0142 # Energy conversion factor - calibrated to reproduce 3.5e+12 mol C yr-1 locb (Lenton et al. 2018)

# Temperature limitation of productivity (limiting high and low temperature prouctivity)
output$Tair_limit <- approx(c(-10, 0, 45, 55), y = c(0, 1, 1, 0), xout = output$Tair, method = "linear", n = 50,
                            yleft = 0, yright = 0, rule = 1, f = 0, ties = mean, na.rm = TRUE)$y


# Aridity response function: Defined to be ~1 in humid zones (budyko: 0-1.2), mid-point at budyko 2 (transition semi-humid/semi-arid) and 0 around budyko 4-6 (arid-hyperarid)
output$aridity_limit <- 1 - 1 / (1 + exp(-3*(output$budyko_aridity - 2)))


# Adaptation scaling
output$Tair_adapt_limit <- exp(-0.1 * (output$Tdiff^2)) # This function defines the niche width: fixed as 10°C for temperature
output$Tair_adapt_limit[which(output$Tair_adapt_limit<0)] <- 0

output$aridity_adapt_limit <- exp(-5 * (output$Adiff^2)) # This function defines the aridity niche width: fixed as 2 budyko aridity units
output$aridity_adapt_limit[which(output$aridity_adapt_limit<0)] <- 0

# Set location with no/extinct species adaptation to 0 
output$Tair_adapt_limit[is.na(output$Tdiff) & output$topo >= 0] <- 0
output$aridity_adapt_limit[is.na(output$Tdiff) & output$topo >= 0] <- 0

# Productivity 
burial_rate <- 0.0007 # constant burial fraction derived from present day rates of NPP and locb (an assumed fraction 0.0007 of global NPP gets buried in sediments)
output$productivity_rate <- NPP_init_adapt * output$radiation * output$Tair_limit * output$aridity_limit * output$Tair_adapt_limit * output$aridity_adapt_limit * 12 # multiplication with 12 for record in g C / m2
output$productivity_adapt <- NPP_init_adapt * output$radiation * output$Tair_limit * output$aridity_limit * output$Tair_adapt_limit * output$aridity_adapt_limit *  output$area * burial_rate 
output$productivity_adapt_sum_global <- sum(output$productivity_adapt, na.rm=TRUE)

# Calculate weathering enhancement based on normalized productivity potential [0 - 1] - radiation saturation at rss 5000 MJ per m2 per year (tropical conditions)
output$weathering_limit_adapt <-  (output$radiation/5000) * output$Tair_limit * output$aridity_limit * output$Tair_adapt_limit * output$aridity_adapt_limit
output$weathering_limit_adapt[output$weathering_limit_adapt > 1] = 1

# For record track adaptation limitation and physiological limitatino of productivity
output$adaptation_limit <-  output$Tair_adapt_limit  * output$aridity_adapt_limit
output$physiological_limit <-  (output$radiation/5000) * output$Tair_limit * output$aridity_limit 
output$physiological_limit[output$physiological_limit > 1] = 1

# For calibration
#output$productivity_adapt_sum_global[1]

######################################################################
################ Silicate weathering calculation 
#####

# Calculate slope 
output_raster <- rasterFromXYZ(output, crs="+proj=longlat +datum=WGS84 +no_defs")
topo_slope <- output_raster$topo
topo_slope[is.na(topo_slope)] <- 0 # assume 0 elevation for oceans, else no slope is calculated for coast
slope <- raster::terrain(topo_slope, opt="slope", unit="tangent", neighbors=8)
slope[is.na(output$topo)] <- NA # Crop to continents
output_raster$slope <- slope

# Runoff
output_raster$runoff <- runoff_raster_current_co2
output <- as.data.frame(output_raster, xy=TRUE)

# Erosion based on slope and runoff
k_e <- 8.2e-3 # scaling parameter to obtain present day total erosion of approximately 16 Gt (Mills et al. 2021)
output$erosion_rate <- k_e * (output$runoff)^0.5 * output$slope # Erosion definition from Maffre et al. (2018)
output$erosion_rate_area_scaled <- output$erosion_rate * output$area
# # Calibration of k_e --> 1.6e+10
#sum(output$erosion_rate_area_scaled, na.rm=TRUE)/1e+10

# All silicate weathering parameters in the following adapted from Mills et al. 2021
# Runoff dependency 
k_w <- 1e-3 # Silicate weathering flow dependence, calibration constant
f_Q <- (1-exp(-k_w*output$runoff)) # silicate weathering runoff dependency

# Temperature dependency
E_a <- 20 # Apparent activation energy (20 kJ mol-1)
R <- 8.314e-3 # ideal gas constant
T_0 <- 286 # reference temperature
f_T <- exp((E_a/(R*T_0))-(E_a/(R*(output$Tair+273.15)))) # temperature dependency

# Erosion dependency 
sigma_factor <- 0.9 # Reaction time parameter = sigma + 1 
Z <- 10 # Silicate weathering zone depth [m]
f_E <- ((Z/output$erosion_rate)^sigma_factor)/(sigma_factor) # Erosion dependence
f_E[is.infinite(f_E)] = 2e+7 # Account for zero-erosion grid cells, does not affect results, will be zero in any case

# kinetic limitation of weathering reaction
f_kinetic = f_Q * f_T * f_E

# Silicate weathering 
X_m <- 0.1 # Silicate cation weight fraction
K <- 6e-5 # Silicate weathering grain size dependence

# Weathering reaction following Mills et al. 2021 / West 2012
W_sil = output$erosion_rate * X_m * (1-exp(-K*f_kinetic))

# Vegetation-mediated weathering enhancement following formulation from Lenton et al. 2018
bio_enhancement = ( 1 - pmin(output$weathering_limit_adapt , 1, na.rm=FALSE) ) * PREPLANT * (RCO2^0.5) + (output$weathering_limit_adapt)
output$bio_enhancement <- bio_enhancement

W_sil = W_sil * bio_enhancement
output$silw <- W_sil

W_sil_area_scaled = W_sil * output$area 

W_sil_tot = sum(W_sil_area_scaled, na.rm=T)

k_silw <- 1.0e+13
if(as.factor(PREPLANT) == "0.25"){
  output$silw_total_mol <- k_silw * (W_sil_tot/339296153) # denominator a calibration constant to get 1e+13 burial with present day climate and topography and PREPLANT (1/4)
}
if(as.factor(PREPLANT) != "0.25"){
  output$silw_total_mol <- k_silw * (W_sil_tot/331388706) # denominator a calibration constant to get 1e+13 burial with present day climate and topography and PREPLANT (1/6)
}



######################################################################
################ Carbonate Weathering - to track isotopic effect, not relevant for atmosphere-ocean carbon mass balance
#####
k_carbw_scale <- 2.42e-4 # scaling parameter, calibrated to reproduce present-day carbw of 8e+12
carbw <- k_carbw_scale * output$runoff * output$bio_enhancement
output$carbw <- carbw
carbw_area_scaled <- carbw * output$area
carbw_tot <- sum(carbw_area_scaled, na.rm=T)
output$carbw_total_mol <- carbw_tot

# caliration  - 8e+12
# carbw_tot


######################################################################
################ Marine productivity calculation
#####

output$topo_complete <- as.data.frame(topo) # record topography in output file
output$Tair_landandocean <- as.data.frame(Tair_raster_current_co2)
output$Tair_ocean <- output$Tair_landandocean
output$Tair_ocean[output$topo_complete >= 0] <- NA

mocb_initial = 3.01 # scaling parameter, calibrated to reproduce present-day marine organic carbon burial of 3.5e+13 (Lenton et al. 2018)
output$Tair_ocean[output$Tair_ocean <= 0 | output$Tair_ocean > 33] <- NA # temperature limits for marine productivity
output$Tlim_marine <- -3.27*(10^-8)*(output$Tair_ocean^7) + 3.4132*(10^-6)*(output$Tair_ocean^6) - 1.348*(10^-4)*(output$Tair_ocean^5) + 2.462*(10^-3)*(output$Tair_ocean^4) - 0.0205*(output$Tair_ocean^3) + 0.0617*(output$Tair_ocean^2) + 0.2749*output$Tair_ocean + 1.2956 
nutrient_input = output$silw_total_mol[1] / k_silw # nutrient input to ocean is proportional to silicate weathering (assuming weathering as the main nutrient input to ocean)
output$mocb <- mocb_initial * nutrient_input * output$Tlim_marine
mocb <- sum(output$mocb* output$area * burial_rate, na.rm=TRUE)
# Calibration
# mocb


######################################################################
################ Seafloor weathering
#####

# Optional for sensitivity test - feedback function is uncertain (seafloor weathering is a deep ocean process)
# # Global function following Mills et al.2021
if(seafloor_weathering_feedback){
  k_sfw <- 1.75e+12
} else {
  k_sfw <- 0
}
f_T_sfw <- exp(0.0608*((GAST+273.15)-288))
sfw <- k_sfw * f_T_sfw


#####################################################################
################ Degassing
##### 

# Define initial degassing rate that matchs carbon sinks - assuming the system to be in steady state at initialization
if(timestep == 0){
  degass_0 <- (output$silw_total_mol[1] + output$productivity_adapt_sum_global[1] + mocb + sfw)
  degass_current <- degass_0
}
# Include random variation in degassing rate around mean
if(timestep > 0){
  degass_current <- degass_0 + rnorm(1, mean = 0, sd = 1e+11)
}

# LIP degassing: mol C per year
if(degass_duration == 2 & timestep %in% c(9:10)){
  degass_current <- degass_current + (degass_mass * (1/12) * (1/DT) * (1/2))
}
if(degass_duration == 1 & timestep == 10){
  degass_current <- degass_current + (degass_mass * (1/12) * (1/DT))
}
output$degassing_rate <- degass_current

########################################################################################################
# 4.) Atmosphere-Ocean Carbon Mass Balance
########################################################################################################

# Change of carbon in atmosphere-ocean pool 
dA <- degass_current - (output$silw_total_mol[1] + output$productivity_adapt_sum_global[1] + mocb + sfw) # per year 
# Integrate over time step 
dA_integrated <- dA * DT
A <- A + dA_integrated 

RCO2 = (A/A0)^2 # Atmospheric CO2 responds quadratic to amount of C in Atmosphere-Ocean System, A0 = carbon content at start of simulation
CO2_ppm <- CO2_ppm_start*RCO2


print(paste("global silicate weathering:",output$silw_total_mol[1], sep=" "))
print(paste("global organic carbon burial:",output$productivity_adapt_sum_global[1], sep=" "))
print(paste("Change in atmosphere-ocean carbon mass balance:", dA, sep=" "))
print(paste("Atmosphere-ocean carbon total:", A, sep=" "))
print(paste("Atmospheric concentration of CO2:", CO2_ppm, sep=" "))


output$timestep <- timestep 
output$cont_config <- t
output$dispersal_par <- dispersal_par
output$thermal_adapt_par <- thermal_adapt_par
output$aridity_adapt_par <- aridity_adapt_par
output$preplant_par <- PREPLANT

# Collect all output
output_collection[[paste(timestep)]] <- output

# these summary outs are recorded with one timestep offset because no zero indexing in R
out[timestep + 1, "time"] = timestep 
out[timestep + 1, "CO2"] = CO2_ppm
out[timestep + 1, "locb"] = output$productivity_adapt_sum_global[1] 
out[timestep + 1, "mocb"] = mocb
out[timestep + 1, "silw"] = output$silw_total_mol[1]
out[timestep + 1, "carbw"] = output$carbw_total_mol[1]
out[timestep + 1, "sfw"] = sfw
out[timestep + 1, "degass"] = output$degassing_rate[1]
out[timestep + 1, "GAST"] = GAST

}



# After main loop has finished, save the model output 
saveRDS(output_collection, "model_run_output_detail.rds")
saveRDS(out, "model_run_output_summary.rds")


# Save a quick overview of modeled fluxes of current model
pdf("Quickview_model_run.pdf", width=8, height=10)
par(mfrow=c(3, 2))
plot(out[,1], out[,2], type="l", lwd = 2, ylab="ppm")
plot(out$time, out$locb, type="l", lwd = 2, ylab="NPP")
plot(out$time, out$mocb, type="l", lwd = 2, ylab="NPPocean")
plot(out$time, out$silw, type="l", lwd = 2, ylab="silw")
plot(out$time, out$GAST, type="l", lwd = 2, ylab="GAST")
plot(out$time, out$sfw, type="l", lwd = 2, ylab="sfw")
dev.off()


print("run completed with no error")
}
