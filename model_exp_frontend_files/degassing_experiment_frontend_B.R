########################################################################################################################################################
# Frontend file for degassing model experiments
# This file defines the biological parameter space that will be explored and sets boundary conditions of paleogeography and LIP degassing scenario
#
# Author: Julian Rogger
# Last update: 30.05.24
########################################################################################################################################################

# Define the biological parameter space
dispersal_rates <- c(seq(0, 1000, 100), seq(1200, 2000, 200), 4000, 6000, 10000) # This should be given in km per time step (here: timestep = 1e+5 years)
temperature_niche_adaptation_rate <- (1/1e+6) *  c(seq(0, 0.1, 0.01), seq(0.2, 1, 0.1), seq(2, 10, 1))  # this should be given in °C per year
aridity_niche_adaptation_rate <- c(seq(0, 0.1, 0.01), seq(0.2, 1, 0.1), 5, 10) # This should be given in Budyko aridity units per time step 
PREPLANT <- c(0.25) # This parameter defines the maximum weathering enhancement effect of plants - current options: 1/4 = maximum 4-fold, 1/6 = maximum 6-fold enhancement

# Get all parameter combinations
pars <- expand.grid(dispersal_rates, temperature_niche_adaptation_rate, aridity_niche_adaptation_rate, PREPLANT)
pars$experiment_ID <- c(1:nrow(pars))
colnames(pars) <- c("dispersal_par", "thermal_par", "aridity_par", "PREPLANT","experiment_ID")
pars_list <- split(pars, seq(nrow(pars)))

# To reduce the time needed for the models and to help debugging, it may be 
# beneficial to split the parameter list into several list and run the following 
# for each list separately (i.e., as several jobs if working with a HPC cluster)

# Initialize parallel mode
library(parallel)
cl <- makeCluster(50, setup_strategy = "sequential", outfile="./out.txt", overwrite=T)

parLapply(cl = cl, pars_list, fun = function(df_entry){
  
  # Check if model has already been run, if yes exit the loop 
  if(file.exists(paste("./experiment_run_",df_entry$experiment_ID,"/model_run_output_summary.rds", sep=""))){
  	 print(paste("This experiment has already been run:", df_entry$experiment_ID, sep=" "))
  } else {
  
  # as computation might have already started: clear experiment file for no overwriting
  if(file.exists(paste("experiment_run_",df_entry$experiment_ID,sep=""))){
  unlink(paste("experiment_run_",df_entry$experiment_ID,sep=""), recursive = TRUE)
  }
	
  print(paste("Starting run for:", df_entry$experiment_ID, sep=" "))

  # Load libraries in workers
  library(raster)
  library(gen3sis)
  library(dplyr)
  library(terra)
  library(ncdf4)
  
  # Set up model directory - copy from source directory for each experiment
  experimental_ID <- df_entry$experiment_ID
  dir.create(paste("experiment_run_",experimental_ID,sep=""))
  file.copy("degassing_experiment_main.R", paste("./experiment_run_",experimental_ID,sep=""), overwrite = T)
  file.copy("biology", paste("./experiment_run_",experimental_ID,sep=""), recursive=T, overwrite = T)
  setwd(paste("./experiment_run_",experimental_ID,sep="")) # working in experimental directory now
  
  # source the function 
  source("degassing_experiment_main.R")

  # Fetch function parameters from list or define here for example runs 
  thermal_adapt_par <- df_entry$thermal_par[1] # in °C per Ma
  aridity_adapt_par <- df_entry$aridity_par[1] # in Budyko aridity units per time step
  dispersal_par <- df_entry$dispersal_par[1] # in km per time step 
  PREPLANT <- df_entry$PREPLANT[1] # maximum weathering enhancement factor, either (1/4) or (1/6)
  
  # Paleogeography and LIP degassing boundary conditions 
  CO2_ppm_start <- 500 # This sets the starting CO2 conentration for the model experiment
  t <- -55e+6 # This sets the geologic time period and paleogeography considered, current options are: -250e+6 for the Permian-Triassic, -200e+6 for the Triassic Jurassic, -55e+6 for the Paleocene-Eocene, 0 for the present day configuration
  degass_mass <- 15000e+15 # Set the amount of C degassing during the LIP Period, this should be given in grams
  degass_duration <- 1 # Number of time steps over which degassing occurs
  seafloor_weathering_feedback <- FALSE # Include seafloor weathering as an additional carbon cycle feedback? Parametrization is associated with some uncertainty and can only be considered as a global flux
  
  # This runs the main model using the parameters and boundary conditions specified here
  degassing_experiment_main_function(thermal_adapt_par, aridity_adapt_par, dispersal_par, PREPLANT, CO2_ppm_start, t, degass_mass, degass_duration, seafloor_weathering_feedback)
  
  # Delete unnecessary files to avoid storage problems
  unlink("degassing_experiment_main.R")
  unlink("biology", recursive = TRUE)
  setwd("..")
  }
})
# Close cluster
print("Successfully done with bunch")
warnings()
stopCluster(cl)
