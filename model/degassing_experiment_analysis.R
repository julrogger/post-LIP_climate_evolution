########################################################################################################################################################
# This is a file that can be used to categorize the degassing model output (i.e., die out, new steady state, return or transient trajectory)
# It is not needed to run the model, only for post-processing
#
# Author: Julian Rogger
# Last update: 30.05.24
########################################################################################################################################################

# libraries 
library(raster)
library(ncdf4)
library(parallel)
library(ggplot2)

# 1.) Summary stats 

num_of_experiments <- length(list.files(path=".", pattern="experiment_run", recursive = F))
split_vec <- c(1:num_of_experiments)
exp_list <- split(split_vec, f=c(1:length(split_vec)))
list_function <- function(exp_ID){

  library(raster)
  library(ncdf4)
  library(parallel)
  library(ggplot2)
  print(exp_ID)

  # Read data
  out_summary <- readRDS(paste("experiment_run_", exp_ID,"/model_run_output_summary.rds", sep=""))
  
  # Convert time steps which go from 1-81 to million years after simulation start
  out_summary$time <- (out_summary$time)*(1e+5/1e+6)
  out_detail <- readRDS(paste("experiment_run_", exp_ID,"/model_run_output_detail.rds", sep=""))

  # Get parameter values - they are stored in detailed output, but the same for every time step
  dispersal_par <- out_detail[[1]]$dispersal_par[1]
  thermal_adapt_par <- out_detail[[1]]$thermal_adapt_par[1]
  aridity_adapt_par <- out_detail[[1]]$aridity_adapt_par[1]
  preplant_par <- out_detail[[1]]$preplant_par[1]
  exp_ID <- exp_ID

  # Get Maximum GAST and other marks of extreme
  GAST_max <- max(out_summary$GAST)
  CO2_ppm_max <- max(out_summary$CO2)
  index_maxtemp <- which(out_summary$GAST == max(out_summary$GAST))[1] # first year of maximum temperature in case several times the same

  # Categorise experiments - return, infinite, new_steady
  initial_temp <- out_summary[as.factor(out_summary$time) == "0","GAST"]
  final_temp <- out_summary[as.factor(out_summary$time) == "8", "GAST"]
  exp_type <- NA
  if(abs(final_temp - initial_temp) <= 0.5){
    exp_type <- "return"
    print(exp_type)
  }
  if(max(out_summary$CO2) > 8000 & is.na(exp_type)){
    exp_type <- "infinite" # should not happen
    print(exp_type)
  }
  delta_endphase <- out_summary[as.factor(out_summary$time) == "8", "GAST"] - out_summary[as.factor(out_summary$time) == "6.1", "GAST"]
  delta_locb_endphase <- out_summary[as.factor(out_summary$time) == "8", "locb"] - out_summary[as.factor(out_summary$time) == "6.1", "locb"]
  print(paste(exp_ID,"has endphase delta", delta_endphase, sep=" "))
  if(is.na(exp_type) & abs(delta_endphase) <= 0.5 & abs(delta_locb_endphase) <= 0.05e+12){
    exp_type <- "new_steady"
    print(exp_type)
  }
  if(is.na(exp_type)){
    exp_type <- "otw"
    print(exp_type)
  }

  df <- data.frame("exp_ID" = exp_ID,
                   "dispersal_par" = dispersal_par,
                   "thermal_adapt_par" = thermal_adapt_par,
                   "aridity_adapt_par" = aridity_adapt_par,
                   "preplant_par" = preplant_par,
                   "GAST_max" = GAST_max,
                   "CO2_ppm_max" = CO2_ppm_max,
                   "exp_type" = exp_type)
  return(df)
}
cl <- makeCluster(50, outfile="out_analysis_1st_step.txt")
out <- parLapply(cl = cl, exp_list, fun = list_function)
out_df <- do.call(rbind.data.frame, out)
saveRDS(out_df, "summary_measures_out.rds")
stopCluster(cl)
rm(list=ls())



# 2.) Put all summary output together to get temporal evolution 

num_of_experiments <- length(list.files(path=".", pattern="experiment_run", recursive = F))
split_vec <- c(1:num_of_experiments)
exp_list <- split(split_vec, f=c(1:length(split_vec)))
list_function <- function(exp_ID){
  print(exp_ID)

  # Read data
  out_summary <- readRDS(paste("experiment_run_", exp_ID,"/model_run_output_summary.rds", sep=""))
  # Convert timesteps which go from 1-81 to million years after simulation start
  out_summary$time <- (out_summary$time)*(1e+5/1e+6)
  out_detail <- readRDS(paste("experiment_run_", exp_ID,"/model_run_output_detail.rds", sep=""))

  # Get parameter values - they are stored in detailed output, but the same for every timestep
  dispersal_par <- out_detail[[1]]$dispersal_par[1]
  thermal_adapt_par <- out_detail[[1]]$thermal_adapt_par[1]
  aridity_adapt_par <- out_detail[[1]]$aridity_adapt_par[1]
  preplant_par <- out_detail[[1]]$preplant_par[1]
  exp_ID <- exp_ID

  out_summary$dispersal_par <- dispersal_par
  out_summary$thermal_adapt_par <- thermal_adapt_par
  out_summary$aridity_adapt_par <- aridity_adapt_par
  out_summary$preplant_par <- preplant_par
  out_summary$exp_ID <- exp_ID

  return(out_summary)
}
cl <- makeCluster(50, outfile="out_analysis_2nd_step.txt")
out <- parLapply(cl = cl, exp_list, fun = list_function)
out_df <- do.call(rbind.data.frame, out)
saveRDS(out_df, "appended_summary_table_out.rds")
stopCluster(cl)


