
######################################
###            METADATA            ###
######################################
# gen3sis configuration
#
# Author: Julian Rogger
#
# Date: 30.05.2024
#
# This is the configuration file for the gen3sis-based eco-evolutionary vegetation model
# It defines the behavior of the modeled floras through a series of functions that are sequentially applied: dispersal, evolution, competition
#
# Note the term species and flora are being used interchangeably in the following. The focus on species is inherited from the original gen3sis model. 
######################################


######################################
###         General settings       ###
######################################

# set the random seed for the simulation - this can be useful for debugging
random_seed = NA

# Transfer parameters from frontend file that are needed in the configuration file 
paste("reading in the parameters")
parameter_list <- readRDS("parameter_list.rds")
total_time_steps <- parameter_list[["total_time_steps"]]
thermal_adapt_par <- parameter_list[["thermal_adapt_par"]]
aridity_adapt_par <- parameter_list[["aridity_adapt_par"]]
dispersal_par <- parameter_list[["dispersal_par"]]
timestep <- parameter_list[["timestep"]]

# set the starting time
start_time = timestep 

# set the end time step - to run only one timestep (needed for climate update), set end time equal to starting time
end_time = start_time

# maximum total number of floras in the simulation before it is aborted - should never be reached as only one flora per grid cell survives
max_number_of_species = 100000

# maximum number of floras within one cell before the simulation is aborted - should never be reached as only one flora per grid cell survives
max_number_of_coexisting_species = 100000

# a list of traits to include with each flora
# a "dispersal" trait is implicitly added in any case
trait_names = c("Topt", "Tdiff", "Aopt", "Adiff")

# ranges to scale the input environments with:
# not listed variable:         no scaling takes place
# listed, set to NA:           the environmental variable will be scaled from [min, max] to [0, 1]
# listed with a given range r: the environmental variable will be scaled from [r1, r2] to [0, 1]
environmental_ranges = NA



######################################
###            Observer            ###
######################################

# a place to inspect the internal state of the simulation, save output and collect additional information if desired
end_of_timestep_observer = function(data, vars, config){
  
  save_landscape()
  save_species()
  landscape_coordinates <- as.data.frame(data$landscape$coordinates)
  landscape_coordinates$Topt <- NA
  landscape_coordinates$Tdiff <- NA
  landscape_coordinates$Aopt <- NA
  landscape_coordinates$Adiff <- NA
  landscape_coordinates$species_id <- NA
  for(i in 1:length(data$all_species)){
    traits <- data$all_species[[i]]$traits
    landscape_coordinates$Topt[which(rownames(landscape_coordinates) %in% rownames(traits))] <- traits[, "Topt"]
    landscape_coordinates$Tdiff[which(rownames(landscape_coordinates) %in% rownames(traits))] <- traits[, "Tdiff"]
    landscape_coordinates$Aopt[which(rownames(landscape_coordinates) %in% rownames(traits))] <- traits[, "Aopt"]
    landscape_coordinates$Adiff[which(rownames(landscape_coordinates) %in% rownames(traits))] <- traits[, "Adiff"]
    landscape_coordinates$species_id[which(rownames(landscape_coordinates) %in% rownames(traits))] <- i
  }
  traits_world <- raster::rasterFromXYZ(landscape_coordinates, res=c(3.75, 3.75) , crs="+proj=longlat +datum=WGS84 +no_defs")
  saveRDS(traits_world, file=paste(config$directories$output,"/world_trait_map/map_t_",vars$ti,".rds", sep=""))
}


######################################
###         Initialization         ###
######################################

# the initial abundance of a newly colonized cell, both during setup and later when 
# colonizing a cell during the dispersal.
initial_abundance = 1

# place species in the landscape:
create_ancestor_species <- function(landscape, config){
  
  if(timestep == 0){
  
  landcells <- length(landscape$environment[, 1]) # as many floras as there are land cells - each optimally adapted
  all_species <- list()
  for(i in 1:landcells){
    new_species <- create_species(rownames(landscape$environment)[i], config)
    new_species$traits[ , "dispersal"] <- 0
    new_species$traits[ , "Topt"] <- landscape$environment[rownames(landscape$environment)[i], "Tair"]
    new_species$traits[ , "Tdiff"] <- 0
    
    new_species$traits[ , "Aopt"] <- landscape$environment[rownames(landscape$environment)[i], "budyko_aridity"]
    new_species$traits[ , "Adiff"] <- 0

    all_species <- append(all_species, list(new_species))
  }
  return(all_species)

  } else {
    
  # Import floras from last time step
  filename <- paste('/species/species_t_', start_time-1, '.rds', sep="")
  species <- readRDS(paste(config$directories$output, filename, sep=""))
  
  # Clean gen3sis species object from non-abundant species to save memory/time
  index_list <- lapply(species, FUN = function(input){
    index <- NA
    if(length(input$abundance) >= 1){index <- 1}
    if(length(input$abundance) == 0){index <- 0}
    return(index)
  })
  index_vec <- do.call(rbind, index_list)
  if(any(index_vec == 0)){
    index <- which(index_vec == 0)
    species <- species[-c(index)]
  }
  
  return(species)  
  }
}



######################################
###             Dispersal          ###
######################################

# the maximum range to consider when calculating the distances from local distance inputs.
max_dispersal <- Inf

# returns n dispersal values.
get_dispersal_values <- function(n, species, landscape, config) {
  
  values <- rweibull(n, shape=2, scale=dispersal_par) # Weibull-based dispersal kernel

  # Dispersal ability scales with primary productivity potential
  # All climatic variables are derived from the carbon-climate model
  occurrence <- names(species$abundance)
  
  Tair_occurrence <- landscape$environment[occurrence, "Tair"]
  Tairlim_occurrence <- approx(c(-10, 0, 45, 55), y = c(0, 1, 1, 0), xout = Tair_occurrence, method = "linear", n = 50,
                               yleft = 0, yright = 0, rule = 1, f = 0, ties = mean, na.rm = TRUE)$y
  
  budyko_aridity_occurrence <- landscape$environment[occurrence, "budyko_aridity"]
  ariditylim_occurrence <- 1 - 1 / (1 + exp(-3*(budyko_aridity_occurrence - 2)))
  
  radiation_occurrence <- landscape$environment[occurrence, "radiation"]
  radiation_lim <- radiation_occurrence/5000 ; radiation_lim[radiation_lim > 1] <- 1 # Assume light saturation at tropical intensity of 5000 MJ year-1 for normalization of productivity potential
  
  # Consider limited degree of adaptation
  Tdiff_species <- species$traits[, "Tdiff"]
  Tdiff_lim <- exp(-0.1 * (Tdiff_species^2)) # niche width around 10°C
  
  Adiff_species <- species$traits[, "Adiff"]
  Adiff_lim <- exp(-5 * (Adiff_species^2))
  
  productivity_potential <- Tairlim_occurrence * ariditylim_occurrence * radiation_lim * Tdiff_lim * Adiff_lim

  # Don't let dispersal go to zero - only max reduction of 50%
  productivity_potential[productivity_potential < 0.5] <- 0.5
  scaled_dispersal <- values*productivity_potential
  
  return(scaled_dispersal)
}


######################################
###          Speciation            ###
######################################

# After dispersal, floras should be treated as two separate floras, therefore "speciation"
divergence_threshold = 0.9
get_divergence_factor <- function(species, cluster_indices, landscape, config) {
  return(1) # as soon as modeled species/plants/floras get separated, they are treated separately in next time step
}


######################################
###      Trait Evolution           ###
######################################

# Species' ability of adapting climatic niche through simplified evolution 
adaptation_increment_temperature <- thermal_adapt_par # Rate defined in frontend parameter space
adaptation_increment_aridity <- aridity_adapt_par # Rate defined in frontend parameter space

# mutate the traits of a species and return the new traits matrix.
apply_evolution <- function(species, cluster_indices, landscape, config){

    traits <- species[["traits"]]
    cells <- rownames(traits)
    
    # evolve Topt
    temperature_adaptation_change <- landscape$environment[cells, "Tair"] - traits[cells, "Topt"]
    temperature_adaptation_change[temperature_adaptation_change >= adaptation_increment_temperature] <-  adaptation_increment_temperature
    temperature_adaptation_change[temperature_adaptation_change <= (-adaptation_increment_temperature)] <- (-adaptation_increment_temperature)
    traits[cells, "Topt"] <- traits[cells, "Topt"] + temperature_adaptation_change 
    traits[cells, "Tdiff"] <- abs(traits[cells, "Topt"] - landscape$environment[cells,"Tair"])
    
    # evolve Aopt
    aridity_adaptation_change <- landscape$environment[cells, "budyko_aridity"] - traits[cells, "Aopt"]
    aridity_adaptation_change[aridity_adaptation_change >= adaptation_increment_aridity] <-  adaptation_increment_aridity
    aridity_adaptation_change[aridity_adaptation_change <= (-adaptation_increment_aridity)] <- (-adaptation_increment_aridity)
    traits[cells, "Aopt"] <- traits[cells, "Aopt"] + aridity_adaptation_change 
    traits[cells, "Adiff"] <- abs(traits[cells, "Aopt"] - landscape$environment[cells,"budyko_aridity"])
  return(traits)
}


######################################
###             Ecology            ###
######################################

# called for every cell with all occurring floras, this function calculates abundances and/or 
# who survives for each grid cell based on the degree of adaptation.
# returns a vector of abundances.
# set the abundance to 0 for every species supposed to die.

apply_ecology <- function(abundance, traits, landscape, config) {
  
  # Competition - best adapted biome survives 
  Tair_diff <- abs( traits[, "Topt"] - landscape[, "Tair"])
  budyko_aridity_diff <- abs( traits[, "Aopt"] - landscape[, "budyko_aridity"])
  
  Tdiff_lim <- exp(-0.1 * (Tair_diff^2)) # niche width around 10°C
  Adiff_lim <- exp(-5 * (budyko_aridity_diff^2))
  
  performance <- Tdiff_lim * Adiff_lim # performance function measuring the degree of local adaptation, only the best adapted flora will survive
  
  # Choose the best adapted flora:
  abundance[names(abundance[performance != max(performance)])] = 0
  abundance[names(abundance[performance == max(performance)])] = 1
  
  # In case there are several equally well adapted floras, random choice: 
  survivor_name <- sample(names(abundance[abundance == 1]), 1)  
  abundance[survivor_name] = 1
  abundance[names(abundance) != survivor_name] = 0
  
  # Hard limits on temperature difference and aridity 
  # This will kill floras, in the next step a new one will take over or space if left empty
  abundance[which(traits[, "Tdiff"] >= 15)] <- 0 
  abundance[which(traits[, "Adiff"] >= 5)] <- 0 

  return(abundance)
}
