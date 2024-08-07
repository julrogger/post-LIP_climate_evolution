This repository includes all files needed to reproduce the climate trajectories from the paper 'Biogeographic climate sensitivity controls Earth system response to Large Igneous Province carbon degassing' by Rogger et al. (https://www.science.org/doi/10.1126/science.adn3450)

The model represents a coupling of an eco-evolutionary vegetation model with a geologic carbon cycle model and a lookup structure of pre-run climate simulations. 
It can be used to study carbon cycling and the climatic evolution following a Large Igneous Province (LIP) carbon degassing event. 

To run the main model a directory containing a 'degassing_experiment_frontend.R' file, a 'degassing_experiment_main.R' file, the 'biology', 'climate_data' and 'paleogeography' folder needs to be created.
The frontend file defines the paleogeographic setting and the LIP degassing scenario (mass and duration) as well as the biological parameter space that should be explored concerning the vegetation's dispersal capacity, temperature and aridity niche adaptation rates. 
After model set up, the main model function 'degassing_experiment_main.R' can be run.

The frontend files necessary to reproduce the degassing scenarios from the paper are collected in 'model_exp_frontend_files'. To run a specific scenario, a directory with the specific frontend file and the above mentioned folders needs to be created. 
Frontend files represent the following scenarios: 
A: present-day geography, 400 ppm starting CO2, 10'000 Gt C degassing
B: Paleocene-Eocene geography, 500 ppm starting CO2, 15'000 Gt C degassing 
C: Triassic-Jurassic geography, 800 ppm starting CO2, 30'000 Gt C degassing
D: Permian-Triassic geography, 600 ppm starting CO2, 40'000 Gt C degassing
S1: same as C, but maximum 6-fold weathering enhancement instead of 4-fold 
S2: same as D, but maximum 6-fold weathering enhancement instead of 4-fold
S3: Permian-Triassic geography, 800 ppm starting CO2, 50'000 Gt C degassing
S4: Paleocene-Eocene geography, 500 ppm starting CO2, 8'000 Gt C degassing

For time reasons, it can be convenient to separate the biological parameter space into smaller batches and run them separately. At the moment the frontend file is coded so that 50 parameter combinations are computed in parallel (on a HPC cluster). The post-processing script 'degassing_experiment_analysis.R' can help to append all the model outputs for plotting and exploration, but it is not needed to run the model. It should be noted that the eco-evolutionary vegetation model produces a large number of files, which can limit the number of model experiments that can be run in parallel.

The 'biology folder' containts the 'biology_configuration_file.R' and empty folders for the model output. The configuration file contains a set of rules that define the behavior of the modeled vegetation floras as explained in the methods of the paper and the comments in the code. 

The model makes use of a lookup table of pre-run climate simulations using the PlaSim climate model. For information on how to use the PlaSim model, please see the user guide and the reference manual provided by the developers. A version of the run folder with the namelists used for the production of the climate data used here is provided ('sample_run_folder'). All input files needed for the different paleogeographic configurations can be found in 'PlaSim_inputs'. All climate model outputs used in the vegetation and carbon-cycle modeling are deposited in the 'climate_data' folder.

The folder 'temperature_proxy_analysis' includes all scripts and data needed to reproduce the temperature reconstructions of the Siberian Traps, Central Atlantic Igenous Province and the North Atlantic Igneous Province LIPs based on several geochemical proxy systems. Please see separate README file for further instructions.

In case of questions, feel free to ask: jul.rogger@gmail.com
