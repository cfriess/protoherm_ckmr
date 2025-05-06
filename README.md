
# Application of Close-kin mark–recapture (CKMR) to Protogynous Fishes

## Description

This repository contains the code and data used for a simulation study applying
CKMR to the protogymous hermaphrodite gag, *Mycteroperca microlepis*. The study
evaluates the performance of a newly developed CKMR estimator for male and
female abundance under uncertainty regarding male contribution to reproductive 
success, variable sampling, and model misspecification.

## Overview

We use renv and snakemake. The full simulation scenarios are defined in the snakefile,
but the simulation can also be run without snakemake, by specifying desired parameters
in the main OM and EM simulation scripts.

All scripts are in the 'code' folder which is further organized into 'admb' (where the
ADMB .tpl file lives), 'main_sim' which houses the main simulation scripts, 'analysis'
where the scripts to analyze result are located, and 'helpers' where many functions
are defined.

The data folder consists of 'inputs' where simulation input data files live, 
'intermediate' where data objects are stored that are created by the simulation and
later reused, 'sim_peds' where OM results needed to run the EM are saved,
'admb_run_directory' where ADMB model run files are saved, and 'results' where EM
results are saved.

### To run the Operating Model:

Run code in code/main_sim: 'ckmr_sim.R'. This script sources the following two
scripts in the same folder:

- sim_inputs.R defines the life history and fishery data inputs for the OM
- ckmr_setup.R initializes the population and generates some additional data
objects for the OM

The OM parameters for running the simulation without snakemake are defined in 
lines 67 to 72 of 'ckmr_sim.R'. Additional options that we ended up not using
after preliminary runs are in lines 79 to 83.

The final section of the 'ckmr_sim.R' script defines which data objects are 
saved in 'results/sim_peds'. These are then read back in in the 
'estimate-SMK.R' script.

### To run the Sampling Model and Estimation Model:

Run code in code/main_sim: estimate-SMK.R.

EM parameters are defined in lines 80 to 114 of the 'estimate-SMK.R' script,
if running the simulation without snakemake.

OM results are first read back in, then samples are generated and written out
to data files for ADMB. ADMB is then called from within R and, if converged, results from 
the ADMB report file are read in, processed and written out to 'data/results'. ADMB
files for the most recent run are stored in the data/admb_run_directory for
diagnostic purposes.

### To analyze results:

Run code in code/analysis: 'process_results_SMK.R' for results generated using
snakemake and harvested with the 'harvest_outputs.R' script in code/helpers

The 'process_results_SMK.R' generates several .rds files that are then used for
generating summary plots. The results files from this study for the full
simulation are saved in data/results/full_sim_summary_results.

### To plot results:

Run code in the code/analysis: 'results_figures.R'

This script reads in the summary results .rds files from the 
data/results/full_sim_summary_results folder and generates the plots in the
paper.
