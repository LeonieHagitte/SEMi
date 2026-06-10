[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19891930.svg)](https://doi.org/10.5281/zenodo.19891930) 
 
# SEMi
A serious threat to the validity of psychological assessments is differential item functioning (DIF; Bauer, 2020), that is, a situation, in which test items function differently across different groups.
The moderated nonlinear factor analysis (MNLFA) approach allows for simultaneous testing of DIF in multiple categorical and continuous covariates if those covariates are specified a priori and if their DIF is linear.
SEM trees and forests (Brandmaier et al., 2013; Brandmaier et al., 2016) could provide a viable alternative to exploring DIF in latent factor models.
SEM trees allow for non-linear DIF effects and provide a way to systematically search for important covariates (and their interactions).
In this project, we design a simulation study to investigate the conditions under which each of the two methods returns favorable results.
   
## Simulation 
Find the simulation in the "future_build" folder.
Then run the contained files in the following order:
  - dataprep_rt.R
  - analysis_rt.R
  - simulation_rt.R (This script contains the actual driver.)
  - measures_rt.R

Please note that this is currently under construction, e.g. the simulation_rt.R script currently contains two design grids, one that is complete, and one test-grid.
Currently the test-grid overwrites the full grid, because we are still in the process of properly setting up the simulation. 
Later on this will be streamlined. 

Please, feel free to report bugs if you encounter problems, to hagitte@mpib-berlin.mpg.de

## SLURM

The script is intended to be run on a SLURM-based high performance computing environment.
To this end, simulation_rt.R can be started with command line arguments to run only slices of
the simulation. For example
`Rscript simulation_rt.R 5 100`
runs the 5-th slice of the simulation when all conditions are sliced into 100 slices.
To run the all slices in parallel using the SLURM array job function, use the submit_jobs.sh file.
Here are a few useful commands on a SLURM server:

- Start the simulation with all conditions
`./submit_jobs.sh`

- Check how many of your jobs are currently running
`squeue -u "$USER" -t RUNNING -h | wc -l`

-