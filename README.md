# soil-plant-model
This repository contains description of a model that couples stomatal conductance, plant hydraulics, ABA production, and soil water balance to quantify plant responses during a designated dry period. Specifically, the responses are quantified using 
1. hydraulic risk
2. net cumulative carbon assimilation. 
The model does not allow for xylem refilling. For more details of model implementation, please refer to: 

**Feng et al. (2018), The ecohydrological context of drought and classification of plant responses, Ecology Letters.**

This code was written in Python 2.7 and requires the Python packages numpy, scipy, matplotlib, and xlrd.
They can be downloaded through an Anaconda Distribution, which also includes the development environment Spyder. 
See: https://www.anaconda.com/download/ and https://docs.anaconda.com/anaconda/packages/pkg-docs.html (for a list of packages). 

Author contact: feng@umn.edu 
Date: 2018/07/14

## Code 
* soil-plant-model.py - contains the model 
* utility_functions.py - contains functions needed to run the simulations
* simulate_plant_trajectories.py - used to simulate soil moisture, minimum plant water potentials, and net cumulative assimilation over a dry period, with designated mean rainfall frequencies and vapor pressure deficit
* params_soil.py - parameters for different soil types
* params_constants.py - contains constants used in the models

## Data
* hydraulic_traits.xls - contains the trait parameters for *Juniperus monosperma* and *Pinus edulis* used to simulate plant trajectories
