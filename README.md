# Alts
An altruism-based simulation

## Description ##

Alts is a simulation that creates a population of individuals following a series of configurable parameters found in the
.ini files. It's developed to work with altruism, but it can be used for any trait just by coding the genetic model for
the desired characteristic. The altruism models can be found in the models/ directory.
The parameters that can be modified are:
* Initial number of groups
* Initial size of the groups
* Maximum size of a group allowed
* Maximum size of the whole population allowed
* Number of descendants per survivor in a generation
* Migration between groups
* Emigration rate
* Immigration rate
* Immigrants’ phenotype
* Life expectancy normal distribution's mean
* Life expectancy normal distribution's standard deviation
* Survival probability normal distribution's mean
* Survival probability normal distribution's standard deviation
* Relatedness to benefit exponential distribution's factor
* Cost to benefit ratio of altruism
* Minimum benefit of altruism given to unrelated individuals
* Maximum cost willing to be sacrificed by altruists
* Number of genealogical generations
* Number of generations per simulation
* Number of simulations parallelized
* Genetic model with its own parameters:
  * Loci's alleles
  * Inheritance pattern
  * Initial frequencies of each allele
  * Size in base pairs of each locus
  * Mutation rate of each locus
  * Recombination rate of each locus

These parameters allow population to be quite realistic and therefore
to create realistic results adapting the parameters to fit the desired population.

## Workflow ##

The simulation is initialized with the given number of groups and the configures size. The initial state of the
population is saved in the output file.
The workflow followed by the simulation is the following:
1. The survival probabilities of the individuals is reset to the initial value that could have been changed by altruism.
2. Altruists help other individuals with lower survival probability than them.
3. The average initial and modified survival probabilities is calculated.
4. Individuals die or survive based on their survival probability, which represents any circumstance that could lead to
the death of an individual, such as foraging or encounters with predators.
5. Migration at the population level. (Note that the immigrants start with an empty genome, so their sequence is
technically the ancestral, do not use if you are going to use haplotypes data.)
6. Reproduction, which include migration between groups.
7. The data of the generation is saved in a HDF5 file.

[//]: # (<img src="https://user-images.githubusercontent.com/96572489/170871211-d75bca92-f345-4022-9f7d-964af75999aa.png" alt="Alts' workflow" width="500"/>)

`relatedness.py` is the script computes relatedness between individuals.  
`simulator.py` is the main script that runs each simulation. The simulation is initialized and populated with the
initial individuals and then the workflow described above is repeated for each generation. It also creates the output
results files.
`alts.py` manages the multiprocessing when running the program with several simulations in parallel. If the simulation
ended because of lack of individuals or altruists the simulation will start over.
`plotter.py` takes the results from `simulator.py` and generates the configured plots in `config.ini`
* Proportion of individuals of each phenotype at the beginning of each generation.
* Proportion of alleles of each locus at the beginning of each generation.
* Proportion of selfish individuals at the beginning of each generation in each simulation.
* Proportion of selfish alleles at the beginning of each generation in each simulation.
* Altruist to selfish ratio in groups at the beginning of each generation.
* Number of total survivors by the end of the generation
* Number of total survivors by the end of the generation based on altruist/selfish group ratio

Inside the directory `models` the scripts with the altruist step for each genetic model are stored
along with their corresponding configuration file.
Inside the directory `tests` the pytest files to check the correct operation of the program can be found.

## Installation
### Linux (Debian-based)
1. Install Python:
```bash
sudo apt-get install python3
```
2. Install dependencies:
    ```bash
    pip install h5py
    pip install matplotlib
    pip install numpy
    pip install scipy
    pip install pytest # If you want to run the tests
    ```
3. Download [project files](https://github.com/NoahJorCal/Alts/archive/refs/heads/main.zip).
4. Unzip the file.

### Windows
1. Install Python through [Microsoft Store](https://www.microsoft.com/store/productId/9PJPW5LDXLZ5).
2. Install dependencies:
    ```bash
    pip install h5py
    pip install matplotlib
    pip install numpy
    pip install scipy
    pip install pytest # If you want to run the tests
    ```
3. Download [project files](https://github.com/NoahJorCal/Alts/archive/refs/heads/main.zip).
4. Unzip the file.

## Use
### Configure
The general configuration can be found in `config.ini` file.

The loci and model configuration can be found in the `.ini` file in the model directory

### Run
To run the simulation:
```bash
python3 alts.py --directoy <output_directory> --output <output_files_name> --cpu <cpu_number> 
```
Use ```--quiet ``` to run the script without printing feedback.
Use ```--seed <seed> ``` to set a seed for the simulations.
Use ```--help ``` for help.

To plot results:
```bash
python3 plotter.py --input <input_file_name>
```
Use ```--show ``` to show the plots as popup windows.
Use ```--no-save ``` to prevent the script to save the plots.
Use ```--help ``` for help.
