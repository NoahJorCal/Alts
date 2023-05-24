# Alts
An altruism-based simulation

## Description ##

Alts is a simulation that creates a population of individuals following a series of configurable parameters found in the .ini files. It's developed to work with altruism, but it can be used for any trait just by coding the genetic model for the desired characteristic. The altruism models can be found in the folder models. The parameters that can be modified are:
* Population size
* Emigration rate
* Immigration rate
* Immigrants’ phenotype
* Lifespan
* Reproduction
* Size of isolated groups
* Survival probability range
* Altruism probability
* Relatedness requirement
* Altruism cost
* Altruism benefit
* Generations per simulation
* Number of simulations
* Genetic model which has its own parameters:
  * Inheritance pattern
  * Initial frequencies of each gene
  * Mutation rate of each gene  

These parameters allow population to be quite realistic and therefore to create realistic results adapting the parameters to fit the desired population.

## Workflow ##

The workflow followed by the simulation is the following:
1. Every generation starts with a fixed size population.
2. Each individual is assigned a random survival probability from a pre-established range.
3. Individuals are grouped into random isolated groups at the indicated group size.
4. The altruistic individuals have a probability of sacrificing part of their survival probability to increase that of another, the recipient will be altruistic or selfish depending on the selected model.
5. The selection event occurs, which represents any circumstance that may face the individuals, such as foraging or encounters with predators.
6. The death or survival of an individuals is simulated based on the survival probability.
7. Immigrants join the population if immigration is configured, then the survivors from the previous generation reproduce until the population size reaches the established number.
8. The individuals whose age have reached the lifespan die.
9. The data of the generation is saved, and a new generation starts.

<img src="https://user-images.githubusercontent.com/96572489/170871211-d75bca92-f345-4022-9f7d-964af75999aa.png" alt="Alts' workflow" width="500"/>

_pedigree.py_ is the module that manages the relationship between individuals and computes relatedness.  
_simulator.py_ is the main script that runs each simulation. The simulation is initialized and populated with the initial individuals and then the workflow described above is repeated for each generation.   
_alts.py_ manages the multithreading when running the program with several simulations at once. The results and the mean of all the simulations are saved and serialized.  
_plotter.py_ takes the results from _alts.py_ and generates three plots:
* Stackplot of the proportion of each phenotype at the beginning of the simulation.
* Line plot of survivors per generation by phenotype.
* Line plot of the proportion of selfish individuals in every simulation.

Inside the folder _models_ the scripts with the selection event for each genetic model are stored, along with their corresponding configuration file.

## Installation
### Linux (Debian-based)
1. Install Python:
```bash
sudo apt-get install python3
```
2. Install dependencies:
```bash
pip install numpy
pip install configparser
pip install matplotlib
pip install scipy
pip install h5py
pip install pandas
pip install multiprocessing
```
3. Download [project files](https://github.com/NoahJorCal/Alts/archive/refs/heads/main.zip).
4. Unzip the file.

### Windows
1. Install Python through [Microsoft Store](https://www.microsoft.com/store/productId/9PJPW5LDXLZ5).
2. Install dependencies:
```bash
pip install numpy
pip install configparser
pip install matplotlib
pip install scipy
pip install h5py
pip install pandas
pip install multiprocessing
```
3. Download [project files](https://github.com/NoahJorCal/Alts/archive/refs/heads/main.zip).
4. Unzip the file.

## Use
### Configure
The general configuration can be found in `config.ini` file.

The gene and model configuration can be found in the `.ini` file in the model folder

### Run
To run the simulation:
```bash
python3 alts.py --outfile <output_file_name> --cpu <cpu_number>
```
To represent results graphically:
```bash
python3 plotter.py --input <input_file_name>
