# Alts
An altruism-based simulation

## Description ##

Alts is a simulation that creates a population of individuals following a series of configurable parameters found in the .ini files. It's developed to work with altruism but it can be used for any trait just by coding the genetic model for the desired characteristic. The altruism models can be found in the folder models. The parameters that can be modified are:
Population size
* Emigration rate
* Immigration rate
* Immigrants’ phenotype
* Lifespan
* Reproduction
* Size of isolated groups
* Survival probability range
* Altruism probability
* Altruism cost
* Altruism benefit
* Generations per simulation
* Number of simulations
* Genetic model
* Inheritance pattern
* Initial frequencies of each gene
* Mutation rate of each gene
These parameters allow population to be quite realistic and therefore to create realistic results adapting the parameters to fit the desired population.

The workflow followed by the simulation is the following:
1. Every generations starts with a fixed size population.
2. Each individual is assigned a random survival probability from a pre-established range.
3. Individuals are  grouped into random isolated groups at the indicated group size.
4. The altruistic individuals have a probability of sacrificing part of their survival probability to increase that of another, the recipient will be altruistic or selfish depending on the selected model.
5. The selection event occurs, which represents any circumstance that may face the individuals, such as foraging or encounters with predators.
6. The death or survival of an individuals is calculated based on the survival probability.
7. The survivors reproduce until the population size reaches the stablished number.
8. The individuals whose age have reached the lifespan die.
9. The data of the generation is saved and a new generation starts.

<img src="https://user-images.githubusercontent.com/96572489/170871211-d75bca92-f345-4022-9f7d-964af75999aa.png" alt="Alts' workflow" width="500"/>

## Instalation ##
This simulation works with Python3. The depedencies needed for this simulation are:
* configparser
* matplotlib (only for plotting)



