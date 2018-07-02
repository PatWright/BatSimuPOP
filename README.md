# BatSimuPOP
# Simulation of population changes on bat colonies over time to assess the use of temporal sampling

# Input file:
1 Genepop file with 2 subpopulations (total number of individuals in genepop is inferior to the number of individuals in the simulation)

subpop 1: local colony that will undergo changes (additional info on age, age category and sex).

subpop 2: Rest of the population. Present to maintain gene flow over time. (additional info on sex and whether juvenile or adult)

Information on sex, age in years (for 1 subpop) or age category (juvenile-adult) is also available (file in fstat format 14 microsatellite loci + additional info)


# Life cycle of species to consider for overlapping generations:
Maxlifespan ~ 17 years old

minMatingAge ~ 2 years old

maxMatingAge ~ 17 years old

Number of young per year ~ 1

1 generation ~ 5 years

# Other settings:

Random mating in subpopulations

Maintain gene flow between both sub-populations using Migrator. This is sex-biased as females remain in the same colony and males disperse.

Female migration rate	0.001

Male migration rate	0.1


# Example scenario:

Simulation of 20 generations (~200 years)

The initial population size is based on Ne estimates from the Genepop file (N genotypes < Actual population size), for example:

N0subpop1=130

N0subpop2=550

Subpop 2 remains constant over time

Subpop 1 undergoes an instant population change (InstantChangeModel) at generation 5 and declines to 80.

Then at generation 12 it undergoes a linear growth at rate 0.2 to a population of 150 (LinearGrowthModel) => Need of MultiStageModel 

# Outputs:

Need Genepop file outputs for each generation (or every 5 years) after breeding to include juveniles. The effective population size estimates will be done separately.

Output info on age and sex of individuals would also help with Ne estimates.

