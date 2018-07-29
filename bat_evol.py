from simuOpt import setOptions
setOptions(quiet=True)
import simuPOP as sim
from simuPOP.utils import importPopulation

pop = importPopulation('GENEPOP', 'Genpop_2subpops.txt')

import random

pop.addInfoFields(['age'])

# age_cat is not imported as they are defined by age 
# 0, 1, 2 as juvenile and 2-above as adult
with open('Fstat_dat(age,agecat,sex).txt') as fs:
    for idx, line in enumerate(fs):
        if idx < 10:
            # skip the first few lines
            continue
        if line.startswith('1 '):
            age, age_cat, sex = line.strip().split()[-3:]
            pop.individual(idx-18).setSex(sim.MALE if sex == '2' else sim.FEMALE )
            pop.individual(idx-18).age = int(age)
            #pop.individual(idx-18).age_cat = int(age_cat)
        elif line.startswith('2 '):
            age_cat, sex = line.strip().split()[-2:]
            pop.individual(idx-18).setSex(sim.MALE if sex == '2' else sim.FEMALE )
            # no age info
            pop.individual(idx-18).age = random.randint(0, 17)
            #pop.individual(idx-18).age_cat = int(age_cat)

pop.setVirtualSplitter(
    sim.CombinedSplitter([
        sim.SexSplitter(),
        sim.InfoSplitter(field='age', cutoff=[3, 17] ),
    ]))

    
# population 0, from 68 female to 68 female + 62 male
pop.resize([130, 338], propagate=True)
# set sex
for idx in range(68, 130):
    pop.individual(idx, 0).setSex(sim.MALE)
# add 110-37=73 male
for idx in range(192, 192+73):
    pop.individual(idx, 1).setSex(sim.MALE)
for idx in range(192+73, 338):
    pop.individual(idx, 1).setSex(sim.FEMALE) 


migr = sim.Migrator(
    rate=[
        [0, 0.1 ],
        [0, 0.001 ],
        [0.1, 0],
        [0.001, 0],
       ],
    mode = sim.BY_PROBABILITY,
    subPops=[(0, 'Male'), (0, 'Female'), (1, 'Male'), (1, 'Female')],
    toSubPops=[0, 1],
)

# let us test the migration
pop.addInfoFields('migrate_to')

def demoModel(gen, pop):
    if gen < 50:
        sz = 130, 550
    elif gen < 120:
        sz = 80, 550
    else:
        sz = min(150, 80*1.2**(gen-120)), 550
    print(f"{gen} young/adult sp0 : {pop.subPopSize([0, 'age < 3'])} / {pop.subPopSize([0, '3 <= age < 17'])}" +
          f" sp1 : {pop.subPopSize([1, 'age < 3'])} / {pop.subPopSize([1, '3 <= age < 17'])}")
    return sz

    
def exportPop(pop):
    ''' Export in genepop format with sex and age'''
    with open(f'pop_{pop.dvars().gen}.txt', 'w') as output:
        #
        # first line: title
        #
        output.write(f'bat population at generation {pop.dvars().gen}\n')
        #
        # second line: allele names
        #
        names = pop.lociNames()
        if names:
            # if names are specified
            output.write(', '.join(names) + ', sex, age, fitness\n')
        else:
            names = []
            for ch in range(pop.numChrom()):
                for loc in range(pop.numLoci(ch)):
                    names.append('ch%d-loc%d' % (ch + 1, loc + 1))
            output.write(', '.join(names) + ', sex, age, fitness\n')
        # sex and age as special names
        # 
        # output genotype
        #
        # progress bar might be wrong with subPops parameter...
        alleleWidth = 3 if max(pop.genotype()) >= 99 else 2
        format_string = '%%0%dd%%0%dd' % (alleleWidth, alleleWidth)
        count = 0
        numLoci = pop.totNumLoci()
        for vsp in [0, 1]:
            # 
            # for each subpopulation, output pop
            #
            output.write('POP\n')
            # the name might contain space etc
            name = ''.join([x for x in pop.subPopName(vsp) if x.isalnum()])
            if not name:
                name = 'SubPop%d' % (vsp if type(vsp) == type(0) else vsp[0])
            #
            for idx, ind in enumerate(pop.individuals(vsp)):
                #
                # label
                #
                output.write('%s-%d, ' % (name, idx + 1))
                #
                # genotype
                #
                geno = ind.genotype()
                output.write(' '.join([format_string % (geno[x], geno[numLoci + x]) for x in range(numLoci)]) +
                        f'\t{"M" if ind.sex() == sim.MALE else "F"}\t{ind.age}\t{ind.fitness}\n')
        return True



pop1 = pop.clone()
pop1.addInfoFields('fitness')
pop1.evolve(
    preOps=[
        migr,
        sim.InfoExec('age += 1'),
        sim.InfoExec('fitness = (50 - ind.age)/50', exposeInd='ind'),
    ],
    matingScheme=sim.HeteroMating([
        # only adult individuals with age >=3 will mate and produce
        # offspring. The age of offspring will be zero.
        sim.RandomMating(ops=[
            sim.MendelianGenoTransmitter()],
            subPops=[(sim.ALL_AVAIL,'3 <= age < 17')],
            weight=-0.1),
        # individuals with age < 17 will be kept, but might be removed due to
        # population size decline
        sim.CloneMating(subPops=[(sim.ALL_AVAIL, 'age < 3'), (sim.ALL_AVAIL, '3 <= age < 17')]),
        ],
        subPopSize=demoModel),
    postOps=[
        sim.Stat(popSize=True),
        sim.PyEval(r'f"{gen} {subPopSize}\n"'),
        sim.PyOperator(func=exportPop, step=10),
    ],
    gen=200
)
