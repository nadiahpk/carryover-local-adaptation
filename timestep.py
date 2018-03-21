import random as random

#import sys
#sys.path.append("../../current_code/")
from carryover import noOffspringFnc
from carryover import gene2bin
from carryover import offspringGenotypeFnc
from carryover import dispFnc
from carryover import compnSimpleFnc

def timestep(parameters, ecosystem, landscape):

    # reproduction and dispersal

    for locn, flock in enumerate(ecosystem):

        habType = landscape[locn]
        adults = flock['adults']
        noOffspring = noOffspringFnc(parameters, adults, habType)

        if noOffspring > 0: # create and disperse each offspring

            # rewrite mum and dads genotypes as binary strings (saves time to do it first) into a dictionary of the form 
            #  {polygenes type: (mum's polygenes as binary string, dad's polygenes as binary str)}
            parGenotypeBin = {
                    geneType: ( gene2bin(parameters, adults[0]['genotype'][geneType],geneType), gene2bin(parameters, adults[1]['genotype'][geneType],geneType) )
                    for geneType in parameters['genetics'] }

            for cnt in range(noOffspring):

                # create offspring
                offspring = {
                        'sex': random.choice( parameters['sexes'] ),
                        'natalHabType': habType,
                        'genotype': offspringGenotypeFnc(parameters, parGenotypeBin)
                        }

                # disperse offspring
                newLocn = dispFnc(parameters, offspring, locn, landscape)
                ecosystem[newLocn]['juveniles'].append(offspring)

    # survival

    for flock in ecosystem:

        flock['adults'] = list() # all adults assumed to die
        # all juveniles live

    # competition

    for flock in ecosystem:

        flock = compnSimpleFnc(parameters, flock) # rearranges each flock according to competition process

    return ecosystem, landscape


