import random
from math import exp
import itertools as it
from bisect import bisect
import numpy as np

def gene2bin(parameters, gene, geneType):
    """
    geneBin = gene2bin(parameters, gene, geneType)

    Turns genes into the binary-string form

    gene:
        integer
    geneType:
        string
    geneBin:
        string, of all ones and zeros

    >>> geneType = 'repn'
    >>> parameters = {'genetics': {geneType: {'noLoci': 20} } }
    >>> gg = gene2bin(parameters, 0, geneType)
    >>> len(gg) # see it's the same length as the number of loci
    20
    >>> gg # and it is the number 0 in binary of length 20
    '00000000000000000000'
    >>> gene2bin(parameters, 50,geneType) # another example
    '00000000000000110010'
    """

    return format(gene, '0' + str(parameters['genetics'][geneType]['noLoci']) + 'b')

def gene2phen(parameters, gene, geneType):
    """
    phen = gene2phen(parameters, gene, geneType)

    Turns the genes into a phenotype value.

    gene:
        integer
    geneType:
        string
    phen:
        value

    >>> geneType = 'repn' # choose the reproduction genes type
    >>> # they have 20 alleles each, with a maximum phenological value of 1, and a minimum of -1
    >>> parameters = {'genetics': {geneType: {'noLoci': 20, 'maxPhen': 1, 'minPhen': -1, 'isInt': False} } }
    >>> gene = 0 # therefore if I take the 0th gene, which has all 0 and no 1 alleles
    >>> gene2phen(parameters, gene, geneType) # then its phenological value will be negative 1
    -1.0
    >>> parameters = { 'sexes': ( 'm', 'f' ), 'genetics': { 'repn': {'noLoci': 20, 'maxPhen': 2,  'minPhen': -2,  'isInt': False, 'pMut': 0.001} }, 'habitats': { 'L': {'rMax': 10, 'phenOpt': -1, 'sd': 1.1}, 'H': {'rMax': 10, 'phenOpt':  1, 'sd': 1.1}, }, 'competition': { 'L': 1, 'H': 10, } }
    >>> parameters = {'genetics': {geneType: {'noLoci': 20, 'maxPhen': 2, 'minPhen': -2, 'isInt': False} } }
    >>> gene2phen(parameters, 31, 'repn')
    -1.0
    """

    maxPhen = parameters['genetics'][geneType]['maxPhen']
    minPhen = parameters['genetics'][geneType]['minPhen']
    noLoci = parameters['genetics'][geneType]['noLoci']

    if parameters['genetics'][geneType]['isInt']:
        phenotype = round(minPhen + (bin(gene).count('1') / noLoci) * ( maxPhen - minPhen ))
    else:
        phenotype = minPhen + (bin(gene).count('1') / noLoci) * ( maxPhen - minPhen )

    return phenotype

def randIdxWeights(weights):
    """
    idx = randIdxWeights(weights)

    Return a randomly chosen index where choice is weighted.
    Special because weights can sum to whatever, don't have to be probabilities.
    If all weights are zero, it will return the index that is 1 + length of weights list.

    weights:
        iterable, containing weightings to apply to each index
    idx:
        integer, the random index chosen according to weights
    """

    cumWeights = list( it.accumulate(weights) )
    idx = bisect( cumWeights, random.random()*cumWeights[-1] )

    return idx

def noOffspringFnc(parameters, adults, habType):
    """
    noOffspring = noOffspringFnc(adults, habType)

    Given adults (e.g. [mum, dad]), and the habitat type they reside on, use the Gaussian fitness function
    to calculate how many offspring they'll have.

    adults:
        two-element list or similar, with each element a dictionary, describing an individual
        e.g. {'genotype': {'disp': 274226, 'dist': 133608, 'repn': 951436}, 'natalHabType': 'L', 'sex': 'f'}
    habType:
        string, describes habitat type the pair reside on
        e.g. 'H' is a habitat type that confers high competitive ability and 'L' confers low
    noOffspring:
        integer, how many offspring they have


    >>> parameters = { 'sexes': ( 'm', 'f' ), 'genetics': { 'repn': {'noLoci': 20, 'maxPhen': 2,  'minPhen': -2,  'isInt': False, 'pMut': 0.001} }, 'habitats': { 'L': {'rMax': 10, 'phenOpt': -1, 'sd': 1.1}, 'H': {'rMax': 10, 'phenOpt':  1, 'sd': 1.1}, }, 'competition': { 'L': 1, 'H': 10, } }
    >>> habType = 'L' # reside on low competitive ability conferring habitat type
    >>> parameters['habitats'][habType]['phenOpt'] # where the optimal phenotype is
    -1
    >>> parameters['genetics']['repn']['minPhen'] # the minimum phenotype when all alleles are 0 is
    -2
    >>> parameters['genetics']['repn']['maxPhen'] # the maximum phenotype when all alleles are 1 is
    2
    >>> parameters['genetics']['repn']['noLoci'] # and the number of alleles is
    20
    >>> # therefore to have the optimal number of offspring, parents need 0.25 of alleles being 1, or 5
    >>> int('0'*(20-5) + '1'*5,2)
    31
    >>> mum = {'genotype': {'repn': 31}} # define a minimal mum and dad w only reproduction genes
    >>> dad = {'genotype': {'repn': 31}} # set to 0 so will have the minimum phenotype
    >>> adults = [mum, dad]
    >>> noOffspringFnc(parameters, adults,habType)
    10
    """

    if len(adults) != 2:

        noOffspring = np.nan

    else:

        # get the reproductive phenotype of pair

        phen0 = gene2phen(parameters, adults[0]['genotype']['repn'], 'repn')
        phen1 = gene2phen(parameters, adults[1]['genotype']['repn'], 'repn')
        phenPair = (phen0+phen1)/2

        # get the attributes of mum and dad's habitat type that relate to reproduction

        rMax = parameters['habitats'][habType]['rMax']
        phenOpt = parameters['habitats'][habType]['phenOpt']
        sd = parameters['habitats'][habType]['sd']

        # calculate the number of offspring the pair will have based on the Gaussian function

        noOffspring = round( rMax * exp( - (phenPair - phenOpt)**2 / (2*sd**2) ) )

    return noOffspring

# parGenotypeBin looks like: {polygenes type: (mum's polygenes as binary string, dad's polygenes as binary str)}
def offspringGenotypeFnc(parameters,parGenotypeBin):
    """
    offGenotypeBin = offspringGenotypeBinFnc(parameters, parGenotypeBin)

    Accepts the parents' genotype in binary format and returns offspring's genotype in integer format

    parGenotypeBin:
        dictionary, keys are genetypes and values are tuples of mum and dad genes as a binary string
        i.e.  {polygenes type: (mum's polygenes as binary string, dad's polygenes as binary str)}
        e.g. {'disp': ('01000010111100110010', '10110110110100100010'),
              'dist': ('00100000100111101000', '10101110100011000110'),
              'repn': ('11101000010010001100', '00101110110001101110')}
    offGenotype:
        dictionary, keys are gene types and values are genes in integer format
        e.g. {'disp': 591984, 'dist': 792214, 'repn': 378862}
    """

    # create new genotype for offspring by randomly choosing genes from mum and dad

    offGenotypeBinList = { geneType: list(map( lambda locus: random.choice(locus), zip(mumGeneBin, dadGeneBin) ))
            for geneType, (mumGeneBin,dadGeneBin) in parGenotypeBin.items() } # retain as list for now so mutations easier

    # find positions to mutate

    mutposns = list()

    for geneType in offGenotypeBinList.keys(): # for each polygenes type

        pMut = parameters['genetics'][geneType]['pMut'] # get its mutations probability
        noLoci = parameters['genetics'][geneType]['noLoci'] # get the number of loci

        p = [ i for i in range(noLoci) if random.random() < pMut ] # create list of locus positions at which to flip alleles

        if p: # if any positions to be modified for these genes type

            mutposns.append( (geneType,p) ) # add them to the list


    # apply mutations

    if mutposns: # if any mutations to apply

        for geneType, p in mutposns: # for each genes type to apply it to

            for pp in p: # for each position to be modified

                # flip the allele at that position
                if offGenotypeBinList[geneType][pp] == '0':

                    offGenotypeBinList[geneType][pp] = '1'

                else:

                    offGenotypeBinList[geneType][pp] = '0'

    # create the genotype of the offspring by joining those lists of alleles and converting to integers

    offGenotype = { geneType: int( ''.join(geneBinList), 2 ) for geneType, geneBinList in offGenotypeBinList.items() }

    return offGenotype

def dispFnc(parameters, offspring, locn, landscape):
    """
    newLocn = dispFnc(parameters, offspring, locn, landscape)

    Accepts an offspring and its location and finds its new location after dispersal

    offspring:
        dictionary, describing an individual e.g. {'genotype': {'disp': 274226, 'dist': 133608, 'repn': 951436}, 'natalHabType': 'L', 'sex': 'f'}
    locn:
        integer, an index to a location in the landscape
    landscape:
        string, describes the habitat types in the landscape e.g. 'LLLLLLLLLLLLLLLLLLLLLLHHHHHHHHHLLLLLLLLLLLLLLLLLLLL'
    """

    # find the maximum dispersal distance of the offspring
    distMax = parameters['distMax']

    if distMax is None:

        # will need to use offspring's genotype to find its dispersal distance
        distMax = gene2phen( parameters, offspring['genotype']['dist'], 'dist' )

    # find the new location it disperses too

    offGenotype = offspring['genotype']

    if ('pref' not in offGenotype) and ('phil' not in offGenotype): # assume random dispersal

        newLocn = ( locn + random.randint(-distMax,distMax) ) % len(landscape)

    else: # has genes controlling habitat type preferences

        # find weighting (weight) this offspring gives to preferred habitat type (prefHabType)

        if 'pref' in offGenotype: # preference by habitat type

            # find the phenotype value expressing the preference
            phen = gene2phen( parameters, offGenotype['pref'], 'pref' )

            # by convention, a positive value phenotype prefers 'H' and negative prefers 'L'

            # identify preferred habitat type
            if phen < 0:
                prefHabType = 'L'
            else:
                prefHabType = 'H'


        else: # prefence by NHPI

            # find the phenotype value expressing the preference
            phen = gene2phen( parameters, offGenotype['phil'], 'phil' )

            # by convention, a positive value phenotype prefers natal habitat type and negative prefers non-natal habitat type

            # identify preferred habitat type
            if phen < 0:
                prefHabType = 'L' if offspring['natalHabType'] == 'H' else 'H'
            else:
                prefHabType = offspring['natalHabType']

        # convert phenotype value to weight
        weight = 1+abs(phen) # > 1 because non-preferred habitat type has weighting 1, "has weight times the probability ..."

        # choose new location using the weight and prefHabType

        # get the locations and habitat types of the neighbourhood around it to which it may disperse given its dispersal distance
        lenLandscape = len(landscape)
        neighbourHabTypes = [ landscape[ i % lenLandscape ] for i in range(locn-distMax, locn+distMax+1) ]
        neighbourLocns    = [ i % lenLandscape for i in range(locn-distMax, locn+distMax+1) ]

        # find out how strongly it weights each neighbouring habitat type
        neighbourWeights = [ weight if habType == prefHabType else 1 for habType in neighbourHabTypes ]

        # use preference weighting of neighbouring locations to choose a new location
        newLocn = neighbourLocns[ randIdxWeights(neighbourWeights) ]

    return newLocn

def compnSimpleFnc(parameters, flock):
    """
    A simple competition function in which one juvenile of each mating-pair sex
    becomes a new adult in the flock and all other juveniles die.
    Winner found by random weighted choice, where weighting determined by natal habitat type
    """

    # find out which positions are open for this flock's mating pair

    openPositions = [ x for x in parameters['sexes'] ]

    for a in flock['adults']:

        openPositions.remove(a)

    # for each position that is open in the mating pair, choose a juvenile to take it

    for sex in openPositions:

        # create competition weights for the juveniles, where a 0 is assigned if the juvenile does not match the specified sex
        compnWeights = [ parameters['competition'][ competitor['natalHabType'] ] if competitor['sex'] == sex else 0 for competitor in flock['juveniles'] ]

        if any( w > 0 for w in compnWeights ): # if any of the competitors can be chosen

            # determine the winner using the competition weights
            winner = flock['juveniles'][ randIdxWeights(compnWeights) ]

            # remove winner from juvenile flock and add to adult list
            flock['juveniles'].remove(winner)
            flock['adults'].append(winner)

    # for this run we assume all remaining juveniles die
    flock['juveniles'] = list()

    return flock

def ecosystemIsEmpty(ecosystem):
    """
    Checks to see that there is at least one mating pair in the ecosystem
    If there is at least one mating pair, returns False, or else returns True
    """

    returnValue = True

    for flock in ecosystem:

        if len(flock['adults']) == 2: # having mating pair

            returnValue = False
            break

    return returnValue

def phenInSpace(parameters, ecosystem):
    """
    Returns a dictionary with keys geneType and values as a list of phenotypes
    Note: will only return value if location has a mating *pair*
    """

    phenDict = { geneType: list() for geneType in parameters['genetics'] }

    for flock in ecosystem:

        adults = flock['adults']

        if len(adults) == 2: # has a mating pair

            adult = random.choice(adults) # choose a random adult

            for geneType in parameters['genetics']: # calculate each phenotype value and append to the list of that phenotypes in space
                phenDict[geneType].append( gene2phen( parameters, adult['genotype'][geneType], geneType ) )

        else: # no mating pair, so append nans

            for geneType in parameters['genetics']:
                phenDict[geneType].append( np.nan )

    return phenDict

def gendiffInSpace(parameters, ecosystem):
    """
    Returns a dictionary with keys geneType and values as a list of the number of alleles that differ between mating partners

    """

    gendiffDict = { geneType: list() for geneType in parameters['genetics'] }
    for flock in ecosystem:

        adults = flock['adults']

        if len(adults) == 2: # has a mating pair

            # calculate the genetic difference between them and append
            for geneType in parameters['genetics']:
                gendiffDict[geneType].append( genediffFnc( parameters, adults[0]['genotype'][geneType], adults[1]['genotype'][geneType], geneType ) )

        else: # no mating pair, so append nans

            for geneType in parameters['genetics']:

                gendiffDict[geneType].append( np.nan )

    return gendiffDict

def pcolormeshCorrectionXY(xV, yV):
    """
    Allows pcolormesh to plot the axes correctly

    e.g. if zM[yi,xi] = f(xV[xi],yV[yi]) then can use
    plt.pcolormesh(x2V,y2V,zM)
    """

    dx = xV[1]-xV[0]
    lx = len(xV)
    xmin = xV[0]; xmax = xV[-1]

    dy = yV[1]-yV[0]
    ly = len(yV)
    ymin = yV[0]; ymax = yV[-1]

    x2V = np.linspace( xmin-dx/2, xmax+dx/2, lx+1 )
    y2V = np.linspace( ymin-dy/2, ymax+dy/2, ly+1 )

    return x2V, y2V

def genediffFnc(parameters, gene0, gene1, geneType):

    geneBin0 = gene2bin( parameters, gene0, geneType ) # each gene as binary string
    geneBin1 = gene2bin( parameters, gene1, geneType )

    return sum(1 for g0, g1 in zip(geneBin0, geneBin1) if g0 != g1) # count differences and return

if __name__ == "__main__":

    import doctest
    doctest.testmod()

    # run doctests like this:
    # python3 -m doctest -v carryover.py




