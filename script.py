# a run in which all individuals are started with the same gene that has the highest reproduction in the majority L habitat type

from simulate import simulate
import random as random

# landscape
landscape = 21*'L' + 9*'H' + 21*'L' # read as 21 LC territories, 9 HC territories, 21 HC territories

idxRun = 0
tf = 600        # number of timesteps
burnInT = 0     # length of burn-in (not stored in pickle result file)
suffix = '_1' # suffix of file name for pickle results file

parameters = {
        # Sex
        'sexes': ( 'h', 'h' ), # sexes reqd for mating pair, could also be set ( 'm', 'f' ) for non-hermaphroditic
        # Genetics
        'genetics': {
               'repn': {'noLoci': 20, 'maxPhen': 2,  'minPhen': -2,  'isInt': False, 'pMut': 0.001},    # reproduction genes
               'neut': {'noLoci': 20, 'maxPhen': 1,  'minPhen': -1,  'isInt': False, 'pMut': 0.001},    # neutral genes
               #'pref': {'noLoci': 20, 'maxPhen': 10, 'minPhen': -10, 'isInt': False, 'pMut': 0.001},   # habitat preference genes
               #'phil': {'noLoci': 20, 'maxPhen': 10, 'minPhen': -10, 'isInt': False, 'pMut': 0.001},   # NHPI genes
               #'dist': {'noLoci': 20, 'maxPhen': 25, 'minPhen':   0, 'isInt': True,  'pMut': 0.001},   # dispersal distance genes
                },
        'distMax': 7, # set to some value if 'dist' genes not specified, else set to None
        # Habitat types and reproduction
        'habitats': {
                 'L': {'rMax': 10, 'phenOpt': -1, 'sd': 1.11},
                 'H': {'rMax': 10, 'phenOpt':  1, 'sd': 1.11},
                },
        # Natal habitats and competition
        'competition': {
                 'L': 1,
                'H': 10,
                #'H': 1,
                }
        }

# specify initial ecosystem

# this reproduction genotype maximises reproduction in LC
g_bin = '000000000000000011111'
g_int = int( g_bin, 2 ) # = 31

# random start for the rest of the genotype
geneFnc = lambda geneType: g_int if geneType == 'repn' else random.randint(0, 2**parameters['genetics'][geneType]['noLoci']-1)

initial_ecosystem = [ {
            'adults': [
            {'sex': 'h', 'natalHabType': habType, 'genotype': {geneType: geneFnc(geneType) for geneType in parameters['genetics'].keys()}},
            {'sex': 'h', 'natalHabType': habType, 'genotype': {geneType: geneFnc(geneType) for geneType in parameters['genetics'].keys()}} ],
            'juveniles': list(),
            }
    for habType in landscape]

simulate(parameters, landscape, burnInT, tf, initial_ecosystem, None, suffix)
