import random as random
import pickle
import os # 
import copy

from timestep import timestep # locally defined timestep function

#import sys
#sys.path.append("../../current_code/")
from carryover import ecosystemIsEmpty # checks if the population has gone extinct

# allows me to construct suffixes for files according to parameter values
def parameters2filesuffix(tf, landscape, parameters):

    w = str( parameters['competition']['H'] )

    distMax = parameters['distMax']
    if distMax is None:
        d = 'x'
    else:
        d = str(distMax)

    H = sum( c == 'H' for c in landscape )
    L = len(landscape)

    sd = round( 10*parameters['habitats']['H']['sd'] )

    pMut = round( 1000*parameters['genetics']['repn']['pMut'] )

    nL = parameters['genetics']['repn']['noLoci']

    geneticsChars = { 'repn': 'r', 'pref': 'h', 'phil': 'p', 'dist': 'd', 'neut': 'n' }
    genetics = ''.join(sorted([ geneticsChars[geneType] for geneType in parameters['genetics'] ]))

    suffix = str(tf) + '_w' + str(w) + '_d' + str(d) + '_H' + str(H) + '_L' + str(L) + '_sd' + str(sd) + '_pMut' + str(pMut) + '_nL' + str(nL) + '_' + genetics

    return suffix

def simulate(parameters, landscape, burnInT, tf, initial_ecosystem = None, idxRun=None, suffix=None):
    '''
    parameters: 
        dictionary, see script.py for example
    landscape: 
        string, defines the landscape as a string of H and L for high and low-quality habitats
    burnInT: 
        integer, some number of timesteps (generations) to not store in the pickle file
    tf:
        integer, number of timesteps
    idxRun: 
        integer, when multiple runs on the same parameter-set/suffix are run, 
        may want to separate pickled results file by additional suffix, 
        like *_run0.pkl, *_run1.pkl, etc.
    suffix:
        string, the string to identify the pickled results file i.e. ecosystems_suffix_run0.pkl
    '''

    # initialise ecosystem

    if initial_ecosystem == None:

        # random start for individuals
        geneRandFnc = lambda geneType: random.randint(0, 2**parameters['genetics'][geneType]['noLoci']-1)
        #sexRandFnc = lambda: random.choice( parameters['sexes'] )
        natalHabRandFnc = lambda: random.choice( landscape )

        initial_ecosystem = [ {
                    'adults': [
                    {'sex': 'h', 'natalHabType': natalHabRandFnc(), 'genotype': {geneType: geneRandFnc(geneType) for geneType in parameters['genetics'].keys()}},
                    {'sex': 'h', 'natalHabType': natalHabRandFnc(), 'genotype': {geneType: geneRandFnc(geneType) for geneType in parameters['genetics'].keys()}} ],
                    'juveniles': list(),
                    }
            for habType in landscape]

    ecosystem = copy.deepcopy( initial_ecosystem ) # make a deep copy so we can store the initial conditions in pickle file

    # simulate ecosystem for burn-in timesteps, but don't record results

    t = 1
    while t <= burnInT and not ecosystemIsEmpty(ecosystem):

        ecosystem, _ = timestep(parameters, ecosystem, landscape)
        t += 1


    # simulate ecosystem for remaining timesteps and record results

    ecosystems = list() # a place to store the ecosystem at each timestep
    while t <= tf and not ecosystemIsEmpty(ecosystem):

        # one timestep of simulation
        ecosystem, _ = timestep(parameters, ecosystem, landscape)

        # store info
        ecosystems.append( copy.deepcopy(ecosystem) )

        t += 1


    # pickle the info

    # build the pickled results file's name
    fName = 'ecosystems'

    if suffix == None:
        suffix = parameters2filesuffix(tf, landscape, parameters)
    fName += suffix

    if idxRun != None:
        fName += '_run' + str(idxRun)

    fName += '.pkl'

    # open the file with the name fName
    f = open(fName, 'wb')

    # a string explaining the pickle file
    path = os.path.dirname(os.path.realpath('simulate.py'))
    ss  = 'Created by simulate.py in ' + path + '.\n'
    ss += 'Contains the following:\n'
    ss += '0. ss, string: this string you are reading now.\n'
    ss += '1. burnInT, integer: the number of unrecorded timesteps of burn-in.\n'
    ss += '2. t, integer: the total number up to which timesteps run (so range(burnInT+1,t)).\n'
    ss += '3. tf, integer: the total maximum number up to which timesteps were attempted (if pop did not go extinct, t = tf+1).\n'
    ss += '4. landscape, string: a string defining the landscape habitat type composition (e.g. LLLLHHHLLLL).\n'
    ss += '5. ecosystems, list of lists of dictionaries: the ecosystem at the end of each recorded timestep.\n'
    ss += '6. path, string: the path in which the run was performed.\n'
    ss += '7. parameters, dictionary: the parameter values with which the run was performed.\n'
    ss += '8. initial_ecosystem, list of dictionaries: the initial ecosystem.\n'

    pickle.dump( ss, f ) # 0.
    pickle.dump( burnInT, f )
    pickle.dump( t, f )
    pickle.dump( tf, f ) # 3.
    pickle.dump( landscape, f )
    pickle.dump( ecosystems, f )
    pickle.dump( path, f ) # 6.
    pickle.dump( parameters, f )
    pickle.dump( initial_ecosystem, f )

    f.close()

