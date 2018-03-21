# run with e.g.: python3 plotFigure1s.py -f ecosystems_1

import pickle
import numpy as np
import matplotlib.pyplot as plt
#import sys
#sys.path.append("../../current_code/")
from carryover import phenInSpace
from carryover import noOffspringFnc
from carryover import pcolormeshCorrectionXY
import sys, getopt

try:

    opts, args = getopt.getopt(sys.argv[1:],'hf:')

except getopt.GetoptError:

    print('plotFigure1s.py -f <pkl file name>')
    sys.exit(2)

for opt, arg in opts:

    if opt == '-h':

        print('plotFigure1s.py -f <pkl file name>')
        sys.exit()

    elif opt == '-f':

        fName = arg


# some parameters for plotting

# order in which to plot them
geneTypeOrder = ['repn', 'pref', 'phil', 'dist','neut']
phenCmaps = { # which colourmap I want used for each phenotype
        'repn': 'coolwarm',
        'pref': 'coolwarm',
        'phil': 'PuOr',
        'dist': 'bone_r',
        'neut': 'bone_r',}

diffCmaps = { # which colourmap I want used for each genetic difference
        'repn': 'bone_r',
        'pref': 'bone_r',
        'phil': 'bone_r',
        'dist': 'bone_r',
        'neut': 'bone_r'}

titles = { # titles to give each column
        'repn': 'reproductive trait',
        'pref': 'habitat-type preference trait',
        'phil': 'NHPI trait',
        'dist': 'dispersal radius',
        'neut': 'neutral trait'}

cbar_labels = { # colourbar labels
        'repn': r'$\longleftarrow$ LC adapted $\: \vert \:$ HC adapted $\longrightarrow$',
        'pref': r'$\longleftarrow$ LC preferred $\vert$ HC preferred $\longrightarrow$',
        'phil': r'$\longleftarrow$ other preferred $\vert$ natal preferred $\longrightarrow$',
        'dist': 'no. territories',
        'neut': 'trait value'}

if 'nodiff' in fName:
    cbar_labels['repn'] = r'$\longleftarrow$ majority adapted $\: \vert \:$ minority adapted $\longrightarrow$'

# get data of run

f = open(fName + '.pkl', 'rb')
ss = pickle.load( f )
burnInT = pickle.load( f )
t = pickle.load( f )
tf = pickle.load( f )
landscape = pickle.load( f )
ecosystems = pickle.load( f )
path = pickle.load( f )
parameters = pickle.load( f )
f.close()


# calculate and store info needed for plotting

# list of our genetypes in order
geneTypes = [ geneType for geneType in geneTypeOrder if geneType in parameters['genetics'] ]

# prepare storage
phenDictTs = { geneType: list() for geneType in geneTypes } # phenotype values
noOffspringTs = list() # no offspring

for ecosystem in ecosystems: # NOTE may want to modify for burn in

    noOffspring = [ np.nan if len(flock['adults']) != 2 else noOffspringFnc( parameters, flock['adults'], habType) for flock, habType in zip(ecosystem,landscape) ] # number of offspring in space

    phenDict = phenInSpace(parameters, ecosystem) # phenotypes in space

    # store info about timestep
    noOffspringTs.append( noOffspring )

    for geneType in geneTypes:
        phenDictTs[geneType].append( phenDict[geneType] )

# ---

# need a correction of the axes for pcolormesh
lV, tV = pcolormeshCorrectionXY( list(range(len(landscape))), np.arange(burnInT+1,t) )

# plot phenotype values and proportion difference between parents' genes in space and time

nrows = 1
ncols = len(geneTypes) + 2 # plus 2 is for occupancy and number of offspring
f, ax = plt.subplots(nrows, ncols, sharex=True, sharey=True, figsize=(4*ncols,4*nrows))


# first column, occupancy

col = 0; aax = ax[col]

# sort out colourmap
base = plt.cm.get_cmap( diffCmaps['neut'] )
color_list = [ base( i ) for i in np.linspace(0, 1, 2+1) ]
cmap_name = base.name + 'occ'
newCmap = base.from_list(cmap_name, color_list, 2)

# plot occupancy
m = np.zeros( np.shape(noOffspringTs) )
m[ np.isnan(noOffspringTs) ] = 1
pp0 = aax.pcolormesh(lV, tV, m, cmap=newCmap, vmin=-1/2, vmax=1.5)
aax.set_xlim( (lV[0],lV[-1]) )
aax.set_ylim( (tV[0], tf+0.5) )
aax.set_ylabel('generation')
aax.set_title( 'territory occupancy' )
cbar = plt.colorbar(pp0, ax=aax, ticks=[0,1])
cbar.ax.set_yticklabels(['occupied','unoccupied'], rotation=90)
aax.set_xlabel('location')


# second column, number of offspring

col = 1; aax = ax[col]
maxVal = max( (parameters['habitats']['L']['rMax'],parameters['habitats']['H']['rMax']) )

# sort out colourmap
base = plt.cm.get_cmap( diffCmaps['neut'] )
color_list = [ base( i ) for i in np.linspace(0, 1, maxVal+1) ]
cmap_name = base.name + str(maxVal+1)
newCmap = base.from_list(cmap_name, color_list, maxVal+1)

m = np.array( noOffspringTs )
m = np.ma.masked_invalid(m)
pp0 = aax.pcolormesh(lV, tV, m, cmap=newCmap, vmin=-1/2, vmax=maxVal+1/2)
aax.set_xlim( (lV[0],lV[-1]) )
aax.set_ylim( (tV[0], tf+0.5) )
aax.set_title( 'reproduction' )
cbar = plt.colorbar(pp0, ax=aax, ticks=range(maxVal+1))
cbar.ax.set_ylabel('no. offspring')
aax.set_xlabel('location')


for col, geneType in enumerate(geneTypes):

    aax = ax[col+2]


    # info about this gene type

    # the phenotype values I want to plot on the z range
    if geneType == 'repn': # range on reproduction to -1 to 1 so clearer

        phens = [-1, -0.8, -0.6, -0.4, -0.2, 0, .2, .4, .6, .8, 1]

    else: # otherwise, full range

        noLoci = parameters['genetics'][geneType]['noLoci']
        minPhen = parameters['genetics'][geneType]['minPhen']
        maxPhen = parameters['genetics'][geneType]['maxPhen']
        phens = np.linspace(minPhen, maxPhen, noLoci+1) # possible phenotype values


    # use info about this gene type to sort out colour map and decorations

    # colourmap
    base = plt.cm.get_cmap( phenCmaps[geneType] )
    color_list = [ base( i ) for i in np.linspace(0, 1, len(phens)) ]
    cmap_name = base.name + 'new'
    newCmap = base.from_list(cmap_name, color_list, len(phens))

    # decorations
    tickPhens = [ phen for i,phen in enumerate(phens) if i%2 == 0 ] # tick every second
    delPhens = phens[1] - phens[0]
    aax.set_xlim( (lV[0],lV[-1]) )
    aax.set_ylim( (tV[0], tf+0.5) )
    aax.set_title( titles[geneType] )
    aax.set_xlabel('location')

    # sort out the data and draw the map

    m = np.array( phenDictTs[geneType] )
    m = np.ma.masked_invalid(m)
    pp0 = aax.pcolormesh(lV, tV, m, cmap=newCmap, vmin=phens[0]-delPhens/2, vmax=phens[-1]+delPhens/2)
    # colourbar
    cbar = plt.colorbar(pp0, ax=aax, ticks=tickPhens)
    cbar.ax.set_ylabel(cbar_labels[geneType])


plt.tight_layout()
plt.savefig(fName + '_Fig1s.png')
plt.close()
#plt.show()

