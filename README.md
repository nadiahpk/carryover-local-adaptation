# Carryover effects and local adaptation

Code to accompany a working manuscript: Kristensen, N.P., Johansson, J., Chisholm, R.A., Smith, H.G., Kokko, H. Carryover effects from natal habitat type upon competitive ability lead to trait divergence or source-sink dynamics

Local adaptation to rare habitat may occur if the habitat has sufficiently higher quality, i.e. per-capita growth rate, to reverse gene-flow asymmetry. We focus on another aspect of quality: improvement of offsprings' phenotype. We model a two-habitat-type scenario where the rarer habitat improves offsprings' later competitive ability; this carryover effect operates on top of local adaptation to the one or the other habitat type. We find that when localised dispersal is modelled so that offspring tend to settle in similar habitat to the natal type, the superior competitive ability of offspring remaining in the rarer habitat hampers immigration from the majority type. This initiates a positive feedback between local adaptation and trait divergence, which can thereafter be reinforced by coevolution with ecotype-habitat matching dispersal traits.

The code runs an individual-based, genetically- and spatially-explicit model of a single population on a landscape of two habitat types. Options include adding sexual reproduction, additional habitat types, and genetically-determined dispersal characteristics (see `script.py` and the dictionary `parameters` for options).

## Quickstart

In Python3:
```
In [1]: %run script.py # runs a script that simulates ecoevolutionary dynamics for 600 generations and stores results in ecosystems_1.pkl
In [2]: %run plotFigure1s.py -f ecosystems_1 # uses ecosystems_1.pkl to create a figure ecosystems_1_Fig1s.png
```

## License

Creative Commons.

