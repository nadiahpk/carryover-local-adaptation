# Carryover effects and local adaptation

Code to accompany a working manuscript: Kristensen, N.P., Johansson, J., Chisholm, R.A., Smith, H.G., Kokko, H. Carryover effects from natal habitat type upon competitive ability lead to trait divergence or source-sink dynamics

Local adaptation to rare habitats is difficult due to gene flow, but can occur if the habitat has higher productivity. Differences in offspring phenotypes have attracted little attention in this context. We model a scenario where the rarer habitat improves offsprings' later competitive ability -- a carryover effect that operates on top of local adaptation to one or the other habitat type. Assuming localised dispersal, so offspring tend to settle in similar habitat to the natal type, the superior competitive ability of offspring remaining in the rarer habitat hampers immigration from the majority habitat. This initiates a positive feedback between local adaptation and trait divergence, which can thereafter be reinforced by coevolution with dispersal traits that match ecotype to habitat type. Rareness strengthens selection on dispersal traits and promotes linkage disequilibrium between locally-adapted traits and ecotype-habitat matching dispersal. We propose that carryover effects may initiate isolation by ecology.

The code runs an individual-based, genetically- and spatially-explicit model of a single population on a landscape of two habitat types. Options include adding sexual reproduction, additional habitat types, and genetically-determined dispersal characteristics (see `script.py` and the dictionary `parameters` for options).

## Quickstart

### Installing

If you have git:

```
$ git clone https://github.com/nadiahpk/carryover-local-adaptation
$ cd carryover-local-adaptation
```
If you don't have git, just download and unpack the latest zip file:
```
[https://github.com/nadiahpk/carryover-local-adaptation/archive/master.zip](https://github.com/nadiahpk/carryover-local-adaptation/archive/master.zip)
```

### Running

In Python3:
```
In [1]: %run script.py # runs a script that simulates ecoevolutionary dynamics for 600 generations and stores results in ecosystems_1.pkl
In [2]: %run plotFigure1s.py -f ecosystems_1 # uses ecosystems_1.pkl to create a figure ecosystems_1_Fig1s.png
```

## License

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or distribute this software, either in source code form or as a compiled binary, for any purpose, commercial or non-commercial, and by any means.

In jurisdictions that recognize copyright laws, the author or authors of this software dedicate any and all copyright interest in the software to the public domain. We make this dedication for the benefit of the public at large and to the detriment of our heirs and successors. We intend this dedication to be an overt act of relinquishment in perpetuity of all present and future rights to this software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>

