# arnold

Arnold is a python package to analyze muscle forces 3D traction force microscopy. It provides interfaces to the open-source finite element mesh generator [`Gmsh`](http://gmsh.info/) and to the network optimizer [`SAENO`](https://github.com/Tschaul/SAENO). 

By assuming a simple cylindircal geometry,  traction forces can be computed for single fibers or a series of fibers by individual simulations. Fibers with different geometry can be  for different strains and materials (linear and non-linear).

A scaling law, derived from such simulation in linear materials, enables computing forces just by cell geometry an strain without the need of individual simulations.

## Getting started


Code example
```python
import arnold as ar
import matplotlib.pyplot as plt

...
```



## Installation
The easiest way to install the latest release version of *arnold* is via `pip`:
```
pip install -e . 
```
Alternatively, a zipped version can be downloaded [here](https://github.com/...). The module is installed by calling `python setup.py install`.

### Development version
The latest development version of *arnold* can be installed from the master branch using pip (requires git):
```
pip install git+https://github.com/...
```
Alternatively, use this [zipped version](https://github.com/.../zipball/master) or clone the repository.

## Intro 


![Loading GIF...](https://raw.githubusercontent.com/davidbhr/arnold/master/docs/GIFs/FDB_contraction(SP-10-25-50-75-100Hz).gif)



## Individual Simulation


<img src="https://raw.githubusercontent.com/davidbhr/arnold/master/docs/PNGs/GMSH_arnold.png" width="600" >



## Scaling law



## Series

starting

series Evaluation


## Dependencies
*arnold* is tested on Python 3.6. It depends on .... All except ... are already included in the [Anaconda distribution](https://www.continuum.io/downloads) of Python. Windows users may also take advantage of pre-compiled binaries for all dependencies, which can be found at [Christoph Gohlke's page](http://www.lfd.uci.edu/~gohlke/pythonlibs/).

## License
[The MIT License (MIT)](https://github.com/.../blob/master/LICENSE)
