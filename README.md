# arnold

Arnold is a python package to analyze muscle forces using 3D traction force microscopy. It provides interfaces to the open-source finite element mesh generator [`Gmsh`](http://gmsh.info/) and to the network optimizer [`SAENO`](https://github.com/Tschaul/SAENO). 

By assuming a simple cylindircal geometry, traction forces can be computed for single fibers or a series of fibers by individual simulations. Fibers with different geometry can be evaluated for different strains and different tissue environments (linear and non-linear materials).


A scaling law, derived from such simulation in linear materials, enables computing forces just by the cell geometry, the strain and the material stiffness without the need of further individual simulations.



## Installation

The current version of *arnold* can be installed by cloning this repository or by downloading the zipped version [here](https://github.com/davidbhr/arnold/zipball/master). By running the following command within the unzipped folder, all other required packages are automatically downloaded and installed.

```
pip install -e . 
```




Then the module can be simply imported to python as following

```python
import arnold as ar
```




We still need to set up the interfaces to the open-source finite element mesh generator [`Gmsh`](http://gmsh.info/) and to the network optimizer [`SAENO`](https://github.com/Tschaul/SAENO) .


GMSH offers a python interface, which are available in the `Gmsh SDK` and can be downloaded [here](http://gmsh.info/#Download),
or by simply running the following command: `pip install --upgrade gmsh-sdk`.


To find the equilibrium configuration to the applied boundary problem *arnold* uses the network optimizer [`SAENO`](https://github.com/Tschaul/SAENO). A precompiled version of `SAENO` for 64bit Windows systems can be downloaded [here](https://github.com/davidbhr/arnold/tree/master/docs/SAENO). To build SAENO on other platforms the project can be found [here](https://github.com/Tschaul/SAENO).


Now we need to tell arnold where both are for one single time. Therefore we simply use

```python
import arnold as ar

ar.set ...

ar. ...
```













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
