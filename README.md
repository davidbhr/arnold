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


Now we need to tell arnold  where `Gmsh` and `SAENO` are stored for one single time. Therefore we simply use

```python
import arnold as ar

ar.set_gmsh_path(r'C:\..\gmsh-4.3.0-Windows64-sdk')
ar.set_saeno_path(r'C:\..\SAENO')
```

Note: If the `Gmsh SDK` has been installed via `pip` the path to `Gmsh` does not have to be set at all.





## Introduction

A typical measurement can look like the following: Muscle fibers are embedded into a certain matrix environment (e.g. matrigel). Then a contraction is induced via electical stimulation. The deformation can be imaged e.g. via confocal, fluorescence or brightfield microscopy using a high framerate. Beads are embedded to check that the fibers are attached to the matrix environment and in order to compare the simulated and measured matrix deformations.

*In the Gif below a single flexor digitorum longus fiber is shown in the relaxed state and during different stimuli (single pulse stimulation and tetanic stimulaitons with frequencies ranging from 10 Hz to 100 Hz)*

![Loading GIF...](https://raw.githubusercontent.com/davidbhr/arnold/master/docs/GIFs/FDB_contraction(SP-10-25-50-75-100Hz).gif)


To measure the exerted muscle forces, we need the following information about the respective contraction: 

- The fiber length and diameter in relaxed state (easily measured from the raw images). 
- The Strain (derived as (Relaxed_Length- Contracted_Length)/Relaxed_Length) 
- The Youngs modulus of the material in Pa. In case of non linear materials more parameters needed (see Julian .. + parameterssee X.(derived from rheometer measurements or literature)



## Individual Simulation

A simulation can be splitted into three parts:

*__First__*, we need to build a mesh model of the geometry using GMSH. Here the fiber in the relaxed state is simulated as an inner cylindric inclusion with length *l_cyl* and diameter *d_cyl*, which is placed into a sphere with radius *r_outer* that simulates the matrix environment. Units are given in *µm* and the *length_factor* allows to tune the fineness of the mesh model (lower values correspond to a finer mesh. The element sizes are increased close to the cylindric inclusion by default)

```python

ar.mesh.cylindrical_inclusion(mesh_file=r'C/../Our_model.msh', d_cyl=30, l_cyl=300,
r_outer=2000, length_factor=0.2)
```

The resulting mesh can be displayed by using the following command:

```python

ar.mesh.show_mesh(r'C/../Our_model.msh')
```


<img src="https://raw.githubusercontent.com/davidbhr/arnold/master/docs/PNGs/GMSH_arnold.png" width="600" >




*__Second__*, we apply corresponding boundary condition to our mesh and simulate the contraction of the cylinder using the network optimizer SAENO. Here *x* is ... strain in .. +zero displacement at outer boundary of the bulk material.


```python
ar.simulation.cylindric_contraction(simulation_folder=r'C/../Simulation', mesh_file=r'C/../Our_model.msh', 
d_cyl=30, l_cyl=300, r_outer=2000, strain=0.1, ar.materials.matrigel10)
```

material  ..   [Steinwachs et al. (2016)](https://www.nature.com/articles/nmeth.3685)
 "pre-configured" material types for collagen gels of three different concentrations (0.6, 1.2, and 2.4mg/ml). Detailed protocols for reproducing these gels can be found in [Steinwachs et al. (2016)](https://www.nature.com/articles/nmeth.3685) and [Condor et al. (2017)](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/cpcb.24). Furthermore, one can define a linear elastic material with a specified `stiffness` (in Pa) with:

```
jf.materials.linear(stiffness)
```
To define a non-linear elastic material, use the `custom` material type:
```
jf.materials.custom(K_0, D_0, L_S, D_S)
```
Non-linear materials are characterized by four parameters:
- `K_0`: the linear stiffness (in Pa)
- `D_0`: the rate of stiffness variation during fiber buckling
- `L_S`: the onset strain for strain stiffening
- `D_S`: the rate of stiffness variation during strain stiffening
A full description of the non-linear material model and the parameters can be found in [Steinwachs et al. (2016)](https://www.nature.com/articles/nmeth.3685)





for displacements set on surface we obtain forces for material we obtain deformation field die eben diese kräfte kompensieren/erklären.
 output files of a simulation can be found the [Wiki of the SAENO project](https://github.com/Tschaul/SAENO/wiki). The file `parameters.txt` contains all parameters used int he simulation.



*__Third__*, we compute the overall contractility from the resultin simulation (by summing up). 

```python
ar.force.reconstruct_contractility(simulation_folder=r'C/../Simulation', d_cyl=30, l_cyl=300, r_outer=2000)
```

We receive a excel document contrctilities. contractility mean x components gves us of contractility of muscefiber only taking x components into account (possible noise errors..) and half of cell (due magdeburger halkugeln..). Also store contractility derived from absolute vectors instead of x components and over whole area instead of half. Additionally residuum forces are calculated to check validity of simulation. For visualization and possible error detection an image of the deformtion and force field is stored. 


<img src="https://raw.githubusercontent.com/davidbhr/arnold/master/docs/PNGs/Force_Displ_Field.png" width="400" >




*__All three abovementioned steps can be performed consecutively by using a single command. For an individual simulation just need to execute the combined function:__*


```python

ar.experiment.cylindrical_inclusion_mesh_simulation_and_contractility(mesh_file=r'C/../Our_model.msh', 
d_cyl=30,l_cyl=300, r_outer=2000, length_factor=0.2, simulation_folder=r'C/../Simulation', 
strain=0.1, ar.materials.matrigel10)
```


## Materials





## Scaling law


## Series

starting

series Evaluation




## Dependencies
*arnold* is tested on Python 3.6. It depends on .... All except ... are already included in the [Anaconda distribution](https://www.continuum.io/downloads) of Python. Windows users may also take advantage of pre-compiled binaries for all dependencies, which can be found at [Christoph Gohlke's page](http://www.lfd.uci.edu/~gohlke/pythonlibs/).

## License
[The MIT License (MIT)](https://github.com/davidbhr/arnold/blob/master/LICENSE)
