# arnold

Arnold is a python package to analyze muscle forces using 3D traction force microscopy. It provides interfaces to the open-source finite element mesh generator [`Gmsh`](http://gmsh.info/) and to the network optimizer [`SAENO`](https://github.com/Tschaul/SAENO). 

By assuming a simple cylindircal geometry, traction forces can be computed for single fibers or a series of fibers by individual simulations. Fibers with different geometry can be evaluated for different contraction strains and different tissue environments (*linear* and *non-linear* materials).


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


GMSH offers a python interface within the `Gmsh SDK` , that can be downloaded [here](http://gmsh.info/#Download)
or by simply running the following command: `pip install --upgrade gmsh-sdk`.


To find the equilibrium configuration to the applied boundary problem *arnold* uses the network optimizer [`SAENO`](https://github.com/Tschaul/SAENO). A precompiled version of `SAENO` for 64bit Windows systems can be downloaded [here](https://github.com/davidbhr/arnold/tree/master/docs/SAENO). To build SAENO on other platforms the project can be found [here](https://github.com/Tschaul/SAENO).


Now we need to tell arnold  where `Gmsh` and `SAENO` are stored for one single time. Therefore we simply use:

```python
import arnold as ar

ar.set_gmsh_path(r'C:\..\gmsh-4.3.0-Windows64-sdk')
ar.set_saeno_path(r'C:\..\SAENO')
```

Note: If the `Gmsh SDK` has been installed via `pip` the path to `Gmsh` does not have to be set at all.





## Introduction

A typical measurement can look like the following: Muscle fibers are embedded into a certain matrix environment (e.g. matrigel). Then a contraction is induced via electical stimulation. The deformation can be imaged e.g. via confocal, fluorescence or brightfield microscopy using a high framerate. Beads are embedded to validate that the fibers are attached to the matrix environment and in order to compare the simulated and measured matrix deformations.

*In the Gif below a single flexor digitorum longus fiber embedded in matrigel is shown in the relaxed state and during different stimuli (single pulse stimulation and tetanic stimulaitons with frequencies ranging from 10 Hz to 100 Hz)*

![Loading GIF...]()


To measure the exerted muscle forces, we need the following information about the respective contraction: 

- The **fiber length** and diameter in relaxed state (can easily be measured from the raw images). 
- The **strain** during the contraction (derived as (Relaxed_Length- Contracted_Length)/Relaxed_Length) 
- The **material properties** of the surrounding matrix (see below for further details)



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





*__Second__*, we apply corresponding boundary condition to our mesh and simulate the contraction of the cylinder using the network optimizer SAENO. We read in the *mesh_file* and apply a symetrical contraction with the given *strain* on the inner cylindrical inclusion with specific length *l_cyl* and diameter *d_cyl*. At the outer boundary of the bulk material we constrain the deformations to zero.


```python
ar.simulation.cylindric_contraction(simulation_folder=r'C/../Simulation', mesh_file=r'C/../Our_model.msh', 
d_cyl=30, l_cyl=300, r_outer=2000, strain=0.1, ar.materials.matrigel_10mg_ml)
```


The material properties can be defined as following: 

 
- To define a non-linear elastic material by it's linear stiffness (for instance 1000 Pa) use:
 ```python
ar.materials.linear_stiffness(1000)
```

- To define a specific Young's modulus (for instance 500 Pa) for a linear-elastic hydrogel (such as Matrigel) with a poission ratio of 0.25 (see  [Steinwachs                 et al. (2016)](https://www.nature.com/articles/nmeth.3685)) use: 
```python
ar.materials.youngs_modulus(500)
````

- To define non-linear materials use the `custom` material type:
```python
ar.materials.custom(K_0, D_0, L_S, D_S)
```

Non-linear materials are characterized by four parameters:
- `K_0`: the linear stiffness (in Pa)
- `D_0`: the rate of stiffness variation during fiber buckling
- `L_S`: the onset strain for strain stiffening
- `D_S`: the rate of stiffness variation during strain stiffening

A full description of the non-linear material model and the parameters can be found in [Steinwachs et al. (2016)](https://www.nature.com/articles/nmeth.3685). Also "pre-configured" material types for Matrigel (10mg/ml) and collagen gels of three different concentrations (0.6, 1.2, and 2.4mg/ml) are available. Detailed protocols for reproducing these gels can be found in [Steinwachs et al. (2016)](https://www.nature.com/articles/nmeth.3685) and [Condor et al. (2017)](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/cpcb.24). 


```python
ar.materials.collagen_06
ar.materials.collagen_12
ar.materials.collagen_24
ar.materials.matrigel_10mg_ml
```







The results are stored in the set *simulation_folder*. For the nodes at the fiber surface we obtain the respective forces, which compensate for the fiber deformation, for the nodes in the bulk material we obtain the deformations, which are thereby generated in the surrounding matrix. More information about the output files of a simulation can be found in the [Wiki of the SAENO project](https://github.com/Tschaul/SAENO/wiki). The file `parameters.txt` contains all set parameters, that are used in the simulation.








*__Third__*, we compute the overall contractility from the resulting simulation by summing up the x-components of all forces on the cell surface from one end of the fiber to its center.

```python
ar.force.reconstruct_contractility(simulation_folder=r'C/../Simulation', d_cyl=30, l_cyl=300, r_outer=2000)
```

In an excel document the resulting contractility is saved under *'Contractility mean x-components'*. Further, the contractility is calculated considering the absolute forces instead of the x-components only and individually over the left half, the right half and total fiber volume. Additionally the residuum forces are calculated, which can be used to verify the simulation (forces at the outer boundary must compensate for the forces within the mesh moddel). For visualization and error detection an image of the deformtion- and force-field is stored. 


<img src="https://raw.githubusercontent.com/davidbhr/arnold/master/docs/PNGs/Force_Displ_Field.png" width="400" >




*__All three abovementioned steps can be performed consecutively by using a single command. For an individual simulation you just need to execute the combined function:__*

```python

ar.experiment.cylindrical_inclusion_mesh_simulation_and_contractility(mesh_file=r'C/../Our_model.msh', 
d_cyl=30,l_cyl=300, r_outer=2000, length_factor=0.2, simulation_folder=r'C/../Simulation', 
strain=0.1, ar.materials.matrigel10)
```
Note: *mesh_file* and *simulation_folder* now only describe the desired directory and are created and used automatically. 




## Scaling law

For linear-elastic materials, the total contractility of muscle fibers scales linearly with matrix elasticity, linearly with fiber strain,  quadratically with the fiber length and linearly with fiber diameter (see REF). Therefore we can estimate the contractility of a muscle fiber from the fiber diameter `d`, fiber length `l`, the fiber strain `s` during contraction, and the Young’s modulus `E` of the surrounding Matrix by comparison with a reference simulation by uing the following simple scaling law (returns the contractility in µN): 

```python

ar.force.scaling_law(d,l,s,E)
```



## Series
by set of simulations..  for nonllineear materials different scaling. (ToDo).

To start a series of simulation preset functions. to simulate fibers with differend diameter ´s for fixed strain etc use

```python

ar.experiment.series (..)
```

Here *n_cores* number of cores to start simulations in paralell (for no entry by default detect cores). In a similar manner series for different lengths `ar.experiment.series (..)`, different strains  `ar.experiment.series (..)` and different stiffness  `ar.experiment.series (..)` can be used.


For all these again excel sheet with evaluated contractility within simulation folder. To compare several simulations in a more convenient way, the following funcion scans a given directiory (and all subdirecorys) for all simulations and creates an overview of the important values (cell geometrie strain etc..) from these.

```python

ar.experiment.evaluate (..)
```


BILD




## Dependencies
*arnold* is tested on Python 3.6. It depends on `numpy`, `pandas`, `matplotlib` and XY. It is recommended to use the Anaconda distribution of Python. Windows users may also take advantage of pre-compiled binaries for all dependencies, which can be found at [Christoph Gohlke's page](http://www.lfd.uci.edu/~gohlke/pythonlibs/).

## License
The `arnold` package itself is licensed under the [MIT License](https://github.com/davidbhr/arnold/blob/master/LICENSE). Everything is provided "as is", without any warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement.
