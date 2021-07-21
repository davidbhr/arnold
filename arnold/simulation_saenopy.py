import os
import numpy as np
from arnold import SAENOPATH
import subprocess
from saenopy import Solver
from saenopy.materials import SemiAffineFiberMaterial


def read_meshfile(meshfile):
    # open mesh file
    with open(meshfile, 'r') as f:
        lines = f.readlines()
 
    # transform nodes and connection in SAENO format
    # nodes
    index_nodes = lines.index('$Nodes\n')
    n_nodes = int(lines[index_nodes + 1])

    coords = np.zeros((n_nodes, 3))
    for i in range(n_nodes):
        coords[i] = np.array([np.float(x) for x in lines[i + index_nodes + 2].split()[1:]])
            
    # connections
    index_elements = lines.index('$Elements\n')
    n_elements = int(lines[index_elements + 1])

    tets = np.zeros((n_elements, 4))
    for i in range(n_elements):
        tets[i] = lines[i + index_elements + 2].split()[-4:]

    # to start with 0 and not 1
    tets -= 1
    return coords, tets




# Cylindric Inclusion --------------------------------------------------------------------------------------------------

def cylindric_contraction(simulation_folder, mesh_file, d_cyl, l_cyl, r_outer, strain, material, logfile = False,  iterations= 300 , step=0.3, conv_crit = 1e-11):
    """
    Simulates an symetric contraction (with constant strain) of the cylindric inclusion inside the spherical bulk.
    
    Args:
        simulation_folder(str): File path to save simulation results
        mesh_file(str): File path to read in mesh file
        d_cyl(float): Diameter of the spherical inclusion in the mesh model (in µm) 
        l_cyl(float): Length of the spherical inclusion in the mesh model (in µm)
        r_outer(float): Outer radius of the bulk mesh in the mesh model (in µm)
        strain(float): Strain as (Length_base - Length_contracted)/Length_base. 
        Deformation is applied in x-direction and split equally to both poles,
        from which defomation-size decreases linearly to the center (symetric contraction with constant strain). 
        material (dict): Material properties in the form {'K_0': X, 'D_0':X, 'L_S': X, 'D_S': X} (see materials)
        logfile(boolean): If True a reduced logfile of the saeno system output is stored. Default: False.
        iterations(float): The maximal number of iterations for the saeno simulation. Default: 300.
        step(float): Step width parameter for saeno regularization. Higher values lead to a faster but less robust convergence. Default: 0.3.
        conv_crit(float): Saeno stops if the relative standard deviation of the residuum is below given threshold. Default: 1e-11.         
    """
    

    
    # read in material parameters
    K_0 = material['K_0']
    D_0 = material['D_0']
    L_S = material['L_S']
    D_S = material['D_S']

    
    """
    Apply Boundary Conditions and start Saeno Simulation
    """
    
    #convert input in um------------------------------------------------------------------------------------------------
    d_cyl *= 1e-6
    l_cyl *= 1e-6
    r_outer *= 1e-6
    deformation = strain*l_cyl
    
    
    # create data folder if they do not exist ------------------------------------------------------------------------------------------------  
    if not os.path.exists(simulation_folder):
        os.makedirs(simulation_folder)
    print('+ Created output folder')
      
    # create coords.dat and tets.dat -----------------------------------------------------------------------------------
    coords, tets = read_meshfile(mesh_file)
    
    # create solver object in saenopy
    M = Solver()
    material_saenopy = SemiAffineFiberMaterial( k=K_0, d0=D_0, lambda_s=L_S, ds=D_S)  
    M.setMaterialModel(material_saenopy)
    M.setNodes(coords)
    M.setTetrahedra(tets)
    

    # create bcond.dat and iconf.dat
    distance = np.sqrt(np.sum(coords**2., axis=1))
    x, y, z = coords.T
    mask_inner = ((y**2. + z**2.) <= (d_cyl/2)**2.) & (x**2. <= (l_cyl/2)**2.)
    mask_outer = distance > r_outer * 0.999  # include close elements to the boundary 
      
    
    # Save Node Density at Surface/outer   
    # Area per outer node
    A_node_outer = (np.pi*4*(r_outer)**2)/np.sum(mask_outer)   
    # simply sqrt as spacing
    outer_spacing = np.sqrt(A_node_outer)     
        
    # Area per inner node
    A_node_inner = (2*np.pi*(d_cyl/2)*((d_cyl/2)+(l_cyl)))/np.sum(mask_inner)
    # simply sqrt as spacing
    inner_spacing = np.sqrt(A_node_inner)
    
    print ('Outer node spacing: '+str(outer_spacing*1e6 )+'µm')
    print ('Inner node spacing: '+str(inner_spacing*1e6 )+'µm')
  
    
    #Set Boundaries-----------------------------------------------------------------------------------------------------
   
    
    # set displacements  
    bcond_displacement = np.zeros((len(coords), 3))*np.nan  # nan in the bulk, there they will be calculated
    bcond_displacement[mask_outer] = 0   # constraint to zero at the outer border
    #fixed displacements in x direction at surface (linear decreasing in size from both poles)
    bcond_displacement[mask_inner, 0] =   (-deformation/2) * (coords[mask_inner][:,0]/ (l_cyl/2))   # fi
    #Set Displacements for volume conservation  (pi r0² h0 = pi r1² h1)
    r1 = (d_cyl/2) * (1/np.sqrt(((l_cyl/2) + (-deformation/2))/(l_cyl/2)))  
    dr = r1 - (d_cyl/2)
    bcond_displacement[mask_inner, 1] = ((coords[mask_inner][:,1])/(d_cyl/2))*dr
    bcond_displacement[mask_inner, 2] = ((coords[mask_inner][:,2])/(d_cyl/2))*dr
   
    # set forces 
    # forces are everywhere 0 in the bulk
    bcond_forces = np.zeros((len(coords), 3))
    # only at outer space and surface calulate fores
    bcond_forces[mask_outer] = np.nan
    # and at the inner border 
    bcond_forces[mask_inner] = np.nan
   
    # Apply zero deformations as intial state 
    # for all fixed nodes (here cell/outer) initial displacements areignored
    
    
    # set initial displacements  
    iconf = np.zeros((len(coords), 3))  # nan in the bulk, there they will be calculated
    iconf[mask_outer] = 0   # constraint to zero at the outer border
    #fixed displacements in x direction at surface (linear decreasing in size from both poles)
    iconf[mask_inner, 0] =   (-deformation/2) * (coords[mask_inner][:,0]/ (l_cyl/2))   # fi
    #Set Displacements for volume conservation  (pi r0² h0 = pi r1² h1)
    r1 = (d_cyl/2) * (1/np.sqrt(((l_cyl/2) + (-deformation/2))/(l_cyl/2)))  
    dr = r1 - (d_cyl/2)
    iconf[mask_inner, 1] = ((coords[mask_inner][:,1])/(d_cyl/2))*dr
    iconf[mask_inner, 2] = ((coords[mask_inner][:,2])/(d_cyl/2))*dr

    
    #M.setInitialDisplacements(np.zeros((len(coords), 3)))
    
    
    # give the boundary conditions to the solver
    M.setBoundaryCondition(bcond_displacement, bcond_forces)
    
    M.setInitialDisplacements(iconf)
    
    
    # create parameters.txt--------------------------------------------------------
    parameters = r"""
    K_0 = {}
    D_0 = {}
    L_S = {}
    D_S = {}
    OUTER_RADIUS = {}
    SURFACE_NODES = {}
    TOTAL_NODES = {}
    Mesh_file = {}
    Output_folder = {}
    d_cyl = {}
    l_cyl = {}
    deformation = {}
    strain = {}
    Inner node spacing = {}
    Outer node spacing = {}
    iterations = {}
    step = {}
    conv_crit = {}
       
    """.format(K_0, D_0, L_S, D_S, str(r_outer*1e6)+' µm', np.sum(mask_inner), len(coords), mesh_file, simulation_folder, 
    str(d_cyl*1e6)+' µm', str(l_cyl*1e6)+' µm', str(deformation*1e6)+' µm',  str(strain), str(inner_spacing*1e6)+' µm',
    str(outer_spacing*1e6)+' µm', iterations, step, conv_crit)
    

    with open(simulation_folder + "/parameters.txt", "w") as f:
        f.write(parameters)
    print('+ Created parameters.txt')
            
    # start SAENOPY ----------------------------------------------------------------------------------------------------------
    # solve the boundary problem
    M.solve_boundarycondition(stepper=step, i_max=iterations, rel_conv_crit=conv_crit, relrecname=simulation_folder + "/relrec.txt", verbose=True)
    M.save(simulation_folder + "/solver.npz")
    

