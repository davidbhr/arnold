import os
import numpy as np
from subprocess import call
from arnold import SAENOPATH



# Cylindric Inclusion --------------------------------------------------------------------------------------------------

def cylindric_contraction(simulation_folder, mesh_file, d_cyl, l_cyl, r_outer, strain, material):
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
        saeno (str): File path to Saeno.exe
        material (dict): Material properties in the form {'K_0': X, 'D_0':X, 'L_S': X, 'D_S': X} (see materials)
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
    deformation = (strain*l_cyl)
    
    
    # create data folder if they do not exist ------------------------------------------------------------------------------------------------  
    if not os.path.exists(simulation_folder):
        os.makedirs(simulation_folder)
    print('+ Created output folder')
    
    
    # create coords.dat and tets.dat -----------------------------------------------------------------------------------
    with open(mesh_file, 'r') as f:
        lines = f.readlines()
    
    print (simulation_folder + '/coords.dat')
    
    # coords.dat 
    index_nodes = lines.index('$Nodes\n')
    number_of_nodes = int(lines[index_nodes+1])
    coords = np.zeros((number_of_nodes, 3))
    for i in range(number_of_nodes):
        coords[i] = np.array([np.float(x) for x in lines[i+index_nodes+2].split()[1:]]) 
    np.savetxt(simulation_folder + '/coords.dat', coords)
       
    
    # tets.dat
    index_elements = lines.index('$Elements\n')
    number_of_elements = int(lines[index_elements+1])
    tets = np.zeros((number_of_elements, 4))
    for i in range(number_of_elements):
        tets[i] = lines[i+index_elements+2].split()[-4:]
    np.savetxt(simulation_folder + '/tets.dat', tets, fmt='%i')
    print('+ Created coords.dat and tets.dat')
    
    
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
  
    #Plot selected nodes -----------------------------------------------------------------------------------------------
    # fig1 = plt.figure()
    # ax1 = Axes3D(fig1)
    # ax1.set_xlabel('X')
    # ax1.set_ylabel('Y')
    # ax1.set_zlabel('Z')
    # ax1.scatter(coords[mask_inner][:,0], coords[mask_inner][:,1], coords[mask_inner][:,2])
    # ax1.scatter(coords[mask_outer][:,0], coords[mask_outer][:,1], coords[mask_outer][:,2])
    # ax1.xaxis.set_major_formatter(FormatStrFormatter('%0.0e'))
    # ax1.yaxis.set_major_formatter(FormatStrFormatter('%0.0e'))
    # ax1.zaxis.set_major_formatter(FormatStrFormatter('%0.0e'))
    # plt.show()
    
    
    
    #Set Boundaries-----------------------------------------------------------------------------------------------------
    bcond = np.zeros((len(coords), 4))
    bcond[:, 3] = 1.
    
    
    #setdispl mode at outer boundary (fixed displacements)
    bcond[mask_outer, 3] = 0
         
    #setdispl mode at surface (fixed displacements)
    bcond[mask_inner, 3] = 0   
         
    #fixed displacements in x direction at surface (linear decreasing in size from both poles)
    bcond[mask_inner, 0] =   (-deformation/2) * (coords[mask_inner][:,0]/ (l_cyl/2))    
    
          
    #Set Displacements for volume conservation  (pi r0² h0 = pi r1² h1)
    r1 = (d_cyl/2) * (1/np.sqrt(((l_cyl/2) + (-deformation/2))/(l_cyl/2)))  
    dr = r1 - (d_cyl/2)
    bcond[mask_inner, 1] = ((coords[mask_inner][:,1])/(d_cyl/2))*dr
    bcond[mask_inner, 2] = ((coords[mask_inner][:,2])/(d_cyl/2))*dr
       
    #save bcond.dat
    np.savetxt(simulation_folder + '/bcond.dat', bcond)
    
    #create iconf
    iconf = np.zeros((len(coords), 3))
    iconf[mask_inner, 0] = (-deformation/2) * (coords[mask_inner][:,0]/(l_cyl/2))
    iconf[mask_inner, 1] = ((coords[mask_inner][:,1])/(d_cyl/2))*dr
    iconf[mask_inner, 2] = ((coords[mask_inner][:,2])/(d_cyl/2))*dr
    np.savetxt(simulation_folder + '/iconf.dat', iconf)
    print('+ Created bcond.dat and iconf.dat')
    
    
    #  Default values:
    #
    # REL_CONV_CRIT = 1e-11
    # REL_ITERATIONS = 300
    # REL_SOLVER_STEP = 0.066
    #
    # ToDo: Iterations and solver as function parameter
    
    
    
    
    # create config.txt -----------------------------------------------------------
    config = r"""
    MODE = relaxation
    BOXMESH = 0
    FIBERPATTERNMATCHING = 0
    REL_CONV_CRIT = 1e-11
    REL_ITERATIONS = 300
    REL_SOLVER_STEP = 0.066
    K_0 = {}
    D_0 = {}
    L_S = {}
    D_S = {}
    CONFIG = {}\config.txt
    DATAOUT = {}
    """.format(K_0, D_0, L_S, D_S, simulation_folder, simulation_folder)
    
    with open(simulation_folder + "/config.txt", "w") as f:
        f.write(config)
    print('+ Created config.txt')
    
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
       
    """.format(K_0, D_0, L_S, D_S, str(r_outer*1e6)+' µm', np.sum(mask_inner), len(coords), mesh_file, simulation_folder, str(d_cyl*1e6)+' µm', str(l_cyl*1e6)+' µm', str(deformation*1e6)+' µm',  str(strain), str(inner_spacing*1e6)+' µm',str(outer_spacing*1e6)+' µm')
    

    with open(simulation_folder + "/parameters.txt", "w") as f:
        f.write(parameters)
    print('+ Created parameters.txt')
            
    # start SAENO ----------------------------------------------------------------------------------------------------------
    print ('SAENO Path: '+str(SAENOPATH))
    print('+ Starting SAENO...')
    
    call(SAENOPATH+"/saeno CONFIG {}/config.txt".format(os.path.abspath(simulation_folder)))
    
    
    




