import os
import numpy as np
from arnold import SAENOPATH
import subprocess


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
       
    
    # create config.txt -----------------------------------------------------------
    config = r"""
    MODE = relaxation
    BOXMESH = 0
    FIBERPATTERNMATCHING = 0
    REL_CONV_CRIT = {}
    REL_ITERATIONS = {}
    REL_SOLVER_STEP = {}
    K_0 = {}
    D_0 = {}
    L_S = {}
    D_S = {}
    CONFIG = {}\config.txt
    DATAOUT = {}
    """.format(conv_crit, iterations, step, K_0, D_0, L_S, D_S, simulation_folder, simulation_folder)
    
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
    iterations = {}
    step = {}
    conv_crit = {}
       
    """.format(K_0, D_0, L_S, D_S, str(r_outer*1e6)+' µm', np.sum(mask_inner), len(coords), mesh_file, simulation_folder, 
    str(d_cyl*1e6)+' µm', str(l_cyl*1e6)+' µm', str(deformation*1e6)+' µm',  str(strain), str(inner_spacing*1e6)+' µm',
    str(outer_spacing*1e6)+' µm', iterations, step, conv_crit)
    

    with open(simulation_folder + "/parameters.txt", "w") as f:
        f.write(parameters)
    print('+ Created parameters.txt')
            
    # start SAENO ----------------------------------------------------------------------------------------------------------
    print ('SAENO Path: '+str(SAENOPATH))
    print('+ Starting SAENO...')
 
    
    # Create log file if activated
    if logfile == True:
        
        # create log file with system output
        logfile = open(simulation_folder + "/saeno_log.txt", 'w')
        cmd = subprocess.Popen(SAENOPATH+" CONFIG {}/config.txt".format(os.path.abspath(simulation_folder)), stdout=subprocess.PIPE , universal_newlines=True, shell=False)
        # print and save a reduced version of saeno log
        for line in cmd.stdout:
            if not '%' in line:
                print (line, end='')
                logfile.write(str(line))
        # close again to avoid loops            
        cmd.stdout.close()            
        
    # if false just show the non reduced system output    
    else:
        #cmd = subprocess.call(SAENOPATH+"\saeno CONFIG {}/config.txt".format(os.path.abspath(simulation_folder))) 
        cmd = subprocess.call(SAENOPATH+" CONFIG {}/config.txt".format(os.path.abspath(simulation_folder))) 
       

  
# two point forces in certain distances --------------------------------------------------------------------------    
                
def point_forces(simulation_folder, mesh_file, r_outer, distance, strain, material, logfile = True,  iterations= 300 , step=0.3, conv_crit = 1e-11):
    """
    Simulates an symetric contraction (with constant strain) of two nodes inside the spherical bulk. 
    Finds nodes matching best to given distance and updates strain accordingly.
    
    Args:
        simulation_folder(str): File path to save simulation results
        mesh_file(str): File path to read in mesh file
        distance(float): Desired distance beetween the two point forces (in µm). Closest possible points are searched and distance is updated
        r_outer(float): Outer radius of the bulk mesh in the mesh model (in µm)
        strain(float):  Strain as (distance_base - dist_contracted)/dist_base. 
        Deformation is applied on the axis between the 2 nodes split equally to both poles.
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
    r_outer *= 1e-6
    distance *= 1e-6
    #deformation = displacement*1e-6
    
    
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
    # now only save tetrahedron elements from mesh file 
    # not necessary for cylinder and spherical inclusion but for pointforces
    tets = [] #np.zeros((number_of_elements, 4))
    for i in range(number_of_elements):
        if lines[i+index_elements+2].split()[1]=='4':  
             tets.append(lines[i+index_elements+2].split()[-4:])
    tets = np.array(tets, dtype='int')
    np.savetxt(simulation_folder + '/tets.dat', tets, fmt='%i')
    print('+ Created coords.dat and tets.dat')
    
    
    # create bcond.dat and iconf.dat
    distance_c = np.sqrt(np.sum(coords**2., axis=1))
    x, y, z = coords.T
    
    
    
    # search for closest grid point on the x-axis for the given distances
    # point on positive x-side
    dist_p = np.sqrt((x-(distance/2))**2+(y**2)+(z**2))   
    min_dist_p = np.min(np.sqrt((x-(distance/2))**2+(y)**2+(z**2))) 
    mask_p1 = (dist_p==min_dist_p)
    
    
    # point on negative x-side
    dist_n = np.sqrt((x+(distance/2))**2+(y**2)+(z**2))
    min_dist_n = np.min(np.sqrt((x+(distance/2))**2+(y)**2+(z**2))) 
    mask_p2 = (dist_n==min_dist_n)
    
    # Coordinates for first (positive) point 
    p1_x = x [mask_p1][0]
    p1_y = y [mask_p1][0]
    p1_z = z [mask_p1][0]
    print ('Set point 1 at: x='+str(p1_x)+', y='+str(p1_y)+', z='+str(p1_z))
    
    # Coordinates for second (negative) point 
    p2_x = x [mask_p2][0]
    p2_y = y [mask_p2][0]
    p2_z = z [mask_p2][0]
    print ('Set point 2 at: x='+str(p2_x)+', y='+str(p2_y)+', z='+str(p2_z))
    
    # update new distance for the found points to have the correct strain 
    distance = np.linalg.norm(np.array([p1_x,p1_y,p1_z])-np.array([p2_x,p2_y,p2_z]))
    print ('Distance between found points is updated for correct strain. Distance: '+str(distance*1e6)+'µm')
    
    
    
    
    mask_outer = distance_c > r_outer * 0.999  # include close elements to the boundary 
    
        
    #Set Boundaries-----------------------------------------------------------------------------------------------------
    bcond = np.zeros((len(coords), 4))
    bcond[:, 3] = 1.
    
    
    #set displ mode at outer boundary (fixed displacements)
    bcond[mask_outer, 3] = 0
        
    
    #set displ mode at both nodes (fixed displacements)
    bcond[mask_p1, 3] = 0   
    bcond[mask_p2, 3] = 0  
    
    # total deformation
    deformation = strain*distance
    
    
    # project (0,0,x_defo) on the disctance vector because points are not perfectly alligned on xaxis
    vec_ideal = np.array([deformation/2,0 ,0])   # ideal deformation in positive direction for negative point 2
    vec_calc = np.array([p1_x,p1_y,p1_z])-np.array([p2_x,p2_y,p2_z]) #pointing in positive direction
            
    
    projection = vec_calc * (vec_ideal @ vec_calc) / (vec_calc @ vec_calc)  # projection for p2 pointing in positive direction

    
    #fixed displacements for both points directed to center 
    bcond[mask_p1, 0] =   -projection[0]
    bcond[mask_p1, 1] =   -projection[1]
    bcond[mask_p1, 2] =   -projection[2]
    
    bcond[mask_p2, 0] =   projection[0]
    bcond[mask_p2, 1] =   projection[1]
    bcond[mask_p2, 2] =   projection[2]

    
    #save bcond.dat
    np.savetxt(simulation_folder + '/bcond.dat', bcond)
    
    #create iconf
    iconf = np.zeros((len(coords), 3))
    iconf[mask_p1, 0] =   -deformation/2    
    iconf[mask_p2, 0] =   +deformation/2
    


    
    np.savetxt(simulation_folder + '/iconf.dat', iconf)
    print('+ Created bcond.dat')
      
    
    # create config.txt -----------------------------------------------------------
    config = r"""
    MODE = relaxation
    BOXMESH = 0
    FIBERPATTERNMATCHING = 0
    REL_CONV_CRIT = {}
    REL_ITERATIONS = {}
    REL_SOLVER_STEP = {}
    K_0 = {}
    D_0 = {}
    L_S = {}
    D_S = {}
    CONFIG = {}\config.txt
    DATAOUT = {}
    """.format(conv_crit, iterations, step, K_0, D_0, L_S, D_S, simulation_folder, simulation_folder)
    
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
    TOTAL_NODES = {}
    Mesh_file = {}
    Output_folder = {}
    distance = {}
    strain = {}
    total_deformation (half applied at each point) = {}
    point_1 =  {}
    point_2 =  {}
    iterations = {}
    step = {}
    conv_crit = {}
      
    """.format(K_0, D_0, L_S, D_S, str(r_outer*1e6)+' µm',  len(coords), mesh_file, simulation_folder, 
    str(distance*1e6)+' µm', strain, deformation,
      'x='+str(p1_x)+', y='+str(p1_y)+', z='+str(p1_z) , 
      'x='+str(p2_x)+', y='+str(p2_y)+', z='+str(p2_z), 
      iterations, step, conv_crit)
    
    with open(simulation_folder + "/parameters.txt", "w") as f:
          f.write(parameters)
    print('+ Created parameters.txt')
           
    # start SAENO ----------------------------------------------------------------------------------------------------------
    print ('SAENO Path: '+str(SAENOPATH))
    print('+ Starting SAENO...')
     
    # Create log file if activated
    if logfile == True:
       
        # create log file with system output
        logfile = open(simulation_folder + "/saeno_log.txt", 'w')
        cmd = subprocess.Popen(SAENOPATH+" CONFIG {}/config.txt".format(os.path.abspath(simulation_folder)), stdout=subprocess.PIPE , universal_newlines=True, shell=False)
        # print and save a reduced version of saeno log
        for line in cmd.stdout:
            if not '%' in line:
                print (line, end='')
                logfile.write(str(line))
        # close again to avoid loops            
        cmd.stdout.close()            
       
    # if false just show the non reduced system output    
    else:
        cmd = subprocess.call(SAENOPATH+" CONFIG {}/config.txt".format(os.path.abspath(simulation_folder))) 
     
  
        
        
def spherical_contraction(meshfile, simulation_folder, pressure, material, r_inner=100, r_outer=20000, logfile = False, half_sphere = False, iterations= 300 , step=0.3, conv_crit = 1e-11):
    """
    Simulates an spherical contraction of the inner inclusion for a given pressure.
    
    Args:
        simulation_folder(str): File path to save simulation results
        mesh_file(str): File path to read in mesh file
        r_inner(float): Radius of inner inclusion
        r_outer(float): Outer radius of the bulk mesh in the mesh model (in µm)
        pressure(float): Pressure that is applied on inner inclusion
        saeno (str): File path to Saeno.exe
        material (dict): Material properties in the form {'K_0': X, 'D_0':X, 'L_S': X, 'D_S': X} (see materials)
        logfile(boolean): If True a reduced logfile of the saeno system output is stored. Default: False.
        iterations(float): The maximal number of iterations for the saeno simulation. Default: 300.
        step(float): Step width parameter for saeno regularization. Higher values lead to a faster but less robust convergence. Default: 0.3.
        conv_crit(float): Saeno stops if the relative standard deviation of the residuum is below given threshold. Default: 1e-11.      
        half_sphere(boolean): If true all forces on spheroid a halved (to compare spheroid with half sphere)  WATCH OUT: SPHEROID POSITION+SHAPE IS NOT FIXED IF TRUE
    """
      
    # open mesh file
    with open(meshfile, 'r') as f:
        lines = f.readlines()
    
    # scale radii to meter
    r_inner *= 10**-6
    r_outer *= 10**-6

    # read in material parameters
    K_0 = material['K_0']
    D_0 = material['D_0']
    L_S = material['L_S']
    D_S = material['D_S']

    # create output folder if it does not exist, print warning otherwise
    if not os.path.exists(simulation_folder):
        os.makedirs(simulation_folder)
    else:
        print('WARNING: Output folder already exists! ({})'.format(simulation_folder))

    # transform nodes and connection in SAENO format
    # nodes
    index_nodes = lines.index('$Nodes\n')
    n_nodes = int(lines[index_nodes + 1])

    coords = np.zeros((n_nodes, 3))
    for i in range(n_nodes):
        coords[i] = np.array([np.float(x) for x in lines[i + index_nodes + 2].split()[1:]])
    np.savetxt(simulation_folder + '/coords.dat', coords)

    # connections
    index_elements = lines.index('$Elements\n')
    n_elements = int(lines[index_elements + 1])

    tets = np.zeros((n_elements, 4))
    for i in range(n_elements):
        tets[i] = lines[i + index_elements + 2].split()[-4:]
    np.savetxt(simulation_folder + '/tets.dat', tets, fmt='%i')

    # define boundary conditions
    distance = np.sqrt(np.sum(coords ** 2., axis=1))
    mask_inner = distance < r_inner * 1.001
    mask_outer = distance > r_outer * 0.999

    # Save Node Density at inner and outer sphere
    # Area per inner node
    A_node_inner = (np.pi*4*(r_inner)**2)/np.sum(mask_inner) 
    # simple sqrt as spacing
    inner_spacing = np.sqrt(A_node_inner)    
    
    # Area per outer node
    A_node_outer = (np.pi*4*(r_outer)**2)/np.sum(mask_outer)   
    # simple sqrt as spacing
    outer_spacing = np.sqrt(A_node_outer)
    
    print ('Inner node spacing: '+str(inner_spacing*1e6)+'µm')
    print ('Outer node spacing: '+str(outer_spacing*1e6)+'µm')


    bcond = np.zeros((len(coords), 4))
    bcond[:, 3] = 1.

    # fixed displacements for outer boundary
    bcond[mask_outer, 3] = 0

    # fixed non-zero force at spheroid surface
    bcond[mask_inner, :3] = coords[mask_inner, :3]
    bcond[mask_inner, :3] /= distance[mask_inner, None]

    # halven forces if half sphere True
    if half_sphere == True:

        A_inner = 4 * np.pi * r_inner ** 2.
        force_per_node = pressure * A_inner / np.sum(mask_inner)
        bcond[mask_inner, :3] *= force_per_node/2
     
    # compute regular forces per node for a full sphere    
    else:
        A_inner = 4 * np.pi * r_inner ** 2.
        force_per_node = pressure * A_inner / np.sum(mask_inner)
        bcond[mask_inner, :3] *= force_per_node


    np.savetxt(simulation_folder + '/bcond.dat', bcond)

    # define initial configuration
    iconf = np.zeros((len(coords), 3))
    np.savetxt(simulation_folder + '/iconf.dat', iconf)

    # create config file for SAENO
    config = r"""MODE = relaxation
    BOXMESH = 0
    FIBERPATTERNMATCHING = 0
    REL_CONV_CRIT = {}
    REL_ITERATIONS = {}
    REL_SOLVER_STEP = {}
    K_0 = {}
    D_0 = {}
    L_S = {}
    D_S = {}
    CONFIG = {}\config.txt
    DATAOUT = {}""".format(conv_crit, iterations, step, K_0, D_0, L_S, D_S, os.path.abspath(simulation_folder), os.path.abspath(simulation_folder))

    with open(simulation_folder + "/config.txt", "w") as f:
        f.write(config)

    # create info file with all relevant parameters of the simulation
    parameters = r"""K_0 = {}
    D_0 = {}
    L_S = {}
    D_S = {}
    PRESSURE = {}
    FORCE_PER_SURFACE_NODE = {}
    INNER_RADIUS = {} µm
    OUTER_RADIUS = {} µm
    INNER_NODE_SPACING = {} µm
    OUTER_NODE_SPACING = {} µm
    SURFACE_NODES = {}
    TOTAL_NODES = {}
    REL_CONV_CRIT = {}
    REL_ITERATIONS = {}
    REL_SOLVER_STEP = {}
    half_sphere = {}""".format(K_0, D_0, L_S, D_S, pressure, force_per_node, r_inner*1e6 , r_outer*1e6 , inner_spacing*1e6, outer_spacing*1e6, np.sum(mask_inner), len(coords), conv_crit, iterations, step, half_sphere)

    with open(simulation_folder + "/parameters.txt", "w") as f:
        f.write(parameters)
        
        
     # Create log file if activated
    if logfile == True:
        
        # create log file with system output
        logfile = open(simulation_folder + "/saeno_log.txt", 'w')
        cmd = subprocess.Popen(SAENOPATH+" CONFIG {}/config.txt".format(os.path.abspath(simulation_folder)), stdout=subprocess.PIPE , universal_newlines=True, shell=False)
        # print and save a reduced version of saeno log
        for line in cmd.stdout:
            if not '%' in line:
                print (line, end='')
                logfile.write(str(line))
        # close again to avoid loops            
        cmd.stdout.close()            
        
    # if false just show the non reduced system output    
    else:
        cmd = subprocess.call(SAENOPATH+" CONFIG {}/config.txt".format(os.path.abspath(simulation_folder)))   



     
        
        
     
def spherical_contraction_strain(meshfile, simulation_folder, strain, material, r_inner=100, r_outer=20000, logfile = False, iterations= 300 , step=0.3, conv_crit = 1e-11):
    """
    Simulates an spherical contraction of the inner inclusion for a given pressure.
    
    Args:
        simulation_folder(str): File path to save simulation results
        mesh_file(str): File path to read in mesh file
        r_inner(float): Radius of inner inclusion
        r_outer(float): Outer radius of the bulk mesh in the mesh model (in µm)
        strain(float): Strain that is applied on inner inclusion
        saeno (str): File path to Saeno.exe
        material (dict): Material properties in the form {'K_0': X, 'D_0':X, 'L_S': X, 'D_S': X} (see materials)
        logfile(boolean): If True a reduced logfile of the saeno system output is stored. Default: False.
        iterations(float): The maximal number of iterations for the saeno simulation. Default: 300.
        step(float): Step width parameter for saeno regularization. Higher values lead to a faster but less robust convergence. Default: 0.3.
        conv_crit(float): Saeno stops if the relative standard deviation of the residuum is below given threshold. Default: 1e-11.      
     """
      
    # open mesh file
    with open(meshfile, 'r') as f:
        lines = f.readlines()
    
    # scale radii to meter
    r_inner *= 10**-6
    r_outer *= 10**-6

    # read in material parameters
    K_0 = material['K_0']
    D_0 = material['D_0']
    L_S = material['L_S']
    D_S = material['D_S']

    # create output folder if it does not exist, print warning otherwise
    if not os.path.exists(simulation_folder):
        os.makedirs(simulation_folder)
    else:
        print('WARNING: Output folder already exists! ({})'.format(simulation_folder))

    # transform nodes and connection in SAENO format
    # nodes
    index_nodes = lines.index('$Nodes\n')
    n_nodes = int(lines[index_nodes + 1])

    coords = np.zeros((n_nodes, 3))
    for i in range(n_nodes):
        coords[i] = np.array([np.float(x) for x in lines[i + index_nodes + 2].split()[1:]])
    np.savetxt(simulation_folder + '/coords.dat', coords)

    # connections
    index_elements = lines.index('$Elements\n')
    n_elements = int(lines[index_elements + 1])

    tets = np.zeros((n_elements, 4))
    for i in range(n_elements):
        tets[i] = lines[i + index_elements + 2].split()[-4:]
    np.savetxt(simulation_folder + '/tets.dat', tets, fmt='%i')

    # define boundary conditions
    distance = np.sqrt(np.sum(coords ** 2., axis=1))
    mask_inner = distance < r_inner * 1.001
    mask_outer = distance > r_outer * 0.999

    # Save Node Density at inner and outer sphere
    # Area per inner node
    A_node_inner = (np.pi*4*(r_inner)**2)/np.sum(mask_inner) 
    # simple sqrt as spacing
    inner_spacing = np.sqrt(A_node_inner)    
    
    # Area per outer node
    A_node_outer = (np.pi*4*(r_outer)**2)/np.sum(mask_outer)   
    # simple sqrt as spacing
    outer_spacing = np.sqrt(A_node_outer)
    
    print ('Inner node spacing: '+str(inner_spacing*1e6)+'µm')
    print ('Outer node spacing: '+str(outer_spacing*1e6)+'µm')


    bcond = np.zeros((len(coords), 4))
    bcond[:, 3] = 1.
    
    # fixed displacements for outer boundary
    bcond[mask_outer, 3] = 0

    # fixed displacements at spheroid surface
    bcond[mask_inner, 3] = 0  # displacement mode
    bcond[mask_inner, :3] = coords[mask_inner, :3]
    bcond[mask_inner, :3] /= distance[mask_inner, None]  # unit vectors
    deformation = strain * r_inner
    bcond[mask_inner, :3] *= -deformation
 

    np.savetxt(simulation_folder + '/bcond.dat', bcond)

    # define initial configuration,  same as bconds   
    iconf = np.zeros((len(coords), 3))
    iconf[mask_inner, :3] = coords[mask_inner, :3]
    iconf[mask_inner, :3] /= distance[mask_inner, None]  # unit vectors
    iconf[mask_inner, :3] *= -deformation
    
    
    np.savetxt(simulation_folder + '/iconf.dat', iconf)

    # create config file for SAENO
    config = r"""MODE = relaxation
    BOXMESH = 0
    FIBERPATTERNMATCHING = 0
    REL_CONV_CRIT = {}
    REL_ITERATIONS = {}
    REL_SOLVER_STEP = {}
    K_0 = {}
    D_0 = {}
    L_S = {}
    D_S = {}
    CONFIG = {}\config.txt
    DATAOUT = {}""".format(conv_crit, iterations, step, K_0, D_0, L_S, D_S, os.path.abspath(simulation_folder), os.path.abspath(simulation_folder))

    with open(simulation_folder + "/config.txt", "w") as f:
        f.write(config)

    # create info file with all relevant parameters of the simulation
    parameters = r"""K_0 = {}
    D_0 = {}
    L_S = {}
    D_S = {}
    STRAIN = {}
    DEFORMATION = {}
    INNER_RADIUS = {} µm
    OUTER_RADIUS = {} µm
    INNER_NODE_SPACING = {} µm
    OUTER_NODE_SPACING = {} µm
    SURFACE_NODES = {}
    TOTAL_NODES = {}
    REL_CONV_CRIT = {}
    REL_ITERATIONS = {}
    REL_SOLVER_STEP = {}
   """.format(K_0, D_0, L_S, D_S, strain, deformation, r_inner*1e6 , r_outer*1e6 , inner_spacing*1e6, outer_spacing*1e6, np.sum(mask_inner), len(coords), conv_crit, iterations, step)

    with open(simulation_folder + "/parameters.txt", "w") as f:
        f.write(parameters)
        
        
     # Create log file if activated
    if logfile == True:
        
        # create log file with system output
        logfile = open(simulation_folder + "/saeno_log.txt", 'w')
        cmd = subprocess.Popen(SAENOPATH+" CONFIG {}/config.txt".format(os.path.abspath(simulation_folder)), stdout=subprocess.PIPE , universal_newlines=True, shell=False)
        # print and save a reduced version of saeno log
        for line in cmd.stdout:
            if not '%' in line:
                print (line, end='')
                logfile.write(str(line))
        # close again to avoid loops            
        cmd.stdout.close()            
        
    # if false just show the non reduced system output    
    else:
        cmd = subprocess.call(SAENOPATH+" CONFIG {}/config.txt".format(os.path.abspath(simulation_folder)))          
       





    
#if __name__ == "__main__":
    
    

    

