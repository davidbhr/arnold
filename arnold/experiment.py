from __future__ import print_function, division
from time import sleep
import arnold as ar
from subprocess import Popen
import os
import numpy as np
import pandas as pd



def mesh_and_show(mesh_file, d_cyl, l_cyl, r_outer, length_factor=0.4):
    
    """
    Creates a spherical bulk mesh with a centered cylindrical inclusion and shows the mesh afterwards 
    
    Args:
        mesh_file(str): File path to save mesh file (add .msh ending)
        d_cyl(float): Diameter of the cylindrical inclusion (in µm)
        l_cyl(float): Length of the cylindrical inclusion (in µm)
        r_outer(float): Outer radius of the bulk mesh (in µm)
        length_factor(float): Mesh element size is determined by curvature and then multipled with this factor
    """
        
    ar.mesh.cylindrical_inclusion(mesh_file, d_cyl, l_cyl, r_outer, length_factor)
        
    ar.mesh.show_mesh(mesh_file)
        
    return
    
    
    
def cylindrical_inclusion_mesh_and_simulation(mesh_file, d_cyl, l_cyl, r_outer, length_factor, simulation_folder, 
                                              strain, material, logfile = False,  iterations= 300 , step=0.3, conv_crit = 1e-11):  
    """
    Creates a spherical bulk mesh with a centered cylindrical inclusion and simulates a symetric contraction 
    (with constant strain) of the cylindric inclusion 
    
    Args:
        mesh_file(str): File path to save mesh file (add .msh ending)
        d_cyl(float): Diameter of the cylindrical inclusion (in µm)
        l_cyl(float): Length of the cylindrical inclusion (in µm)
        r_outer(float): Outer radius of the bulk mesh (in µm)
        length_factor(float): Mesh element size is determined by curvature and then multipled with this factor
        simulation_folder(str): File path to save simulation results
        strain(float): Strain as (Length_base - Length_contracted)/Length_base. 
        Deformation is applied in x-direction and split equally to both poles,
        from which defomation-size decreases linearly to the center (symetric contraction with constant strain). 
        Deformation can be determed as (Length_base - Length_contracted) 
        material (dict): Material properties in the form {'K_0': X, 'D_0':X, 'L_S': X, 'D_S': X}, see materials)
        logfile(boolean): If True a reduced logfile of the saeno system output is stored. Default: False.
        iterations(float): The maximal number of iterations for the saeno simulation. Default: 300.
        step(float): Step width parameter for saeno regularization. Higher values lead to a faster but less robust convergence. Default: 0.3.
        conv_crit(float): Saeno stops if the relative standard deviation of the residuum is below given threshold. Default: 1e-11.      
    """
   
    
    
    ar.mesh.cylindrical_inclusion(mesh_file, d_cyl, l_cyl, r_outer, length_factor)    
    
    # wait a bit longer until mesh is stored 
    sleep(5)  
    
    ar.simulation.cylindric_contraction(simulation_folder, mesh_file, d_cyl, l_cyl, r_outer, strain, material, 
                                        logfile = logfile,  iterations= iterations , step=step, conv_crit = conv_crit)
   
    
    return
        
        
def cylindrical_inclusion_mesh_simulation_and_contractility(mesh_file, d_cyl, l_cyl, r_outer, length_factor, simulation_folder, 
                                              strain, material, logfile = False,  iterations= 300 , step=0.3, conv_crit = 1e-11):  
    """
    Creates a spherical bulk mesh with a centered cylindrical inclusion and simulates a symetric contraction 
    (with constant strain) of the cylindric inclusion and computes the contractile forces.
    
    Args:
        mesh_file(str): File path to save mesh file (add .msh ending)
        d_cyl(float): Diameter of the cylindrical inclusion (in µm)
        l_cyl(float): Length of the cylindrical inclusion (in µm)
        r_outer(float): Outer radius of the bulk mesh (in µm)
        length_factor(float): Mesh element size is determined by curvature and then multipled with this factor
        simulation_folder(str): File path to save simulation results
        strain(float): Strain as (Length_base - Length_contracted)/Length_base. 
        Deformation is applied in x-direction and split equally to both poles,
        from which defomation-size decreases linearly to the center (symetric contraction with constant strain). 
        material (dict): Material properties in the form {'K_0': X, 'D_0':X, 'L_S': X, 'D_S': X}, see materials)
        logfile(boolean): If True a reduced logfile of the saeno system output is stored. Default: False.
        iterations(float): The maximal number of iterations for the saeno simulation. Default: 300.
        step(float): Step width parameter for saeno regularization. Higher values lead to a faster but less robust convergence. Default: 0.3.
        conv_crit(float): Saeno stops if the relative standard deviation of the residuum is below given threshold. Default: 1e-11.      
    """


    # build mesh
    ar.mesh.cylindrical_inclusion(mesh_file, d_cyl, l_cyl, r_outer, length_factor)   
    
    # wait a bit longer until mesh is stored 
    sleep(5)  
    
    # start saeno simulation
    ar.simulation.cylindric_contraction(simulation_folder, mesh_file, d_cyl, l_cyl, r_outer, strain, material,
                                        logfile = logfile,  iterations= iterations , step=step, conv_crit = conv_crit)
    
    # wait a bit  
    sleep(1)  
   
    # compute individual contractilities
    ar.force.reconstruct_contractility(simulation_folder, d_cyl, l_cyl, r_outer)
    
    return
        

def cylindrical_inclusion_simulation_and_contractility(mesh_file, d_cyl, l_cyl, r_outer, simulation_folder, 
                                              strain, material, logfile = False,  iterations= 300 , step=0.3, conv_crit = 1e-11):  
    """
    Simulates a symetric contraction (with constant strain) of the cylindric inclusion and
    computes the contractile forces.
    
    Args:
        mesh_file(str): File path to save mesh file (add .msh ending)
        d_cyl(float): Diameter of the cylindrical inclusion (in µm)
        l_cyl(float): Length of the cylindrical inclusion (in µm)
        r_outer(float): Outer radius of the bulk mesh (in µm)
        simulation_folder(str): File path to save simulation results
        strain(float): Strain as (Length_base - Length_contracted)/Length_base. 
        Deformation is applied in x-direction and split equally to both poles,
        from which defomation-size decreases linearly to the center (symetric contraction with constant strain). 
        material (dict): Material properties in the form {'K_0': X, 'D_0':X, 'L_S': X, 'D_S': X}, see materials)
        logfile(boolean): If True a reduced logfile of the saeno system output is stored. Default: False.
        iterations(float): The maximal number of iterations for the saeno simulation. Default: 300.
        step(float): Step width parameter for saeno regularization. Higher values lead to a faster but less robust convergence. Default: 0.3.
        conv_crit(float): Saeno stops if the relative standard deviation of the residuum is below given threshold. Default: 1e-11.      
    """


    ar.simulation.cylindric_contraction(simulation_folder, mesh_file, d_cyl, l_cyl, r_outer, strain, material,
                                        logfile = logfile,  iterations= iterations , step=step, conv_crit = conv_crit)
   
    ar.force.reconstruct_contractility(simulation_folder, d_cyl, l_cyl, r_outer)
    
    return
        


def simulation_series_lengths(d_cyl, l_cyl_min, l_cyl_max, n, r_outer, length_factor, simulation_folder, strain, material,  
                              log_scaling=True, n_cores=None, dec=10, logfile = True,  iterations= 300 , step=0.3, conv_crit = 1e-11):
    """
    Starts a series of simulation for different fiber lengths and evaluates fiber contractility
    
     Args:
        d_cyl(float): Diameter of the cylindrical inclusion (in µm)
        l_cyl_min(float): Minmal length of the cylindrical inclusion (in µm) 
        l_cyl_max(float): Maximal length of the cylindrical inclusion (in µm) 
        n(float): Number of simulates to be made between minimal and maximal length 
        r_outer(float): Outer radius of the bulk mesh (in µm)
        length_factor(float): Mesh element size is determined by curvature and then multipled with this factor
        simulation_folder(str): File path to save simulation results
        strain(float): Strain as (Length_base - Length_contracted)/Length_base. 
        Deformation is applied in x-direction and split equally to both poles,
        from which defomation-size decreases linearly to the center (symetric contraction with constant strain). 
        log_scaling(boolean): Logarithmic or Linear scaling between cell lengths
        n_cores(float): Amount of simultanious processes to be started, default detects the amount of CPU cores
        dec(int): Decimal value to round the lengths,  default is 10    
        material (dict): Material properties in the form {'K_0': X, 'D_0':X, 'L_S': X, 'D_S': X}, see materials)
        logfile(boolean): If True a reduced logfile of the saeno system output is stored. Default: True.
        iterations(float): The maximal number of iterations for the saeno simulation. Default: 300.
        step(float): Step width parameter for saeno regularization. Higher values lead to a faster but less robust convergence. Default: 0.3.
        conv_crit(float): Saeno stops if the relative standard deviation of the residuum is below given threshold. Default: 1e-11.      
    """
              
    
    # detect number of cores for paralell computing 
    if n_cores is None:
        n_cores = os.cpu_count()
        print(str(n_cores)+' cores detected')

    # List of lengths in logarithmic steps to start simulation
    if log_scaling:
        L = np.logspace(np.log10(l_cyl_min), np.log10(l_cyl_max), num=n, endpoint=True)
   
    # List of lengths in linear steps to start simulation   
    else:
        L = np.linspace(np.log10(l_cyl_min), np.log10(l_cyl_max), num=n, endpoint=True)
        
    # Limit decimal     
    L = np.around(L, decimals=dec)    
      
    print('Lengths: '+str(L))

    L = list(L)
    
    # create output folder if it doesn't exist
    if not os.path.exists(simulation_folder):
        os.makedirs(simulation_folder)
 
    
    # Start single simulations in paralell
    processes = []
    while True:
        if len(L) == 0:
            break
        
        if len(processes) < n_cores:

            command =  '''python -c "import arnold as ar; ar.experiment.cylindrical_inclusion_mesh_simulation_and_contractility(r'{}',{},{},{},{},r'{}',{},{},{},{},{},{})"'''.format(simulation_folder+'\\'+str(L[0])+'.msh', d_cyl, 
                                                                                                                                L[0], r_outer, length_factor, simulation_folder+'\\'+str(L[0]), strain, material,
                                                                                                                                logfile, iterations , step, conv_crit )
                        
            processes.append(Popen(command))

            del L[0]
        
        sleep(1.)
        
        processes = [p for p in processes if p.poll() is None]

        
 
def simulation_series_diameter(d_cyl_min, d_cyl_max, l_cyl, n, r_outer, length_factor, simulation_folder, strain, material,  
                               log_scaling=True, n_cores=None, dec=2 , logfile = True,  iterations= 300 , step=0.3, conv_crit = 1e-11):
    """
    Starts a series of simulation for different fiber diameters and evaluates the fiber contractility
    
     Args:
        d_cyl_min(float): Minimal diameter of the cylindrical inclusion (in µm)
        d_cyl_max(float): Maximal diameter of the cylindrical inclusion (in µm) 
        l_cyl(float): Length of the cylindrical inclusion (in µm)
        n(float): Number of simulates to be made between minimal and maximal length 
        r_outer(float): Outer radius of the bulk mesh (in µm)
        length_factor(float): Mesh element size is determined by curvature and then multipled with this factor
        simulation_folder(str): File path to save simulation results
        Strain(float): Strain applied on the muscle fiber (Length_base - Length_contracted/Length_base) 
        # Total deformation of the cylindric inclusion (in µm). Deformation is applied in x-direction 
        # and split equally to both poles, from which defomation-size decreases linearly to the center
        # (symetric contraction with constant strain). Positive values denote a contraction. 
        # Deformation can be determed as (Length_base - Length_contracted) 
        log_scaling(boolean): Logarithmic or Linear scaling between cell lengths
        n_cores(float): Amount of simultanious processes to be started, default detects the amount of CPU cores
        dec(int): Decimal value to round the lengths,  default is 10    
        material (dict): Material properties in the form {'K_0': X, 'D_0':X, 'L_S': X, 'D_S': X}, see materials)
        logfile(boolean): If True a reduced logfile of the saeno system output is stored. Default: True.
        iterations(float): The maximal number of iterations for the saeno simulation. Default: 300.
        step(float): Step width parameter for saeno regularization. Higher values lead to a faster but less robust convergence. Default: 0.3.
        conv_crit(float): Saeno stops if the relative standard deviation of the residuum is below given threshold. Default: 1e-11.      
    """
              
    
    # detect number of cores for paralell computing 
    if n_cores is None:
        n_cores = os.cpu_count()
        print(str(n_cores)+' cores detected')

    # List of lengths in logarithmic steps to start simulation
    if log_scaling:
        d = np.logspace(np.log10(d_cyl_min), np.log10(d_cyl_max), num=n, endpoint=True)
   
    # List of lengths in linear steps to start simulation   
    else:
        d = np.linspace(np.log10(d_cyl_min), np.log10(d_cyl_max), num=n, endpoint=True)
        
    # Limit decimal     
    d = np.around(d, decimals=dec)    
      
    print('Diameters: '+str(d))

    d = list(d)
    
    # create output folder if it doesn't exist
    if not os.path.exists(simulation_folder):
        os.makedirs(simulation_folder)
 
    
    # Start single simulations in paralell
    processes = []
    while True:
        if len(d) == 0:
            break
        
        if len(processes) < n_cores:

            command =  '''python -c "import arnold as ar; ar.experiment.cylindrical_inclusion_mesh_simulation_and_contractility(r'{}',{},{},{},{},r'{}',{},{},{},{},{},{})"'''.format(simulation_folder+'\\'+str(d[0])+'.msh', d[0],  l_cyl, r_outer,
                                                                                                                                length_factor, simulation_folder+'\\'+str(d[0]), strain, material,
                                                                                                                                logfile, iterations , step, conv_crit )
                        
            processes.append(Popen(command))

            del d[0]
        
        sleep(1.)
        
        processes = [p for p in processes if p.poll() is None]
    

    
    
def simulation_series_strain(d_cyl, l_cyl, n, r_outer, length_factor, simulation_folder, strain_min, strain_max, material,  
                             log_scaling=True, n_cores=None, dec=2 , logfile = True,  iterations= 300 , step=0.3, conv_crit = 1e-11):
    """
    Starts a series of simulation for different cell lengths and evaluates cell contractility
    
     Args:
        d_cyl(float): Diameter of the cylindrical inclusion (in µm) 
        l_cyl(float): Length of the cylindrical inclusion (in µm)
        n(float): Number of simulates to be made between minimal and maximal Diameter 
        r_outer(float): Outer radius of the bulk mesh (in µm)
        length_factor(float): Mesh element size is determined by curvature and then multipled with this factor
        simulation_folder(str): File path to save simulation results
        strain_min(float): Minimal strain applied on the muscle fiber (Length_base - Length_contracted/Length_base)
        strain_max(float): Maximal strain applied on the muscle fiber (Length_base - Length_contracted/Length_base)
        # Total deformation of the cylindric inclusion (in µm). Deformation is applied in x-direction 
        # and split equally to both poles, from which defomation-size decreases linearly to the center
        # (symetric contraction with constant strain). Positive values denote a contraction. 
        # Deformation can be determed as (Length_base - Length_contracted) 
        log_scaling(boolean): Logarithmic or Linear scaling between cell lengths
        n_cores(float): Amount of simultanious processes to be started, default detects the amount of CPU cores
        dec(int): Decimal value to round the lengths,  default is 10    
        material (dict): Material properties in the form {'K_0': X, 'D_0':X, 'L_S': X, 'D_S': X}, see materials)
        logfile(boolean): If True a reduced logfile of the saeno system output is stored. Default: True.
        iterations(float): The maximal number of iterations for the saeno simulation. Default: 300.
        step(float): Step width parameter for saeno regularization. Higher values lead to a faster but less robust convergence. Default: 0.3.
        conv_crit(float): Saeno stops if the relative standard deviation of the residuum is below given threshold. Default: 1e-11.      
    """
     

    # detect number of cores for paralell computing 
    if n_cores is None:
        n_cores = os.cpu_count()
        print(str(n_cores)+' cores detected')

    # List of lengths in logarithmic steps to start simulation
    if log_scaling:
        e = np.logspace(np.log10(strain_min), np.log10(strain_max), num=n, endpoint=True)
   
    # List of lengths in linear steps to start simulation   
    else:
        e = np.linspace(np.log10(strain_min), np.log10(strain_max), num=n, endpoint=True)
        
    # Limit decimal     
    e = np.around(e, decimals=dec)    
      
    print('Strains: '+str(e))   
    
    # create output folder if it doesn't exist  DOPPELT
    if not os.path.exists(simulation_folder):
        os.makedirs(simulation_folder)
 
    e = list(e)
    # Start single simulations in paralell
    processes = []
    while True:
        if len(e) == 0:
            break
        
        if len(processes) < n_cores:

            command =  '''python -c "import arnold as ar; ar.experiment.cylindrical_inclusion_mesh_simulation_and_contractility(r'{}',{},{},{},{},r'{}',{},{},{},{},{},{})"'''.format(simulation_folder+'\\'+str(e[0])+'.msh', d_cyl,  l_cyl,
                                                                                                                                r_outer, length_factor, simulation_folder+'\\'+str(e[0]), e[0], material,
                                                                                                                                logfile, iterations , step, conv_crit )
                        
            processes.append(Popen(command))

            del e[0]
        
        sleep(1.)
    
        processes = [p for p in processes if p.poll() is None]
    
    
   
def simulation_series_stiffness(d_cyl, l_cyl, n, r_outer, length_factor, simulation_folder, strain, material_min, material_max, 
                                log_scaling=True, n_cores=None, dec=2 , logfile = True,  iterations= 300 , step=0.3, conv_crit = 1e-11):
    """
    Starts a series of simulation for different cell lengths and evaluates cell contractility
    
     Args:
        d_cyl(float): Diameter of the cylindrical inclusion (in µm) 
        l_cyl(float): Length of the cylindrical inclusion (in µm)
        n(float): Number of simulates to be made between minimal and maximal Diameter 
        r_outer(float): Outer radius of the bulk mesh (in µm)
        length_factor(float): Mesh element size is determined by curvature and then multipled with this factor
        simulation_folder(str): File path to save simulation results
        strain_(float): Strain applied on the muscle fiber (Length_base - Length_contracted/Length_base)
        # Total deformation of the cylindric inclusion (in µm). Deformation is applied in x-direction 
        # and split equally to both poles, from which defomation-size decreases linearly to the center
        # (symetric contraction with constant strain). Positive values denote a contraction. 
        # Deformation can be determed as (Length_base - Length_contracted) 
        log_scaling(boolean): Logarithmic or Linear scaling between cell lengths
        n_cores(float): Amount of simultanious processes to be started, default detects the amount of CPU cores
        dec(int): Decimal value to round the lengths,  default is 10    
        material_min (dict): Material properties in the form {'K_0': X, 'D_0':X, 'L_S': X, 'D_S': X}, see materials) ranging from min K_0 stiffness to max K_0 stiffness
        material_max (dict): Material properties in the form {'K_0': X, 'D_0':X, 'L_S': X, 'D_S': X}, see materials) ranging from min K_0 stiffness to max K_0 stiffness
        logfile(boolean): If True a reduced logfile of the saeno system output is stored. Default: True.
        iterations(float): The maximal number of iterations for the saeno simulation. Default: 300.
        step(float): Step width parameter for saeno regularization. Higher values lead to a faster but less robust convergence. Default: 0.3.
        conv_crit(float): Saeno stops if the relative standard deviation of the residuum is below given threshold. Default: 1e-11.      
    """
     
    


    # detect number of cores for paralell computing 
    if n_cores is None:
        n_cores = os.cpu_count()
        print(str(n_cores)+' cores detected')

    # List of lengths in logarithmic steps to start simulation   
    if log_scaling:
        k = np.logspace(np.log10(material_min['K_0']), np.log10(material_max['K_0']), num=n, endpoint=True)
   
    # List of lengths in linear steps to start simulation   
    else:
        k = np.linspace(np.log10(material_min['K_0']), np.log10(material_max['K_0']), num=n, endpoint=True)
        
    # Limit decimal     
    k = np.around(k, decimals=dec)    
      
    print('Strains: '+str(k))   
    
    # create output folder if it doesn't exist
    if not os.path.exists(simulation_folder):
        os.makedirs(simulation_folder)
 
    
    k = list(k)
    
    
    # Start single simulations in paralell
    processes = []
    while True:
        if len(k) == 0:
            break
        
        if len(processes) < n_cores:
            
            
            material_min['K_0'] = str(k[0])
            command =  '''python -c "import arnold as ar; ar.experiment.cylindrical_inclusion_mesh_simulation_and_contractility(r'{}',{},{},{},{},r'{}',{},{},{},{},{},{})"'''.format(simulation_folder+'\\'+str(k[0])+'.msh', d_cyl, 
                                                                                                                                l_cyl, r_outer, length_factor, simulation_folder+'\\'+str(k[0]), strain, material_min,
                                                                                                                                logfile, iterations , step, conv_crit )
                        
            processes.append(Popen(command))

            del k[0]
        
        sleep(1.)
    
        processes = [p for p in processes if p.poll() is None]
    



    
def evaluate_series(path, comment=''):
    """
    Evaluates several simulations and stores a overview in an excel file
    
     Args:
         path(str): File path to search subfolders for simulations
         comment(str): Comment which is added to the evaluation file containing all simulations of all found subfolders 
    """
   
    
    # initialize summary dictionary
    results = {'Diameter [µm]': [], 'Length [µm]': [], 'Strain [%]': [],'Contractility Absolute (mean) [N]': [], 'Contractility x-components (mean) [N]': [],
                'Residuum Forces [N]': [], 'K_0': [], 'D_0': [],'L_S': [],'D_S': [], 'Deformation [µm]': [], 'Simulation Folder': []}
    
    
    
    print ('Files found:')
    for root, dirs, files in os.walk(path):
        for name in files:
            if name.endswith(("Contractilities.xlsx")):
                       
                print (str(root)+"\\parameters.txt") 
                print (str(root)+"\\Contractilities.xlsx")
                
                pd.set_option('display.max_colwidth', -1) # use full length strings in pandas                
                            
                Parameters = pd.read_csv(root+"\\parameters.txt", delimiter='='  ,  encoding = "ISO-8859-1", header=None)
                Parameters.columns = ['name','value']
                
                
                # Append values from parameter.txt
                results['K_0'].append(float(Parameters[Parameters['name'].str.contains('K_0')]['value']))    # K0 value
                results['D_0'].append(float(Parameters[Parameters['name'].str.contains('D_0')]['value']))             
                results['L_S'].append(float(Parameters[Parameters['name'].str.contains('L_S')]['value']))              
                results['D_S'].append(float(Parameters[Parameters['name'].str.contains('D_S')]['value']))              
                results['Diameter [µm]'].append(float(str(Parameters[Parameters['name'].str.contains('d_cyl')]['value']).split()[1]))   # one more split needed due to um string inside..
                length = float(str(Parameters[Parameters['name'].str.contains('l_cyl')]['value']).split()[1])
                results['Length [µm]'].append(length)   
                deformation = float(str(Parameters[Parameters['name'].str.contains('deformation')]['value']).split()[1])
                results['Deformation [µm]'].append(deformation)   # one more split needed due to um string inside..
                results['Strain [%]'].append((deformation/length)*100) 
                #results['Simulation Folder'].append(str(Parameters[Parameters['name'].str.contains('Output_folder')]['value']).split()[1])    #simulation folder from paramteres
                results['Simulation Folder'].append(str(root))   # in case the name is manually changed after simulation
                
                
                # Append values from Contractilities.xlsx
                Contractilities = pd.read_excel(root+"\\Contractilities.xlsx")
                results['Contractility Absolute (mean) [N]'].append(float(Contractilities['Contractility Absolute (mean) [N]']))
                results['Contractility x-components (mean) [N]'].append(float(Contractilities['Contractility x-components (mean) [N]']))
                results['Residuum Forces [N]'].append(float(Contractilities['Residuum Forces [N]']))
    
           
            
    df = pd.DataFrame.from_dict(results)
    df.columns = ['Diameter [µm]',
     'Length [µm]',
     'Strain [%]',
     'Contractility Absolute (mean) [N]', 
     'Contractility x-components (mean) [N]',
     'Residuum Forces [N]',
     'K_0', 
     'D_0',
     'L_S',
     'D_S',
     'Deformation [µm]',
     'Simulation Folder']
     
    df.to_excel(path+"\\series_evaluation_"+str(comment)+".xlsx")
    
    return df        



  
        


