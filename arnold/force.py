from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FormatStrFormatter
import pandas as pd




def scaling_law(d,l,s,E):
    """
    Reconstruct the contractility of a cylindric inclusion for different geometry, strain and material properties
    
    Args:
        d(float): Diameter of the cyldindric inclusion (muscle fiber) (in µm)
        l(float): Length of the cyldindric inclusion (muscle fiber) (in µm)
        s(float): Strain of the cylindric inclusion (e.g. muscle fiber) in percent 
        (Derived as: [(length_relaxed - length_contracted)/(len_relaxed)])
        E(float): Young's Modulus of the surrounding Material (in Pascal)
      
    """   
         
    d0 = 30.
    l0 = 300.
    s0 = 0.1 # equals 10%
    E0 = 200.  
            
    Contractility = (  4.553581636199e-07 * (E/E0) * (l/l0) * ((0.5 * l/l0) + ((0.5 *d/d0))) * (s/s0)) * 1e6          
    
    print ('Contractility = '+str(Contractility)+'µN')
    return (Contractility)




def reconstruct_contractility(simulation_folder, d_cyl, l_cyl, r_outer,   scalef = 1000, scaleu = 1, scaleb = 1):
    """
    Reconstruct the contractility of a given cylindric simulation. Also calculates residuum forces of Inclusion surface and matrix, 
    which shoould be equal (use to certify simulation and detect possible errors due incorrect lengths as input).
    
    Saves the contractility calculated from absolute values, the contractility calculated only from 
    x-components and the residuum forces to .txt files and a quiver plot is as .png file.
       
    Args:
        simulation_folder(str): File path to the simulation
        d_cyl(float): Diameter of the spherical inclusion in the mesh model (in µm) 
        l_cyl(float): Length of the spherical inclusion in the mesh model (in µm)
        r_outer(float): Outer radius of the bulk mesh in the mesh model (in µm)   
        scalef ,scalu, scaleb: To scale the arrows for deformation , force and 
        boundary cond. in quiver plot - only visually no impact on valeus
         
    """ 

    # Read in results                   
    coords = np.genfromtxt(simulation_folder+'/coords.dat')                            
    U = np.genfromtxt(simulation_folder+'/U.dat')  
    force =  np.genfromtxt(simulation_folder+'/F.dat')
    bcond =  np.genfromtxt(simulation_folder+'/bcond.dat')  
    
    # Mask Boundaries on which the tractions are summed up
    x, y, z = coords.T
    mask_inner = ((y**2. + z**2.) <= (d_cyl*1e-6/2)**2.) & (x**2. <= (l_cyl*1e-6/2)**2.)
    mask_outer = ~mask_inner
    # Left and right half of cell 
    mask_inner_l = ((y**2. + z**2.) <= (d_cyl*1e-6/2)**2.) & (x**2. <= (l_cyl*1e-6/2)**2.) & (x<=0)
    mask_inner_r = ((y**2. + z**2.) <= (d_cyl*1e-6/2)**2.) & (x**2. <= (l_cyl*1e-6/2)**2.) & (x>=0)
    
    "Display Results in 3D Plot"
    fig1 = plt.figure()
    ax1 = Axes3D(fig1)
    ax1.set_xlabel('X [\u03bcm]', labelpad=6, fontsize = 13)
    ax1.set_ylabel('Y [\u03bcm]', labelpad=6, fontsize = 13)
    ax1.set_zlabel('Z [\u03bcm]', labelpad=6, fontsize = 13)
    # Display Deformation, Forces and selected Boundary
    u_x = np.array(U[:,0])
    u_y = np.array(U[:,1])
    u_z = np.array(U[:,2])
    f_x = np.array(force[:,0])
    f_y = np.array(force[:,1])
    f_z = np.array(force[:,2])
    bcond_x = np.array(bcond[:,0])
    bcond_y = np.array(bcond[:,1])
    bcond_z = np.array(bcond[:,2])
    #Scale quiver lengths 
    scale_um = 1e6   # for axis in um
    # Deformationfield
    ax1.quiver(coords[:,0]*scale_um,coords[:,1]*scale_um,coords[:,2]*scale_um,
                u_x*scaleu*scale_um,u_y*scaleu*scale_um,u_z*scaleu*scale_um, 
                lw=2,length = 1, colors = 'cadetblue' , normalize=False, label='Displacement')
    # Forcefield
    ax1.quiver(coords[:,0]*scale_um,coords[:,1]*scale_um,coords[:,2]*scale_um,
                f_x*scalef*scale_um,f_y*scalef*scale_um,f_z*scalef*scale_um, 
                lw=2,length = 1, colors = 'orange' ,normalize=False,label='Force')

       
    ax1.set_xlim(-(l_cyl/2)*1.1, (l_cyl/2)*1.1)   
    ax1.set_ylim(-(l_cyl/2)*1.1, (l_cyl/2)*1.1)          
    ax1.set_zlim(-(l_cyl/2)*1.1, (l_cyl/2)*1.1)
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
    ax1.zaxis.set_major_formatter(FormatStrFormatter('%3.0f')) 
    plt.legend(fontsize = 15) 
    plt.savefig(simulation_folder+"\Force_Displ_Field.png")
    plt.close()
    
           
    # Plot 2: Left / Right Forces and boundary conditions
    fig2 = plt.figure()
    ax2 = Axes3D(fig2)
    ax2.set_xlabel('X [\u03bcm]', labelpad=14, fontsize = 15)
    ax2.set_ylabel('Y [\u03bcm]', labelpad=14, fontsize = 15)
    ax2.set_zlabel('Z [\u03bcm]', labelpad=14, fontsize = 15)
    
    # Left and right forces
    ax2.quiver(coords[mask_inner_l][:,0]*scale_um,coords[mask_inner_l][:,1]*scale_um,coords[mask_inner_l][:,2]*scale_um,
      f_x[mask_inner_l]*scalef*scale_um,f_y[mask_inner_l]*scalef*scale_um,f_z[mask_inner_l]*scalef*scale_um, 
               lw=2,length = 1, colors = 'r' ,normalize=False,label='Force (Left)')
    
    ax2.quiver(coords[mask_inner_r][:,0]*scale_um,coords[mask_inner_r][:,1]*scale_um,coords[mask_inner_r][:,2]*scale_um,
      f_x[mask_inner_r]*scalef*scale_um,f_y[mask_inner_r]*scalef*scale_um,f_z[mask_inner_r]*scalef*scale_um, 
               lw=2,length = 1, colors = 'b' ,normalize=False,label='Force (Right)')
    
    # Boundary Cond.
    ax2.quiver(coords[:,0]*scale_um,coords[:,1]*scale_um,coords[:,2]*scale_um,
                bcond_x*scaleb*scale_um,bcond_y*scaleb*scale_um,bcond_z*scaleb*scale_um, 
                lw=2,length = 1, colors = 'grey' ,normalize=False,label='Boundary Cond.' , alpha= 0.1 , zorder= -1)
             

    
    ax2.set_xlim(-(l_cyl/2)*1.1, (l_cyl/2)*1.1)   
    ax2.set_ylim(-(l_cyl/2)*1.1, (l_cyl/2)*1.1)          
    ax2.set_zlim(-(l_cyl/2)*1.1, (l_cyl/2)*1.1) 
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
    ax2.zaxis.set_major_formatter(FormatStrFormatter('%3.0f')) 
    
    plt.legend(fontsize = 15) 
    plt.savefig(simulation_folder+"\SplitForce_bcond.png")
    plt.close()




    # Plot 3: Distance Deformation Curve --linscale
    fig3 = plt.figure()

    u_plot = np.sqrt(np.sum(coords ** 2., axis=1)) * 1e6
    v_plot = np.sqrt(np.sum(U ** 2., axis=1)) * 1e6
    plt.scatter(u_plot, v_plot, lw=0, alpha=1, s=10, c='C2')

    plt.xlabel('Distance [µm]')
    plt.ylabel('Deformation [µm]')

    plt.savefig(simulation_folder+'\Deformation_Distance_Lin.png', format='png', dpi=700)
    plt.close()


    # Plot 3: Distance Deformation Curve --logscale
    fig4 = plt.figure()

    print (np.min(v_plot))

    plt.scatter(u_plot, v_plot, lw=0, alpha=1, s=10, c='C2')

    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(10**-10, np.max(v_plot)+10)  # fixed, as autoscale does not work with log scale
    plt.xlabel('Distance [µm]')
    plt.ylabel('Deformation [µm]')
    plt.savefig(simulation_folder + '\Deformation_Distance_Log.png', format='png', dpi=700)
    plt.close()


    "Compute Contractilities"
    # initialize result dictionary
    results = {'Contractility Absolute (left)': [], 'Contractility Absolute (right)': [],
               'Contractility Absolute (mean)': [], 'Contractility x-components (left)': [],
               'Contractility x-components (right)': [], 'Contractility x-components (mean)': [],
               'Contractility Absolute (total length)': [], 'Residuum Forces': [],
               'Residuum Inner': [], 'Residuum Outer': []}
    
    # Residuum Forces: Outer and inner mask
    ResOut = np.sum(np.sum(force[mask_outer]))
    ResIn =  np.sum(np.sum(force[mask_inner]))
    # Absolute Contractility computed for left and right site 
    Contr_abs_l = np.sum(np.sqrt(np.sum(force[mask_inner_l]**2., axis=1)))
    Contr_abs_r = np.sum(np.sqrt(np.sum(force[mask_inner_r]**2., axis=1)))
    # Contractility computed by x-components for left and right site 
    Contr_x_l = np.sum(np.abs(force[mask_inner_l][:,0]))
    Contr_x_r = np.sum(np.abs(force[mask_inner_r][:,0]))
    
    results['Contractility Absolute (left)'].append(Contr_abs_l)
    results['Contractility Absolute (right)'].append(Contr_abs_r)
    results['Contractility Absolute (mean)'].append(0.5*(Contr_abs_l+Contr_abs_r))
    results['Contractility x-components (left)'].append(Contr_x_l)
    results['Contractility x-components (right)'].append(Contr_x_r)  
    results['Contractility x-components (mean)'].append(0.5*(Contr_x_l+Contr_x_r)) 
    results['Contractility Absolute (total length)'].append(np.sum(np.abs(force[mask_inner][:,0])))  
    results['Residuum Inner'].append(ResIn)  
    results['Residuum Outer'].append(ResOut)  
    results['Residuum Forces'].append(ResOut+ResIn) 
    
    df = pd.DataFrame.from_dict(results)
    df.columns = ['Contractility Absolute (left) [N]',
                  'Contractility Absolute (right) [N]',
                  'Contractility Absolute (mean) [N]',
                  'Contractility x-components (left) [N]',
                  'Contractility x-components (right) [N]',
                  'Contractility x-components (mean) [N]',
                  'Contractility Absolute (total length) [N]',
                  'Residuum Inner [N]',
                  'Residuum Outer [N]',
                  'Residuum Forces [N]'        ]
 
    df.to_excel(simulation_folder+"\Contractilities.xlsx")
    
    return df

    
    






def reconstruct_pointforces(simulation_folder,  scalef = 2000, scaleu = 8):
    """
    Evaluation of a point force simulations      
    Only input is simulatio folder. Scaling only visual for quiver plot scale f for force arrows and scaleu for deformation arrows.
    """ 
    import csv
    # Read in results                   
    coords = np.genfromtxt(simulation_folder+'/coords.dat')                            
    u = np.genfromtxt(simulation_folder+'/U.dat')  
    f =  np.genfromtxt(simulation_folder+'/F.dat')
    f_abs = np.sqrt(np.sum(f**2., axis=1))*1e6
    # rather complicated way to read in data from parameters.txt
    txt = open(simulation_folder+'/parameters.txt') 
    reader = csv.reader(txt)
    filelist=[]
    for line in reader:
        for i in line:
            filelist.append(i.replace(' ','').split('='))   
    table = pd.DataFrame(filelist) 
    distance = float(np.array(table[1][table[0] == 'distance'])[0][:-2]) # in um

    # get force maxima for two pointnodes and third as an control
    f1 = f_abs[f_abs==np.unique(f_abs)[-1]]
    f2 = f_abs[f_abs==np.unique(f_abs)[-2]]
    f3 = f_abs[f_abs==np.unique(f_abs)[-3]]
    mean_force = np.mean([f1,f2])
    print ('Force maxium  1+2 (µN) = point forces: ' +  str(f1) + ',  ' +  str(f2) )
    print ('Mean Force (µN): ' + str(mean_force) )   
    print ('Control - Third maximum (µN): '+ str(f3) )
    np.savetxt(simulation_folder+'\MeanForce-uN.txt', [mean_force])
    np.savetxt(simulation_folder+'\max1-pointforce-uN.txt', [f1])
    np.savetxt(simulation_folder+'\max2-pointforce-uN.txt', [f2])
    np.savetxt(simulation_folder+'\max3-control-uN.txt', [f3])


    # Plot Deformation distance curves
    # logscale
    plt.figure()
    plt.subplot(121)
    plt.xscale('log')
    plt.yscale('log')
    u1 = np.sqrt(np.sum(coords**2., axis=1))*1e6
    v1 = np.sqrt(np.sum(u**2., axis=1))*1e6
    plt.scatter(u1, v1, lw=0, alpha=1, s=7, label='Displ.') 
    plt.xlabel('Distance r [µm]')
    plt.ylabel('Deformation u [µm]')
    plt.xlim([10**0, 10**3.5])
    plt.ylim([10**-12.5 *1e6, 10])  
    # linscale
    plt.subplot(122)
    # Deformation distance curves
    u1 = np.sqrt(np.sum(coords**2., axis=1))*1e6
    v1 = np.sqrt(np.sum(u**2., axis=1))*1e6
    plt.scatter(u1, v1, lw=0, alpha=1, s=10, label='Displ.') 
    plt.xlabel('Distance r [µm]')
    plt.ylabel('Deformation u [µm]')    
    # save plot
    plt.tight_layout()
    plt.savefig(simulation_folder+'/deformations.png', dpi=300)
    plt.close()


    # Plot Forcce over distance
    plt.figure() 
    plt.scatter(u1, f_abs, lw=0, alpha=1, s=10) 
    plt.xlabel('Distance r [µm]')
    plt.ylabel('Force [µN]') 
    plt.tight_layout()
    plt.savefig(simulation_folder+'/force-distance.png', dpi=300)
    plt.close()

   
    # Plot 3D Force Disp Field
    fig1 = plt.figure()
    ax1 = Axes3D(fig1)
    ax1.set_xlabel('X [\u03bcm]', labelpad=6, fontsize = 13)
    ax1.set_ylabel('Y [\u03bcm]', labelpad=6, fontsize = 13)
    ax1.set_zlabel('Z [\u03bcm]', labelpad=6, fontsize = 13)
    # # Display Deformation, Forces and selected Boundary
    u_x = np.array(u[:,0])
    u_y = np.array(u[:,1])
    u_z = np.array(u[:,2])
    f_x = np.array(f[:,0])
    f_y = np.array(f[:,1])
    f_z = np.array(f[:,2])
    #Scale quiver lengths
    scale_um = 1e6   # for axis in um
    shown_dist = distance 
    # Deformationfield
    ax1.quiver(coords[:,0]*scale_um,coords[:,1]*scale_um,coords[:,2]*scale_um,
                u_x*scaleu*scale_um,u_y*scaleu*scale_um,u_z*scaleu*scale_um, 
                lw=2,length = 1, colors = 'cadetblue' , normalize=False, label='Displacement')
    # Forcefield
    ax1.quiver(coords[:,0]*scale_um,coords[:,1]*scale_um,coords[:,2]*scale_um,
                f_x*scalef*scale_um,f_y*scalef*scale_um,f_z*scalef*scale_um, 
                lw=2,length = 1, colors = 'orange' ,normalize=False,label='Force')

    ax1.set_xlim(-shown_dist,shown_dist)   
    ax1.set_ylim(-shown_dist,shown_dist)          
    ax1.set_zlim(-shown_dist,shown_dist)   
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
    ax1.zaxis.set_major_formatter(FormatStrFormatter('%3.0f')) 
    plt.legend(fontsize = 15) 
    plt.savefig(simulation_folder+'/Force_disp_field.png', dpi=300)
    
    
    return mean_force
 
