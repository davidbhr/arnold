from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
import glob as glob
import os
from PIL import Image
from tqdm import tqdm
from arnold import SAENOPATH
import subprocess

def convert_to_png(input_folder, output_folder, img_format = "*.tif"):
    """
    Converts all images in input_folder with ending "img_format" to png and 
    saves them to output folder.
    
    On some systems saeno has troubles to read in .tif-files. In that case images
    can be converted to .png using this function.
    """
    # scan images in folder with img_format extension
    img_list  =  glob.glob(os.path.join(input_folder,img_format))
   
    # create outputfolder if not exists
    if not os.path.exists(output_folder):
     os.makedirs(output_folder)
     
    for i in tqdm(img_list):
        im = Image.open(i)
        # new output file 
        newfile = os.path.join(  output_folder ,  os.path.splitext(os.path.basename(i))[0] )+".png"
        # save conversion 
        im.save(newfile, quality=100)

    return img_list



def regularization_stacks(simulation_folder,    # path to output folder
                          stack_deformed, 
                          stack_relaxed,   # eg  C/../Series001_z*_ch00.tif  
                          K_0,              # material properties for simulation
                          L_S, 
                          D_0,
                          D_S, 
                          voxelsize_x,     # voxelsize ( in m - e.g use formt like 0.3611e-6 )
                          voxelsize_y,
                          voxelsize_z, 
                          AUTOMATIC_Drift_Values = True,  # if True these are chosen automatically
                          drift_range = None,     # ~20 voxelsize                                    - ignored if automatic
                          drift_step = None,       # ~ 2 voxelsize                                   - ignored if automatic
                          AUTOMATIC_BM_Values = True,    # if True these are chosen automatically
                          BM_GRAIN = None,    # ~ 10x voxelsize                                      - ignored if automatic
                          BM_MULOUT=0,  # can be use to make mesh more coarse further away           - ignored if automatic
                          BM_N = None,       #  BM_GRAIN * BM_N must cover total image               - ignored if automatic
                          BM_RIN=0,    # can be use to make mesh more coarse further away            - ignored if automatic
                          FM_RMAX = None,    # define maximal distance where deformation values are   - ignored if automatic
                                      # used for force reconstruction  
                          
                          REG_CONV_CRIT=0.01, # regularization
                          REG_ITERATIONS=50,
                          REG_LAPLACEGRAIN=15.0,
                          REG_SOLVER_PRECISION=1e-18,
                          REG_SOLVER_STEP=0.33,
                          REL_CONV_CRIT=0.01, # relaxation
                          REL_ITERATIONS= 300,
                          REL_SOLVER_STEP = 0.066, 
                          SUBTRACTMEDIANDISPL=0,
                          VB_MINMATCH=0.7,
                          VB_REGPARA=0.01,
                          VB_REGPARAREF=0.1,
                          VB_SX=12,
                          VB_SY=12, 
                          VB_SZ=12  ):
    """
    Starts a full 3D saeno regularization for the given relaxed and alive input  image stacks.
    
    Simulation folder is created. The deformed and relaxed image stacks must be 
    specified to read in a list (analog to  : C/../Series001_z*_ch00.tif ). 
    
    Voxelsize must be given in m. Parameters are described above and in the SAENO
    documentation. Material properties specified with K_0, L_S,D_0,D_S,
    

    
    USEFULL: for 512x512 input images we use BM_Grain which is approx.  10 x voxelsize (a bit larger)
    and BM_N of 54 to capture the total image.
    If  AUTOMATIC_BM_Values = True the mesh fineness is applied automatically
    analog to the upper example for the corresponding image dimensions. Use
    AUTOMATIC_Drift_Values = True to change the lengths in the same manner
    according to the image dimensions.   Only works for quadratic image shapes.
    """
    # create outputfolder if not exists---------------------------------------
    if not os.path.exists(simulation_folder):
     os.makedirs(simulation_folder)
 
     
    # automatic parameter adjustment if chosen--------------------------------- 
     
    if AUTOMATIC_Drift_Values: # x y images should have equal length 
        example_image_for_shape = Image.open(glob.glob(stack_deformed)[0])
        dimension = example_image_for_shape.size[0]  #
        print("Image dimension:" + str(example_image_for_shape.size))
        drift_range = 20 * voxelsize_x     # ~20 voxelsize                                  
        drift_step =  2 * voxelsize_x     # ~ 2 voxelsiz
        print("Automatic value for drift_range (20x voxelsize):" + str(drift_range))
        print("Automatic value for drift_step (2x voxelsize):" + str(drift_step))
        
        
    if AUTOMATIC_BM_Values:  # x y images should have equal length 
        example_image_for_shape = Image.open(glob.glob(stack_deformed)[0])
        dimension = example_image_for_shape.size[0]  
        print("Image dimension:" + str(example_image_for_shape.size)+"(check if x-y has equal length for automatic parameters)")
       # fineness
        BM_GRAIN =  10 * voxelsize_x    # ~ 10x voxelsize                                      
        # Tolat amount  -  BM_GRAIN * BM_N must cover all pixels  
        # -> BM_N = (dimension * voxelsize_x) / BM_GRAIN  add some buffer and convert to integer          
        BM_N = int (dimension * voxelsize_x  / BM_GRAIN) + 4
        # maximal distance (to center of image) for deformations that are used for force reconstruction
        FM_RMAX = 0.5 * (dimension * voxelsize_x)  # corresponds to total image except without the "buffer" nodes
        
        print("Automatic value for BM_GRAIN (10x voxelsize):" + str(BM_GRAIN))
        print("Automatic value for BM_N (scaled to include total image + small buffer):" + str(BM_N))
        print("Automatic value for FM_RMAX (total image except the 'buffer' nodes):" + str(FM_RMAX))
        
     
    # create config.txt -------------------------------------------------------
    config = r"""
    MODE = regularization

    DATAOUT = {}
    STACKA = {}
    STACKR = {}
    
    K_0 = {}
    L_S = {}
    D_0 = {}
    D_S = {}
    
    VOXELSIZEX = {}
    VOXELSIZEY = {}
    VOXELSIZEZ = {}
    
    ALLIGNSTACKS = 1
    DRIFTCORRECTION = 1
    DRIFT_RANGE = {}
    DRIFT_STEP = {}
    
    ALPHA = 3.0e9
    BEAMS = 300
    
    BOXMESH = 1
    BM_GRAIN = {}
    BM_MULOUT = {}
    BM_N = {}
    BM_RIN = {}
    
    EPSMAX = 4.0
    EPSSTEP = 0.000001
    FIBERPATTERNMATCHING = 1
    FM_RMAX = {}
    
    REFINEDISPLACEMENTS = 0
    REGMETHOD = robust
    REG_CONV_CRIT = {}
    REG_ITERATIONS = {}
    REG_LAPLACEGRAIN = {}
    REG_SIGMAZ = 1.0
    REG_SOLVER_PRECISION = {}
    REG_SOLVER_STEP = {}
    
    REL_CONV_CRIT = {}
    REL_ITERATIONS = {}
    REL_SOLVER_STEP = {}
    
    ROBUSTMETHOD = huber
    SAVEALLIGNEDSTACK = 0
    SAVEEPSILON = 1
    SUBPIXEL = 0.0005
    SUBTRACTMEDIANDISPL =  {}
    USESPRINTF = 0
    
    VB_MINMATCH = {}
    VB_N = 1
    VB_REGPARA = {}
    VB_REGPARAREF = {}
    VB_SX = {}
    VB_SY = {}
    VB_SZ = {}
    """.format(simulation_folder, stack_deformed, stack_relaxed, 
            K_0, L_S, D_0, D_S, voxelsize_x,voxelsize_y,voxelsize_z, 
            drift_range, drift_step,BM_GRAIN, BM_MULOUT, BM_N,BM_RIN,
             FM_RMAX,REG_CONV_CRIT, REG_ITERATIONS, REG_LAPLACEGRAIN,
             REG_SOLVER_PRECISION, REG_SOLVER_STEP,
             REL_CONV_CRIT, REL_ITERATIONS,REL_SOLVER_STEP, 
             SUBTRACTMEDIANDISPL,VB_MINMATCH, VB_REGPARA, VB_REGPARAREF,
             VB_SX,VB_SY, VB_SZ )

    with open(simulation_folder + "/config.txt", "w") as f:
        f.write(config)
    print('+ Created config.txt')
 
    #sytem call to start simulation
    subprocess.call(SAENOPATH+"/saeno CONFIG {}/config.txt".format(os.path.abspath(simulation_folder))) 
    return






def display_displacements(simulation_folder, lower_percentile=0, upper_percentile=100, save_plot = True):
    """
    show the total displacement field or the masked range of displacements by size
    simulation_folder: path to simulation folder
    lower_percentile and upper_percentile give range of deformationsize which is used
    as a mask 
    save_plot: option to save the plot to the simulation folder
    
    To increase plotting speed you might increase the lower percentile (e.g 30 instead 0)

    """
    # load in deformations and coordinates
    r =  np.genfromtxt(os.path.join(simulation_folder, "R.dat"))       # positions
    u =  np.genfromtxt(os.path.join(simulation_folder, "Ufound.dat"))     #deformations  
    uabs = np.sqrt(np.sum(u ** 2., axis=1))   # absolute values for filtering
    #S = np.genfromtxt('../Sfound.dat')          # accuracy of deformations

    # filter displacements by absolute size
    mask = (uabs > np.percentile(uabs, lower_percentile)) & (uabs < np.percentile(uabs, upper_percentile))
    r2 = r[mask]
    u2 = u[mask]
    uabs2 = uabs[mask]
    
    # plot the masked deformation
    fig = plt.figure()
    ax1 = fig.gca(projection='3d', label='fitted-displ', rasterized=True)

    color_bounds1 =  np.array([np.percentile(uabs2, 0), np.percentile(uabs2, 100)]) * 10 ** 6  #[0, 2.5]  #
    
    u_scale = 1.0e-06/np.median(uabs2)
    
    np.random.seed(1234)
    for r2i, u2i, uabs2i in tqdm(zip(r2 * 10 ** 6, u2 * 10 ** 6, uabs2 * (10 ** 6)*u_scale)):
        # if np.random.uniform(0, 1) < 0.2:  # uabs2i/u_upper:
        color = plt.cm.jet(((uabs2i - color_bounds1[0]) / (color_bounds1[1] - color_bounds1[0])))
        alpha = 1. - (r2i[0] - r2i[1]) / (270. * 0.5)
    
        if alpha > 1:
            alpha = 1.
        if alpha < 0:
            alpha = 0.
    
        plt.quiver(r2i[0], r2i[1], r2i[2], u2i[0], u2i[1], u2i[2], length=uabs2i * 4,
                   color=color, arrow_length_ratio=0, alpha=alpha, pivot='tip', linewidth=0.5)
    
    #ax1.set_xlim([-135, 135])
    #ax1.set_ylim([-135, 135])
    #ax1.set_zlim([-135, 135])
    #ax1.set_xticks([-100, -50, 0, 50, 100])
    #ax1.set_yticks([-100, -50, 0, 50, 100])
    #ax1.set_zticks([-100, -50, 0, 50, 100])
    #ax1.set_xticklabels(['']*5)
    #ax1.set_yticklabels(['']*5)
    #ax1.set_zticklabels(['']*5)
    ax1.w_xaxis.set_pane_color((0.2, 0.2, 0.2, 1.0))
    ax1.w_yaxis.set_pane_color((0.2, 0.2, 0.2, 1.0))
    ax1.w_zaxis.set_pane_color((0.2, 0.2, 0.2, 1.0))
    
    if save_plot:
        plt.savefig( os.path.join(simulation_folder,'deformations_plot_lower_{}_upper_{}.png'.format(lower_percentile, upper_percentile)), dpi=600 )
    
     
    return
    

    

def filter_displacements():
    """
    filters out deformations by size and stores the resulting files in output
    """ 
    
# TODO  add paths and stuff
    
    r =  np.genfromtxt('../R.dat')           # position
    u =  np.genfromtxt('../Ufound.dat')     #deformations  
    uabs = np.sqrt(np.sum(u ** 2., axis=1))   # absolute values for filtering
    S = np.genfromtxt('../Sfound.dat')          # accuracy of deformations
    
    # mask deformation between 60 and 99.7 percent
    mask = (uabs > np.percentile(uabs, 60)) & (uabs < np.percentile(uabs, 99.7))
    
    # set masked values to zero
    u2 = u.copy()
    u2[~mask] = 0
    S2 = S.copy()
    S2[~mask] = 0
    
    # save filtered data for new simulation
    np.savetxt('Rfound.dat',r) # die 3 musst du dann für die neue Simulationin den Ordner der Config legen und Config wie besprochen ändern
    np.savetxt('Ufound.dat',u2)
    np.savetxt('Sfound.dat',S2)
        
    
    
    
    
    return






def regularization_files():
    """
    starts a full 3D saeno regularization with deformations from .dat-files
    (in particular to restart the filtered deformations from filter_displacements() )
    """
    return

def force_displacement_fields():
    """
    visualizes and saves the force and displacement field of a simulation
    """
    return


if __name__ == "__main__":    
    
    
    # Test display stacks
    #display_displacements(r"\\131.188.117.96\biophysDS\dboehringer\Platte_4\Eval-data\Example-Data-test\nk-testing\2019-10-15-nk21-series88\out-arnold-test", 30,100)    
    display_displacements(r"\\131.188.117.96\biophysDS\dboehringer\Platte_4\Eval-data\Example-Data-test\nk-testing\2019-10-15-nk21-series88\out-arnold-test", 0,99.7) 
    
    
    
    
    
    # Test regularization_stacks
    
    # stacka =  r"\\131.188.117.96\biophysDS\dboehringer\Platte_4\Measurements_NK_TFM\2019-10-15_Konfokal-Primarynk\nk21-1_zoom6_eval\88\Series088_t14_z*_ch00.tif"
    # stackb = r"\\131.188.117.96\biophysDS\dboehringer\Platte_4\Measurements_NK_TFM\2019-10-15_Konfokal-Primarynk\nk21-1_zoom6_eval\89_rel\Series089_t12_z*_ch00.tif"
    # out= r"\\131.188.117.96\biophysDS\dboehringer\Platte_4\Eval-data\Example-Data-test\nk-testing\2019-10-15-nk21-series88\out-arnold-test"
    
    # regularization_stacks(out,stacka,stackb,    
    #                     K_0=1645,              # material properties for simulation
    #                     L_S=0.0075, 
    #                     D_0=0.0008,
    #                     D_S=0.033, 
    #                     voxelsize_x=0.241e-6,     # voxelsize ( in m - e.g use formt like 0.3611e-6 )
    #                     voxelsize_y=0.241e-6,
    #                     voxelsize_z=1.007e-6, 
    #                     AUTOMATIC_Drift_Values = True,  # if True these are chosen automatically                  
    #                     AUTOMATIC_BM_Values = True,    # if True these are chosen automatically
#                                 )
    
    
    
    
    
    
    
   















#   stack a  \\131.188.117.96\biophysDS\dboehringer\Platte_4\Measurements_NK_TFM\2019-10-15_Konfokal-Primarynk\nk21-1_zoom6_eval\88\Series088_t14_z*_ch00.tif
#   stack r  \\131.188.117.96\biophysDS\dboehringer\Platte_4\Measurements_NK_TFM\2019-10-15_Konfokal-Primarynk\nk21-1_zoom6_eval\89_rel\Series089_t12_z*_ch00.tif

#stacka= r'\\131.188.117.96\biophysDS\dboehringer\Platte_4\Measurements_NK_TFM\2019-10-15_Konfokal-Primarynk\nk21-1_zoom6_eval\88\Series088_t14_z*_ch00.tif'
