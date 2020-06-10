import os
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import imageio
from PIL import Image
import os
import roipoly
import glob as glob
from tqdm import tqdm
from scipy.ndimage import gaussian_filter



# maximum projections
def projections(files, output, name, zrange=None, mode="max", percentile_mode =99.5, gauss_filter = None):
    """
    Compute and store maximum intensity projections of a given stack to combine z information within a single image
    
    files: list of file paths for the max projection (glob.glob() might be usefull to create such list )
    output: path to the folder to store projections
    name: name of output projection image
    zrange: Defines the number of imgas around the z-plane which are combined within the  output
        default is None, which combines all input data (if overexposed might use as smaller z_range)
    projection: different mode for "max", "mean", "median" and "percentile"
    percentile_mode: percentile for projection -  only valid if mode = percentile
    """
    
    stacksize = len(files)
    if zrange == None:
        zrange = len(files)
   
    zmin=int(stacksize/2-(zrange/2))
    zmax=int(stacksize/2+(zrange/2))
             
             
    # creates folder if it doesn't exist
    if not os.path.exists(output):
        os.makedirs(output)

    # iterate over files, read img (.tif), append to list
    im_list = []
    for im in files:
        im_array = np.array(imageio.imread(im))
        im_list.append(im_array)
    # make list to array
    images_array = np.array(im_list)

    # limit z range
    images_array = images_array[zmin:zmax,:,:]
    
    # make the chosen projection
    if mode == "max":
        proj = np.max(images_array,axis=0)
    if mode == "percentile":
        proj = np.percentile(images_array, percentile_mode, axis=0)
    if mode == "mean":    
        proj = np.mean(images_array,axis=0)  
    if mode == "median":       
        proj = np.median(images_array,axis=0)
    # create names    
    if mode == "max":
       filename = 'maxproj_' + str(name) 
    if mode == "percentile":
        filename = 'percentile_'+ str(percentile_mode)+"_" + str(name)
    if mode == "mean":    
        filename = 'meanproj_' + str(name) 
    if mode == "median":       
        filename = 'mdianproj_' + str(name) 
    
    # apply gaussian if active
    if gauss_filter:
        proj =  gaussian_filter(proj, sigma=gauss_filter)
        filename = filename+'_gauss_'+str(gauss_filter)+'.tif'
    else:
        filename = filename+'.tif'

       
    im = Image.fromarray(proj)    
       
    savepath = os.path.join(output, filename)
    im.save(savepath)

    return  proj




def piv_analysis( contr,relax, outfolder,scale, winsize_um = 10 , overlap_um =None, searchsize_um = None, drift_correction =True, threshold=1.2, scale_quiver=None):
    """
    Computes deformations between 2 images by crosscorrelation using openpiv (must be installed - see Readme).
    Saves several quiver plots and the maximal found deformation for later analysis.
    
    contr: Path to active/contracted image file to calculate deformations between contracted and relaxed state
    relax: Path to relaxed image file to calculate deformations between contracted and relaxed state
    scale: resolution of image in µm per pixel
    winsize_um: size of the search window to be applied in µm 
    drift_correction: Applies a drift correction before piv analysis
    threshold: filters displacement
    scale_quiver: can be used to scale the arrows in quiver plot (only for visualization)
                Default is None meaning automatically scaling ,  scale is inverse
  
    """ 
    from openpiv import tools, process, validation, filters, scaling 


    winsize = int(winsize_um/scale)  # for pixel 
    
     # if not specified use defaults
    if not overlap_um:
        overlap = winsize/2
    else:
        overlap = int(overlap_um/scale)
        
    if not searchsize_um:
        searchsize = 2*winsize
    else:
        searchsize = int(searchsize_um/scale)
    
    print (winsize,overlap,searchsize)
    # odd winsize raise problems due to the overlap, thefore decrease by one pixel if odd and give warning
    if not (winsize % 2) == 0:
        print('Odd pixelnumbers raise problems due to the overlap: Winsize changed to '+str( (winsize-1) * scale) +' um')
        winsize -= 1
    
    # creates folder if it doesn't exist
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
        
    #read in images  
    relax = tools.imread(relax)
    contr = tools.imread(contr)
    
    # convert for openpiv
    relax =  relax.astype(np.int32)
    contr =  contr.astype(np.int32)
  
    dt = 1 # sec

    u0, v0, sig2noise = process.extended_search_area_piv( relax , contr, window_size=winsize, 
                                                         overlap=overlap, dt=dt, search_area_size=searchsize, sig2noise_method='peak2peak' )

    x, y = process.get_coordinates( image_size=relax.shape, window_size=winsize, overlap=overlap )
    
    # drift correction
    if drift_correction:
        u0 -= np.nanmean(u0)
        v0 -= np.nanmean(v0)
    
    # filtered deformations
    u1, v1, mask = validation.sig2noise_val( u0, v0, sig2noise, threshold = threshold)
    
    # save deformations + coords
    np.save(outfolder+"/u.npy", u1) 
    np.save(outfolder+"/v.npy", v1) 
    np.save(outfolder+"/x.npy", x) 
    np.save(outfolder+"/y.npy", y) 
    
    # show filtered+unfiltered data    
    plt.figure(figsize=(6,3))
    plt.subplot(121)
    plt.quiver( x, y[::-1], u0, v0 ,   alpha=1,  facecolor='orange', scale_units='xy', scale=scale_quiver) 
    plt.title('unfiltered')
    plt.subplot(122)
    plt.title('filtered')
    plt.quiver( x, y[::-1], u1, v1 ,   alpha=1,   facecolor='b', scale_units='xy', scale=scale_quiver)
    plt.savefig(outfolder+'/filtered+unfilterd.png',bbox_inches='tight',  pad_inches=0)
    plt.close()
    # save overlay 
    # different color channels
    plt.figure()   
    overlay = np.zeros((contr.shape[0],contr.shape[1],3), 'uint8')
    overlay[..., 1] = contr   #1
    overlay[..., 2] = relax   #2
    plt.imshow(overlay, origin='upper')   #extent=[0, contr.shape[0], contr.shape[1], 0])  # turn Y
    plt.quiver( x, y, u1, v1,  facecolor='orange', alpha=1, scale_units='xy', scale=scale_quiver)  
    plt.axis('off')
    plt.savefig(outfolder+'/overlay-filtered.png',  bbox_inches='tight',  pad_inches=0)
    # difference image
    plt.figure()   
    plt.imshow(np.abs(contr-relax), cmap='viridis', origin='upper' )  
    plt.quiver( x, y, u1, v1,  facecolor='orange', alpha=1, scale_units='xy', scale=scale_quiver)  
    plt.axis('off')
    plt.savefig(outfolder+'/difference-img-filtered.png',  bbox_inches='tight',  pad_inches=0)
    plt.close()

    # save deformation and image data
    deformation_unfiltered = np.sqrt(u0**2+v0**2)
    deformation_filtered = np.sqrt(u1**2+v1**2)
    maxdefo_um_unfiltered = np.nanmax(deformation_unfiltered)*scale
    maxdefo_um_filtered = np.nanmax(deformation_filtered)*scale
    np.savetxt(outfolder+'/maxdefo_um_unfiltered.txt',[maxdefo_um_unfiltered])
    np.savetxt(outfolder+'/maxdefo_um_filtered.txt',[maxdefo_um_filtered])
    print('Maximal deformation unfiltered: ',maxdefo_um_unfiltered,'Maximal deformation filtered: ',maxdefo_um_filtered)
    plt.imsave(outfolder+'/raw_contr.png', contr, cmap='gray')
    plt.imsave(outfolder+'/raw_relax.png', relax, cmap='gray')
    plt.close()        
 
    return



def click_length(img_path, out= None, scale = None):
    """
    img: path to image to determine length
    scale: scale of input image in µm per pixel - if None length is stored in px-values
    out: path to store length within a txt file ( µm or px if  scale=None ) 
            if None no text file is saved
    return: length (if scale given in um otherwiese in px)
    """
    # if outout folder is given creates that folder if it doesn't exist
    if out:
        if not os.path.exists(out):
            os.makedirs(out)
    
    # read in image
    a = plt.imread(img_path)
    # determine distance in raw image
    fig = plt.figure()
    plt.imshow(a, cmap='viridis')  
    plt.text(0.5, 1.05,'Click two-point line with left click, finish with right click',  fontsize=12,
        horizontalalignment='center',
        verticalalignment='center',c='darkred', transform= plt.gca().transAxes)
    # click length
    roi1 = roipoly.RoiPoly(color='r', fig=fig)
    # calculate and save length + plot - if no scale given pixel values instead of um
    if scale:
        length = np.round(np.sqrt((roi1.x[0]-roi1.x[1])**2+(roi1.y[0]-roi1.y[1])**2) * scale)  # convert to um with pxsize
        length_px = np.round(np.sqrt((roi1.x[0]-roi1.x[1])**2+(roi1.y[0]-roi1.y[1])**2) * 1 )  
        center= (0.5*(roi1.x[0]+roi1.x[1]), 0.5*(roi1.y[0]+roi1.y[1]))
        print (center)
        print(length)
        if out:
            np.savetxt(out+'\length_um.txt',[length])
            np.savetxt(out+'\length_px.txt',[length_px])
            np.savetxt(out+'\center_xy_px.txt',[center])
            np.savetxt(out+'\scale.txt',[scale])
            fig.savefig(out+'\length_clicked.png',dpi=150, bbox_inches='tight', pad_inches=0)
    else:
        length = np.round(np.sqrt((roi1.x[0]-roi1.x[1])**2+(roi1.y[0]-roi1.y[1])**2) * 1 )  
        center= (0.5*(roi1.x[0]+roi1.x[1]), 0.5*(roi1.y[0]+roi1.y[1]))
        print(length)
        if out:
            np.savetxt(out+'\length_px.txt',[length])
            fig.savefig(out+'\length_clicked.png',dpi=150, bbox_inches='tight', pad_inches=0)
            np.savetxt(out+'\center_xy_px.txt',[center])
    plt.close()  

    return length




def leica_projections_experiment(experiment_list, output_folder, zrange=None, mode="max", percentile_mode=99,gauss_filter=None):
    """
    specially tailored to leica software data structure - creates individual maximum projections for several  experiments 
    (stacks of several positions for several time steps stored in several sub-series)
    Fixed for 2 channels (ch00 + ch01) here
    
    
    experiment_list: list of all mark_and_find experiments to be evaluated (can contain several subseries)
                     e.g.    [r"..\data\Sample_1",  r"..\data\Sample_2"] where Sample_1 contains all tiff images
    
    output_folder: to store output projections - subfolder are created automatically
    gauss_filter: gauss filter projecteions if not none 
    zrange: range of images around stack center for the maximum projection
    """
   # create outputfolder if not exists
    if not os.path.exists(output_folder):
     os.makedirs(output_folder)                                         
    # loop through all experiments
    for exp in (experiment_list):
            print ("Current:"+str(exp))
            data_ch00 = glob.glob(exp+"\*ch00.tif")  # read in imagedata also containing all subseries 
            data_ch01 = glob.glob(exp+"\*ch01.tif")  # read in imagedata also containing all subseries 
            # find the number of subseries for this expereiment
            number_series = np.unique([int(os.path.basename(data_ch00[i])[6:9]) for i in range(len(data_ch00)) ])  
            print("Subseries found: "+str(number_series))
            # loop through subseries (counting from 1)
            for n in tqdm(number_series): 
                print ("Suberies number: "+str(n))
                # mask all data of subseries
                mask_subseries =  "Series{}".format(str(n).zfill(3))  
                subseries_ch00 =  [x for x in data_ch00 if mask_subseries in x] 
                subseries_ch01 =  [x for x in data_ch01 if mask_subseries in x]
                # find maximum timestep of the eries
                t_list = np.unique([ ((os.path.basename(subseries_ch00[i]).split('_'))[1][1:]) for i in range(len(subseries_ch00)) ])
                max_t = np.max([int(tt) for tt in t_list])  
                print ("Timesteps found: "+ str(max_t))
                # create maximumprojection for each time step
                for t in tqdm(t_list):
                    # stacks ch00 and ch01 ot make the current maxproj
                    stack_data_ch00 = [x for x in subseries_ch00 if "t{}".format(t) in x]
                    stack_data_ch01 = [x for x in subseries_ch01 if "t{}".format(t) in x]
                    # nameing for current stacks
                    current_experiment = os.path.basename(exp)   
                    name =  os.path.split(stack_data_ch00[0])[1][:-14] 
                    #current subfolder
                    exp_folder = os.path.join(output_folder,current_experiment)
                    # make the projection
                    projections(stack_data_ch00, os.path.join(exp_folder,"ch00_zrange_{}".format(zrange)), name, zrange=zrange, mode=mode, percentile_mode =percentile_mode, gauss_filter = gauss_filter)
                    projections(stack_data_ch01, os.path.join(exp_folder,"ch01_zrange_{}".format(zrange)), name, zrange=zrange, mode=mode,  percentile_mode =percentile_mode, gauss_filter = gauss_filter)              
    return

def leica_projections_markfind(experiment_list, output_folder, zrange=None, mode="max", percentile_mode=99):
    """
    specially tailored to leica software data structure - creates individual maximum projections for several mark and find experiments 
    (stacks of several positions for several time steps stored in several sub-series)
    Fixed for 2 channels (ch00 + ch01) here
    
    
    experiment_list: list of all mark_and_find experiments to be evaluated (can contain several subseries)
                     e.g.    [r"..\data\Sample_1",  r"..\data\Sample_2"] where Sample_1 contains Mark_and_Find_001, Mark_and_Find_002 etc..
    
    output_folder: to store output projections - subfolder are created automatically
    
    zrange: range of images around stack center for the maximum projection
    """
                     
    # create outputfolder if not exists
    if not os.path.exists(output_folder):
     os.makedirs(output_folder)                                         
    # loop through all experiments
    for exp in (experiment_list):
            print ("Current:"+str(exp))
            data_ch00 = glob.glob(exp+"\*\*ch00.tif")  # read in imagedata also containing all subseries 
            data_ch01 = glob.glob(exp+"\*\*ch01.tif")  # read in imagedata also containing all subseries 
            # find the number of subseries for this expereiment
            number_series = np.max([int(os.path.split(os.path.split(data_ch00[i])[0])[1][-3:]) for i in range(len(data_ch00))])  # second las entry mark and find number
            print("Subseries found: "+str(number_series))
            # loop through subseries (counting from 1)
            for n in tqdm(range(1,number_series+1)): 
                print ("Suberies number: "+str(n))
                # mask all data of subseries
                mask_subseries =  "Mark_and_Find_{}".format(str(n).zfill(3))  
                subseries_ch00 =  [x for x in data_ch00 if mask_subseries in x] 
                subseries_ch01 =  [x for x in data_ch01 if mask_subseries in x]
                # find maximum position of the series
                position_list = set([int(str.split((os.path.split(os.path.split(subseries_ch00[i])[1])[1]), "_")[0][3:]) for i in range(len(subseries_ch00))])
                print ("Positions found: "+ str(position_list))
                # loop through positions
                for z in (position_list):
                    # reduce subseries to individual positions postion
                    mask_position =  "Pos{}".format(str(z).zfill(3))  #starts from 1
                    subseries_pos_ch00 =  [x for x in subseries_ch00 if mask_position in x] 
                    subseries_pos_ch01 =  [x for x in subseries_ch01 if mask_position in x]
                    # find maximum timestep of the eries
                    t_list = [int(str.split((os.path.split(os.path.split(subseries_pos_ch00[i])[1])[1]), "_")[2][1:]) for i in range(len(subseries_pos_ch00))]
                    max_t = np.max(t_list)  
                    print ("Timesteps found: "+ str(max_t))
                    # create maximumprojection for each time step
                    for t in (range(max_t+1)):
                        # stacks ch00 and ch01 ot make the current maxproj
                        if max_t<10:
                            stack_data_ch00 = [x for x in subseries_pos_ch00 if "t{}".format(str(t).zfill(0)) in x]
                            stack_data_ch01 = [x for x in subseries_pos_ch01 if "t{}".format(str(t).zfill(0)) in x]
                        if max_t>=10:
                            stack_data_ch00 = [x for x in subseries_pos_ch00 if "t{}".format(str(t).zfill(2)) in x]
                            stack_data_ch01 = [x for x in subseries_pos_ch01 if "t{}".format(str(t).zfill(2)) in x]
                            
                        # nameing for current stacks
                        current_experiment = os.path.basename(exp)
                        upper = os.path.split(os.path.split(stack_data_ch00[0])[0])[1]
                        name =  os.path.split(stack_data_ch00[0])[1][:-14] 
                        
                        #current subfolder
                        exp_folder = os.path.join(output_folder,current_experiment)
                        Pos_folder = os.path.join(exp_folder,"Pos{}".format(str(z).zfill(3)))
       
                        # make the projection
                        projections(stack_data_ch00, os.path.join(Pos_folder,"ch00_zrange_{}".format(zrange)), upper+"_"+name, zrange=zrange, mode=mode, percentile_mode =percentile_mode)
                        projections(stack_data_ch01, os.path.join(Pos_folder,"ch01_zrange_{}".format(zrange)), upper+"_"+name, zrange=zrange, mode=mode,  percentile_mode =percentile_mode)
    return


def video_maker(files,output, fps=15, no_dir=False):
    """
    creates a video for given files (path in the "glob"-structure  for example r"C/User/cell1/*.tiff")
    
    files: list of imagepaths for the video (see above)
    output: directiory to store the .mp4 video
    fps: framerate for the video
    no_dir: if True directories are not created automatically
    """
    # create outputfolder if not exists
    if not no_dir:
        if not os.path.exists(output):
         os.makedirs(output)      
    # read in files with glob 
    file_list = glob.glob(files)
    # create video
    writer = imageio.get_writer(output+'.mp4', fps=fps)
    for im in file_list:
        writer.append_data(imageio.imread(im))
    writer.close()    
    return
    

def leica_video_maker(leicamaxpath, output=None, fps = 12):

    """
    creates videos of all cells from the output folder of leica_projections.  
    just uses video_maker  for the special data structure
       

    leicamaxpath: output folder of leica_projections containing all subfodlers
    output: fodler containing all videos - by default created in input folder
    fps: framerate for videos        
  
    """
    if not output:
        output = os.path.join(leicamaxpath, "videos")

    # all maxproj series from subfolders    
    vid_list =  glob.glob(leicamaxpath+"\*\Pos*\*")    
    # loop through series
    for v in tqdm(vid_list):
        # create subfolder for experiment
        out = os.path.join(output,v.split(os.path.sep)[-3])  
        # create if not exiting
        if not os.path.exists(out):
            os.makedirs(out)  
        # and add right filename with full information
        out_vid = os.path.join(out, v.split(os.path.sep)[-2]+"_"+v.split(os.path.sep)[-1])
        # make video from tiff images
        video_maker(v+"\*.tif" ,out_vid, fps=fps, no_dir=True)

      
    return  
        
        
        
        
        