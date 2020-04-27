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

test = r'\\131.188.117.96\biophysDS\dboehringer\Platte_4\Measurements_NK_TFM\2019-10-15_Konfokal-Primarynk\nk21-1_zoom6_eval\88'
out = r"\\131.188.117.96\biophysDS\dboehringer\Platte_4\Eval-data\Example-Data-test\nk-testing\2019-10-15-nk21-series88\out-arnold-test"
a = convert_to_png(test, out)

def regularization_stacks():
    """
    starts a full 3D saeno regularization for the given relaxed and alive input stacks.
    """
    
    
    return


def display_displacements():
    """
    show the total displacement field or the masked range of displacements by size
    """
    return

def filter_displacements():
    """
    filters out deformations by size and stores the resulting files in output
    """ 
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


#   stack a  \\131.188.117.96\biophysDS\dboehringer\Platte_4\Measurements_NK_TFM\2019-10-15_Konfokal-Primarynk\nk21-1_zoom6_eval\88\Series088_t14_z*_ch00.tif
#   stack r  \\131.188.117.96\biophysDS\dboehringer\Platte_4\Measurements_NK_TFM\2019-10-15_Konfokal-Primarynk\nk21-1_zoom6_eval\89_rel\Series089_t12_z*_ch00.tif

stacka= r'\\131.188.117.96\biophysDS\dboehringer\Platte_4\Measurements_NK_TFM\2019-10-15_Konfokal-Primarynk\nk21-1_zoom6_eval\88\Series088_t14_z*_ch00.tif'
