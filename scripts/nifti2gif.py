"""
Noée Ducros-Chabot 
23 Juillet 2021
Script turn two nifti files into intertwined gif
"""
import os, sys
import shutil
import argparse
import warnings
import copy
import time

def pngfolder(niftifile):
    head, tail = os.path.split(niftifile)
    imName = tail.split('.')[0]
    pngfolder = os.path.abspath(os.path.join(niftifile,'..', imName+'_png'))
    return pngfolder

def createLegendPNG(png_path):
    images = sorted([os.path.abspath(os.path.join(png_path, p)) for p in os.listdir(png_path)])

    split_path = png_path.split('/')
    method = split_path[-2]
    folder = split_path[-1]

    out_folder = png_path + '_legend'
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    
    for image in images:
        ofile = copy.copy(image)
        ofile = ofile.replace(folder, folder +'_legend')
        os.system(f'convert {image} -pointsize 18 -annotate +20+20 {method} {ofile}')

    return method, out_folder

def gifCreation(path1, path2, ofile):
    set1 = sorted([os.path.abspath(os.path.join(path1, p)) for p in os.listdir(path1)])
    set2 = sorted([os.path.abspath(os.path.join(path2, p)) for p in os.listdir(path2)])
    if len(set1)>200 and len(set2)>200:
        set1 = set1[40:]
        set2 = set2[40:]

    paths = [item for pair in zip(set1, set2) for item in pair]
    
    os.system(f'convert -delay 50 -loop 0 {" ".join(x for x in paths)} {ofile}')

def main(im1_path, im2_path, outpath, legend = False, delPNG = False):
    # im1_path : path to nifti file of first image
    # im2_path : path to nifti file of the second image
    # out_path : output nifti filename
    
    pngpath1 = pngfolder(im1_path)
    pngpath2 = pngfolder(im2_path)
    
    if os.path.exists(pngpath1):
        warnings.warn(f"{pngpath1} already exists. Deleting folder and its content.")
        time.sleep(2)
        shutil.rmtree(pngpath1)
    else: 
        os.makedirs(pngpath1)
        
    if os.path.exists(pngpath2):
        warnings.warn(f'{pngpath2} already exists. Deleting folder and its content.')
        time.sleep(2)
        shutil.rmtree(pngpath2)
    else: 
        os.makedirs(pngpath2)  

    #Creation of png images 
    fpng = os.getcwd()
    fpng = os.path.join(fpng, 'nii2png.py')
    os.system(f'python {fpng} -i {im1_path} -o {pngpath1}')
    os.system(f'python {fpng} -i {im2_path} -o {pngpath2}')

    
    if legend: 
        #Create png with legend
        print("Adding legend to png files")
        method1, legPath1 = createLegendPNG(pngpath1)
        method2, legPath2 = createLegendPNG(pngpath2)

        #Create gif 
        print("Creating GIF")
        gifCreation(legPath1, legPath2, outpath)

        if delPNG: 
            print("Deleting PNG folders")
            pngfolders = [legPath1, legPath2, pngpath1, pngpath2]
            for folder in pngfolders: 
                shutil.rmtree(folder)
    else: 
        #Create gif 
        print("Creating GIF")
        gifCreation(pngpath1, pngpath2, outpath)

        if delPNG:
            print("Deleting PNG folders")
            pngfolders = [pngpath1, pngpath2]
            for folder in pngfolders:
                shutil.rmtree(folder)

if __name__ == '__main__':
    # Example of command line input : /home/magic-chusj-2/Desktop/Registration/dhcp-volumetric-atlas-groupwise/mean/ga_40/average_t1_masked.nii /home/magic-chusj-2/Desktop/Registration/results/dipy/HIE_baby_016/registered_mag.nii /home/magic-chusj-2/Desktop/Registration/results/dipy/HIE_baby_016/reg.gif
    parser = argparse.ArgumentParser(description = 'Script creatin alternating gif from two nifti.')
    parser.add_argument('img1', type = str, help='first niftifile')
    parser.add_argument('img2', type = str, help='second niftifile')
    parser.add_argument('gif_file', type = str, help='filename of output gif')
    parser.add_argument('-l', '--legend', type = bool, default = False, help='option to add legend to nifti images')
    parser.add_argument('-d', '--delete', type = bool, default = False, help='option to delete intermediate png images')
    args = parser.parse_args()

    main(args.img1, args.img2, args.gif_file, args.legend, args.delete)
