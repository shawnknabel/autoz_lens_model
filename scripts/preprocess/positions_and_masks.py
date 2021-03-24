#!/usr/bin/env python
# coding: utf-8

# # 3/10/21 - This notebook creates cutouts and background noise images from KiDS coadds and weight images, creates a noise map, converts all images to eps, creates a psf, and saves the files to either fits or csv.

# In[1]:


### libraries
#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
#from autoconf import conf
import autolens as al
import autolens.plot as aplt
#import autofit as af
import pandas as pd
import numpy as np
import astropy.io.fits as fits
#from astropy.visualization import astropy_mpl_style
#plt.style.use(astropy_mpl_style)
#from astropy.stats import sigma_clip as clip
#from astropy.coordinates import SkyCoord
#import astropy.units as u
#from astropy.nddata import Cutout2D
#from astropy.wcs import WCS
import time
import os

from pyprojroot import here
print('We are here: ', str(here()))


# In[2]:


# set datetime variable
date = time.strftime("%d%m%y")

# paths
autoz_path = f'{str(here())}/'
file_path = f'{autoz_path}files/'
csv_path = f'{file_path}csv/'
fits_path = f'{file_path}fits/'
png_path = f'{autoz_path}visuals/png/'
pdf_path = f'{autoz_path}visuals/pdf/'




# In[38]:


# write functions

def effective_radius(gama_id, band):
    # effective radius and error from DR3 sersic catalog
    hdul = fits.open(f'{fits_path}SersicCatSDSS.fits')
    hdul.verify('fix')
    data = hdul[1].data
    candidate = data[data.CATAID == gama_id]
    if band == 'r':
        re = candidate.GALRE_r[0]
        re_err = candidate.GALREERR_r[0]
    elif band == 'g':
        re = candidate.GALRE_g[0]
        re_err = candidate.GALREERR_g[0]
    else:
        print('Get your image bands straight!')
    return(re, re_err)

# function to identify positions to create a json (place after the fits file has been generated)
def position_suspicion (image, gama_id, links_id, band, pixel_scale):
    visuals_2d = aplt.Visuals2D()
    correct = 'n'
    while correct == 'n':
        print('Plotting those suspicious positions...')
        array_plotter = aplt.Array2DPlotter(array=image, visuals_2d=visuals_2d)
        array_plotter.figure()
        print('Are you satisfied? y/n/skip')
        correct = str(input())
        if correct == 'y':
            print('Coooool. Outputting to csv.')
            positions_array = np.array(positions_grid)
            print(positions_array)
            np.savetxt(f'{fits_path}G{gama_id}_{links_id}/{links_id}_{band}_positions_grid.csv', positions_array, delimiter = ',')
            break
        elif correct == 'n':
            print('How many positions?')
            n = int(input())
            positions_grid = []
            for i in range(n):
                print(f'Position {i+1}... Where shall we probe? (y, x)')
                y, x = [float(input()), float(input())]
                positions_grid.append((y, x))
            print(f'Showing positions at {positions_grid}.')
            positions = al.Grid2DIrregular(grid=positions_grid)
            visuals_2d = aplt.Visuals2D(positions)
        else:
            print('Skipping')
            break
    
# function to draw masks
def masquerade (image, gama_id, links_id, band, pixel_scale, re, re_err):
    correct = 'n'
    print('Plotting those suspicious masks...')
    while correct == 'n':
        # overall mask
        print('Overall mask? y/n')
        answer = input()
        if answer == 'y':
            satisfied = 'n'
            radius=3.0
            centre=(0.0, 0.0)
            while satisfied == 'n':
                mask = al.Mask2D.circular(
                        shape_native=image.shape_native,
                        pixel_scales=image.pixel_scales,
                        sub_size=2,
                        radius=radius,
                        centre=centre,
                            )
                print('Plotting mask.')
                array_plotter = aplt.Array2DPlotter(array=image, visuals_2d=aplt.Visuals2D(mask=mask))
                array_plotter.figure()
                print('Are you satisfied? y/n/skip')
                satisfied = str(input())                
                if satisfied == 'y':
                    print('Coooool. Outputting to fits.')
                    mask.output_to_fits(file_path=f'{fits_path}G{gama_id}_{links_id}/{links_id}_{band}_mask.fits', overwrite=True)
                    break
                elif satisfied == 'n':
                    print('Coooool. We can change it.')
                    print('New radius? y/n')
                    answer = str(input())
                    if answer == 'y':
                          print('New radius value?')
                          radius = float(input())
                    else:
                          print('Radius unchanged.')
                    print('New center?')
                    answer = str(input())
                    if answer == 'y':
                          print('New center values? (y, x)')
                          centre = (float(input()), float(input()))
                    else:
                          print('Center unchanged.')
                else:
                    print('Skipping')
                    break                    
        else:
            print('Skipping overall mask')
        
        # lens mask
        print('Lens mask? y/n')
        answer = input()
        if answer == 'y':
            satisfied = 'n'
            radius=re+re_err
            centre=(0.0, 0.0)
            while satisfied == 'n':
                lens_mask = al.Mask2D.circular(
                    shape_native=image.shape_native,
                    pixel_scales=image.pixel_scales,
                    sub_size=2,
                    radius=radius,
                    centre=centre)
                print('Plotting lens mask')
                array_plotter = aplt.Array2DPlotter(array=image, visuals_2d=aplt.Visuals2D(mask=lens_mask))
                array_plotter.figure()
                print('Are you satisfied? y/n/skip')
                satisfied = str(input())                
                if satisfied == 'y':
                    print('Coooool. Outputting to fits.')
                    lens_mask.output_to_fits(file_path=f'{fits_path}G{gama_id}_{links_id}/{links_id}_{band}_lens_mask.fits', overwrite=True)
                    break
                elif satisfied == 'n':
                    print('Coooool. We can change it.')
                    print('New radius? y/n')
                    answer = str(input())
                    if answer == 'y':
                          print('New radius value?')
                          radius = float(input())
                    else:
                          print('Radius unchanged.')
                    print('New center?')
                    answer = str(input())
                    if answer == 'y':
                          print('New center values? (y, x)')
                          centre = (float(input()), float(input()))
                    else:
                          print('Center unchanged.')
                else:
                    print('Skipping')
                    break
        else:
            print('Skipping lens mask')
            
        # source mask
        print('Source mask? y/n')
        answer = input()
        if answer == 'y':
            satisfied = 'n'
            inner_radius=re+re_err
            outer_radius=3.5
            centre=(0.0, 0.0)
            while satisfied == 'n':
                source_mask = al.Mask2D.circular_annular(
                    shape_native=image.shape_native,
                    pixel_scales=image.pixel_scales,
                    sub_size=2,
                    inner_radius=inner_radius,
                    outer_radius=outer_radius,
                    centre=centre)  
                print('Plotting source mask')
                array_plotter = aplt.Array2DPlotter(array=image, visuals_2d=aplt.Visuals2D(mask=source_mask))
                array_plotter.figure()
                print('Are you satisfied? y/n/skip')
                satisfied = str(input())                
                if satisfied == 'y':
                    print('Coooool. Outputting to fits.')
                    source_mask.output_to_fits(file_path=f'{fits_path}G{gama_id}_{links_id}/{links_id}_{band}_source_mask.fits', overwrite=True)
                    break
                elif satisfied == 'n':
                    print('Coooool. We can change it.')
                    print("Perhaps you'd like an elliptical annular mask? y/n")
                    answer = str(input())
                    if answer == 'y':
                          satisfied = 'y'
                          print('Do it later.')
                          break
                    else:
                        print('New radii? y/n')
                        answer = str(input())
                        if answer == 'y':
                              print('New radiu values? (inner, outer)')
                              inner_radius, outer_radius = [float(input()), float(input())]
                        else:
                            print('Radius unchanged.')
                        print('New center?')
                        answer = str(input())
                        if answer == 'y':
                            print('New center values? (y, x)')
                            centre = (float(input()), float(input()))
                        else:
                            print('Center unchanged.')
                else:
                    print('Skipping')
                    break
            print('Is everything correct? y/n')
            correct = str(input())
            if correct == 'y':
                break
            else:
                print('Then why did you say you were satisfied with the masks??? Do it later.')
                break
        else:
            print('Skipping source mask')
            break        
    print('Masks off!')    

## Ask for input (for script)
print('Thank you for using ScriptBotPlus Lensing Ed. copyright Idaho Potatoes, Ltd.')
print('Please follow the prompts without failure, or you will be punished severely.')
print('\n Enjoy your visit!')

print('GAMA ID?')
gama_id = int(input())
print('LiNKS ID?')
links_id = int(input())
print('Pixel scale?')
pixel_scale = float(input())
print('Enter any phun phrase to get us started...')
input()

        
print('Great, you have ordered the dirty turnips with basil pesto. That will be $45.87. Pull up to the next window.')
print()
print(':)')

# load/show images
r_image = al.Array2D.from_fits(file_path= f"{fits_path}G{gama_id}_{links_id}/{links_id}_r_image.fits", pixel_scales=pixel_scale)
g_image = al.Array2D.from_fits(file_path= f"{fits_path}G{gama_id}_{links_id}/{links_id}_g_image.fits", pixel_scales=pixel_scale)
print('Plotting r image')
r_array_plotter = aplt.Array2DPlotter(array=r_image)
r_array_plotter.figure()
print('Plotting g image')
g_array_plotter = aplt.Array2DPlotter(array=g_image)
g_array_plotter.figure()

# positions
print('Do you want to indicate positions? y/n')
answer = input()
if answer == 'y':
    print('r band? y/n')
    answer = input()
    if answer == 'y':
        position_suspicion(r_image, gama_id, links_id, 'r', pixel_scale)
    else:
        print('Skipping r band')
    print('g band? y/n')
    answer = input()
    if answer == 'y':
        position_suspicion(g_image, gama_id, links_id, 'g', pixel_scale)
    else:
        print('Skipping g band')
elif answer == 'n':
    print('No positions today.')
else:
    print('This was a yes or no question...')
    


# masks
print('Do you want to create masks? y/n')
answer = input()
if answer == 'y':
    print('r band?')
    answer = input()
    if answer == 'y':
        band = 'r'
        re, re_err = effective_radius(gama_id, band)
        masquerade(r_image, gama_id, links_id, band, pixel_scale, re, re_err)
    else:
        print('Skipping r band')
    print('g band?')    
    answer = input()
    if answer == 'y':
        band = 'g'
        re, re_err = effective_radius(gama_id, band)
        masquerade(g_image, gama_id, links_id, band, pixel_scale, re, re_err)
    else:
        print('Skipping g band')
elif answer == 'n':
    print('No masks today.')
else:
    print('This was a yes or no question...')
    
print('\n\n\n Work complete!')
    
    
    





