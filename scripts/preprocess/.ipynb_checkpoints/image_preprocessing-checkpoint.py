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
from astropy.stats import sigma_clip as clip
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
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


# In[3]:


# load links data
links = pd.read_csv(f'{csv_path}latest/links_sample_latest.csv')
#ra = links[links.GAMA_ID == 250289].RA.iloc[0]
#dec = links[links.GAMA_ID == 250289].DEC.iloc[0]
#print(ra, dec)


# In[38]:


# write functions
# open fits file
def open_sesame (gama_id, links_id, band, weight=False):
    if weight == True:
        band=f'{band}_weight'
    print(f'Opening {gama_id}_{links_id} {band} image')
    hdul = fits.open(f'{fits_path}G{gama_id}_{links_id}/{links_id}_{band}.fits')
    header = hdul[0].header
    image = hdul[0].data
    hdul.close()
    return(image, header)

# create the cutout
def cut_it_out (gama_id, image, header, pixel_scale):
    print(f'Making cutout of {gama_id}')
    # attached real-world coordinates to pixel locations
    wcs = WCS(header) # let wcs pull info from header
    #print(f'wcs: {wcs}')
    
    # get ra and dec from links data
    candidate = links[links.GAMA_ID==gama_id]
    ra = candidate.RA.iloc[0]
    dec = candidate.DEC.iloc[0]
    
    answer = 'n'
    
    while answer == 'n':
            
        coord=SkyCoord(ra=ra, dec=dec, unit='deg', frame='icrs') # international celestial reference frame
        position = wcs.world_to_pixel(coord)
        print(f'Coordinates: {coord}')
        #print(wcs.wcs.crval)
        print(f' Position: {position}')
        size = u.Quantity(101, u.pixel)

        cutout = Cutout2D(data=image, position=position, size=size, wcs=wcs, mode='trim')
        cutout_image = cutout.data
    
        # plot cutout
        print('Cutout')
        plt.figure()
        plt.imshow(cutout_image, origin='lower', cmap='gray')  
        plt.scatter(50, 50, color='r')
        plt.show()
        print(f'Cutout shape: {cutout_image.shape}')
        print(f'Is this centered correctly? (Remember to do the same adjustments to the image and weight map!) y/n')
        answer = str(input())
        if answer == 'y':
            break
            print('Lovely, let us continue.')
        elif answer == 'n':
            print(f'Give adjustments in pixels. Pixel scale is {pixel_scale}. Up, Down, Left, Right')
            up, down, left, right = [float(input()), float(input()), float(input()), float(input())]
            ra = ra+(right-left)*pixel_scale/3600
            dec = dec+(down-up)*pixel_scale/3600
        else:
            print('Please answer y or n')
            answer = str(input())
    return(cutout_image)
        
def count_chocula (image, header, band, noise=False):
    gain = header['GAIN']
    print(f'Gain: {gain}')
    if band == 'r':
        exp_time = 1800
    elif band == 'g':
        exp_time = 900
    elif band == 'i':
        exp_time = 1200
    else:
        print('Error - filter not g, r, or i')
    print(f'Exposure time is {exp_time}')
    if noise == True:
        image = 1/np.sqrt(image)
        print('Calculating rms noise...')
    #print(np.mean(rms_noise))
    image_counts = image*gain#*exp_time
    if noise == True:
        image_counts = image_counts**2
    image_eps = image_counts/exp_time
    print(f'Mean/Min/Max image counts: {np.mean(image_counts), np.min(image_counts), np.max(image_counts)}')
    print(f'Mean/Min/Max image eps: {np.mean(image_eps), np.max(image_eps), np.max(image_eps)}')
    return(image_counts, exp_time)

def plot_image (image, units):
    print(f'Image ', {units})
    plt.figure()
    plt.imshow(image, cmap='gray') # show image in grayscale
    plt.colorbar(label="pixel value", orientation="vertical")
    plt.show()

def reconstruct_image (image, weight):
    added_image = image + weight
    if np.min(added_image) < 0:    
        added = added_image - np.min(added_image)
    print(f'Reconstructed image min pixel value: {added_image} (should be >= 0)')
    return(added_image)

def get_noisey (reconstructed_image):
    print('Gettin noisey!')
    noise_map = np.sqrt(reconstructed_image)
    return(noise_map)

def divide_the_time (image_counts, exp_time):
    image_eps = image_counts/exp_time
    return(image_eps)
    
# resize
def resize_image(image, new_size):
    print(f'Resizing image to {new_size}.')
    size=image.shape[0] # 101
    center=int(image.shape[0]/2) # 50
    new_center=(new_size-1)/2
    lower = int(center-new_center) 
    upper = int(center+new_center)
    resized_image=image[lower:upper+1,lower:upper+1]
    print(f'Middle pixel at index {center}. New image created from indices {lower} to {upper} in axes 0 and 1.')
    print(f'New shape: {resized_image.shape}')
    return(resized_image)

# generate psf
def point_to_the_spread(image, header, pixel_scale, new_size):

    # define psf values
    avg_psf = header['PSF_RAD'] # arcsec
    avg_psf_pxl = avg_psf/pixel_scale # pixels
    sigma_psf = avg_psf_pxl/2
    size = int(np.around(image.shape[0]/2)) # gives a grid of 101 (50 on either side of the center)
    
    # set psf for 101, 101 image
    y, x = np.mgrid[-size:size+1, -size:size+1]
    psf = np.exp(-(x**2 + y**2)/(2.0*sigma_psf**2))
    psf /= np.sum(psf)
    print(f'A psf of {psf} with size {size} has been generated')
    
    # resize psf
    psf_resized = resize_image(psf, new_size) # cut to 21x21
    print(f'New psf shape: {psf_resized.shape}')
    
    # good vibes
    print('This has been fun, right? Very fun! :)')
    
    return(psf_resized)

def if_i_fits_i_sits (image, gama_id, links_id, band, noise=False, psf=False):
    if noise == True:
        band=f'{band}_noise_map'
    if psf == True:
        band=f'{band}_psf'
    print(f'Saving {gama_id}_{links_id} {band} image')
    hdu = fits.PrimaryHDU(image)
    hdu.writeto(f'{fits_path}G{gama_id}_{links_id}/{links_id}_{band}_image.fits', overwrite=True)
    print(f'Image sent to {fits_path}G{gama_id}_{links_id}/{links_id}_{band}_image.fits')
                
                #stopped here to play soccer!
                
def if_i_csv_i_cannosee (image, gama_id, links_id, band, noise=False, psf=False):
    if noise == True:
        band=f'{band}_noise_map'
    if psf == True:
        band=f'{band}_psf'
    print(f'Saving {gama_id}_{links_id} {band} image')
    np.savetxt(f'{csv_path}G{gama_id}_{links_id}/{links_id}_{band}_image.csv', image, delimiter=',')
    print(f'Image sent to {csv_path}G{gama_id}_{links_id}/{links_id}_{band}_image.csv')



# In[39]:


def one_ring_to_rule_them_all (gama_id, links_id, band, pixel_scale, psf_kernel_size, filetype):
    print('One ring to rule them all... \n One ring to find them... \n One ring to bring them all... \n and in the darkness bind them!')
    
    #load image
    print('\n Loading coadd image.')
    image, image_header = open_sesame(gama_id, links_id, band)
    
    #cut out image
    print('\n Producing cutout of coadd image.')
    cutout_image = cut_it_out(gama_id, image, image_header, pixel_scale)
    
    #convert to counts
    print('\n Converting cutout image to counts.')
    image_counts, exp_time = count_chocula(cutout_image, image_header, band)

    #plot counts
    plot_image(image_counts, 'counts')
    
    #load weight
    print('\n Loading weight image.')
    weight, weight_header = open_sesame(gama_id, links_id, band, weight=True)
    
    #cut out weight
    print('\n Producing cutout of weight image.')
    cutout_weight = cut_it_out(gama_id, weight, weight_header, pixel_scale)
    
    print('\n All this work makes me hungry...')
    
    #convert to counts
    print('\n Converting weight image to noise image in counts.')
    background_counts, exp_time = count_chocula(cutout_weight, image_header, band, noise=True)
    
    #plot counts
    plot_image(background_counts, 'counts')

    #reconstruct the image
    print('\n Reconstructing image with background.')
    reconstructed_image = reconstruct_image(image_counts, background_counts)
    
    print('\n I am hungry for human food.')
    
    #create noise
    print('\n Creating noise map in counts.')
    print('Gettin noisey!')
    noise_map_counts = get_noisey(reconstructed_image)
    
    # convert image to eps
    print('\n Converting image to eps.')
    image_eps = divide_the_time(image_counts, exp_time)

    # plot eps
    plot_image(image_eps, 'eps')
    
    #convert noise to eps
    print('\n Converting noise map to eps.')
    noise_map_eps = divide_the_time(noise_map_counts, exp_time)

    plot_image(noise_map_eps, 'eps')
    

    #create psf
    print('\n Creating psf.')
    psf = point_to_the_spread(image_eps, image_header, pixel_scale, psf_kernel_size)
    
    print('\n I mean foods that humans eat, not humans as food.')
    
    if filetype=='fits':
        os.makedirs(f'{fits_path}G{gama_id}_{links_id}', exist_ok=True)
        #save image to hdu
        print('\n If I fits, I sits.')
        if_i_fits_i_sits(image_eps, gama_id, links_id, band)
    
        #save noise to hdu
        if_i_fits_i_sits(noise_map_eps, gama_id, links_id, band, noise=True)
    
        #save psf to hdu
        if_i_fits_i_sits(psf, gama_id, links_id, band, psf=True)
    
        
    elif filetype=='csv':
        os.makedirs(f'{csv_path}G{gama_id}_{links_id}', exist_ok=True)
        #save image to hdu
        print('\n If I csv, I cannosee.')
        if_i_csv_i_cannosee(image_eps, gama_id, links_id, band)
    
        #save noise to hdu
        if_i_csv_i_cannosee(noise_map_eps, gama_id, links_id, band, noise=True)
    
        #save psf to hdu
        if_i_csv_i_cannosee(psf, gama_id, links_id, band, psf=True)
    else:
        print('\n Get your filetypes straight!')
    
    
    print('\n\n\n Work complete!')
    
    
    


# In[40]:


## Ask for input (for script)
print('Thank you for using ScriptBotPlus Lensing Ed. copyright Idaho Potatoes, Ltd.')
print('Please follow the prompts without failure, or you will be punished severely.')
print('\n Enjoy your visit!')

print('GAMA ID?')
gama_id = int(input())
print('LiNKS ID?')
links_id = int(input())
print('Band?')
band = str(input())
print('Pixel scale?')
pixel_scale = float(input())
print('PSF size?')
psf_kernel_size = float(input())
print('Filetype?')
filetype = input()
print('Enter any phun phrase to get us started...')
input()

        
print('Great, you have ordered the dirty turnips with basil pesto. That will be $45.87. Pull up to the next window.')
print()
print(':)')

one_ring_to_rule_them_all(gama_id, links_id, band, pixel_scale, psf_kernel_size, filetype)


# In[ ]:




