#!/usr/bin/env python
# coding: utf-8

# # 3/31/21

# ## Experiment 9 (White): phase1_nlive = 200, phase1_tolerance = 0.5, phase2_nlive = 300, phase2_tolerance = 0.5, phase3_nlive = 500, phase3_tolerance = 0.25, 4 Positions, Positions threshold = 1.0, use effective radius from GAMA DR3 Sersic Photometry catalog with sigma the observed error and upper limit (effective radius + error 1.3435+0.1146=1.4581), smaller mask (effective radius) centered at (0.075, -0.075), fix lens light bulge profile in second phase. Source R_E (0->10*lens.inensity)
# 
# ## Don't take the instance for the ellipticity and effective radius. Set same max effective radius


# In[5]:


### libraries
#get_ipython().run_line_magic('matplotlib', 'inline')
# set workspace path
from pyprojroot import here
workspace_path = str(here())#'/data/sknabel/autolens_workspace'
#get_ipython().run_line_magic('cd', '$workspace_path')
print(f"Working Directory has been set to `{workspace_path}`")

import matplotlib.pyplot as plt
from autoconf import conf
import autolens as al
import autolens.plot as aplt
import autofit as af
import pandas as pd
import numpy as np
from astropy.io import fits
#from astropy.visualization import astropy_mpl_style
#plt.style.use(astropy_mpl_style)
#from astropy.stats import sigma_clip as clip
from os import path
import time

# set datetime variable
datetime = time.strftime("%d%m%y")
print(f'Today is {datetime}')
# paths
autoz_path = f'{workspace_path}/'
config_path = f'{here()}/config'
conf.instance.push(new_path=config_path)
file_path = f'{autoz_path}files/'
csv_path = f'{file_path}csv/'
fits_path = f'{file_path}fits/'
png_path = f'{autoz_path}visuals/png/'
pdf_path = f'{autoz_path}visuals/pdf/'


# In[6]:


#gama_id = 250289
#links_id = 2730

# set up for input
print('GAMA ID?')
gama_id = int(input())
print('LinKS ID?')
links_id = int(input())
print('Experiment #?')
experiment_number = int(input())
print('Shall we have a go?')
input()

object_folder = f'{fits_path}G{gama_id}_{links_id}/'
output_folder = f'{autoz_path}output/G{gama_id}_{links_id}/'

# load object data table
links = pd.read_csv(f'{csv_path}/latest/links_sample_latest.csv')
lens_galaxy_data = links[links.GAMA_ID == gama_id]
zlens=lens_galaxy_data.zlens.values
zsource=lens_galaxy_data.zsource.values
einstein_radius=np.mean([lens_galaxy_data.theta_e_pm.values, lens_galaxy_data.theta_e_sis.values]) # take average of einstein radius estimates for prior
print(f'Lens and source redshifts at {zlens} and {zsource}.')
#print(f'Einstein radius prior: {einstein_radius}')

# effective radius and error from DR3 sersic catalog
hdul = fits.open(f'{fits_path}SersicCatSDSS.fits')
hdul.verify('fix')
data = hdul[1].data
#print(data.columns)
#print(data.CATAID)
candidate = data[data.CATAID == gama_id]
re_r = candidate.GALRE_r[0]
re_r_err = candidate.GALREERR_r[0]
#re_g = candidate.GALRE_g[0]
#re_g_err = candidate.GALREERR_g[0]

# load r-band imaging
imaging_r = al.Imaging.from_fits(image_path=path.join(object_folder, f'{links_id}_r_image.fits'),
                             noise_map_path=path.join(object_folder, f'{links_id}_r_noise_map_image.fits'),
                             psf_path=path.join(object_folder, f'{links_id}_r_psf_image.fits'),
                              pixel_scales=0.2)
#imaging_g = al.Imaging.from_fits(image_path=path.join(object_folder, f'{links_id}_g_image.fits'),
#                             noise_map_path=path.join(object_folder, f'{links_id}_g_noise_map_image.fits'),
#                             psf_path=path.join(object_folder, f'{links_id}_g_psf_image.fits'),
#                              pixel_scales=0.2)
imaging_white = al.Imaging.from_fits(image_path=path.join(object_folder, f'{links_id}_white_image.fits'),
                             noise_map_path=path.join(object_folder, f'{links_id}_white_noise_map_image.fits'),
                             psf_path=path.join(object_folder, f'{links_id}_white_psf_image.fits'),
                              pixel_scales=0.2)

# set up masks
mask = al.Mask2D.from_fits(f'{object_folder}{links_id}_white_mask.fits', pixel_scales=imaging_r.pixel_scales)
lens_mask = al.Mask2D.from_fits(f'{object_folder}{links_id}_r_lens_mask.fits', pixel_scales=imaging_r.pixel_scales)
source_mask = al.Mask2D.elliptical_annular(
    shape_native=imaging_white.shape_native, pixel_scales=imaging_white.pixel_scales, sub_size=2,# outer_radius=2.5, inner_radius=0.75,
    #centre=phase1_result.instance.galaxies.lens.bulge.centre
    centre=(0.2, -0.3),
    inner_major_axis_radius=1.2,
    inner_axis_ratio=0.6,
    inner_phi=60.0,
    outer_major_axis_radius=2.25,
    outer_axis_ratio=0.9,
    outer_phi=160.0,
)

# set white positions
pos_white = np.genfromtxt(f'{object_folder}{links_id}_white_positions_grid.csv', delimiter=',', skip_header=0)
imaging_white.positions = al.Grid2DIrregular(
    [(pos_white[0]), (pos_white[1]), (pos_white[2]), (pos_white[3])]#, (pos[4]),]
)

#visuals_2d = aplt.Visuals2D(mask=mask)

# plot subplots for first view
#print('Plotting r-band')
#imaging_plotter_r = aplt.ImagingPlotter( # this is where the noise is coming up weird
#   imaging=imaging_r, visuals_2d=aplt.Visuals2D(mask=mask)
#)
#imaging_plotter_r.subplot_imaging()

# load g-band imaging
#print('Plotting g-band')
#imaging_plotter_g = aplt.ImagingPlotter( # this is where the noise is coming up weird
#   imaging=imaging_g, visuals_2d=aplt.Visuals2D(mask=mask)
#)
#imaging_plotter_g.subplot_imaging()

# set up grid and settings
settings_masked_imaging = al.SettingsMaskedImaging(grid_class=al.Grid2D)#, psf_shape_2d=imaging.psf.shape_2d)
settings = al.SettingsPhaseImaging(settings_masked_imaging=settings_masked_imaging)

#set up lens light profile
lens_start = al.GalaxyModel(
   redshift=zlens, bulge=al.lp.EllipticalSersic#, mass=al.mp.EllipticalIsothermal
)

# set priors
# lens position
lens_start.bulge.centre_0 = af.UniformPrior(lower_limit=-0.1, upper_limit=0.1)
lens_start.bulge.centre_1 = af.UniformPrior(lower_limit=-0.1, upper_limit=0.1)
# effective radius
#lens_start.bulge.effective_radius = af.UniformPrior(lower_limit=0.0, upper_limit=5.0) # why have I chosen 3 here? because the mask is 3...


# In[ ]:


# Start a pandas dataframe
#performance_log = pd.DataFrame(columns=['Experiment', 
#                                        'phase1_time', 
#                                        'phase1_likelihood', 
#                                        'phase2_time', 
#                                        'phase2_likelihood' 
#                                        ])

# load performance log from csv
performance_log = pd.read_csv(f'{csv_path}G{gama_id}_{links_id}_performance_log.csv')
#print(performance_log)
print(performance_log)


# ### Phase 1

# In[ ]:


# set experiment number
#experiment_number = '4'

lens_start.bulge.effective_radius = af.GaussianPrior(mean=re_r, sigma=re_r_err, lower_limit=0.0, upper_limit=re_r+3*re_r_err)
#lens_start.bulge.centre_0 = af.UniformPrior(lower_limit=0.0, upper_limit=0.1)
#lens_start.bulge.centre_1 = af.UniformPrior(lower_limit=-0.1, upper_limit=0.0)

# set up mask
#lens_mask = al.Mask2D.circular(
#    shape_native=imaging_r.shape_native, pixel_scales=imaging_r.pixel_scales, sub_size=2, radius=re_r, #centre=(0.075, -0.075) 
#)
#visuals_2d = aplt.Visuals2D(mask=lens_mask)

# plot subplots for first view
#print('Plotting r-band')
#imaging_plotter_r = aplt.ImagingPlotter(
#    imaging=imaging_r, visuals_2d=visuals_2d
#)
#imaging_plotter_r.subplot_imaging()

# set up phase
phase1 = al.PhaseImaging(
    search=af.DynestyStatic(
        path_prefix=f"{output_folder}", name=f"experiment_{experiment_number}_phase1_{datetime}", n_live_points=200,
        evidence_tolerance=0.5, walks = 10
    ),
    settings=settings,
    galaxies=af.CollectionPriorModel(lens=lens_start)#, source=source)
)

print(lens_start)


# In[ ]:


# run the phase
print('Phase running...')
tick = time.perf_counter()
phase1_result = phase1.run(dataset=imaging_r, mask=lens_mask)
tock = time.perf_counter()
print(f'Work complete! Took us {tock-tick} seconds or {(tock-tick)/60} minutes.')


# In[ ]:


# write the results to the log
print(f'Log likelihood: {phase1_result.log_likelihood}')
print(f'Model: {phase1_result.model}')

phase1_time = tock-tick

# plot it!
#fit_imaging_plotter = aplt.FitImagingPlotter(fit=phase1_result.max_log_likelihood_fit)
#fit_imaging_plotter.subplot_fit_imaging()


# ### Phase 2

# In[ ]:


# now phase 2!

print('Phasors set to stun.')

#set up lens and source

# set stellar mass profile
mass = af.PriorModel(al.mp.EllipticalIsothermal)
mass.take_attributes(source=phase1_result.model.galaxies.lens.bulge)

dark = af.PriorModel(al.mp.SphericalNFWMCRLudlow)
dark.mass_at_200 = af.LogUniformPrior(lower_limit=1e8, upper_limit=1e15)
dark.redshift_object = zlens
dark.redshift_source = zsource

lens = al.GalaxyModel(
    redshift=zlens, mass=mass, dark=dark
)

source = al.GalaxyModel(
    redshift=zsource, bulge=al.lp.SphericalExponential)


# fix mass center
lens.mass.centre = phase1_result.instance.galaxies.lens.bulge.centre
lens.dark.centre = phase1_result.instance.galaxies.lens.bulge.centre

# fix lens elliptical comps
lens.mass.elliptical_comps = phase1_result.instance.galaxies.lens.bulge.elliptical_comps

# einstein radius
lens.mass.einstein_radius = af.GaussianPrior(mean=einstein_radius, sigma=0.5*einstein_radius, 
                                             lower_limit=0.0, upper_limit=2.5) # take sigma to be 50% of mean # hmmm

# source position
source.bulge.centre_0 = af.UniformPrior(lower_limit=-5, upper_limit=5)
source.bulge.centre_1 = af.UniformPrior(lower_limit=-5, upper_limit=5)
source.bulge.effective_radius = af.UniformPrior(lower_limit=0.0, upper_limit=3.0)
#source.bulge.intensity = af.UniformPrior(lower_limit=0.0, upper_limit=10*lens.bulge.intensity)

print(f'Lens: {lens}')
print(f'Source: {source}')


# In[25]:



# Set up the positions... (GUI is not working...)

# plot the r-band image to see it
#print('Plotting r-band')
#imaging_plotter_r = aplt.ImagingPlotter( # this is where the noise is coming up weird
#    imaging=imaging_r, visuals_2d=aplt.Visuals2D(mask=source_mask)
#)
#imaging_plotter_r.figures(image=True)

# set positions
#imaging_r.positions = al.Grid2DIrregular(
#    [(0.75, -1.7), (-0.7, 1.0), (-0.95, 0.5), (0.25, 1.0)]
#)

# plot the image
#visuals_2d = aplt.Visuals2D(mask=source_mask, positions=imaging_r.positions)
#imaging_plotter_r = aplt.ImagingPlotter(imaging=imaging_r, visuals_2d=visuals_2d)
#imaging_plotter_r.figures(image=True)

# set the settings to include the positions
settings_lens = al.SettingsLens(positions_threshold=1.0)

settings = al.SettingsPhaseImaging(
    settings_masked_imaging=settings_masked_imaging, settings_lens=settings_lens
)


# In[ ]:


# set up phase
phase2 = al.PhaseImaging(
    search=af.DynestyStatic(
        path_prefix=f'{output_folder}', name=f"experiment_{experiment_number}_phase2_fit_{datetime}", n_live_points=300,
        evidence_tolerance=0.25, walks=10, facc=0.3
    ),
    settings=settings,
    galaxies=af.CollectionPriorModel(lens=lens, source=source)#, source=source)
)


# In[ ]:


# run the phase
print('Phase running...')
tick = time.perf_counter()
phase2_result = phase2.run(dataset=imaging_white, mask=source_mask)
tock = time.perf_counter()
print(f'Work complete! Took us {tock-tick} seconds or {(tock-tick)/60} minutes.')


# In[ ]:


# write the results to the log
print(f'Log likelihood: {phase2_result.log_likelihood}')
print(f'Model: {phase2_result.model}')

phase2_time=tock-tick

# plot it!
#fit_imaging_plotter = aplt.FitImagingPlotter(fit=phase2_result.max_log_likelihood_fit)
#fit_imaging_plotter.subplot_fit_imaging()


# ### Phase 3

# In[ ]:


# now phase 3!

print('Phasors set to kill.')

#update mask to be centered on lens
#mask.centre = phase1_result.model.galaxies.lens.bulge.centre

#set up lens and source

# set stellar mass/light profile
bulge = af.PriorModel(al.lmp.EllipticalSersic)
bulge.take_attributes(source=phase1_result.model.galaxies.lens.bulge)

# set dark matter profile, take prior from previous
dark = af.PriorModel(al.mp.SphericalNFWMCRLudlow)
dark.take_attributes(source=phase2_result.model.galaxies.lens.dark)
                    
lens = al.GalaxyModel(
    redshift=zlens, bulge=bulge, dark=dark
)

source_bulge = af.PriorModel(al.lp.SphericalExponential)
source_bulge.take_attributes(source=phase2_result.model.galaxies.source.bulge)

source = al.GalaxyModel(
    redshift=zsource, bulge=source_bulge)


# make dark matter centered at stellar mass center
lens.bulge.centre = phase1_result.instance.galaxies.lens.bulge.centre
lens.dark.centre = lens.bulge.centre

# make lens effective radius upper limit
#lens.bulge.effective_radius.upper_limit=re_r+3*re_r_err

# einstein radius
#lens.mass.einstein_radius = af.GaussianPrior(mean=einstein_radius, sigma=0.3*einstein_radius) # take sigma to be 30% of mean # hmmm

#source.bulge.intensity = af.UniformPrior(lower_limit=0.0, upper_limit=10*lens.bulge.intensity)

print(f'Lens: {lens}')
print(f'Source: {source}')


# In[ ]:


# set the settings to include the positions
settings_lens = al.SettingsLens(positions_threshold=1.0)

settings = al.SettingsPhaseImaging(
    settings_masked_imaging=settings_masked_imaging, settings_lens=settings_lens
)


# In[ ]:


# set up phase
phase3 = al.PhaseImaging(
    search=af.DynestyStatic(
        path_prefix=f'{output_folder}', name=f"experiment_{experiment_number}_phase3_fit_{datetime}", n_live_points=500,
        evidence_tolerance=0.25, walks=10, facc=0.3
    ),
    settings=settings,
    galaxies=af.CollectionPriorModel(lens=lens, source=source)#, source=source)
)


# In[ ]:


# run phase
# run the phase
print('Phase running...')
tick = time.perf_counter()
phase3_result = phase3.run(dataset=imaging_white, mask=mask)
tock = time.perf_counter()
print(f'Work complete! Took us {tock-tick} seconds or {(tock-tick)/60} minutes.')


# In[ ]:


# write the results to the log
print(f'Log likelihood: {phase3_result.log_likelihood}')
print(f'Model: {phase3_result.model}')

phase3_time=tock-tick

log = open(f"{output_folder}experiment_log", 'a') # append the log
lines = [f'Experiment {experiment_number} \n',
         f'Phase 1 \n',
         f'Time to convergence: {phase1_time} seconds \n', 
         f'Log likelihood: {phase1_result.log_likelihood} \n',
         f'Model: {str(phase1_result.model)} \n',
         f'Phase 2 \n',
         f'Time to convergence: {phase2_time} seconds \n', 
         f'Log likelihood: {phase2_result.log_likelihood} \n', 
         f'Model: {str(phase2_result.model)} \n',
         f'Phase 3 \n',
         f'Time to convergence: {phase3_time} seconds \n', 
         f'Log likelihood: {phase3_result.log_likelihood} \n', 
         f'Model: {str(phase3_result.model)} \n','\n'] # set lines to write the model result
log.writelines(lines) # write lines
log.close()


# In[44]:


# append performance log
data_list = [[experiment_number, 
            phase1_time,
            phase1_result.log_likelihood,
            phase2_time, 
            phase2_result.log_likelihood,
            phase3_time, 
            phase3_result.log_likelihood,
             ]]
new_dataframe_entry = pd.DataFrame(data_list,
                                  columns=['Experiment', 
                                           'phase1_time', 
                                           'phase1_likelihood', 
                                           'phase2_time', 
                                           'phase2_likelihood',
                                           'phase3_time', 
                                           'phase3_likelihood'
                                           ])
performance_log = pd.concat([performance_log, new_dataframe_entry])
print(performance_log)
performance_log.to_csv(f'{csv_path}G{gama_id}_{links_id}_performance_log.csv')
#print(phase_result.model)
# get samples to see how it ran?
#log_likelihoods = phase_result.samples.log_likelihoods
#x = np.arange(0, len(log_likelihoods), 1)
#plt.plot(x[500:], log_likelihoods[500:])
#plt.xlim()

# plot it!
#fit_imaging_plotter = aplt.FitImagingPlotter(fit=phase3_result.max_log_likelihood_fit)
#fit_imaging_plotter.subplot_fit_imaging()

