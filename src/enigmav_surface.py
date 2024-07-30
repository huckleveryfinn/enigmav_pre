#!/usr/bin/env python
import os
# os.environ['ETS_TOOLKIT'] = 'qt4'
import sys
import mayavi
import numpy as np
import nibabel as nib
#from surfer import Brain
from enigmatoolbox.datasets import load_summary_stats
from enigmatoolbox.utils.parcellation import parcel_to_surface
from enigmatoolbox.plotting import plot_cortical

# Load summary statistics for ENIGMA-22q
# sum_stats = load_summary_stats('22q')
# Get case-control cortical thickness and surface area tables
# CT = sum_stats['CortThick_case_vs_controls']
#SA = sum_stats['CortSurf_case_vs_controls']
# Extract Cohen's d values
#CT_d = CT['d_icv']
#SA_d = SA['d_icv']
# Map parcellated data to the surface
#CT_d_fsa5 = parcel_to_surface(CT_d, 'aparc_fsa5')
# Project the results on the surface brain
# plot_cortical(array_name=CT_d_fsa5, surface_name="fsa5", size=(800, 400),
#               cmap='RdBu_r', color_bar=True, color_range=(-0.5, 0.5))

surf_inds = [
    'curvind',
    'foldind',
    'gauscurv',
    'meancurv',
    'surfavg',
    'thickavg'
    ]

from numpy import genfromtxt

for s in surf_inds:
 mydata = genfromtxt('../output/tables/enigmav_' + s + '.csv')

 # Project the results on the surface brain
 r = parcel_to_surface(mydata, 'aparc_fsa5')

 # Find absolute maximum of z values
 clim=max(abs(min(mydata)), max(mydata))

 # Create brain figure
 fn='../output/figures/enigmav_fig_' + s + '.png'
 plot_cortical(array_name=r, surface_name="fsa5", size=(2000, 500),
               cmap='RdBu_r', color_bar=True, 
               color_range=(-clim, clim),
               screenshot=True,
               filename=fn)
 print(fn)


