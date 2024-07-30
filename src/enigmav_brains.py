#!/usr/bin/env python
import os
os.environ['SUBJECTS_DIR'] = '/usr/local/opt/freesurfer/subjects'
# os.environ['ETS_TOOLKIT'] = 'qt5'
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

"""IMPORTANT: To run cortical function, you need A NumPy version >=1.16.5 and <1.23.0"""

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

import sys
import mayavi
import numpy as np
import nibabel as nib
from surfer import Brain
from numpy import genfromtxt
from enigmatoolbox.datasets import load_example_data
from enigmatoolbox.utils.useful import zscore_matrix
from enigmatoolbox.plotting import plot_subcortical
from enigmatoolbox.utils.useful import reorder_sctx
import pandas as pd
import matplotlib.pylab as plt
#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = ['Arial']
plt.rc('font', family='Arial')
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.pyplot import figure
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap


def get_hex_color(val_lst):
    """
    Converts a z-value to a hex color using a given colormap.
    
    Parameters:
    z_value (float): The z-value to convert.
    cmap (str): The name of the colormap to use (default: 'viridis').
    
    Returns:
    str: The hex color string.
    """
    # Normalize the z-value to [0, 1]
    hex_lst = []
    cmap = 'RdBu_r'
    z_min = min(val_lst)
    z_max = max(val_lst)

    norm = colors.Normalize(vmin=z_min, vmax=z_max)

    for val in val_lst:
        norm_value = norm(val)
        
        # Get the RGBA color tuple from the colormap
        cmap = cm.get_cmap(cmap)
        rgba_color = cmap(norm_value)
        
        # Convert the RGBA color tuple to a hex color string
        hex_color = colors.rgb2hex(rgba_color)
        hex_lst.append(hex_color)

    
    return hex_lst


def add_colorbar(grid, fig, cmap_name,
                 y_min=0, y_max=1, cbar_min=0,
                 cbarlabelsize=18,
                 cbarminlabel=0, cbarmaxlabel=1,
                 cbar_max=1, vert=False, label='Percentile',
                 show_ticks=False, pad=-18.5):
    '''
    Add a colorbar to the big_fig in the location defined by grid 
    
    grid       :  grid spec location to add colormap
    fig        :  figure to which colorbar will be added
    cmap_name  :  name of the colormap
    x_min      :  the minimum value to plot this colorbar between
    x_max      :  the maximum value to plot this colorbar between
    cbar_min   :  minimum value for the colormap (default 0)
    cbar_max   :  maximum value for the colormap (default 1)
    vert       :  whether the colorbar should be vertical 
    label      :  the label for the colorbar (default: None)
    ticks      :  whether to put the tick values on the colorbar 
    pad        :  how much to shift the colorbar label by (default: 0)
    '''


    # Add an axis to the big_fig
    ax_cbar = plt.Subplot(fig, grid)
    fig.add_subplot(ax_cbar)
    
    # Normalise the colorbar so you have the correct upper and
    # lower limits and define the three ticks you want to show
    norm = mpl.colors.Normalize(vmin=cbar_min, vmax=cbar_max)

    if show_ticks:
        #ticks = [y_min, np.average([y_min, y_max]), y_max]
        ticks = [y_min, y_max]
    else:
        ticks=[]
            
    # Figure out the orientation
    if vert:
        orientation='vertical'
        rotation=270
    else:
        orientation='horizontal'
        rotation=0
        
    # Add in your colorbar:
    cb = mpl.colorbar.ColorbarBase(ax_cbar, cmap=cmap_name, norm=norm,
                                   orientation=orientation,
                                   ticks=ticks,
                                   boundaries=np.linspace(y_min,
                                                          y_max, 100))
         
    cb.outline.set_visible(False)                              
    cb.ax.set_xticklabels([cbarminlabel, cbarmaxlabel])

    cb.ax.tick_params(labelsize=cbarlabelsize, length=0) 
    #cb.set_family("Arial")
    #ax = cb.ax
    #text = ax.xaxis.label
    #font = plt.font_manager.FontProperties(family='Arial',
    #text.set_font_properties(font)
    #show()

    if label:
        cb.set_label(label, rotation=rotation, labelpad=pad,
                     size=cbarlabelsize, weight='bold', family='Arial')
    return fig


def plot_dk_parc(subid='fsaverage', atlas='dk'):
    import os
    from os.path import join as pjoin
    from surfer import Brain

    subid = 'fsaverage'
    hemi = 'lh'
    surf = 'pial'
    view = 'lateral'

    if atlas == 'dk':
        annot = 'aparc'
    else:
        annot = 'aparc.a2009s'
        

    """
    Bring up the visualization
    """
    brain = Brain(subid, hemi, surf, views=view,
                  cortex="bone", background="white")
    brain.add_annotation(annot, borders=False)
    prefix = 'enigmav_freesurfer_dkparc_lh'
    brain.save_imageset(prefix=prefix,
                        views=['med', 'lat'], filetype='png')


def cortical():
    # Load percentile scores for each region
    # surf_inds = [
    #     'CurvInd',
    #     'FoldInd',
    #     'GausCurv',
    #     'MeanCurv',
    #     'SurfAvg',
    #     'ThickAvg'
    #     ]
    
    surf_inds = [
        'curvind',
        'foldind',
        'gauscurv',
        'meancurv',
        'surfavg',
        'thickavg'
        ]

    hemis = ['lh', 'rh']
    subid = "fsaverage"
    figpath = "../output/figures/"

    for s in surf_inds:
        #print(os.getcwd())
        mydata = pd.read_csv('../output/tables/enigmav_' + s + '.csv')
        #x = np.array([1])
        #mydata = [x[1] for x in mydata]
        mydata = mydata.zval.values
        x = mydata
        ml = len(mydata)
        rhdata = np.insert(mydata[34:68], 3, x)
        rhdata = np.insert(rhdata, 0, x)
        lhdata = np.insert(mydata[0:34], 3, x)
        lhdata = np.insert(lhdata, 0, x)

        cmap="RdBu_r"
        y_min=0
        y_max=1
        clim=round(max(abs(min(mydata)), max(mydata)), 2)
        cbarminlabel=-clim
        cbarmaxlabel=clim
        cbarlabel="z"
        cbarlabelsize=12
        hemidata = [lhdata, rhdata] 


        # Save brain images for both hemispheres
        for h in range(0, len(hemis)):
            aparc_file = os.path.join(os.environ["SUBJECTS_DIR"],
                                        "fsaverage", "label",
                                        hemis[h] + ".aparc.annot")
            labels, ctab, names = nib.freesurfer.read_annot(aparc_file)
            vtx_data = hemidata[h][labels]
            vtx_data[labels==-1] = -99 
            brain = Brain(subject_id='fsaverage', hemi=hemis[h], surf='pial',
                            background='white',
                            cortex='white')
            #brain = Brain(subid, 'split', surf, views=['lat', 'med'])
            brain.add_data(vtx_data, min=-clim, max=clim, thresh=-clim,
                            colormap=cmap, alpha=.8)
            prefix = figpath + subid + "_" + s + "_" + hemis[h]
            #print prefix
            brain.save_imageset(prefix=prefix, colorbar=None,
                                views=['med', 'lat'], filetype='png')
            
        figsize = (7, 2)
        fig = plt.figure(figsize=figsize, facecolor='white')

        grid = gridspec.GridSpec(1, 4)
        grid.update(left=0.01, right=0.99, top=0.99, bottom=0.2, wspace=0,
                    hspace=0)

        f_list = [ '_'.join([figpath + subid, s, "lh", 'lat.png']),
                '_'.join([figpath + subid, s, "lh", 'med.png']),
                '_'.join([figpath + subid, s, "rh", 'med.png']),
                '_'.join([figpath + subid, s, "rh", 'lat.png']) ]
        # Plot each figure in turn
        for g_loc, f in zip(grid, f_list):
            ax = plt.Subplot(fig, g_loc)
            fig.add_subplot(ax)
            img = mpimg.imread(f)

            # Crop the figures appropriately NOTE: this can change depending
            # on which system you've made the images on originally - it's a
            # bug that needs to be sorted out!
            #if 'lat' in f:
            #    img_cropped = img[130:(-130), 5:(-40), :]
            #else:
            #    img_cropped = img[120:(-100), :, :]
            #ax.imshow(img_cropped, interpolation='none')
            ax.imshow(img, interpolation='none')
            ax.set_axis_off()

        # Add the bottom of one of the images as the color bar
        # at the bottom of the combo figure
        grid_cbar = gridspec.GridSpec(1, 1)
        grid_cbar.update(left=0.4, right=0.6, top=0.26, bottom=0.18,
                        wspace=0, hspace=0)

        # Save the figure
        fname=figpath + subid + "_" + s + "_" + "combined.png"
        fig = add_colorbar(grid_cbar[0], fig, cmap, show_ticks=True,
                        cbarminlabel=cbarminlabel,
                        cbarmaxlabel=cbarmaxlabel,
                        cbarlabelsize=cbarlabelsize,
                        pad=-10, label=cbarlabel, y_max=y_max, y_min=y_min)
        fig.savefig(fname, bbox_inches=0, dpi=300)
        #plt.show()
        plt.close()
        del(brain)
    return fig

def subcortical():
    """
    Extracts labels from enigima toolbox,
    compares them to labels in data file,
    assigns labels a zval,
    plots subcortical zvals
    """
    # Load all example data from an individual site
    cov, metr1_SubVol, metr2_CortThick, metr3_CortSurf = load_example_data()

    # Re-order the subcortical data alphabetically and by hemisphere
    metr1_SubVol_r = reorder_sctx(metr1_SubVol)

    # Z-score patients' data relative to controls (lower z-score = more atrophy)
    group = cov['Dx'].to_list()
    controlCode = 0
    SV_z = zscore_matrix(metr1_SubVol_r.iloc[:, 1:-1], group, controlCode)

    # Mean z-score values across individuals with from a specific group (e.g., left TLE, that is SDx == 3)
    SV_z_mean = SV_z.iloc[cov[cov['SDx'] == 3].index, :].mean(axis=0)

    # get our subcortical labels
    sub_data = pd.read_csv('../output/tables/enigmav_landrvolumes.csv')

    # rename ventricles to match enigma atlas
    sub_data.roi = sub_data.roi.replace('Right_Lateral_Ventricle', 'RLatVent')
    sub_data.roi = sub_data.roi.replace('Left_Lateral_Ventricle', 'LLatVent')

    #find matching labels
    idx_eng_lbl = np.in1d(sub_data.roi.values, SV_z_mean.keys()) 
    df_match = sub_data[idx_eng_lbl]

    # customize order according to enigmatoolbox
    custom_order = SV_z_mean.keys().tolist()
    df_match['enigmaorder'] = df_match.roi.apply(lambda x: custom_order.index(x))
    df_match = df_match.sort_values(by='enigmaorder')
    df_match = df_match.drop(columns=['enigmaorder'])
    
    sub_zmean = df_match.zval
    clim=round(max(abs(min(sub_zmean)), max(sub_zmean)), 2)
    cbarminlabel=-clim
    cbarmaxlabel=clim

    # Project the results on the subcortical surface
    plot_subcortical(array_name=sub_zmean, size=(1500, 500),background=(1, 1, 1), embed_nb=True, screenshot=True, filename='../output/figures/subcortical_both_screenshot.png',
                 cmap='RdBu_r', color_bar=False, transparent_bg=False, color_range=(cbarminlabel, cbarmaxlabel))
    
    # generate fig from picture 
    wm3dp = '../output/figures/subcortical_both_screenshot.png'
    im = mpimg.imread(wm3dp)

    mydata = sub_zmean
    my_dpi = 144

    cmap="RdBu_r"
    y_min=0
    y_max=1
    clim=round(max(abs(min(mydata)), max(mydata)), 2)
    cbarminlabel=-clim
    cbarmaxlabel=clim

    f = plt.figure(figsize=(2672/my_dpi, 866/my_dpi), dpi=my_dpi)
    f.add_subplot(1,2, 1)
    plt.imshow(im)
    plt.axis('off')

    #im_ratio = (600/my_dpi)/(1800/my_dpi)
    #plt.colorbar(label='z', cbarminlabel, cbarmaxlabel,fraction=0.047*im_ratio,orientation='horizontal')

    cax = f.add_axes([0.25, 0.33, 0.1, 0.04]) #has to be as a list - starts with x, y coordinates for start and then width and height in % of figure width
    norm = mpl.colors.Normalize(vmin = cbarminlabel, vmax = cbarmaxlabel)     
    cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm = norm, orientation = 'horizontal')#, boundaries=np.linspace(cbarminlabel,cbarmaxlabel, 100))
            
    cb.outline.set_visible(False)                      
    cb.ax.set_xticklabels(['',cbarminlabel,'',cbarmaxlabel])
    #cb.ax.set_xticklabels([cbarminlabel,cbarmaxlabel])
    #cb.ax.tick_params(labelsize=10, length=0) 
    cb.set_label('z', rotation=0, labelpad=-10,
                        size=10, weight='bold', family='Arial')   

    #plt.show(block=True)
    plt.savefig('../output/figures/enigmav_subcort_plot.png')
    plt.close()
    return f


def white_matter():
    """
    enigma_wm_right/left_masked.png images were generated with enigmav_wmtracts.R script,
    that script adjusts the color depending on hex_code in output/tables/enigmav_zval_wm_20_3d.csv,
    hex-code can be generated running get_hex_color(val_lst) in this script.
    images are imported as figures and plotted together with colorbar
    """
    wm_3d_path = '../output/tables/enigmav_zval_wm_20_3d.csv'
    df_wm_3d = pd.read_csv(wm_3d_path)
    
    # generate fig from picture 
    #wm3dp_rh = '../output/figures/enigmav_wm_right_masked.png'
    #wm3dp_lh = '../output/figures/enigmav_wm_left_masked.png'
    wm3dp_rh = '../output/figures/enigmav_wm_right.png'
    wm3dp_lh = '../output/figures/enigmav_wm_left.png'
    im_rh = mpimg.imread(wm3dp_rh)
    im_lh = mpimg.imread(wm3dp_lh)
    mydata = df_wm_3d.zval.values
    im = [im_rh,im_lh]
    my_dpi = 300

    cmap="RdBu_r"
    y_min=0
    y_max=1
    clim=round(max(abs(min(mydata)), max(mydata)), 2)
    cbarminlabel=-clim
    cbarmaxlabel=clim

    # generate subplots and merge pictures into one
    f = plt.figure(figsize=(1800/my_dpi, 600/my_dpi), dpi=my_dpi)
    f.add_subplot(1,2, 1)
    plt.imshow(im_rh)
    plt.axis('off')
    f.add_subplot(1,2, 2)
    plt.imshow(im_lh)
    plt.axis('off')

    # add colorbar horizontally in the middle of both pictures
    cax = f.add_axes([0.41, 0.2, 0.2, 0.1]) #has to be as a list - starts with x, y coordinates for start and then width and height in % of figure width
    norm = mpl.colors.Normalize(vmin = cbarminlabel, vmax = cbarmaxlabel)     
    cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm = norm, orientation = 'horizontal')#, boundaries=np.linspace(cbarminlabel,cbarmaxlabel, 100))
            
    cb.outline.set_visible(False)                      
    cb.ax.set_xticklabels(['',cbarminlabel,'',cbarmaxlabel])
    #cb.ax.set_xticklabels([cbarminlabel,cbarmaxlabel])
    #cb.ax.tick_params(labelsize=10, length=0) 
    cb.set_label('z', rotation=0, labelpad=-10,
                        size=10, weight='bold', family='Arial')   

    #plt.show(block=True)
    plt.savefig('../output/figures/enigmav_wm_fa_plot.png')
    plt.close()
    return f

# def main(measure):
#     m_dict = {'cortical':cortical(),'subcortical':subcortical(),'fa':white_matter()}
#     fig = m_dict[measure]
#     return fig

if __name__ == "__main__":
    m_dict = {'cortical':cortical(),'subcortical':subcortical(),'fa':white_matter()}
    keys_lst = ['cortical','subcortical','fa']
    fig_lst = []
    for con in range(len(m_dict)):
        fig = m_dict[keys_lst[con]]
        fig_lst.append(fig)
