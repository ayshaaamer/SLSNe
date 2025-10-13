"""
This file contains general utilities to plot and use the spectra of SLSNe.

Written by Aysha Aamer, 2025.
"""

# Import necessary packages
import matplotlib.pyplot as plt
import pandas as pd
from astropy import table
import numpy as np
import glob
import os
from google.oauth2 import service_account
import gspread
import json
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
import random 
import scipy as sp
import math
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


# Get directory with reference data
current_file_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(current_file_dir, 'ref_data')

# Import reference data and supernovae names
data_table = table.Table.read(os.path.join(data_dir, 'sne_data.txt'), format='ascii')
use_names = np.genfromtxt(os.path.join(data_dir, 'use_names.txt'), dtype='str')

# Only use names that are in the reference data directory
directories = glob.glob(os.path.join(data_dir, 'supernovae', '*'))
exists = [i.split('/')[-1] for i in directories]
# Only select the items in data_table that are in use_names and in good
data_table = data_table[np.isin(data_table['Name'], use_names) & np.isin(data_table['Name'], exists)]

# Creating custom colourmap for plots
# Define the colours for the custom colormap
colors = ["#8ea0cb", "#efc2d9", "#ffffb4"]  # Pastel purple, pastel pink, pale yellow
# Create the custom colormap
custom_cm = LinearSegmentedColormap.from_list("pastel_cmap", colors)


def get_plot_data():
    """This function extracts the redshifts and phases of SLSNe spectra from folders to make plots.
    
    Returns
    -------
    phases : list of float
        Rest-frame phases (in days) of each spectrum relative to peak light.
    phases_exp : list of float
        Rest-frame phases (in days) of each spectrum relative to explosion.
    objects : list of str
        Names of the SLSNe corresponding to each spectrum. 
    """

    # Reading in SN names and parameter file
    all_sne = glob.glob(f'supernovae_final/*')
    sne_params = pd.read_csv('../sne_data.txt', delim_whitespace=True)

    # Initialising empty arrays
    phases = []
    phases_exp = []
    objects = []

    # Looping through each folder containing the spectra
    for sn in all_sne:
        spectra = glob.glob(sn + '/*')

        # Getting the SN name from the filename and reading file containing parameters
        obj = sn.split('/')[-1]
        sn_params = sne_params[sne_params['Name'] == obj]
        exp_date = float(sn_params['Explosion'].values[0])      # Time of explosion (MJD)
        max_date = float(sn_params['Peak'].values[0])           # Time of maximum light (MJD)

        # Looping through each spectrum within the file
        for filename in spectra:
            data = table.Table.read(filename, format='ascii')
            header = table.Table.read(data.meta['comments'], delimiter='=', format='ascii.no_header',
                                        names=['key', 'val'])
            MJD = float(data.meta['comments'][np.where(header['key'] == 'MJD')[0][0]].split('=')[1].strip())

            # Calculating phase from explosion and peak, and converting to rest frame
            phase = round((MJD-max_date)/(1+z), 3)
            phase_exp = round((MJD-exp_date)/(1+z), 3)

            # Only accounting for spectra with phases less than 200 days post peak
            if phase < 200: 
                phases.append(phase)
                phases_exp.append(phase_exp)
                objects.append(obj)


    return phases, phases_exp, objects


def make_spec_phase_distribution_plots(output_dir='.', binwidth=10):
    """Generate and save histograms showing the distribution of SLSN spectra versus phase. Each bin is colour 
    coded according to the number of unique objects in the bin.

    This function creates two histograms:
      1. Number of spectra versus phase relative to maximum light.
      2. Number of spectra versus phase relative to explosion.

    Parameters
    ----------
    output_dir : str, optional
        Directory where the output plots will be saved (default is current directory).
    binwidth : int, optional
        Width of the histogram bins in days (default is 10).

    Returns
    -------
    None
        The function saves two PDF plots:
          - `phase_histogram_peak_<binwidth>.pdf`
          - `phase_histogram_exp_<binwidth>.pdf`"""

    # Setting bin width and edges of bins
    min_bin = math.floor(min(phases)/binwidth)*binwidth
    max_bin = math.ceil(max(phases)/binwidth)*binwidth
    min_bin_exp = math.floor(min(phases_exp)/binwidth)*binwidth
    max_bin_exp = math.ceil(max(phases_exp)/binwidth)*binwidth

    # Getting phase information of spectra
    phases, phases_exp, objects = get_plot_data()

    # Computing the number of unique SNe contributing to each phase bin (for colour weighting)
    weights = []
    df = pd.DataFrame({'Phase':phases, 'Object':objects})
    for min_phase in range(min_bin, max_bin+binwidth, binwidth):
        # print(min_phase)
        df_subset = df[(df['Phase'] >= min_phase) & (df['Phase'] < min_phase+binwidth)]
        weights.append(len(df_subset['Object'].unique()))

    weights_exp= []
    df_exp = pd.DataFrame({'Phase':phases_exp, 'Object':objects})
    for min_phase in range(min_bin_exp, max_bin_exp+binwidth, binwidth):
        df_subset = df_exp[(df_exp['Phase'] >= min_phase) & (df_exp['Phase'] < min_phase+binwidth)]
        weights_exp.append(len(df_subset['Object'].unique()))
    

    # Normalise colour scale for each histogram based on number of unique objects
    print(range(min_bin, max_bin+binwidth, binwidth))
    hist, bins = np.histogram(phases, bins=range(min_bin, max_bin+binwidth, binwidth))
    norm = Normalize(vmin=min(weights), vmax=max(weights))
    colors = custom_cm(norm(weights))

    hist_exp, bins_exp = np.histogram(phases_exp, bins=range(min_bin_exp, max_bin_exp+binwidth, binwidth))
    norm_exp = Normalize(vmin=min(weights_exp), vmax=max(weights_exp))
    colors_exp = custom_cm(norm(weights_exp))


    # Plotting a histogram of the phases from peak
    plt.clf()
    fig, ax = plt.subplots(figsize=(6,4))
    bars = ax.bar(bins[:-1], hist, width=np.diff(bins), edgecolor='black', color=colors, align='edge')
    # Addding a colour bar to the plot
    sm = plt.cm.ScalarMappable(cmap=custom_cm, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Number of Objects')
    ax.set_xlabel('Phase relative to maximum light (days)')
    ax.set_ylabel('Number of Spectra')
    plot_name = 'phase_histogram_peak_{binwidth}.pdf'
    plot_dir = os.path.join(output_dir, plot_name)
    plt.savefig(plot_dir, bbox_inches='tight')

    # Plotting a histogram of the phases from explosion
    plt.clf()
    fig, ax = plt.subplots(figsize=(6,4))
    bars = ax.bar(bins_exp[:-1], hist_exp, width=np.diff(bins_exp), edgecolor='black', color=colors_exp, align='edge')
    # Adding a colour bar
    sm = plt.cm.ScalarMappable(cmap=custom_cm, norm=norm_exp)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Number of Objects')
    ax.set_xlabel('Phase relative to explosion (days)')
    ax.set_ylabel('Number of Spectra')
    plot_name = 'phase_histogram_exp_{binwidth}.pdf'
    plot_dir = os.path.join(output_dir, plot_name)
    plt.savefig(plot_dir, bbox_inches='tight')



def redshit_plots():
    # Looping through all supernovae
    all_sne = glob.glob(f'supernovae_final/*')
    
    # Initialising arrays for redshift
    redshifts = []
    redshifts_all = []
    len_spec = []
    objects = []
    cm = plt.get_cmap('Set3')

    # Getting redshift information of spectra
    phases, phases_exp, redshifts, objects = get_plot_data()

    # Each SN in new folder so looping through them all
    for sn in all_sne:
        spectra = glob.glob(sn + '/*')
        print('Processing', sn)
        # Extracting object name and redshift from filename
        obj = sn.split('/')[-1]
                z = float(sn_params['Redshift'].values)                 # Redshift 
        z = get_redshift(obj)
        redshifts_all.append(z)
        len_spec.append(len(spectra))

        if len(spectra) >= 3:
            redshifts.append(z)
            objects.append(obj)


    # bin width and edges
    binwidth = 0.1
    max_bin = 2
    colors = [cm(9), cm(5)]
    

    plt.clf()
    plt.hist([redshifts_all, redshifts], bins=np.arange(0, max_bin+binwidth, binwidth), alpha=0.5, 
             label=['All objects', 'Objects with $\geq$3 spectra'], color=colors)
    plt.xlabel('Redshift')
    plt.ylabel('Number of Objects')
    plt.legend()
    plt.xlim(-0.05, 2.05)
    
    
    # plt.savefig('Figures/redshift_distribution_combined.pdf', bbox_inches='tight')



def no_spec_plots(master_spreadsheet):

    cm = plt.get_cmap('Set3')

    all_sne = glob.glob(f'supernovae_final/*')
    no_spectra = []
    objects = []
    redshifts = []
    redshifts_all = []
    peak_mags = []
    peak_mags_all = []
    counter1 = 0

    # Each SN in new folder so looping through them all
    for sn in all_sne:
        # spectra = glob.glob(sn + "/*.txt")
        spectra = glob.glob(sn + '/*')
        # print('Processing', sn)
        # print(len(spectra))

        obj = sn.split('/')[-1]

        no_spectra.append(len(spectra))
        objects.append(obj)
        z = get_redshift(obj)
        redshifts_all.append(z)
        peak_mag = str((master_spreadsheet[master_spreadsheet['Name'] == obj]['Peak Abs. Mag']).values[0])

        # if peak mag available, adding to list
        if peak_mag.split(' ')[0] != '':
            peak_mags_all.append(float(peak_mag.split(' ')[0]))
        
        # if only 1 spectra available
        if len(spectra) == 1:
            counter1 +=1
            redshifts.append(z)

            if peak_mag.split(' ')[0] != '':
                peak_mags.append(float(peak_mag.split(' ')[0]))

        # if over 50 spectra for the object
        if len(spectra) >= 20:
            print('this object has lots of spectra')
            print(obj)
            print('This object has ', len(spectra), ' spectra')
            print('redshift is: ', get_redshift(obj))
            print('peak observed mag is: ', peak_mag)
            # redshifts_high.append(get_redshift(obj))
            print()

        # # if redshift over 1.75
        # if z >= 1.75:
        #     print('this object has a high redshift')
        #     print(obj)
        #     print('This object has ', len(spectra), ' spectra')
        #     print('redshift is: ', get_redshift(obj))
        #     print('peak observed mag is: ', peak_mag)
        #     print()


    print(len(no_spectra))
    print('number of objects with 1 spectra: ', counter1)
    # print('redshifts: ', redshifts)
    # print('peak mag: ', peak_mags)
    # print('min redshift: ', min(redshifts_all))



    # bin width and edges
    binwidth = 1
    max_bin = max(no_spectra)
    
    plt.figure(1, figsize=(6,4))
    plt.clf()
    plt.hist(no_spectra, bins=np.arange(1, max_bin+binwidth, binwidth), color=cm(0))
    plt.yscale('log')
    plt.xlabel('Number of Spectra')
    plt.ylabel('Number of Objects') 
    # plt.xlim(-200, 750)
    plt.savefig('Figures/number_spectra_objects_log.pdf', bbox_inches='tight')
    # plt.show()

    # plt.figure(2, figsize=(6,3))
    # plt.clf()
    # plt.hist(redshifts)
    # plt.xlabel('Redshift')
    # plt.ylabel('Number of Objects')
    # # plt.xlim(-200, 750)
    # plt.savefig('Figures/redshift_1spec.pdf', bbox_inches='tight')
    # plt.show()

    # plt.figure(3, figsize=(6,3))
    # plt.clf()
    # plt.hist(peak_mags)
    # plt.xlabel('Peak Apparent Magnitude')
    # plt.ylabel('Number of Objects')
    # # plt.xlim(-200, 750)
    # plt.savefig('Figures/peak_mag_1spec.pdf', bbox_inches='tight')
    # plt.show()

    # plt.figure(4, figsize=(6,3))
    # plt.clf()
    # plt.hist(peak_mags_all)
    # plt.xlabel('Peak Apparent Magnitude')
    # plt.ylabel('Number of Objects')
    # # plt.xlim(-200, 750)
    # plt.savefig('Figures/peak_mag.pdf', bbox_inches='tight')
    # plt.show()

    # plt.figure(5, figsize=(6,3))
    # plt.clf()
    # plt.hist(peak_mags_all)
    # plt.xlabel('Peak Absolute Magnitude')
    # plt.ylabel('Number of Objects')
    # # plt.xlim(-200, 750)
    # plt.savefig('Figures/peak_abs_mag.pdf', bbox_inches='tight')
    # plt.show()



