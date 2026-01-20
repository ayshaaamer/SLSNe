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
    sne_params = pd.read_csv('ref_data/sne_data.txt', delim_whitespace=True)

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


def redshit_plots(output_dir='.', binwidth=0.1, max_bin=2):
    """
    This function a histogram of the redshift distribution of the SLSNe spectra. It overlays the distribution of all events with those that over 3 or more spectra.

    Parameters
    ----------
    output_dir : str, optional
        Directory in which to save the output plot (default is the current directory).
    binwidth : float, optional
        Width of the histogram bins in redshift (default is 0.1).
    max_bin : float, optional
        Upper limit of the redshift range for the histogram (default is 2).

    Returns
    -------
    None
        The function saves a PDF plot named
        `redshift_distribution_{binwidth}.pdf` in the specified output directory.
    """


    # Looping through all supernovae
    all_sne = glob.glob(f'supernovae_final/*')
    sne_params = pd.read_csv('ref_data/sne_data.txt', delim_whitespace=True)

    # Colours of bins
    cm = plt.get_cmap('Set3')
    colors = [cm(9), cm(5)]
    
    # Initialising arrays for redshift
    redshifts = []
    redshifts_all = []
    len_spec = []

    # Each SN in new folder so looping through them all
    for sn in all_sne:
        spectra = glob.glob(sn + '/*')
        print('Processing', sn)
        # Extracting object name and redshift from filename
        obj = sn.split('/')[-1]
        sn_params = sne_params[sne_params['Name'] == obj]
        z = float(sn_params['Redshift'].values) 
        # Appending redshfit and number of spectra to lists
        redshifts_all.append(z)
        len_spec.append(len(spectra))

        # If over 2 spectra, also appending to a separate list
        if len(spectra) >= 3:
            redshifts.append(z)

    # Plotting histogram
    plt.clf()
    plt.hist([redshifts_all, redshifts], bins=np.arange(0, max_bin+binwidth, binwidth), alpha=0.5, 
             label=['All objects', 'Objects with $\geq$3 spectra'], color=colors)
    plt.xlabel('Redshift')
    plt.ylabel('Number of Objects')
    plt.legend()
    plt.xlim(-0.05, 2.05)
    
    plot_name = 'redshift_distribution_{binwidth}.pdf'
    plot_dir = os.path.join(output_dir, plot_name)
    plt.savefig(plot_dir, bbox_inches='tight')


def load_and_plot_spectra(sn_name, flux_type='processed', output_dir='.', plot=True):
    """
    This function loads and plots spectra for a given SLSNe. The flux type can be specified to 
    return either 'raw' or 'processed flux.

    Parameters
    ----------
    sn_name : str
        Name of the supernova.
    flux_type : {'raw', 'processed'}, optional
        Which flux column to use when loading spectra (default is 'processed').
    output_dir : str, optional
        Directory in which to save the output plot (default is current directory).
    plot : bool, optional
        Whether to display and save the overlaid spectra plot (default is True).

    Returns
    -------
    spectra_list : list of dict
        A list of dictionaries, each containing:
            {
                'filename': str
                'mjd': float,
                'wavelength': numpy.ndarray,
                'flux': numpy.ndarray,
                'error': numpy.ndarray or None,
            }

    """

    # Locate all spectra files
    spectra_files = sorted(glob.glob(os.path.join(sn_name, 'spectra_raw_processed', '*.txt')))
    if not spectra_files:
        raise FileNotFoundError(f'No spectra found in {sn_name}/spectra_raw_processed/')

    spectra_list = []

    for spec_file in spectra_files:
        data = table.Table.read(spec_file, format='ascii')
        header = table.Table.read(data.meta['comments'], delimiter='=', format='ascii.no_header',
                                    names=['key', 'val'])
        MJD = float(data.meta['comments'][np.where(header['key'] == 'MJD')[0][0]].split('=')[1].strip())
        
        # Expecting columns: wavelength, raw flux, (optional error), processed flux
        if data.shape[1] < 3:
            raise ValueError(f'Unexpected column format in {spec_file}')

        wavelength = np.array(data['Wavelength'].value).astype(float)
        raw_flux = np.array(data['Raw_Flux'].value).astype(float)
        processed_flux = np.array(data['Processed_Flux'].value).astype(float)
        error = np.array(data['Error'].value).astype(float)

        flux = raw_flux if flux_type.lower() == 'raw' else processed_flux

        spectra_list.append({
            'mjd': MJD,
            'wavelength': wavelength,
            'flux': flux,
            'error': error,
            'filename': os.path.basename(spec_file)
        })

    # Sort by MJD
    spectra_list.sort(key=lambda x: x['mjd'])

    # Calculate a global median flux for scaling the offsets
    all_flux_values = np.concatenate([s['flux'] for s in spectra_list])
    global_median_flux = np.median(all_flux_values)
    if global_median_flux == 0:
        global_median_flux = 1.0  # avoid division by zero

    # Plot if requested
    if plot:
        plt.clf()
        fig, ax = plt.subplots(figsize=(7, 5))

        # Create a colour gradient based on MJD
        mjds = [s['mjd'] for s in spectra_list]
        cmap = plt.get_cmap('plasma')
        norm = plt.Normalize(min(mjds), max(mjds))

        for i, s in enumerate(spectra_list):
            ax.plot(s['wavelength'], s['flux']+i*global_median_flux,
                    color=cmap(norm(s['mjd'])),
                    lw=1,
                    label=f"MJD {s['mjd']:.1f}")

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label('MJD')

        ax.set_xlabel('Wavelength (Å)')
        ax.set_ylabel('Flux (arbitrary units)')
        ax.set_title(f'{sn_name} spectra ({flux_type} flux)')
        plt.tight_layout()

        # Save plot
        plot_name = f'{sn_name}_spectra_{flux_type}.pdf'
        plot_path = os.path.join(output_dir, plot_name)
        plt.savefig(plot_path, bbox_inches='tight')

    return spectra_list





