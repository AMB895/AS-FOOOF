#! /usr/bin/env python3
import fooof
import matplotlib
import matplotlib.pyplot as plt
# Import MNE, as well as the MNE sample dataset
import mne
from matplotlib.pyplot import title
from mne import io
from mne.datasets import sample
from mne.viz import plot_topomap
# import pandas
# FOOOF imports
from fooof import FOOOF
from fooof import FOOOFGroup
from fooof.sim.gen import gen_power_spectrum
from fooof.sim.utils import set_random_seed
from fooof.plts.spectra import plot_spectrum
from fooof.plts.annotate import plot_annotated_model
from fooof.bands import Bands
from fooof.analysis import get_band_peak_fg
from fooof.plts.spectra import plot_spectrum
import os.path as op
import glob
import re
import specparam
#x = glob.glob("*")
#print(x[1])
#print(f"Hello {x[1]} World")
#for a in x :
#    print(a)

from fooof import FOOOF, FOOOFGroup
from fooof.plts.spectra import plot_spectra
from fooof.plts.annotate import plot_annotated_model
from fooof.plts.aperiodic import plot_aperiodic_params
from fooof.sim.params import Stepper, param_iter
from fooof.sim import gen_power_spectrum, gen_group_power_spectra
from fooof.utils.params import compute_time_constant, compute_knee_frequency
from numpy.ma.core import size

# files = glob.glob(/Volumes/Hera/Abby/Lindsay/preprocessed_data/anti/AfterWhole/ICAwholeclean_homogenize/*.set)
# set up something where it extracts the luna id and eeg date
# Example power spectrum
#fname = '/Volumes/Hera/Abby/Lindsay/preprocessed_data/anti/AfterWhole/ICAwholeclean/10129_20180919_anti_Rem_renamedtrials_rerefwhole_ICA_icapru.set'
#raw = mne.io.read_raw_eeglab(fname, preload=True)
#raw = raw.pick_types(meg=False, eeg=True, eog=False)
#raw.plot()

#spectrum = raw.compute_psd(method='welch', fmin=1, fmax=50, tmin=0, tmax=None,
                                                      #     picks='all', n_fft=512, n_overlap=128,
                                                       #   window='hamming')
#freqs, powers = gen_power_spectrum([3,40],[1, 1],
 #                                    [[10, 0.2, 1.25], [30, 0.15, 2]])
#psds, freqs = spectrum.get_data(return_freqs=True)  # grab frequency values corresponding to spectrum
#spectra = spectrum._data  # grab spectra values
#spectraslice = spectra[1,:]
#plot_spectra(freqs,psds)
#plt.show()
# Initialize a FOOOFGroup object, with desired settings
# should either do this across people or across trials
#fg = FOOOF(peak_width_limits=[0.5, 12], min_peak_height=0,
 #                   peak_threshold=2, aperiodic_mode='fixed', max_n_peaks=4, verbose=False)

# Define the frequency range to fit
#freq_range = [1, 50]
    # Fit the power spectrum model across all channels
    # spectraSlice = spectra[i,:,:]
#fg.report(freqs, spectraslice, freq_range)
    # Check the overall results of the group fits
#fg.plot()
#plt.show()
    # # Report: fit the model, print the resulting parameters, and plot the reconstruction
#fm = fg.get_fooof(regenerate=True)
    # # # # # Print results and plot extracted model fit
#fg.print_results()
#fg.plot()
    # plt.ylim([-2.5, 2.5])
#plt.show()
    # plt.savefig(savefile)
    # plt.close()


## STRUCTURE FOR SCRIPT
# load in all AS eeg files (full runs)
files = glob.glob("/Volumes/Hera/Abby/Lindsay/preprocessed_data/anti/AfterWhole/ICAwholeclean/*.set")
for fname in files :
    # get lunaid and scandate
    iddate_pattern = r'(\d{5})_(\d{8})'
    iddate = re.search(iddate_pattern, fname)
    lunaid = iddate.group(1)
    eegdate = iddate.group(2)
    print(lunaid)
    print(eegdate)
    # load in raw data
    raw = mne.io.read_raw_eeglab(fname, preload=True)
    raw = raw.pick_types(meg=False, eeg=True, eog=False)
    # compute power spectrum
    spectrum = raw.compute_psd(method='welch', fmin=1, fmax=50, tmin=0, tmax=None,
                               picks='all', n_fft=512, n_overlap=128,
                               window='hamming')
    psds, freqs = spectrum.get_data(return_freqs=True)  # grab frequency values corresponding to spectrum
    spectra = spectrum._data  # grab spectra values
    # Loop over channels
    n_channels = spectra.shape[0]
    for c in range(n_channels) :
        # get spectra slice
        spectra_slice = spectra[c, :]
        # initiate  FOOOF object
        fg = FOOOF(peak_width_limits=[0.5, 12], min_peak_height=0,
                       peak_threshold=2, aperiodic_mode='fixed', max_n_peaks=4, verbose=False)
        # Define the frequency range to fit
        freq_range = [1, 50]
        fg.report(freqs, spectra_slice, freq_range)
        plt.show()
        # save each channel's aperiodic and periodic components then outside of channel loop save that to an entire subject's file