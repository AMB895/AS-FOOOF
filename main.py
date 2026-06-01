#!/usr/bin/env python3
import fooof
import matplotlib
import matplotlib.pyplot as plt
import mne
from matplotlib.pyplot import title
from mne import io
from mne.datasets import sample
from mne.viz import plot_topomap
from fooof import FOOOF, fit_fooof_3d
import scipy.io
import os
from fooof import FOOOFGroup
from fooof.sim.gen import gen_power_spectrum
from fooof.sim.utils import set_random_seed
from fooof.plts.spectra import plot_spectrum
from fooof.plts.annotate import plot_annotated_model
from fooof.bands import Bands
from fooof.analysis import get_band_peak_fg
from fooof.plts.spectra import plot_spectrum
import os.path as op
import glob2
import re
import specparam
from fooof import FOOOF, FOOOFGroup
from fooof.plts.spectra import plot_spectra
from fooof.plts.annotate import plot_annotated_model
from fooof.plts.aperiodic import plot_aperiodic_params
from fooof.sim.params import Stepper, param_iter
from fooof.sim import gen_power_spectrum, gen_group_power_spectra
from fooof.utils.params import compute_time_constant, compute_knee_frequency
from numpy.ma.core import size
import numpy as np
import pandas as pd



print(fooof.__version__)
files = glob2.glob('/Volumes/Hera/Abby/Lindsay/preprocessed_data/anti/AfterWhole/kept_epoch/*.set', recursive=True)
print(len(files))
for fname in files:
   # extract id and date from file name   
    pattern = r'(\d{5})_(\d{8})'
    match = re.search(pattern, fname)
    id = match.group(1)
    date = match.group(2)
    filename = f"{match}"
    print(id,date)

    # final file path and name
    save_path = '/Volumes/Hera/Abby/FOOOF/SubjectFiles/'
    sub_name = "{id}_{date}.mat".format(id=id, date=date)
    full_file = save_path + sub_name
    print(full_file)

    # skip if ran subject
    if os.path.exists(full_file):
        print(f"skipping {id} {date} - output file already exists")
        continue

    # catch files that fail to load (e.g. fewer than 2 epochs)
    try:
        raw = mne.io.read_epochs_eeglab(fname)
    except Exception as e:
        print(f"Skipping {id} {date} - failed to load: {e}")
        continue

    raw = raw.pick_types(meg=False, eeg=True, eog=False)
   
    pattern = r'(\d{5})_(\d{8})'
    match = re.search(pattern, fname)
    id = match.group(1)
    date = match.group(2)
    filename = f"{match}"
    print(id,date)
    # save all data frames
    save_path = '/Volumes/Hera/Abby/FOOOF/SubjectFiles/'
    sub_name = "{id}_{date}.mat".format(id=id, date=date)
    full_file = save_path + sub_name
    print(full_file)

    if os.path.exists(full_file):
        print(f"skipping {id} {date} - output file already exists")
        continue    # raw.plot()
    # plt.show()
    spectrum = raw.compute_psd(method='welch', fmin=1, fmax=60, tmin=0, tmax=None,
                               picks='all', n_fft=512, n_overlap=128, window='hamming')
    psds, freqs = spectrum.get_data(return_freqs=True)  # grab frequency values corresponding to spectrum
    spectra = spectrum._data  # grab spectra values
    #print(*spectra.shape)
    fm = FOOOFGroup(peak_width_limits=[0.5, 12], max_n_peaks=4, min_peak_height=0, verbose=False,
                    aperiodic_mode='fixed')
    # fm.fit(freqs, spectra)
    fms = fit_fooof_3d(fm,freqs,spectra)
    ntrials = len(fms)
    #print(len(fms))

    # aperiodic parameters
    ap_params = fm.get_params('aperiodic_params')
    flat_idx = np.arange(len(ap_params))
    epochs_cols = (flat_idx // 64).reshape(-1, 1)
    channels_cols = (flat_idx % 64).reshape(-1, 1)
    ap_params_full = np.hstack([ap_params, epochs_cols, channels_cols])
    ap_params_df = pd.DataFrame(ap_params_full, columns=['offset','exponent','epoch','channel'])
    #print(ap_params_df)

    # periodic parameters
    peak_params = fm.get_params('peak_params')
    flat_idx = peak_params[:,-1].astype(int)
    epochs_cols = (flat_idx // 64).reshape(-1,1)
    channels_cols = (flat_idx % 64).reshape(-1,1)
    peak_params_full = np.hstack([peak_params,epochs_cols,channels_cols])
    peak_params_df = pd.DataFrame(peak_params_full, columns=['CF','PW','BW','spectrum_idx','epoch','channel'])
   # print(peak_params_df)

    # gaussian parameters
    gauss_params = fm.get_params('gaussian_params')
    flat_idx = gauss_params[:,-1].astype(int)
    epochs_cols = (flat_idx // 64).reshape(-1,1)
    channels_cols = (flat_idx % 64).reshape(-1,1)
    gauss_params_full = np.hstack([gauss_params,epochs_cols,channels_cols])
    gauss_params_df = pd.DataFrame(gauss_params_full, columns=['mean', 'height', 'SD','spectrum_idx','epoch','channel'])
    #print(gauss_params_df)

    # fit params
    errors = fm.get_params('error')
    r2s = fm.get_params('r_squared')
    flat_idx = np.arange(len(r2s))
    epochs_cols = (flat_idx // 64)
    channels_cols = (flat_idx % 64)
    fit_df = pd.DataFrame({
        'r_squared': r2s,
        'error': errors,
        'epoch': epochs_cols,
        'channel': channels_cols
    })

   
    scipy.io.savemat(full_file,{
        'ap_params': ap_params_df.values,
        'peak_params': peak_params_df.values,
        'gauss_params': gauss_params_df.values,
        'fit_params': fit_df.values,
        'ap_cols': ap_params_df.columns.tolist(),
        'peak_cols': peak_params_df.columns.tolist(),
        'gauss_cols': gauss_params_df.columns.tolist(),
        'fit_cols': fit_df.columns.tolist()
    })


