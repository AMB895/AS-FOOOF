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
import random
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
# files = glob.glob('/Volumes/Hera/Abby/Lindsay/preprocessed_data/anti/AfterWhole/kept_epoch/*kept_prep.set')
task = input("Enter task: ")
print(f"running fooof on {task}")

# files = glob.glob('/Volumes/Hera/Abby/Lindsay/auto/anti/AfterWhole/kept_epoch/*kept_2_scored.set')
files = glob.glob(f"/Volumes/Hera/Abby/Lindsay/auto/{task}/AfterWhole/kept_epoch/10129_*kept_2_scored.set")

print(len(files))
for fname in files:

   # extract id and date from file name   
    pattern = r'(\d{5})_(\d{8})'
    match = re.search(pattern, fname)
    lunaid = match.group(1)
    date = match.group(2)
    filename = f"{match}"
    print(lunaid, date)

    # final file path and name
    save_path = f"/Volumes/Hera/Abby/FOOOF/{task}/SubjectFiles/"
    sub_name = "{id}_{date}_{task}.mat".format(id=lunaid, date=date,task=task)
    full_file = save_path + sub_name
    # print(full_file)

    # skip if ran subject
    if os.path.exists(full_file):
        print(f"skipping {lunaid} {date} {task} - output file already exists")
        continue

    # catch files that fail to load (e.g. fewer than 2 epochs)
    try:
        raw = mne.io.read_epochs_eeglab(fname)
    except Exception as e:
        print(f"Skipping {lunaid} {date} {task} - failed to load: {e}")
        continue

    raw = raw.pick_types(meg=False, eeg=True, eog=False) # data type is EEG
    # get trial score and rename
    eventid_to_name = {v: k for k, v in raw.event_id.items()}
    epoch_labels = [eventid_to_name[e] for e in raw.events[:,2]]
    # print(raw.event_id)
    # print(epoch_labels)
    prep_epoch_labels = [re.search(r'\b(19|20|21|22)\b',e).group() for e in epoch_labels]
    # print(prep_epoch_labels)
    label_map = {'19':-1, '20':0, '21':1, '22':2}
    # data frame that stores epoch number and score for that trial
    prep_epoch_df = pd.DataFrame({
        'epoch': range(len(prep_epoch_labels)),
        'label': [label_map[e] for e in prep_epoch_labels]
    })

    spectrum = raw.compute_psd(method='welch', fmin=1, fmax=60, tmin=0, tmax=None,
                               picks='all', n_fft=512, n_overlap=128, window='hamming')
    psds, freqs = spectrum.get_data(return_freqs=True)  # grab frequency and psds values corresponding to spectrum
    # check the shape of psds
    print(psds.shape)
    # MNE gives PSDs in V^2/Hz; convert V^2/Hz to uV^2/Hz
    psds_uv = psds * 1e12

    # define FOOOF group parameters
    fm = FOOOFGroup(peak_width_limits=[2, 12], max_n_peaks=4, min_peak_height=0, verbose=False, aperiodic_mode='fixed')
    # ^^ peak_width_limits should be set to twice the frequency resolution so single points are not fit as peaks
    # 3d fooof fit from FOOOFGroup
    fms = fit_fooof_3d(fm, freqs, psds_uv)

    # aperiodic parameters for each trial and channel
    ap_3d = np.array([fg.get_params('aperiodic_params') for fg in fms])
    print(ap_3d.shape)
    n_epochs, n_channels, _ = ap_3d.shape
    epoch_idx, chan_idx = np.meshgrid(np.arange(n_epochs), np.arange(n_channels), indexing = 'ij')
    ap_params_df = pd.DataFrame({
        'offset': ap_3d[:, :, 0].ravel(),
        'exponent': ap_3d[:, :, 1].ravel(),
        'epoch': epoch_idx.ravel(),
        'channel': chan_idx.ravel()
    })

    # periodic parameters for each trial and channel
    peak_rows = []
    for i, fg in enumerate(fms):
        pp = fg.get_params('peak_params')
        if pp is not None and len(pp) > 0:
            for row in pp:
                peak_rows.append([row[0], row[1], row[2], row[3], i, int(row[3]) % n_channels])

    peak_params_df = pd.DataFrame(peak_rows, columns=['CF', 'PW', 'BW', 'spectrum_idx', 'epoch', 'channel'])

   # gaussian parameters for each trial and channel
    gauss_rows = []
    for i, fg in enumerate(fms):
        gp = fg.get_params('gaussian_params')
        if gp is not None and len(gp) > 0:
            for row in gp:
                gauss_rows.append([row[0], row[1], row[2], row[3], i, int(row[3]) % n_channels])

    gauss_params_df = pd.DataFrame(gauss_rows, columns=['mean', 'height', 'SD', 'spectrum_idx', 'epoch', 'channel'])

   # fit paramters for each trial and channel
    r2s = np.array([fg.get_params('r_squared') for fg in fms])  # (n_epochs, n_channels)
    errors = np.array([fg.get_params('error') for fg in fms])
    fit_df = pd.DataFrame({
        'r_squared': r2s.ravel(),
        'error': errors.ravel(),
        'epoch': epoch_idx.ravel(),
        'channel': chan_idx.ravel()
    })

    scipy.io.savemat(full_file,{
        'ap_params': ap_params_df.values,
        'peak_params': peak_params_df.values,
        'gauss_params': gauss_params_df.values,
        'fit_params': fit_df.values,
        'trial_score': prep_epoch_df.values,
        'ap_cols': ap_params_df.columns.tolist(),
        'peak_cols': peak_params_df.columns.tolist(),
        'gauss_cols': gauss_params_df.columns.tolist(),
        'fit_cols': fit_df.columns.tolist(),
        'trial_score_cols': prep_epoch_df.columns.tolist()
    })

    # plot random combination of fits and save
    combination = random.sample([(i, j) for i in range(n_epochs) for j in range(n_channels)], 3)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, (epoch, chan) in zip(axes, combination):
        r2 = r2s[epoch, chan]
        err = errors[epoch, chan]
        fms[epoch].get_fooof(ind=chan).plot(ax=ax)
        ax.set_title(f'Epoch {epoch}, Chan {chan} | R^2 = {r2:.3f}, error = {err:.3f}')

    plt.tight_layout()
    plt.savefig(f'/Volumes/Hera/Abby/FOOOF/{task}/SubjectFiles/fits_{lunaid}_{date}.png')
    plt.close()

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))



