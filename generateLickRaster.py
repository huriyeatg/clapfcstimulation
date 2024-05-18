# This is the code for plotting daily licking pattern

# Make sure you have these libraries installed
import numpy as np
import tkinter # pip install tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import os
from datetime import datetime
import easygui


def paq_read(file_path=None, plot=False, save_path=None):
    """
    Read PAQ file (from PackIO) into python
    Lloyd Russell 2015
    Parameters
    ==========
    file_path : str, optional
        full path to file to read in. if none is supplied a load file dialog
        is opened, buggy on mac osx - Tk/matplotlib. Default: None.
    plot : bool, optional
        plot the data after reading? Default: False.
    Returns
    =======
    data : ndarray
        the data as a m-by-n array where m is the number of channels and n is
        the number of datapoints
    chan_names : list of str
        the names of the channels provided in PackIO
    hw_chans : list of str
        the hardware lines corresponding to each channel
    units : list of str
        the units of measurement for each channel
    rate : int
        the acquisition sample rate, in Hz
    """

    # file load gui
    if file_path is None:
        print('No file path')

        file_path = easygui.fileopenbox(default="*.paq", filetypes=["*.paq"])

    # open file
    fid = open(file_path, 'rb')
    print('Reading file: ' + file_path)
    # get sample rate
    rate = int(np.fromfile(fid, dtype='>f', count=1))
    # get number of channels
    num_chans = int(np.fromfile(fid, dtype='>f', count=1))
    # get channel names
    chan_names = []
    for i in range(num_chans):
        num_chars = int(np.fromfile(fid, dtype='>f', count=1))
        chan_name = ''
        for j in range(num_chars):
            chan_name = chan_name + chr(int(np.fromfile(fid, dtype='>f', count=1)))
        chan_names.append(chan_name)

    # get channel hardware lines
    hw_chans = []
    for i in range(num_chans):
        num_chars = int(np.fromfile(fid, dtype='>f', count=1))
        hw_chan = ''
        for j in range(num_chars):
            hw_chan = hw_chan + chr(int(np.fromfile(fid, dtype='>f', count=1)))
        hw_chans.append(hw_chan)

    # get acquisition units
    units = []
    for i in range(num_chans):
        num_chars = int(np.fromfile(fid, dtype='>f', count=1))
        unit = ''
        for j in range(num_chars):
            unit = unit + chr(int(np.fromfile(fid, dtype='>f', count=1)))
        units.append(unit)

    # get data
    temp_data = np.fromfile(fid, dtype='>f', count=-1)
    num_datapoints = int(len(temp_data)/num_chans)
    data = np.reshape(temp_data, [num_datapoints, num_chans]).transpose()

    # close file
    fid.close()

    return {"data": data,
            "chan_names": chan_names,
            "hw_chans": hw_chans,
            "units": units,
            "rate": rate}

# Functions for reading in data from .paq files
def paq_data(paq, chan_name, threshold=1, threshold_ttl = False, plot=False):
    '''
    Do not include any exclusion of data in this function
    returns the data in paq (from paq_read) from channel: chan_names
    if threshold_tll: returns sample that trigger occured on
    '''
    chan_idx = paq['chan_names'].index(chan_name)
    data = paq['data'][chan_idx, :]
    if threshold_ttl == False:
        data = data
    elif threshold_ttl == 'Mix':
        data = threshold_detect(data,threshold,cutoff = 5)
    elif threshold_ttl == 'Lick':
        threshold = 4.5# to clean the reward signal from licking - there is some cross talk across channels
        data = threshold_detect(data,threshold,cutoff = 5)
    else:
        data = threshold_detect(data,threshold, cutoff=False)

    if plot:     
        if threshold_ttl:
            plt.plot(data, np.ones(len(data)), '.')
        else:
            plt.plot(data)
    return data

def threshold_detect(signal, threshold, cutoff=True):
    '''lloyd russell, cutoff is added by HA'''
    if cutoff:
        thresh_signal = (signal > threshold) & (signal < cutoff)
        thresh_signal[1:][thresh_signal[:-1] & thresh_signal[1:]] = False
        times = np.where(thresh_signal)
    else:
        thresh_signal = signal > threshold
        thresh_signal[1:][thresh_signal[:-1] & thresh_signal[1:]] = False
        times = np.where(thresh_signal)

    return times[0]

def lick_binner(paqData,trial_start, stChanName, stimulation=True ):
    ''' makes new easytest binned lick variable in run object '''
    licks = paq_data (paqData, stChanName, 1, threshold_ttl='Lick')

    binned_licks = []

    for i, t_start in enumerate(trial_start):
        if i == len(trial_start) - 1:
            t_end = np.inf
        else:
            t_end = trial_start[i+1]

        trial_idx = np.where((licks >= t_start) & (licks <= t_end))[0]

        trial_licks = licks[trial_idx] - t_start

        binned_licks.append(trial_licks)

    licks = licks
    # attribute already exists called 'binned_licks' and cannot overwrite it
    binned_licks_easytest = binned_licks

    return licks, binned_licks


##### MAIN FUNCTION #####

# Get the paq file 
file_path = None # 'Z:\\Data\\2023-08-27\\2023-08-27_2303_trainingDay-001.paq'
paqData = paq_read( file_path=file_path, plot=True, save_path=None)
fRate = 20000
ymax = 100
stChanName_lick = 'lickDetection'
stChanName_reward = 'waterDelivery' #'reward'
savepathname = 'C:\\Users\\Huriye\\'
    
# Create behaviour sessions plots

# Get the stim start times 
trialStartTimes = paq_data (paqData, stChanName_reward, 1, threshold_ttl=True)
trialStartTimes = trialStartTimes - (2*fRate)

# Get the lick times
licks, animal_lick = lick_binner(paqData, trialStartTimes,stChanName_lick, stimulation = False)

# Plot the figure
fig = plt.figure(figsize=(10, 5))
for i, array in enumerate(animal_lick):
    plt.plot(array, np.ones_like(array)+i, 'k.',markersize = 1)
             
plt.xlim(0, 8*fRate)
plt. ylim(0, ymax)
plt.yticks(range(0,ymax, 150), range(ymax,0, -150))
plt.xticks (range(0,(8*fRate)+1,fRate), range(-2,7,1))
plt.ylabel('Trials')
plt.xlabel('Time (sec)') 
savefilename = savepathname + datetime.now().strftime('%Y-%m-%d') + '.png'
plt.savefig(savefilename, bbox_inches='tight', transparent=False, dpi = 300)
plt.close(fig)