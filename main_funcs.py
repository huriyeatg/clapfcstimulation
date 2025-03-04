# Functions used more globally

import platform
import os
import pandas as pd
import numpy as np
import tifffile
import utils_funcs as utils
import plot_funcs as pfun
import paq2py
import re
import pickle
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import entropy
from scipy import stats
from sklearn.feature_selection import mutual_info_classif
from sklearn.metrics import mutual_info_score
from sklearn.preprocessing import KBinsDiscretizer
if os.environ['CONDA_DEFAULT_ENV'] == 'clapfcstimulation':
    print('Env: ' + os.environ['CONDA_DEFAULT_ENV'])
    import pims
    import cv2

class analysis:
    
    def __init__(self):
        
        if platform.platform() == 'macOS-12.6-arm64-arm-64bit':
            # for Windows - Huriye PC
            print("Computer: Huriye MAC")
            self.recordingListPath = "/Users/Huriye/Documents/Code/clapfcstimulation/"
            self.rawPath           = "Z:\\Data\\"
            self.rootPath          = "/Users/Huriye/Documents/Code/clapfcstimulation/"
        elif platform.platform() == 'Windows-10-10.0.19045-SP0':
            print("Computer: Huriye Windows")
            # for Windows - Huriye PC
            self.recordingListPath = "C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\"
            self.rawPath           = "Z:\\Data\\"
            self.rootPath          = "C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\"
        else:
            print("Problem detected")
            print(platform.platform())

        self.analysisPath =  os.path.join(self.rootPath, 'analysis')
        self.figsPath     =  os.path.join(self.rootPath, 'figs')
        self.DLCconfigPath = os.path.join(self.rootPath, 'pupil-extraction', 'pupil-Eren CAN-2021-11-21')
        self.DLCconfigPath = self.DLCconfigPath + '\\config.yaml'
        
        self.recordingList = pd.read_csv(self.recordingListPath + "clapfcstimulation - imaging.csv")#error_bad_lines=False
        # Add main filepathname
        self.recordingList['analysispathname'] = np.nan
        for ind, recordingDate in enumerate(self.recordingList.recordingDate):
            filepathname = (self.rawPath + self.recordingList.recordingDate[ind] +
                            '\\'+ self.recordingList.recordingDate[ind]+ '_'+
                            str(self.recordingList.animalID[ind]) )
            self.recordingList.loc[ind,'filepathname'] = filepathname  

            analysispathname = (self.analysisPath +
                                '\\' + self.recordingList.recordingDate[ind] + '_' + 
                                str(self.recordingList.animalID[ind]) + '_' + 
                                format(int(self.recordingList.recordingID[ind]), '03d'))   
            self.recordingList.loc[ind,'analysispathname'] = analysispathname +'\\' 

def update_dff_traces(dff_traces, new_data, key, shape):
    """Update dF/F trace arrays with new data or placeholder NaN arrays."""
    if key in new_data:
        new_trace = new_data[key]
    else:
        new_trace = np.full(shape, np.nan)  # Placeholder for missing data
    
    return np.vstack((dff_traces, new_trace)) if dff_traces is not None else new_trace

def update_pupil_traces(extData, len_cellID, indPupil, pupilID, pupil_traceBoth, pupil_traceVis, pupil_traceOpto):
    """Organizes pupil traces from extData based on quality and recording order.
    
    Args:
        extData (dict): Extracted data containing pupil traces and quality.
        len_cellID (int): The number of cells in the current recording.
        indPupil (int): Indicator if it's the first recording (0) or subsequent ones (1+).
        pupilID (list): List to extend with pupil trace IDs based on quality.
        pupil_traceBoth, pupil_traceVis, pupil_traceOpto (np.ndarray): Existing pupil trace arrays to update.
        
    Returns:
        tuple: Updated indPupil, pupilID, pupil_traceBoth, pupil_traceVis, pupil_traceOpto.
    """
    pupil_trace = extData[17]
    pupilQuality = extData[18]
    trace_shapes = (len_cellID, 240)  # Assuming fixed shape for simplicity

    def get_or_fake(trace_key):
        if trace_key in pupil_trace and pupilQuality > 0.75:
            return np.tile(pupil_trace[trace_key], (len_cellID, 1))
        return np.full(trace_shapes, np.nan)

    if indPupil == 0 and (len(pupil_trace) > 0) & (pupilQuality > 0.75):
        indPupil = 1  # Update for subsequent recordings
    pupilID.extend([int((len(pupil_trace) > 0) & (pupilQuality > 0.75))] * len_cellID)

    # Fetch or fake traces
    new_traceBoth = get_or_fake('Both')
    new_traceVis = get_or_fake('onlyVis')
    new_traceOpto = get_or_fake('onlyOpto')

    # Update existing traces or set new ones
    pupil_traceBoth = new_traceBoth if indPupil == 0 else np.vstack((pupil_traceBoth, new_traceBoth))
    pupil_traceVis = new_traceVis if indPupil == 0 else np.vstack((pupil_traceVis, new_traceVis))
    pupil_traceOpto = new_traceOpto if indPupil == 0 else np.vstack((pupil_traceOpto, new_traceOpto))

    return indPupil, pupilID, pupil_traceBoth, pupil_traceVis, pupil_traceOpto

def normalize_to_baseline(dff_trace, pre_frames):
    ac = np.mean(dff_trace[:, :pre_frames], axis=1, keepdims=True)
    return dff_trace - ac

def convert_tiff2avi (imagename, outputsavename, fps=30.0):
    # path = 'Z:\Data\\2022-05-09\\2022-05-09_22107_p-001\\'
    # filename = path + 'test.tif'
    # outputsavename = 'Z:\\Data\\2022-05-09\\2022-05-09_22107_p-001\\output.avi'
    # fps = 30.0

     #Load the stack tiff & set the params
     imageList = pims.TiffStack(imagename)
     nframe, height, width = imageList.shape
     size = width,height
     duration = imageList.shape[0]/fps
     # Create the video & save each frame
     fourcc = cv2.VideoWriter_fourcc(*"MJPG")
     out = cv2.VideoWriter(outputsavename, fourcc, fps,size,0)
     for frame in imageList:
        # frame = cv2.resize(frame, (500,500))
         out.write(frame)
     out.release()
     
def get_file_names_with_strings(pathIn, str_list):
    full_list = os.listdir(pathIn)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]

    return final_list

def fdr(p_vals):
    #http://www.biostathandbook.com/multiplecomparisons.html
    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1

    return fdr

def calculateDFF (tiff_folderpath, frameClockfromPAQ):
    s2p_path = tiff_folderpath +'\\suite2p\\plane0\\'
    # from Vape - catcher file: 
    flu_raw, _, _ = utils.s2p_loader(s2p_path, subtract_neuropil=False) 

    flu_raw_subtracted, spks, stat = utils.s2p_loader(s2p_path)
    flu = utils.dfof2(flu_raw_subtracted)

   # _, n_frames = tiff_metadata(tiff_folderpath)
    tseries_lens = [len(flu)] #n_frames

    # deal with the extra frames 
    frameClockfromPAQ = frameClockfromPAQ[:tseries_lens[0]] # get rid of foxy bonus frames

    # correspond to analysed tseries
    paqio_frames = utils.tseries_finder(tseries_lens, frameClockfromPAQ, paq_rate=20000)
    paqio_frames = paqio_frames

    if len(paqio_frames) == sum(tseries_lens):
        print('Dff extraction is completed: ' +tiff_folderpath)
        imagingDataQuality = True
       # print('All tseries chunks found in frame clock')
    else:
        imagingDataQuality = False
        print('WARNING: Could not find all tseries chunks in '
              'frame clock, check this')
        print('Total number of frames detected in clock is {}'
               .format(len(paqio_frames)))
        print('These are the lengths of the tseries from '
               'spreadsheet {}'.format(tseries_lens))
        print('The total length of the tseries spreasheets is {}'
               .format(sum(tseries_lens)))
        missing_frames = sum(tseries_lens) - len(paqio_frames)
        print('The missing chunk is {} long'.format(missing_frames))
        try:
            print('A single tseries in the spreadsheet list is '
                  'missing, number {}'.format(tseries_lens.index
                                             (missing_frames) + 1))
        except ValueError:
            print('Missing chunk cannot be attributed to a single '
                   'tseries')
    return {"imagingDataQuality": imagingDataQuality,
            "frame-clock": frameClockfromPAQ,
            "paqio_frames":paqio_frames,
            "n_frames":tseries_lens,
            "flu": flu,
            "spks": spks,
            "stat": stat,
            "flu_raw": flu_raw}

def calculatePupil (filename, frameClockfromPAQ):
    dataPupilCSV = pd.read_csv(filename [0], header = 1)
    dataPupilCSV.head()

    verticalTop_x =np.array(dataPupilCSV['Xmax'] [1:], dtype = float)
    verticalTop_y =np.array(dataPupilCSV['Xmax.1'][1:], dtype = float)
    verticalBottom_x =np.array(dataPupilCSV['Xmin'][1:], dtype = float)
    verticalBottom_y =np.array(dataPupilCSV['Xmin.1'][1:], dtype = float)
    verticallikelihood =np.mean(np.array(dataPupilCSV['Xmax.2'][1:], dtype = float))

    verticalDis = np.array(np.sqrt((verticalTop_x - verticalBottom_x)**2 + (verticalTop_y - verticalBottom_y)**2))

    horizontalTop_x =np.array(dataPupilCSV['Ymax'] [1:], dtype = float)
    horizontalTop_y =np.array(dataPupilCSV['Ymax.1'][1:], dtype = float)
    horizontalBottom_x =np.array(dataPupilCSV['Ymin'][1:], dtype = float)
    horizontalBottom_y =np.array(dataPupilCSV['Ymin.1'][1:], dtype = float)
    horizontallikelihood =np.mean(np.array(dataPupilCSV['Ymax.2'][1:], dtype = float))

    horizontalDis = np.array(np.sqrt((horizontalTop_x - horizontalBottom_x)**2 + (horizontalTop_y - horizontalBottom_y)**2))

    lengthCheck = len(frameClockfromPAQ)==len(horizontalDis)

    return {"verticalTop_x": verticalTop_x,
            "verticalTop_y": verticalTop_y,
            "verticalBottom_x": verticalBottom_x,
            "verticalBottom_y": verticalBottom_y,
            "verticallikelihood": verticallikelihood,
            "verticalDis": verticalDis,
            "horizontalTop_x": horizontalTop_x,
            "horizontalTop_y": horizontalTop_y,
            "horizontalBottom_x": horizontalBottom_x,
            "horizontalBottom_y": horizontalBottom_y,
            "horizontallikelihood": horizontallikelihood,
            "horizontalDis": horizontalDis,
            "lengthCheck": lengthCheck,
            "frameClockfromPAQ": frameClockfromPAQ}

def tiff_metadata(folderTIFF):

    ''' takes input of list of tiff folders and returns 
        number of frames in the first of each tiff folder '''
    
    # First check if tiff file is good and correct
    tiff_list = []
    tseries_nframes = []
    tiffs = utils.get_tiffs(folderTIFF)
    if not tiffs:
        raise print('cannot find tiff in '
                                    'folder {}'.format(tseries_nframes))
    elif len(tiffs) == 1:
        assert tiffs[0][-7:] == 'Ch3.tif', 'channel not understood '\
                                            'for tiff {}'.format(tiffs)
        tiff_list.append(tiffs[0])
    elif len(tiffs) == 2:  # two channels recorded (red is too dim)
        print('There are more than one tiff file - check: '+ folderTIFF)

    with tifffile.TiffFile(tiffs[0]) as tif:
        tif_tags = {}
        for tag in tif.pages[0].tags.values():
            name, value = tag.name, tag.value
            tif_tags[name] = value

    x_px = tif_tags['ImageWidth']
    y_px = tif_tags['ImageLength']
    image_dims = [x_px, y_px]

    n_frames = re.search('(?<=\[)(.*?)(?=\,)', 
                         tif_tags['ImageDescription'])

    n_frames = int(n_frames.group(0))
    tseries_nframes.append(n_frames)

    return image_dims, tseries_nframes

def getIndexForInterestedcellsID ( s_recDate, s_animalID, s_recID, s_cellID ):
    infoPath = 'C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\analysis\\infoForAnalysis-readyForSelectingInterestedCells.pkl'    
    animalID, stimuliFamilarity, dataQuality,recData, recID, cellID, pvalsBoth, pvalsVis, pvalsOpto,dff_meanVisValue, dff_meanBothValue, dff_meanOptoValue, pupilID = pd.read_pickle(infoPath) 
    ind = np.where((np.array(animalID) == s_animalID) & (np.array(recID) == s_recID) & (np.array(cellID) == s_cellID) & (np.array(recData) == s_recDate))
    return ind

def selectInterestedcells ( aGroup, stimType, responsive = True, plotValues = False, pupil = True ):
    infoPath = 'C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\analysis\\infoForAnalysis-readyForSelectingInterestedCells.pkl'    
    animalID, stimuliFamilarity, dataQuality,recData, recID, cellID, pvalsBoth, pvalsVis, pvalsOpto, dff_meanVisValue, dff_meanBothValue, dff_meanOptoValue, pupilID = pd.read_pickle(infoPath) 
    animalID = np.where(animalID == 23040, 2304, animalID) # typo in animal ID
    #animalID = np.where(animalID == 22104, 22105, animalID) # typo in animal ID ( 22104 & 22105 is same animal)
    infoPath = 'C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\analysis\\infoForAnalysis-readyForPlotting_normalisedtoPre.pkl'
    dff_traceBoth, dff_traceVis, dff_traceOpto = pd.read_pickle(infoPath) 
    Chrimson_animals = [ 21104, 21107, 21108, 21109,22101,22102,22103,22105,22106,22107,22108, 2303, 2304]
    NAAP_animals    = [21101, 21102, 21103, 21105, 21106]  
    control_animals = [23040, 23036, 23037]
    OPN_animals     = [2306, 2307,2308,2309,2310, 2311]
    if stimType == 'Trained':
        # exclude 22102 and 22108 as they did not learn
        Chrimson_animals = [22101,22102, 22103,22105,22106,22107, 2303, 2304]
        OPN_animals      = [2306, 2307,2308,2309,2310, 2311] # 2312

    if aGroup == 'Chrimson':
        s = set(Chrimson_animals)
        selectedAnimals = np.array([i in s for i in animalID])
    elif aGroup == 'NAAP':
        s = set(NAAP_animals)
        selectedAnimals = np.array([i in s for i in animalID])
    elif aGroup == 'Control': #pupil only
        s = set(control_animals)
        selectedAnimals = np.array([i in s for i in animalID])
    elif aGroup == 'OPN3':
        s = set(OPN_animals)
        selectedAnimals = np.array([i in s for i in animalID])
    else:
        s =set([aGroup])
        selectedAnimals = np.array([i in s for i in animalID])

    # exclude trained stimuli
    if stimType == 'Trained':
        includeType = [2,9]
        s = set(includeType)
        selectedFamilarity = np.array([i in s for i in stimuliFamilarity])
    elif stimType  == 'Naive':
        includeType = [0, 1, 8]
        s = set(includeType)
        selectedFamilarity = np.array([i in s for i in stimuliFamilarity])
    elif stimType  == 'Pupil-control-coveredMicroscope':
        includeType = [6]
        s = set(includeType)
        selectedFamilarity = np.array([i in s for i in stimuliFamilarity])
    elif stimType  == 'Pupil-control-not-coveredMicroscope':
        includeType = [7]
        s = set(includeType)
        selectedFamilarity = np.array([i in s for i in stimuliFamilarity])
    else:
        print(stimType + ' is not defined')

    # select only pupil
    if pupil:
        includeType = [1]
        s = set(includeType)
        selectedPupil = np.array([i in s for i in pupilID])
    else:
        includeType = [0, 1]
        s = set(includeType)
        selectedPupil = np.array([i in s for i in pupilID])
    
    selectedExpGroup = selectedAnimals & selectedFamilarity & selectedPupil 

    if responsive==False:
        selectedCellIndex = selectedExpGroup
    else:

        temp = fdrcorrection(pvalsVis, alpha=0.05, method='n', is_sorted=False)
        responsiveVis  = temp[0]
        responsiveNoVis  = ~responsiveVis
    
        temp = fdrcorrection(pvalsOpto, alpha=0.05, method='n', is_sorted=False)
        responsiveOpto = temp[0]
        responsiveNoOpto = ~responsiveOpto
        
        temp = fdrcorrection(pvalsBoth, alpha=0.05, method='n', is_sorted=False)
        responsiveBoth = temp[0]
        responsiveNoBoth  = ~responsiveBoth

        responsiveVis  = selectedExpGroup & responsiveVis   # np.logical_and(selectedExpGroup,responsiveVis)
        responsiveOpto = selectedExpGroup & responsiveOpto  # np.logical_and(selectedExpGroup,responsiveOpto)
        responsiveBoth = selectedExpGroup & responsiveBoth  # np.logical_and(selectedExpGroup,responsiveBoth) 

        responsiveNoVis  = selectedExpGroup & responsiveNoVis   # np.logical_and(selectedExpGroup,responsiveVis)
        responsiveNoOpto = selectedExpGroup & responsiveNoOpto  # np.logical_and(selectedExpGroup,responsiveOpto)
        responsiveNoBoth = selectedExpGroup & responsiveNoBoth  # np.logical_and(selectedExpGroup,responsiveBoth) 

        responsiveAll = responsiveVis | responsiveOpto | responsiveBoth
        nonResponsiveAll = responsiveNoVis & responsiveNoOpto & responsiveNoBoth

        responsiveOnlyVis   = responsiveVis & responsiveNoOpto & responsiveNoBoth
        responsiveOnlyOpto   = responsiveOpto & responsiveNoVis & responsiveNoBoth
        responsiveOnlyBoth  = responsiveBoth & responsiveNoVis & responsiveNoOpto

        if plotValues:
            print('All cell number:'+ str(np.sum(selectedExpGroup)))
            # # All responsive cells
            responsiveAll = responsiveVis | responsiveOpto | responsiveBoth
            print('Any responsive cell number:'+ str(np.sum(responsiveAll)))
            VisOptoNoBoth = responsiveVis & responsiveOpto & responsiveNoBoth
            print('Responsice to  visual and opto but not both cell number:'+ str(np.sum(VisOptoNoBoth)))

            # visual cue responsive cells
            excDff = (np.array(dff_meanVisValue)> 0)
            inhDff = (np.array(dff_meanVisValue)<0)
            excOnlyV = excDff & responsiveVis
            inhOnlyV = inhDff & responsiveVis
            
            print('Visual cue - all visual responsive cells: '+ str(np.sum(responsiveVis)))
            print('Visual cue - only visual responsive: '+ str(np.sum(responsiveOnlyVis)))
            print('Visual only cue - EXC opto responsive: '+ str(np.sum(excOnlyV)/np.sum(responsiveVis)))
            print('Visual only cue - INH opto responsive: '+ str(np.sum(inhOnlyV)/np.sum(responsiveVis)))

            # opto cue responsive cells
            excDff = (np.array(dff_meanBothValue) > 0)
            inhDff = (np.array(dff_meanBothValue) < 0)
            excOnlyO = excDff & responsiveOpto
            inhOnlyO = inhDff & responsiveOpto
            print('Opto stimulation - all opto responsive cells: '+ str(np.sum(responsiveOpto)))
            print('Opto stimulation - only opto responsive: '+ str(np.sum(responsiveOnlyOpto)))
            print('Opto stimulation - EXC opto responsive: '+ str(np.sum(excOnlyO)/np.sum(responsiveOpto)))
            print('Opto stimulation - INH opto responsive: '+ str(np.sum(inhOnlyO)/np.sum(responsiveOpto)))

            # Both cue responsive cells
           
            responsiveOnlyBoth  = responsiveBoth & responsiveNoVis & responsiveNoOpto   
            excDff = (np.array(dff_meanBothValue) > 0)
            inhDff = (np.array(dff_meanBothValue) < 0)
            excOnlyB = excDff & responsiveOnlyBoth
            inhOnlyB = inhDff & responsiveOnlyBoth
            print('Both - all both responsive cells:'+ str(np.sum(responsiveBoth)))
            print('Both - only both responsive: '+ str(np.sum(responsiveOnlyBoth)))
            print('Both - EXC opto responsive: '+ str(np.sum(excOnlyB)/np.sum(responsiveOnlyBoth)))
            print('Both - INH opto responsive: '+ str(np.sum(inhOnlyB)/np.sum(responsiveOnlyBoth)))

        else:
            # visual cue responsive cells
            excDff = (np.array(dff_meanVisValue)> 0)
            inhDff = (np.array(dff_meanVisValue)<0)
            excOnlyV = excDff & responsiveVis
            inhOnlyV = inhDff & responsiveVis
            # opto cue responsive cells
            excDff = (np.array(dff_meanBothValue) > 0)
            inhDff = (np.array(dff_meanBothValue) < 0)
            # Both cue responsive cells
            excOnlyO = excDff & responsiveOpto
            inhOnlyO = inhDff & responsiveOpto
            responsiveOnlyBoth  = responsiveBoth & responsiveNoVis & responsiveNoOpto   
            excDff = (np.array(dff_meanBothValue) > 0)
            inhDff = (np.array(dff_meanBothValue) < 0)
            excOnlyB = excDff & responsiveOnlyBoth
            inhOnlyB = inhDff & responsiveOnlyBoth

    if responsive =='Sensory':
        selectedCellIndex = responsiveVis
        excOnly = excOnlyV
        inhOnly = inhOnlyV
        
    elif responsive =='Opto':
        selectedCellIndex = responsiveOpto
        excOnly = excOnlyO
        inhOnly = inhOnlyO

    elif responsive =='Opto-boosted':
        selectedCellIndex = responsiveOnlyBoth
        excOnly = excOnlyB
        inhOnly = inhOnlyB

    elif responsive =='All':
        selectedCellIndex = responsiveAll
        excOnly = responsiveAll
        inhOnly = responsiveAll

    elif responsive=='None':
        selectedCellIndex = nonResponsiveAll
        excOnly = nonResponsiveAll
        inhOnly = nonResponsiveAll
    else:
        selectedCellIndex = 'Select responsive type'

    if plotValues:
            output = {
                "All": np.sum(selectedExpGroup),
                "Any responsive": np.sum(responsiveAll)*100/np.sum(selectedExpGroup),
                "all": np.sum(selectedCellIndex)*100/np.sum(selectedExpGroup),
                "only": np.sum(selectedCellIndex),
                "EXC": np.sum(excOnly)*100/np.sum(selectedCellIndex),
                "INH": np.sum(inhOnly)*100/np.sum(selectedCellIndex),
            }
            return output, selectedCellIndex
    else:
        return selectedCellIndex

def norm_to_zero_one(row):
    min_val = np.nanmin(row)
    max_val = np.nanmax(row)
    normalized_row = (row - min_val) / (max_val - min_val)
    return normalized_row

def ecdf(data):
    """Compute ECDF for a one-dimensional array of measurements."""
    # Number of data points
    n = len(data)
    
    # x-data for the ECDF: sorted data
    x = np.sort(data)
    
    # y-data for the ECDF: evenly spaced sequence between 0 and 1
    y = np.arange(1, n + 1) / n
    
    return x, y          

# Function to calculate SNR
def calculate_SNR(flu, pre_frames, post_frames):
    # Flu is Cell x Time x Trial 
    vars_ = []
    for t in range(flu.shape[0]):
        trial = flu[t, :, :]
        mean_signal = np.nanmax(np.mean(trial[post_frames,:], axis=1))# - np.mean(trial[pre_frames,:], axis = 1))
        std_noise   = np.nanstd(trial[pre_frames,:]) #np.mean(trial[pre_frames,:], axis = 1)
        snr = np.where(std_noise != 0, mean_signal / std_noise, np.inf)
        vars_.append(snr)
    return np.array(vars_)

# Function to calculate absolute magnitude
def calculate_absMagnitude(flu, pre_frames, post_frames):
    # Flu is Cell x Time x Trial 
    vars_ = []
    for t in range(flu.shape[0]):
        trial = flu[t, :, :]
        magnitudeCalc = np.abs(np.nanmean(np.nanmean(trial[post_frames,:], axis=1) - np.nanmean(trial[pre_frames,:], axis = 1)))
        std_noise   = np.nanstd(trial[pre_frames,:]) #np.mean(trial[pre_frames,:], axis = 1)
        magnitude = np.where(std_noise != 0, magnitudeCalc, np.inf) # This line of code does nothing, but it is here to keep the same structure as the other functions
        vars_.append(magnitude)
    return np.array(vars_)

# Function to calculate Coefficient of Variation (CV)
def calculate_CV(signal):
    # Assuming isis is a list/array of inter-spike intervals for a given neuron
    # cv = calculate_CV(isis)
    std_signal = np.nanstd(signal, axis=1)
    mean_signal = np.nanmean(signal, axis=1)
    cv = np.where(mean_signal != 0, std_signal / mean_signal, np.inf)
    return cv

# Function to calculate mutual information
def calculate_MI(flu, pre_frames, post_frames):
    # Calculate Mutual Information
    # Assuming stimulus_conditions is a list/array of different stimulus conditions (e.g., visual cues),
    # and neural_responses is a list/array of the corresponding neural responses
    #mutual_information = calculate_mutual_information(stimulus_conditions, neural_responses)
    # Flu is Cell x Time x Trial 
    vars_ = []
    for t in range(flu.shape[0]):
        trial = flu[t, :, :]
        response_post  = np.nanmean(trial[post_frames,:], axis = 0)
        response_pre   = np.nanmean(trial[pre_frames,:], axis = 0)

        # Use the custom calculate_SNR function
        num_trials = response_pre.shape[0]
        stimuli = np.concatenate([
                    np.ones(num_trials)*2,  # Labels for visual condition
                    np.ones(num_trials),   # Labels for pre-visual (no) condition
                ])
        responses = np.concatenate([response_post, response_pre])
        mi_value     = mutual_info_score(responses, stimuli)
        mi_value     = mi_value.ravel()
        vars_.append(mi_value)
    vars_ = np.array(vars_)
    vars_ = vars_.flatten()
    return vars_

def variance_cell_rates(flu, frames):
    # Flu is Cell x Time x Trial - From Jimmy's code
    # vars_ = []
    # for t in range(flu.shape[0]):
    #     trial = flu[t, :, :]
    #     trial = trial[frames,:]
    #     vars_.append(np.var(np.nanmean(trial, axis=1)))
    
    # return np.array(np.array(vars_))
    # Flu is Cell x Time x Trial
    covs_ = []  # List to hold the coefficient of variation for each trial
    for t in range(flu.shape[0]):
        trial = flu[t, :, :]
        trial = trial[frames, :]  # Assuming 'frames' is an index or slice for the desired time frame
        
        # Calculate the mean and standard deviation across trials for each time point
        mean_trial = np.abs(np.nanmean(trial, axis=1))
        std_trial = np.abs(np.nanstd(trial, axis=1))
        
        # Calculate the CoV for each time point and then take the mean CoV across all time points
        # Make sure to handle the case where the mean is zero
        with np.errstate(divide='ignore', invalid='ignore'):
            cov_trial = np.where(mean_trial != 0, std_trial / mean_trial, 0)
        
        # The CoV for the trial is the mean of the CoVs across time points
        covs_.append(np.nanmean(cov_trial))

    # Return an array of the mean coefficient of variation for each trial
    return np.array(covs_)

def mean_cross_correlation(flu, frames):

    ''' Takes the mean of the absolute off-diagonal
        values of the crosscorrelation matrix


    Parameters
    ----------
    flu : fluoresence array [n_cells x n_trials x n_frames]
    frames : indexing array, frames across which to compute correlation

    Returns
    -------
    trial_corr : vector of len n_trials ->
                 mean of correlation coefficient matrix matrix on each trial.

    '''

    trial_corr = []

    for t in range(flu.shape[0]):
        trial = flu[t, :, :]
        trial = trial[frames,:].transpose()
        matrix = np.corrcoef(trial)
        matrix = matrix[~np.eye(matrix.shape[0], dtype=bool)]
        matrix = np.abs(matrix)
        mean_cov = np.mean(matrix)
        trial_corr.append(mean_cov)

    return np.array(trial_corr)


## END
# def run_suite2p (self, data_path, filename): # Not tested - 05/03/2022 HA
#     from suite2p.run_s2p import run_s2p
#     ops = {
#         # main settings
#         'nplanes' : 1, # each tiff has these many planes in sequence
#         'nchannels' : 1, # each tiff has these many channels per plane
#         'functional_chan' : 1, # this channel is used to extract functional ROIs (1-based)
#         'diameter': 12, # this is the main parameter for cell detection, 2-dimensional if Y and X are different (e.g. [6 12])
#         'tau':  1.26, # this is the main parameter for deconvolution (1.25-1.5 for gcamp6s)
#         'fs': 30.,  # sampling rate (total across planes)
#         # output settings
#         'delete_bin': False, # whether to delete binary file after processing
#         'save_mat': True, # whether to save output as matlab files
#         'combined': True, # combine multiple planes into a single result /single canvas for GUI
#         # parallel settings
#         'num_workers': 0, # 0 to select num_cores, -1 to disable parallelism, N to enforce value
#         'num_workers_roi': 0, # 0 to select number of planes, -1 to disable parallelism, N to enforce value
#         # registration settings
#         'batch_size': 500, # reduce if running out of RAM
#         'do_registration': True, # whether to register data
#         'nimg_init': 300, # subsampled frames for finding reference image
#         'maxregshift': 0.1, # max allowed registration shift, as a fraction of frame max(width and height)
#         'align_by_chan' : 1, # when multi-channel, you can align by non-functional channel (1-based)
#         'reg_tif': False, # whether to save registered tiffs
#         'subpixel' : 10, # precision of subpixel registration (1/subpixel steps)
#         # cell detection settings
#         'connected': True, # whether or not to keep ROIs fully connected (set to 0 for dendrites)
#         'navg_frames_svd': 5000, # max number of binned frames for the SVD
#         'nsvd_for_roi': 1000, # max number of SVD components to keep for ROI detection
#         'max_iterations': 20, # maximum number of iterations to do cell detection
#         'ratio_neuropil': 6., # ratio between neuropil basis size and cell radius
#         'ratio_neuropil_to_cell': 3, # minimum ratio between neuropil radius and cell radius
#         'tile_factor': 1., # use finer (>1) or coarser (<1) tiles for neuropil estimation during cell detection
#         'threshold_scaling': 0.8, # adjust the automatically determined threshold by this scalar multiplier
#         'max_overlap': 0.75, # cells with more overlap than this get removed during triage, before refinement
#         'inner_neuropil_radius': 2, # number of pixels to keep between ROI and neuropil donut
#         'outer_neuropil_radius': np.inf, # maximum neuropil radius
#         'min_neuropil_pixels': 350, # minimum number of pixels in the neuropil
#         # deconvolution settings
#         'prctile_baseline': 8.,# optional (whether to use a percentile baseline)
#         'baseline': 'maximin', # baselining mode
#         'win_baseline': 60., # window for maximin
#         'sig_baseline': 10., # smoothing constant for gaussian filter
#         'neucoeff': .7,  # neuropil coefficient
#     }
#     db = {
#     'data_path': data_path,
#     'tiff_list': data_path + filename, 
#     }
#     opsEnd = run_s2p (ops=ops,db=db)