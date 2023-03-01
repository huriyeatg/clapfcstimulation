# Functions used more globally",

import platform
import os
import pandas as pd

class analysis:
    
    def __init__(self):
        
        if platform.platform() == 'macOS-12.6-arm64-arm-64bit':
            # for Windows - Huriye PC
            print("Computer: Huriye MAC")
            self.recordingListPath = "/Users/Huriye/Documents/Code/clapfcstimulation/"
            self.rawPath           = "Z:\\Data\\"
            self.rootPath          = "/Users/Huriye/Documents/Code/clapfcstimulation/"
        elif platform.platform() == 'Windows-10-10.0.19041-SP0':
            print("Computer: Huriye Windows")
            # for Windows - Huriye PC
            self.recordingListPath = "C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\"
            self.rawPath           = "Z:\\Data\\"
            self.rootPath          = "C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\"
        self.analysisPath =  os.path.join(self.rootPath, 'analysis')
        self.figsPath     =  os.path.join(self.rootPath, 'figs')
        
        self.recordingList = pd.read_csv(self.recordingListPath + "clapfcstimulation - imaging.csv")
        
        
    def run_suite2p(self, data_path, filename) 
        ops = {
            # main settings
            'nplanes' : 1, # each tiff has these many planes in sequence
            'nchannels' : 1, # each tiff has these many channels per plane
            'functional_chan' : 1, # this channel is used to extract functional ROIs (1-based)
            'diameter': 12, # this is the main parameter for cell detection, 2-dimensional if Y and X are different (e.g. [6 12])
            'tau':  1.26, # this is the main parameter for deconvolution (1.25-1.5 for gcamp6s)
            'fs': 30.,  # sampling rate (total across planes)
            # output settings
            'delete_bin': False, # whether to delete binary file after processing
            'save_mat': True, # whether to save output as matlab files
            'combined': True, # combine multiple planes into a single result /single canvas for GUI
            # parallel settings
            'num_workers': 0, # 0 to select num_cores, -1 to disable parallelism, N to enforce value
            'num_workers_roi': 0, # 0 to select number of planes, -1 to disable parallelism, N to enforce value
            # registration settings
            'batch_size': 500, # reduce if running out of RAM
            'do_registration': True, # whether to register data
            'nimg_init': 300, # subsampled frames for finding reference image
            'maxregshift': 0.1, # max allowed registration shift, as a fraction of frame max(width and height)
            'align_by_chan' : 1, # when multi-channel, you can align by non-functional channel (1-based)
            'reg_tif': False, # whether to save registered tiffs
            'subpixel' : 10, # precision of subpixel registration (1/subpixel steps)
            # cell detection settings
            'connected': True, # whether or not to keep ROIs fully connected (set to 0 for dendrites)
            'navg_frames_svd': 5000, # max number of binned frames for the SVD
            'nsvd_for_roi': 1000, # max number of SVD components to keep for ROI detection
            'max_iterations': 20, # maximum number of iterations to do cell detection
            'ratio_neuropil': 6., # ratio between neuropil basis size and cell radius
            'ratio_neuropil_to_cell': 3, # minimum ratio between neuropil radius and cell radius
            'tile_factor': 1., # use finer (>1) or coarser (<1) tiles for neuropil estimation during cell detection
            'threshold_scaling': 0.8, # adjust the automatically determined threshold by this scalar multiplier
            'max_overlap': 0.75, # cells with more overlap than this get removed during triage, before refinement
            'inner_neuropil_radius': 2, # number of pixels to keep between ROI and neuropil donut
            'outer_neuropil_radius': np.inf, # maximum neuropil radius
            'min_neuropil_pixels': 350, # minimum number of pixels in the neuropil
            # deconvolution settings
            'prctile_baseline': 8.,# optional (whether to use a percentile baseline)
            'baseline': 'maximin', # baselining mode
            'win_baseline': 60., # window for maximin
            'sig_baseline': 10., # smoothing constant for gaussian filter
            'neucoeff': .7,  # neuropil coefficient
          }
        db = {
         'data_path': data_path,
         'tiff_list': data_path + filename, 
         }
        opsEnd = run_s2p (ops=ops,db=db)
                    
    def run_DLCpupil(self)