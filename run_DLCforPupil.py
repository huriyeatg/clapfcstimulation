import main_funcs as mfun
import deeplabcut
import glob
info = mfun.analysis()

# Get the config file - the model trained by Eren.
path_config_file = info.DLCconfigPath

# Create csv files
for ind, recordingDate in enumerate(info.recordingList.recordingDate):
    if info.recordingList.pupilImaging[ind] ==1: # if pupil recording is done & converted to .avi file
        try:
            filepathname = (info.recordingList.filepathname[ind] +'_p-' + 
                            format(int(info.recordingList.recordingID[ind]), '03d') )
            videofilename  = [f for f in glob.glob(filepathname + "\\*.avi")]
            print('Pupil size extraction from video:' + filepathname)
            deeplabcut.analyze_videos(path_config_file, videofilename, videotype='.avi',save_as_csv=True)
        except:
            print(str(ind) + ': Failed')