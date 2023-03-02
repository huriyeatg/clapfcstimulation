# Pupil preprocessing - uses Deep Lab Cut a custom model for pupil extraction

import deeplabcut
import glob

# create the configpath
configPath = "C:\\Users\\Huriye\\Documents\\code\\DeepLabCut-master\\pupilExtraction\\Updated 5 Dot Training Model-Eren CAN-2021-11-21\\"
path_config_file = configPath + "config.yaml" # Enter the path of the config file that was just created from the above step (check the folder)
print(path_config_file)


# get the list of videos
for folder in fList:
    folderDateName = (folder[0:10])
   # print(folderDateName)
    fpath     = os.path.join(input_path,str(folderDateName))
    fpath     = fpath.replace(os.sep, '/')
    fpath     = Path(fpath)
   # print(fpath)
    save_path  = os.path.join(output_path,folder)
    save_path  = save_path.replace(os.sep, '/')
   # print(save_path)
videofile_path = ['D:\DLC_output'] #Enter a folder OR a list of videos to analyze.
deeplabcut.analyze_videos(path_config_file,videofile_path, videotype='.avi',save_as_csv=True)

# new way:
configfile, path_train_config = deeplabcut.create_pretrained_project(
    Task,
    YourName,
    video,
    model=MODEL_NAME,
    videotype="avi",
    analyzevideo=True,
    createlabeledvideo=True,
    copy_videos=False,
