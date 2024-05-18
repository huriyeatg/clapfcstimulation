# This code  has plotting functions
import numpy as np
import pandas as pd
import seaborn as sns
import main_funcs as mfun
import utils_funcs as utils
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from scipy import stats
import math
import glob
from scipy.cluster.hierarchy import linkage, leaves_list

import os


def set_figure( size = 'double'):
    from matplotlib import rcParams
    if size == 'single':
        rcParams['figure.figsize'] = [3.5, 7.2]
    elif size == 'double':
        rcParams['figure.figsize'] = [7.2, 7.2]

        # set the plotting values
    rcParams['figure.figsize'] = [7.1, 3]
    rcParams['font.size'] =12
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial'] 
    rcParams['axes.linewidth'] = 1#0.5
    rcParams['xtick.major.pad'] = '0'
    rcParams['xtick.minor.pad'] = '0'
    rcParams['ytick.major.pad'] = '0'
    rcParams['ytick.minor.pad'] = '0'
    rcParams['xtick.major.size'] = 2
    rcParams['xtick.minor.size'] = 1
    rcParams['ytick.major.size'] = 2
    rcParams['ytick.minor.size'] = 1
    rcParams['xtick.major.width'] = 0.1
    rcParams['xtick.minor.width'] = 0


    rcParams['axes.spines.right']  = False
    rcParams['axes.spines.top']    = False
    rcParams['axes.spines.left']   = True
    rcParams['axes.spines.bottom'] = True

    params = {'axes.labelsize': 'large',
            'axes.titlesize':'medium',
            'xtick.labelsize':'small',
            'ytick.labelsize':'small',
            'legend.fontsize': 'large'}
    
    rcParams.update(params)

def plot_cell_trace(ax, analysispath,tTypes, s_recDate, s_animalID, s_recID, s_cellID, cell_title, y_min, y_max):
    pathname = analysispath + s_recDate + '_' + str(s_animalID) + '_' + str(s_recID).zfill(3) + '\\'
    cellID = s_cellID

    plot_dff_mean_traces(pathname, cellID, tTypes, ax)
    ax.set_ylim(ymin=y_min, ymax=y_max)
    ax.set_ylabel('DF/F')
    # ax.set_xticklabels([])
    ax.set_title(cell_title)
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # ax.spines['bottom'].set_visible(False)
    # ax.spines['left'].set_visible(False)

    # ax.axis('off')
    ax.legend(fontsize=9, frameon=False)
    # current_lines = plt.gca().get_lines()
    # for line in current_lines:
    #     line.set_linewidth(0.5)  # Adjust the value as needed

def plot_dff_mean_traces (pathname, cellID, tTypes, axis):
        ## Parameters
    fRate = 1000/30
    responsiveness_test_duration = 1000.0 #in ms 
    simulationDur_ms = 350.0 # in ms
    simulationDur = int(np.ceil(simulationDur_ms/fRate))
    pre_frames    = 2000.0# in ms
    pre_frames    = int(np.ceil(pre_frames/fRate))
    post_frames   = 6000.0 # in ms
    post_frames   = int(np.ceil(post_frames/fRate))
    analysisWindowDur = 1500 # in ms
    analysisWindowDur = int(np.ceil(analysisWindowDur/fRate))
    shutterLength     = int(np.round(simulationDur_ms/fRate))
    #tTypes = [ 'onlyVis', 'Both', 'onlyOpto']

    ########## Organise stimuli times 
    paqData = pd.read_pickle (pathname+'paq-data.pkl')
    paqRate = paqData['rate']
    # Get the stim start times 
    frame_clock    = utils.paq_data (paqData, 'prairieFrame', threshold_ttl=True, plot=False)
    optoStimTimes  = utils.paq_data (paqData, 'optoLoopback', threshold_ttl=True, plot=False)

    # the frame_clock is slightly longer in paq as there are some up to a sec delay from
    # microscope to PAQI/O software.  
    optoStimTimes = utils.stim_start_frame (paq=paqData, stim_chan_name='optoLoopback',
                                        frame_clock=None,stim_times=None, plane=0, n_planes=1)
    visStimTimes = utils.stim_start_frame (paq=paqData, stim_chan_name='maskerLED',
                                        frame_clock=None,stim_times=None, plane=0, n_planes=1)
    shutterTimes = utils.shutter_start_frame (paq=paqData, stim_chan_name='shutterLoopback',
                                        frame_clock=None,stim_times=None, plane=0, n_planes=1)

    # Lets organise it more for analysis friendly format
    trialStartTimes = np.unique(np.concatenate((optoStimTimes,visStimTimes),0))
    trialTypes = []
    for t in trialStartTimes:
        optoexist =  np.any(optoStimTimes== t)
        visexist  =  np.any( visStimTimes == t)
        if  optoexist  & visexist: 
            trialTypes += ['Visual + Opto']
        elif optoexist &~ visexist:
            trialTypes += ['Opto']
        elif ~optoexist & visexist:
            trialTypes += ['Visual']
        else:
            trialTypes += ['CHECK']
    trialStartTimes = shutterTimes
    #t = [idx for idx, t_type in enumerate(trialTypes) if t_type=='Both']

    ########## Organise calcium imaging traces 
    imData = pd.read_pickle (pathname +'imaging-data.pkl')
    fluR      = imData['flu']
    n_frames  = imData['n_frames']
    flu_raw   = imData['flu_raw']

    # Lets put nans for stimulated frames
    frameTimes = np.full((1,fluR.shape[1] ), False) # create a full false array
    for sT in shutterTimes:
        frameTimes[:,sT:(sT+shutterLength)] = True
    fluR[:, frameTimes[0,:]] = np.nan

    # clean detrended traces
    flu = utils.clean_traces(fluR)

    ### Get dff values for 4 trial types
    dffTrace ={} 
    dffTrace_mean ={}
    dffAfterStim1500ms_median ={}
    for indx, t in enumerate(tTypes) :
        if t =='All':
            trialInd = np.transpose(list(range(len(trialStartTimes))))
        else:
            trialInd = [idx for idx, t_type in enumerate(trialTypes) if t_type==t]
        
        if len(trialInd)>1:
            dffTrace[t]      = utils.flu_splitter(flu, trialStartTimes[trialInd], pre_frames, post_frames) # Cell x time x trial

    #create dff for all cells
    for indx, t in enumerate(tTypes):
        data = dffTrace[t][cellID]
        ac = np.nanmean(data[:(2*30),:], axis=0, keepdims=True)
        data = data - ac
        lineplot_withSEM (data=data, colorInd = indx, label=t, axis = axis)

def lineplot_withSEM_pupil (data, colorInd, label, axis):
    #lineplot_matrix(data=pupil_arr[session.outcome=='hit'], x_axis=x_axis, color=COLORS[0], label='hit')
    x_axis = np.linspace(-2, 6, data.shape[0])
    color =  sns.color_palette("bright")
    color  = [  (1,0,0), (0,0,0),color[2]]
    df = pd.DataFrame(data).melt()
    df['Time (sec)'] = np.tile(x_axis, data.shape[1])
    sns.lineplot(x='Time (sec)', y='value', data=df, color=color[colorInd],
                label=label,  ax = axis)
    #ylim_min = np.floor(np.min(np.nanmean(data,1)))*1.5
    #ylim_max = np.ceil(np.max(np.nanmean(data,1)))*1.5
    plt.ylabel('Pupil radius change (cm) ')  
    
def lineplot_withSEM (data, colorInd, label, axis=None):
    #lineplot_matrix(data=pupil_arr[session.outcome=='hit'], x_axis=x_axis, color=COLORS[0], label='hit')
    x_axis = np.linspace(-2, 6, data.shape[0])
    #color =  sns.color_palette("Paired")
    #color  = [ color[5], (0,0,0),color[2]]
    #print(data)
    #ac = np.nanmean(data[:, :(-2*30)], axis=1, keepdims=True)
    #data = data - ac
    palette1 = sns.color_palette('cividis') # for BOTH - blue
    palette2 = sns.color_palette('viridis') # for OPTO - greeen
    color =['grey','red', palette1[0] ,palette2[4]] # visual, opto, both
    df = pd.DataFrame(data).melt()
    df['Time (seconds)'] = np.tile(x_axis, data.shape[1])

    if axis == None:
        axis = sns.lineplot(x='Time (seconds)', y='value', data=df, color=color[colorInd],
                    label=label)
    else: 
        sns.lineplot(x='Time (seconds)', y='value', data=df, color=color[colorInd],
                    label=label, ax = axis )
    #ylim_min = np.floor(np.min(np.nanmean(data,1)))*1.5
    #ylim_max = np.ceil(np.max(np.nanmean(data,1)))*1.5
    ylim_min = np.floor(np.nanmin(np.nanmean(data,1)))*1.5
    ylim_max = np.ceil (np.nanmax(np.nanmax(data,1)))*1.5
    #if np.isnan(ylim_min): ylim_min = -1
    #if np.isnan(ylim_max): ylim_max = 5
    ylength = np.absolute(ylim_max - ylim_min)
    xlength = 0.25
    #add rectangle to plot
 
    axis.add_patch (Rectangle ((0, ylim_min), xlength, ylength, alpha = 1, facecolor="grey",zorder=10))
    plt.ylabel('DFF')
    axis.set_xlim(-1, 6)
   # axis.set_ylim(ylim_min, ylim_max)
    
def save_figure(name, base_path):
    plt.savefig(os.path.join(base_path, f'{name}.png'), 
                bbox_inches='tight', transparent=False, dpi = 300)
   # plt.savefig(os.path.join(base_path, f'{name}.svg'), 
   #             bbox_inches='tight', transparent=True)

def set_analysisParams ():
    #parameters
    fRate = 1000/30.0
    pre_frames    = 2000.0# in ms
    pre_frames    = int(np.ceil(pre_frames/fRate))
    post_frames   = 6000.0 # in ms
    post_frames   = int(np.ceil(post_frames/fRate))
    analysis_time = 1500.0 # in ms
    analysis_time = int(np.ceil(analysis_time/fRate))
    simulationDur_ms = 350.0 # in ms
    simulationDur = int(np.ceil(simulationDur_ms/fRate))
    return pre_frames, post_frames, analysis_time, simulationDur

def lineplot_withSEMWithParams (data, colorInd, label, params, axis=None):
    #lineplot_matrix(data=pupil_arr[session.outcome=='hit'], x_axis=x_axis, color=COLORS[0], label='hit')
    x_axis = np.linspace(-1*params['preStimSec'], params['postStimSec'], data.shape[0])
    color = params['plotColor'] # sns.color_palette('cividis', 3) 

    df = pd.DataFrame(data).melt()
    df['Time (seconds)'] = np.tile(x_axis, data.shape[1])

    #Normalise the data to the first 2 seconds
    basePreStim = np.nanmean(data[:(params['preStimSec']*params['fRate_imaging']),:], axis=0, keepdims=True)
    data = data - basePreStim

    if axis == None:
        axis = sns.lineplot(x='Time (seconds)', y='value', data=df, color=color[colorInd],
                    label=label)
    else: 
        sns.lineplot(x='Time (seconds)', y='value', data=df, color=color[colorInd],
                    label=label, ax = axis )
    ylim_min = np.floor(np.nanmin(np.nanmean(data,1)))*1.5
    ylim_max = np.ceil (np.nanmax(np.nanmax(data,1)))*1.5

    ylength = np.absolute(ylim_max - ylim_min)
    xlength = 0.25
    # add a  vertical line at zero label as stimulus onset
    axis.axvline(x=0, color='black', linestyle='--')
    # add a grey rectangular transparent box to indicate the stimulus duration
    axis.axvspan(0, params['visualStimSec'], color='grey', alpha=0.1)
   # axis.set_ylim(ylim_min, ylim_max)
    axis.set_xlim(-1*params['preStimSec'], params['postStimSec'])
    axis.set_xlabel('Time (seconds)')
    axis.set_ylabel('DFF')

def createTrialvsTraceMatrix (stimType, sortType=None, cohort = None, 
                              trainedLevel=None, condition=None ):# Chrimson, Naive, Visual
    # HA, 25/01/2024

    # Set params
    pre_frames, post_frames, analysis_time, simulationDur = set_analysisParams ()

    # Load data
    infoPath = r'C:\Users\Huriye\Documents\code\clapfcstimulation\analysis\infoForAnalysis-readyForPlotting_normalisedtoPre.pkl'
    dff_traceBoth, dff_traceVis, dff_traceOpto  = pd.read_pickle(infoPath)

    # Select responsive cells for the interested co
    responsiveCells = mfun.selectInterestedcells ( cohort, trainedLevel, responsive = condition, 
                                                  plotValues = False, pupil = False ) # Chrimson, Naive, Visual
    
    # Get the animal ID for each responsive cell
    infoPath = 'C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\analysis\\infoForAnalysis-readyForSelectingInterestedCells.pkl' 
    animalID, stimuliFamilarity, dataQuality,recData, recID, cellID, pvalsBoth, pvalsVis, pvalsOpto,dff_meanVisValue, dff_meanBothValue, dff_meanOptoValue, pupilID = pd.read_pickle(infoPath)    

    # load data
    dff_trace_dict = {'Visual': dff_traceVis, 'Visual + Opto': dff_traceBoth, 'Opto': dff_traceOpto}
    dff_traces = [dff_trace_dict[stim] for stim in stimType]

    plotDatas = []
    animalLists = []

    for dff_trace in dff_traces:
        if sortType:
            if sortType in dff_trace_dict:
                dff_traceMean = np.nanmean(dff_trace_dict[sortType][:, pre_frames:(pre_frames + simulationDur + analysis_time)], axis=1)
                sortedInd = np.argsort(dff_traceMean[responsiveCells])[::-1]
                plotData = dff_trace[responsiveCells][sortedInd]
                sorted_animalList = np.array([animalID[i] for i, responsive in enumerate(responsiveCells) if responsive])
                sorted_animalList = sorted_animalList[sortedInd]
            else:
                print('sortType is not defined in data dictionary')
                plotData = dff_trace[responsiveCells]
                sorted_animalList = np.array([animalID[i] for i, responsive in enumerate(responsiveCells) if responsive])
        else:
            print('sortType is not defined')
            plotData = dff_trace[responsiveCells]
            sorted_animalList = np.array([animalID[i] for i, responsive in enumerate(responsiveCells) if responsive])

        plotDatas.append(plotData)
        animalLists.append(sorted_animalList)

    return [plotDatas,animalLists[0]] if len(plotDatas) > 1 else [plotDatas[0],animalLists[0]]

def get_moreStatsValues(labels, cohort=None, trainedLevel=None, condition=None, perAnimal=False):
    # HA, 01/02/2024
    responsiveCells = mfun.selectInterestedcells(cohort, trainedLevel, responsive=condition, plotValues=False, pupil=False)
    infoPath = 'C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\analysis\\infoForAnalysis-readyForPlotting_moreStats.pkl'
    variance_dict_pre, variance_dict_post, snr_dict, mi_dict, crosscorr_dict_pre, crosscorr_dict_post, abs_dict = pd.read_pickle(infoPath)
    comparison_mapping = {
        'Visual': ('onlyVis',),
        'Opto': ('onlyOpto',),
        'Visual + Opto': ('Both',)
    }

    datasets, snrs, mis, abs = [], [], [], []
    for label in labels:
        key = comparison_mapping.get(label)[0]
        variancePre, variancePost = map(np.array, [crosscorr_dict_pre[key], crosscorr_dict_post[key]])
        datasets.append(np.array(variance_dict_post[key])) # (variancePost - variancePre))) #variance_dict_post[key]))    #
        snrs.append(np.array(snr_dict[key]))
        mis.append(np.array(mi_dict[key]))
        abs.append(np.array(abs_dict[key]))

    cv_values = [{'Condition': label, 'value': value} for label, dataset in zip(labels, datasets) for value in dataset[responsiveCells]]
    snr_values = [{'Condition': label, 'value': value} for label, snr in zip(labels, snrs) for value in snr[responsiveCells]]
    mi_values = [{'Condition': label, 'value': value} for label, mi in zip(labels, mis) for value in mi[responsiveCells]]
    abs_values = [{'Condition': label, 'value': value} for label, ab in zip(labels, abs) for value in ab[responsiveCells]]

    if perAnimal == True:
        infoPath = 'C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\analysis\\infoForAnalysis-readyForSelectingInterestedCells.pkl'    
        animalID, stimuliFamilarity, dataQuality,recData, recID, cellID, pvalsBoth, pvalsVis, pvalsOpto,dff_meanVisValue, dff_meanBothValue, dff_meanOptoValue, pupilID = pd.read_pickle(infoPath) 
        animalList = [animalID[i] for i in range(len(animalID)) if responsiveCells[i]]
        return cv_values, snr_values, mi_values, abs_values, animalList
    else:
        return cv_values, snr_values, mi_values,abs_values

def heatmap_comparison(compare1, compare2, sortType=None, cohort = None, 
                       trainedLevel=None, condition=None, colormapSelection = None,
                       axis=None, cbar_ax=None,savefigname=None, savefigpath=None, colorbarlimits=None):
    # HA, 25/01/2024
    plotData, _ = createTrialvsTraceMatrix([compare1, compare2], sortType, cohort, trainedLevel, condition)
    plotData1 = plotData[0]
    plotData2 = plotData[1]
    yminValue, ymaxValue = (-1, 1) if colorbarlimits is None else colorbarlimits
    # Parameters for plot
    step = 30
    set_figure()
    x_labels = np.linspace(-2, 6, plotData1.shape[1], dtype=int)
    xticks = np.arange(0, len(x_labels), step)
    xticklabels = x_labels[::step]
    colormap = None if colormapSelection == 'none' else None if compare1 == 'Visual' else 'viridis' if compare1 == 'Opto' else 'seagreen'

    # Plot heatmap for each condition ( compare1, compare2, compare1 - compare2)
    for i, (plotData, title) in enumerate(zip([plotData1, plotData2, plotData2 - plotData1],
                                              [compare1, compare2, 'Opto modulation'])):
        ax = plt.subplot(1, 3, i+1) if axis is None else axis[i]
        cbar_arg = {} if i != 2 else {'cbar_ax': cbar_ax, 'cbar_kws': {'label': 'DF/F'}}
        sns.heatmap(plotData, vmin=yminValue, vmax=ymaxValue, cbar=i==2, yticklabels=False,
                    cmap=colormap, ax=ax, **cbar_arg)
        ax.set_xticks(ticks=xticks, labels=xticklabels)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
        ax.set_xlim(30, 240)
        ax.set_xlabel('Time (sec)')
        ax.set_title(title)
        if i == 0:
            ax.set_ylabel(f'Sensory responsive neurons (DF/F)\n(sorted; n= {np.shape(plotData)[0]})')
            yticksValues = np.around(np.linspace(0, np.shape(plotData)[0], 10), -2)
            yticksValues = yticksValues.astype(int)
            ax.set_yticks(ticks=yticksValues)
            ax.set_yticklabels(ax.get_yticks(), rotation=0)

    if axis is not None:
        pos = cbar_ax.get_position()
        cbar_ax.set_position([pos.x0 - 0.01, pos.y0, pos.width, pos.height])

    if savefigname:
        save_figure(savefigname, savefigpath)

def population_plots(compare1, compare2, sortType=None, cohort = None, 
                       trainedLevel=None, condition=None, plotParams = None,
                       axisAll=None, savefigname=None, savefigpath=None):


    # Set params
    pre_frames, post_frames, analysis_time, simulationDur = set_analysisParams ()
    if plotParams is None:
        plotParams = {
            'ylimitsforhist': [0, 750],
            'xlimitsforhist': [-0.75, 0.75],
            'analysis_time': 1500,  # in ms
            'colorbarlimitsForHeatMap': [-1, 1],
            'scatterplotlimits': [-4.5, 4.5],
            'ylimitsforECDF': [0.6, 1.1],
            'xlimitsforABS': [-0.1, 2],
            'ylimitsforCV': [0.1, 0.15],
            'faceColors': ['black','red'],
                    }
    common_props = {'binwidth': 0.035}
        
    plotDatas,_ = createTrialvsTraceMatrix([compare1, compare2], sortType, cohort, trainedLevel, condition)
    plotData1 = plotDatas[0]
    plotData2 = plotDatas[1]
    plot_dataDIFF = plotData2 - plotData1

    # Subplot 1-2: Histogtams of Compare 1 & Compare 2
    titles = [compare1, compare2]

    data_means = [np.nanmean(data[:, pre_frames:(pre_frames + simulationDur + analysis_time)], axis=1) for data in [plotData1, plotData2, plot_dataDIFF]]
    for idx, (title, color, data_mean) in enumerate(zip(titles, plotParams['faceColors'], data_means)):
        ax = axisAll[idx]
        sns.histplot(data_mean, color=color, ax=ax, **common_props)
        ax.set_ylim(plotParams['ylimitsforhist'])
        ax.set_title(title)
        ax.set_xlabel('DF\F')
        ax.set_xlim(plotParams['xlimitsforhist'])
        ax.set_ylabel('# cells') #if idx == 0 else None

    plot_dataVmean = data_means[0] # np.nanmean(plotData1 [:, pre_frames:(pre_frames + simulationDur + analysis_time)], axis = 1)
    plot_dataBmean = data_means[1] # np.nanmean(plotData2 [:, pre_frames:(pre_frames + simulationDur + analysis_time)], axis = 1)
    stat, p_value = stats.ks_2samp(plot_dataVmean, plot_dataBmean)
    print('Stats for ' + condition)
    print(f'KS test for {compare1} and {compare2} is {stat} and p: {p_value*len(data_means[0])}') # Correct for the number of observation
    # Subplot 3:  Scatter plot
    ax = axisAll[2]
    plot_data = pd.DataFrame( {'Visual (DF/F)' :data_means[0], 
                            'Opto modulation (DF/F)':  data_means[2]})
    sns.scatterplot (y = 'Visual (DF/F)', x = 'Opto modulation (DF/F)', data = plot_data, 
                     color='black', ax=ax, linewidth = 0.5, markers='.', s=7)
    correlation_coefficient, p_value = stats.pearsonr(data_means[0], data_means[2])
    print(f'Pearson correlation coefficient for {compare1} and {compare2} is {correlation_coefficient} and p: {p_value*len(data_means[0])} ')

    ax.axhline(y = 0, color = 'black', linestyle = '--', linewidth = 0.2)
    ax.axvline(color = 'black', linestyle = '--', linewidth = 0.2)
    ax.set_ylim(plotParams['scatterplotlimits'])
    ax.set_xlim(plotParams['scatterplotlimits'])
    
    # INSET IN subplot 3:  Violin plots
    xsmall = ['V', 'V + O']
    inset_ax = inset_axes(ax, width="10%", height="25%", loc="upper right")
    s_index = np.where((20>data_means[0]) &(data_means[0]>0))[0] # Select increased response (few extreme outliers excluded )
    plot_data = pd.DataFrame( {'Mean DFF' :np.concatenate((data_means[0][s_index], data_means[1][s_index])), 
                                'Type':  np.concatenate((np.repeat('Visual', len(data_means[0][s_index])), np.repeat('Visual+Opto', len(data_means[1][s_index]))))})

    sns.barplot(x = 'Type', y = 'Mean DFF', data = plot_data, palette=plotParams['faceColors'], ax=inset_ax, linewidth = 0.1)
    #plt.axhline(y = 0, color = 'black', linestyle = '--', linewidth = 0.1)
    plt.annotate('p = {:.3f}'.format(stats.ttest_rel(plot_dataVmean[s_index], plot_dataBmean[s_index], alternative = 'greater')[1], 3), 
                xy=(0.5, 1), xytext=(0, 0), textcoords='offset points', 
                ha = 'center', va = 'top', fontsize = 5)
    inset_ax.set_xticklabels(xsmall, rotation = 45)
    inset_ax.set_ylabel('')
    inset_ax.set_xlabel('')
    inset_ax.set_ylim(-0.05,0.7)
    for item in ([inset_ax.title, inset_ax.xaxis.label, inset_ax.yaxis.label] +
             inset_ax.get_xticklabels() + inset_ax.get_yticklabels()):
        item.set_fontsize(5)

    inset_ax = inset_axes(ax, width="10%", height="25%", loc="lower left",
                          bbox_to_anchor=(0.1, 0.1, 1, 1),
                            bbox_transform=ax.transAxes)
    s_index = np.where((0>data_means[0]) &(data_means[0]>-20))[0] # Select decreased response (few extreme outliers excluded )
    plot_data = pd.DataFrame( {'Mean DFF' :np.concatenate((data_means[0][s_index], data_means[1][s_index])), 
                                    'Type':  np.concatenate((np.repeat('Visual', len(data_means[0][s_index])), np.repeat('Visual+Opto', len(data_means[1][s_index]))))})

    sns.barplot(x = 'Type', y = 'Mean DFF', data = plot_data, palette=plotParams['faceColors'], ax=inset_ax,  linewidth = 0.1)
    inset_ax.annotate('p = {:.3f}'.format(stats.ttest_rel(plot_dataVmean[s_index], plot_dataBmean[s_index], alternative = 'greater')[1], 3), 
                xy=(0.5, 1), xytext=(0, 0), textcoords='offset points', 
                ha = 'center', va = 'top', fontsize = 5)
    inset_ax.set_ylim(-0.4,0.1)
    inset_ax.set_ylabel('')
    inset_ax.set_xlabel('')
    inset_ax.set_xticklabels(xsmall,rotation=45)
    for item in ([inset_ax.title, inset_ax.xaxis.label, inset_ax.yaxis.label] +
            inset_ax.get_xticklabels() + inset_ax.get_yticklabels()):
        item.set_fontsize(5)


    ###### Subplot 4: Cumulative ECDF
    ax = axisAll[4]
    responses_vis = np.nanmean(plotData1[:, pre_frames:(pre_frames + simulationDur + analysis_time)], axis=1)
    x, y = mfun.ecdf(responses_vis)
    ax.plot(x, y, marker='.', markersize=1, color= plotParams['faceColors'][0], linestyle='none',label='Visual')

    responses_both =  np.nanmean(plotData2[:, pre_frames:(pre_frames + simulationDur + analysis_time)], axis=1)
    x, y = mfun.ecdf(responses_both)
    ax.plot(x, y, marker='.', markersize=1, linestyle='none', color = plotParams['faceColors'][1], label='Visual + Opto')
 
    ax.set_xlabel('Neural response (DF/F)')
    ax.set_ylabel('Cumulative density')
    ax.set_ylim(plotParams['ylimitsforECDF'])
    ax.set_xscale('log') 
    #ax.legend(fontsize='small',loc='upper left', frameon=False)
    
    cv_values, snr_values, mi_values, abs_values = get_moreStatsValues([compare1, compare2], cohort = cohort ,trainedLevel=trainedLevel, condition=condition)
    
    ###### Subplot 5: Absolute magnitude
    ax2 = axisAll[3]
    abs_df = pd.DataFrame(abs_values)
    sns.histplot(data=abs_df, x='value', hue="Condition", element="step",
                 palette=plotParams['faceColors'], ax=ax2, **common_props)
    ax2.set_xlabel('Absolute magnitude')
    ax2.set_xlim(plotParams['xlimitsforABS'])
    ax2.set_ylabel('# cells')
    ax2.legend(labels=['Visual + Opto', 'Visual'],fontsize='small', loc='upper right', frameon=False)
    #data1 = abs_values['value'][abs_values['Condition'] == compare1]
    #data2 = abs_values['value'][abs_values['Condition'] == compare2]
    #stat, p_value = stats.ks_2samp(data1, data2)
    #print(f'ABS: KS test: {len(data1)} {stat} and p: {p_value*len(data_means[0])}') # Correct for the number of observation

    ###### Subplot 6: Coefficiency of variation differences between two conditions
    #print(cv_values)
    ax8 = axisAll[5]
    plotParams['faceColors'][0] = 'grey'
    data_df = pd.DataFrame(cv_values)
    common_props = {'ci': 'sem'}
    sns.boxplot(data=data_df, x="Condition", y="value", showfliers = False, width=0.6, saturation=0.75, #ci =' sem',  #outliers= False,  #element="step",
                 palette=plotParams['faceColors'], ax=ax8)
    ax8.set_ylabel ('Coefficient of Variation')
    stat, p_value = stats.wilcoxon(data_df[data_df['Condition'] == 'Visual']['value'],
                                    data_df[data_df['Condition'] == 'Visual + Opto']['value'])
    
    print(p_value)
    y_max = data_df.groupby('Condition')['value'].quantile(0.75).max() + (data_df.groupby('Condition')['value'].quantile(0.75) - data_df.groupby('Condition')['value'].quantile(0.25)).max() * 1.5
    y_min = 0#data_df.groupby('Condition')['value'].quantile(0.25).min() - (data_df.groupby('Condition')['value'].quantile(0.25) - data_df.groupby('Condition')['value'].quantile(0.25)).min() * 1.5

    if p_value < 0.05:
            ax8.annotate('p < {:.3f}'.format(max(p_value, 0.001), 3),
                         xy=(0.5, y_max*1.05), ha='center', color='black')
    ax8.set_ylim(y_min*1.2, y_max*1.2)
    
    # save figure
    if savefigname:
        save_figure(savefigname, savefigpath)

def plot_magnitude(compare1, compare2, cohort = None, trainedLevel=None, condition=None, 
                   plotParams = None,axisAll=None):
    cv_values, snr_values, mi_values, abs_values = get_moreStatsValues([compare1, compare2], cohort = cohort ,trainedLevel=trainedLevel, condition=condition)
    if plotParams is None:
        plotParams = {
            'ylimitsforhist': [0, 750],
            'xlimitsforhist': [-0.75, 0.75],
            'analysis_time': 1500,  # in ms
            'colorbarlimitsForHeatMap': [-1, 1],
            'scatterplotlimits': [-4.5, 4.5],
            'ylimitsforECDF': [0.6, 1.1],
            'xlimitsforABS': [-0.1, 1.2],
            'xlimitsforCV': [-3, 3],
            'faceColors': ['black','red'],
                    }
    common_props = {'binwidth': 0.035}
    # Plot Absolute magnitude
    ax = axisAll
    snr_df = pd.DataFrame(abs_values)
    sns.histplot(data=snr_df, x='value', hue="Condition", element="step",
                 palette=plotParams['faceColors'], ax=ax, **common_props,
                 legend = True)
    ax.set_xlabel('Absolute magnitude')
    ax.set_xlim(plotParams['xlimitsforABS'])
    ax.set_ylabel('# cells')
    ax.set_title(condition)
    ax.legend(labels=['Visual + Opto', 'Visual'], fontsize='small', loc='lower right', bbox_to_anchor=(1.1, 0), frameon=False)

    # Add insert plot for ECDF
    pre_frames, post_frames, analysis_time, simulationDur = set_analysisParams ()
    sortType = 'Visual'
    plotDatas,_ = createTrialvsTraceMatrix([compare1, compare2], sortType, cohort, trainedLevel, condition)
    plotData1 = plotDatas[0]
    plotData2 = plotDatas[1]
    inset_ax = inset_axes(ax, width="30%", height="50%", loc="upper right")
    responses_vis = np.nanmean(plotData1[:, pre_frames:(pre_frames + simulationDur + analysis_time)], axis=1)
    x, y = mfun.ecdf(responses_vis)
    inset_ax.plot(x, y, marker='.', markersize=1, color= plotParams['faceColors'][0], linestyle='none',label='Visual')

    responses_both =  np.nanmean(plotData2[:, pre_frames:(pre_frames + simulationDur + analysis_time)], axis=1)
    x, y = mfun.ecdf(responses_both)
    inset_ax.plot(x, y, marker='.', markersize=1, linestyle='none', color = plotParams['faceColors'][1], label='Visual + Opto')
 
    inset_ax.set_xlabel('DF/F',fontsize='small')
    inset_ax.set_ylabel('ECDF',fontsize='small')
    inset_ax.set_ylim(plotParams['ylimitsforECDF'])
    inset_ax.set_xscale('log') 
    #inset_ax.legend(fontsize='small',loc='upper left', frameon=False)

def plot_paramsDiffPerAnimal(params, cohorts, trainedLevel, ax=None, savefigname=None, 
                             savefigpath=None,ComparePlot = False):
    # HA, 01/02/2024
    colorpalet = [ 'red', 'green', 'blue', 'grey'] #'black',
    xlabel =  ['Sensory', 'Opto-boosted', 'Opto', 'None']#, 'All', 'Not responsive']

    # Create an empty DataFrame to store all differences
    diff_df = pd.DataFrame()

    for cond in trainedLevel:
        for cohort in cohorts:
            for idx, cellresponsiveness in enumerate(['Sensory', 'Opto-boosted', 'Opto', 'None']):
                cv_values_visual, snr_values_visual, mi_values_visual, abs_values_visual, animalID = get_moreStatsValues(['Visual'], cohort = cohort ,trainedLevel=cond, condition=cellresponsiveness, perAnimal = True)
                cv_values_visual_opto, snr_values_visual_opto,  mi_values_visual_opto, abs_values_opto, animalID= get_moreStatsValues(['Visual + Opto'], cohort = cohort ,trainedLevel=cond, condition=cellresponsiveness, perAnimal = True)
                if params == 'SNR':
                    values_visual = snr_values_visual
                    values_visual_opto = snr_values_visual_opto
                    ylabel = '$\Delta$ SNR'
                elif params == 'CV':
                    values_visual = cv_values_visual
                    values_visual_opto = cv_values_visual_opto
                    ylabel = '$\Delta$ CV'
                elif params == 'MI':
                    values_visual = mi_values_visual
                    values_visual_opto = mi_values_visual_opto
                    ylabel = 'MI Difference\n(MI$_{Visual + Opto}$  - MI$_{Visual}$)'
                elif params == 'ABS':
                    values_visual = abs_values_visual
                    values_visual_opto = abs_values_opto
                    ylabel = '$\Delta$ Absolute magnitude'
                elif params == 'ECDF Shift':
                    pre_frames, post_frames, analysis_time, simulationDur = set_analysisParams ()
                    plotDatas,animalListMatrixIndices = createTrialvsTraceMatrix(['Visual','Visual + Opto'], 'Visual', cohort, cond, cellresponsiveness)
                    plotData1 = plotDatas[0]
                    plotData2 = plotDatas[1]
                    ylabel = '$\Delta$ ECDF Shift'
                elif params == 'ECDF Stepness':
                    pre_frames, post_frames, analysis_time, simulationDur = set_analysisParams ()
                    plotDatas, animalListMatrixIndices= createTrialvsTraceMatrix(['Visual','Visual + Opto'], 'Visual', cohort, cond, cellresponsiveness)
                    plotData1 = plotDatas[0]
                    plotData2 = plotDatas[1]
                    ylabel = '$\Delta$ ECDF Stepness'

                if 'Trained' in trainedLevel:
                    _, _, _,_, animalID = get_moreStatsValues(['Visual'], cohort = cohort ,trainedLevel='Trained', condition=cellresponsiveness, perAnimal = True)
                    animalList = np.unique(animalID)
                else:
                    animalList = np.unique(animalID)
                
                for animal in animalList:
                    # if params includes ECDF regardless of Leftwardness or Stepness
                    if params == 'ECDF Stepness':
                        # get index for this animal
                        indices = np.where(animalListMatrixIndices == animal)
                        responses_vis = np.nanmean(plotData1[indices[0], pre_frames:(pre_frames + simulationDur + analysis_time)], axis=1)
                        responses_both = np.nanmean(plotData2[indices[0], pre_frames:(pre_frames + simulationDur + analysis_time)], axis=1)
                        # Let's say we want to calculate the steepness between the 25th and 75th percentiles
                        low_percentile = 25
                        high_percentile = 75

                        # Find the response levels (x-values) at these percentiles
                        x_low_vis = np.percentile(responses_vis, low_percentile)
                        x_high_vis = np.percentile(responses_vis, high_percentile)

                        x_low_both = np.percentile(responses_both, low_percentile)
                        x_high_both = np.percentile(responses_both, high_percentile)

                        # Find the cumulative probability (y-values) at these response levels
                        y_low_vis = np.searchsorted(np.sort(responses_vis), x_low_vis, side='right') / len(responses_vis)
                        y_high_vis = np.searchsorted(np.sort(responses_vis), x_high_vis, side='right') / len(responses_vis)

                        y_low_both = np.searchsorted(np.sort(responses_both), x_low_both, side='right') / len(responses_both)
                        y_high_both = np.searchsorted(np.sort(responses_both), x_high_both, side='right') / len(responses_both)

                        # Calculate the steepness of the ECDF between these points
                        steepness_vis = (y_high_vis - y_low_vis) / (x_high_vis - x_low_vis)
                        steepness_both = (y_high_both - y_low_both) / (x_high_both - x_low_both)
                        diff = steepness_both - steepness_vis
                    elif params == 'ECDF Shift':
                        # get index for this animal
                        indices = np.where(animalListMatrixIndices == animal)
                        responses_vis = np.nanmean(plotData1[indices[0], pre_frames:(pre_frames + simulationDur + analysis_time)], axis=1)
                        responses_both = np.nanmean(plotData2[indices[0], pre_frames:(pre_frames + simulationDur + analysis_time)], axis=1)
                        # Find the median response level
                        median_vis1 = np.percentile(responses_vis, 5)
                        median_vis3 = np.percentile(responses_vis, 95)
                        median_both1 = np.percentile(responses_both,5)
                        median_both3= np.percentile(responses_both, 95)

                        median_vis = np.abs(median_vis3 - median_vis1)
                        median_both = np.abs(median_both3 - median_both1)

                        # Calculate the leftward shift (if the result is negative, it's actually a rightward shift)
                        diff = median_both - median_vis
                    else:
                        # get index for this animal 
                        indices = np.where(animalID == animal)
                        vis = [values_visual[i] for i in indices[0]]
                        visOpt = [values_visual_opto[i] for i in indices[0]]
                        visual_df = pd.DataFrame(vis)
                        visual_opto_df = pd.DataFrame(visOpt)
                        diff = np.nanmedian((visual_opto_df['value'] - visual_df['value']))
                    
                    temp_df = pd.DataFrame({
                        'Type' : cond,
                        'Condition': [xlabel[idx]],
                        'Difference': diff,
                        })
                    diff_df = pd.concat([diff_df, temp_df], ignore_index=True)
    
    y_max = diff_df['Difference'].max() + diff_df['Difference'].std()
    y_min = diff_df['Difference'].min() - diff_df['Difference'].std()  # This will be the y coordinate for our stars

    if ComparePlot:
        colorpalet = [ 'black', 'magenta']
        conditions = ['Sensory', 'Opto-boosted', 'Opto'] 
        sns.swarmplot(data=diff_df, x='Condition', y='Difference', hue = 'Type', palette= colorpalet,
                       ax = ax, dodge= True)
        ax.legend (loc = 'lower right', fontsize = 10, frameon = False)
        p_values = []
        for condition in conditions:
            stat, p = stats.wilcoxon(diff_df[(diff_df['Type'] =='Trained') &(diff_df['Condition'] ==condition)]['Difference'],
                                    diff_df[(diff_df['Type'] =='Naive') &(diff_df['Condition'] ==condition)]['Difference'])
            p_values.append(p)
        print(p_values)
        x_coord = 0
        for p_val in p_values:
            if p_val <= 0.0001:
                ax.annotate('p < 0.001'.format(p_val, 3), xy=(x_coord, y_max), 
                            ha='center', color='black', fontsize=10)
            elif (p_val > 0.0001) & (p_val <= 0.05):
                ax.annotate('p = {:.3f}'.format(p_val, 3), xy=(x_coord, y_max), 
                            ha='center', color='black', fontsize=10)
            x_coord = x_coord+1  # The x coordinate for the star
    else:
        # Calculate the stats - THIS IS NOT CORRECT - IGNORED FOR NOW!
        conditions = ['Sensory', 'Opto-boosted',] # Opto has less animal
        p_values = []
        for condition in conditions:
            stat, p = stats.ttest_rel(diff_df[diff_df['Condition'] == condition]['Difference'],
                                    diff_df[diff_df['Condition'] == 'None']['Difference'])
            p_values.append(p)
        print(ylabel + ' - Total animal number for ' + str(len(np.unique(animalID))))
        print(p_values)
        sns.swarmplot(data=diff_df, x='Condition', y='Difference', palette= colorpalet, ax = ax)
        ax.legend().set_visible(False)
   
    ax.set_ylim(y_min, y_max*1.2)
    ax.set_xlim(-0.5, 2.5)
    ax.axhline(0, color='grey', linestyle='--')
    ax.tick_params(axis='x', labelsize='x-small') # You can specify a numerical value instead of 'small'
    ax.set_ylabel(ylabel)
    ax.set_xlabel('Neural populations')
    
    # save figure
    if savefigname:
        save_figure(savefigname, savefigpath)

def plot_cellRatiosPerAnimal(params, cohorts, trainedLevel, ax=None, savefigname=None, 
                             savefigpath=None, ComparePlot = False):
    # HA, 01/02/2024
    infoPath = 'C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\analysis\\infoForAnalysis-readyForSelectingInterestedCells.pkl' 
    animalID, stimuliFamilarity, dataQuality,recData, recID, cellID, pvalsBoth, pvalsVis, pvalsOpto,dff_meanVisValue, dff_meanBothValue, dff_meanOptoValue, pupilID = pd.read_pickle(infoPath)    
    palette1 = sns.color_palette('cividis') # for OPTO - blue
    palette2 = sns.color_palette('viridis') # for BOTH - green
    xlabel =   ['Sensory','Opto-boosted','Opto','All']
    if params == 'All':
        cellresponsivenessList = ['Sensory', 'Opto-boosted', 'Opto','All']
        colorpalet = ['red',palette2[4], palette1[0] , 'grey'] # all, visual, opto, both
    elif params == 'Sensory':
        cellresponsivenessList = ['Sensory']
        faceColors = 'red'
    elif params == 'Opto':
        cellresponsivenessList = ['Opto']
        faceColors = palette1[0]
    elif params == 'Opto-boosted':
        cellresponsivenessList = ['Opto-boosted']
        faceColors = palette2[4]

    # Create an empty DataFrame to store all differences
    diff_df = pd.DataFrame()

    for tLevel in trainedLevel:
        for cohort in cohorts:
            for idx, cellresponsiveness in enumerate(cellresponsivenessList):
                if 'Trained' in trainedLevel:
                    _, _, _,_, animalID = get_moreStatsValues(['Visual'], cohort = cohort ,trainedLevel='Trained', condition=cellresponsiveness, perAnimal = True)
                    animalList = np.unique(animalID)
                else:
                    _, _, _,_, animalID = get_moreStatsValues(['Visual'], cohort = cohort ,trainedLevel=tLevel, condition=cellresponsiveness, perAnimal = True)
                    animalList = np.unique(animalID)
                animalList = np.unique(animalList)
                print(cellresponsiveness + ' - Total animal number for ' + str(len(animalList)))
                for animal in animalList:
                    output, _= mfun.selectInterestedcells(animal, tLevel, responsive=cellresponsiveness, plotValues=True, pupil=False)
                    temp_df = pd.DataFrame({
                        'Type' : tLevel,
                        'Condition': [xlabel[idx]],
                        'All': output['all'],
                        'EXC': output['EXC'],
                        'INH': output['INH'],
                    })
                    diff_df = pd.concat([diff_df, temp_df], ignore_index=True)
    
    if ComparePlot:
        sns.violinplot(data=diff_df, x='Type', y=params, hue = 'Condition', palette= colorpalet, ax = ax)
        ax.legend (loc = 'upper left', fontsize = 10, frameon = False)
        ax.set_ylabel('Percentage of neurons (%)')
        ax.set_xlabel('Neural population based on responsiveness')

    else:
        if params== 'All':
            sns.violinplot(data=diff_df, x='Condition', y=params,palette= colorpalet, width = 1.5, ax = ax)
            ax.set_ylabel('Percentage of neurons (%)')
            ax.set_xlabel('Neural population based on responsiveness')
            ax.legend().set_visible(False)
            ax.set_xlim(-1, 3.7)
            ax.set_ylim(0,90)
        else:
            ax.axhline(y = 50, color ='grey', linestyle = '--', linewidth = 0.2)
            sns.violinplot(x = 'variable', y='value', data=diff_df.melt(value_vars=['EXC', 'INH']), 
                        ax = ax, color = faceColors, width=0.4)
            statVal, pval = stats.wilcoxon(diff_df['EXC'], diff_df['INH'])
            print(params + ': p ='+ str(pval) + ' Wilcoxon: ' + str(statVal))
            
            if pval<0.05:
                ax.annotate('p = {:.3f}'.format(max(pval, 0.001), 3),xy=(0.5, 90 ), 
                            ha = 'center', fontsize = 10)
            ax.set_xticklabels(['Exc', 'Inh'])  # Set custom x-labels
            ax.set_ylabel('% of neurons', fontsize = 10)
            ax.set_ylim(0,100)
            ax.set_xlabel('')
    # save figure
    if savefigname:
        save_figure(savefigname, savefigpath)

def plot_exampleTrainingBehaviour(s_animalID,duration, ax): 
    preStimDur = duration[0]
    allStimDur = duration[1] + preStimDur
    visualStimDur = 2 # This is a fix value
    # Load the info file
    infoPath = 'C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\analysis\\infoForAnalysis-extracted.pkl'
    info = pd.read_pickle(infoPath) 

    fRate = 20000
    s_stimuliID = 5 # For training sessions 
    animalID = info.recordingList.animalID
    stimuliID = info.recordingList.stimuliFamiliarity


    # Create truncated trials from all sessions for selected animal and stimuli
    indList = np.where((np.array(animalID) == s_animalID) & (np.array(stimuliID) == s_stimuliID ))[0]

    animal_lick  =[]
    animal_water =[]
    for ind in (indList):
        savepathname = info.recordingList.analysispathname[ind]
        pathname = [f for f in glob.glob(savepathname + 'training-paq-data.pkl')]
        paqData = pd.read_pickle (pathname[0])
        # Get the stim start times 
        trialStartTimes = utils.paq_data (paqData, 'maskerLED', 1, threshold_ttl=True)

        trialStartTimes = trialStartTimes - (visualStimDur*fRate)
        if len(trialStartTimes)<500:
            licks, trial_licks = utils.lick_binner(savepathname, trialStartTimes,'lickDetection', stimulation = False)
            waterPoints, trial_water = utils.lick_binner(savepathname,trialStartTimes, 'waterDelivery', stimulation = False)

            animal_lick = animal_lick + trial_licks
            animal_water = animal_water + trial_water

    for i, array in enumerate(animal_lick):
        ax.plot(array, np.ones_like(array)+i, 'k.',markersize = 1)
        #plt.plot(animal_water[i], np.ones_like(animal_water[i])+i+0.3, 'bo',markersize = 2)
    ymax = math.floor(len(animal_lick) / 150.0) * 150
    ax.set_xlim(0, allStimDur*fRate)
    ax.set_ylim(0,ymax)
    ax.set_yticks(range(0,ymax, 150), range(ymax,0, -150))
    ax.set_xticks (range(0,(allStimDur*fRate)+1,fRate), range(preStimDur*-1,(allStimDur-1),1))
    ax.set_ylabel('Trials')
    ax.set_xlabel('Time (sec)')
    # No box on the top and right side
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # add a vertical line at 0 in x axis
    ax.axvline(x=preStimDur*fRate, color='magenta', linestyle='-', linewidth=2)
    ax.text((preStimDur-0.45)*fRate, len(animal_lick)-30, 'Stim', fontsize=8, ha='center')
    # add a dashed line for the end of the stimulation
    ax.axvline(x=(preStimDur+visualStimDur) *fRate, color='k', linestyle='--', linewidth=2)
    ax.text((preStimDur+visualStimDur-0.45)*fRate, len(animal_lick)-30, 'Water', fontsize=8, ha='center')
    ax.set_title('Sample Animal Training (Truncated Sessions)')
    
def plot_lickDensityTraces(cohort,duration, axTrace, axSummary): 
    # Plots the average lick density for each animal & average of all animals
    # Requires 3 axis
    infoPath = 'C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\analysis\\infoForAnalysis-extracted.pkl'
    info = pd.read_pickle(infoPath) 

    animalID = info.recordingList.animalID
    stimuliID = info.recordingList.stimuliFamiliarity
        
    if cohort == 'Chrimson':
        animalList =  [22101,22102, 22103,22105,22107, 2303, 2304]
    elif cohort == 'OPN3':
        animalList = [2306, 2307,2308,2309,2310, 2311, 2312]

    fRate = 20000
    s_stimuliID = 5 # For training sessions
    types = ['Last Sessions', 'First Sessions']
    preStimDur = duration[0]
    allStimDur = duration[1] + preStimDur
    visualStimDur = 2 # This is a fix value
    bin_width = 10  # Adjust this value to change the width of each bin
    num_bins = range(0, allStimDur*fRate, int(fRate/bin_width))

    # Create behaviour sessions plots
    ax = axTrace[0]
    print(animalList)
    for s_animalID in animalList:
        indListAll = np.where((np.array(animalID) == s_animalID) & (np.array(stimuliID) == s_stimuliID ))[0]
        chName = 'lickDetection'
        if s_animalID == 22101:
            indListAll = indListAll[:-2]
        elif s_animalID == 22102:
            indListAll = indListAll[:-3]
        elif s_animalID == 22103:
            indListAll = indListAll[3:]
        elif s_animalID == 2303 or s_animalID == 2304:
            indListAll = indListAll[:-1]
            indListAll = indListAll[::-1]
            #chName = 'waterDelivery' # Channel names mixed in these animals

        for type in types:
            if type == 'Last Sessions':
                indList = indListAll[:1]
                colorV = 'magenta'
            elif type == 'First Sessions':
                indList = indListAll[-1:]
                colorV = 'k'
            
            animal_lick  =[]
            animal_water =[]
            for ind in (indList):
                savepathname = info.recordingList.analysispathname[ind]
                pathname = [f for f in glob.glob(savepathname + 'training-paq-data.pkl')]
                paqData = pd.read_pickle (pathname[0])
                # Get the stim start times 
                trialStartTimes = utils.paq_data (paqData, 'maskerLED', 1, threshold_ttl=True)
                trialStartTimes = trialStartTimes - (visualStimDur*fRate)

                if len(trialStartTimes)<500:
                    licks, trial_licks = utils.lick_binner(savepathname, trialStartTimes,chName, stimulation=False)
                    animal_lick = animal_lick + trial_licks
            
            animal_hist = np.zeros(len(num_bins)-1)
            for i, array in enumerate(animal_lick):
                hist, bins = np.histogram(array, bins=num_bins, range=(0, allStimDur*fRate))
                animal_hist = animal_hist + hist
            
            animal_hist = (animal_hist/ (i+1))*bin_width

            if animal_hist[(4*bin_width)]> 20: #PackI/O crosstalk between channel creates a large TTL ehen reward is given, this is not lick - only exist in two animal
                animal_hist[(4*bin_width)] = np.nan
            ax.plot( bins[1:],animal_hist, colorV, linewidth =1, alpha= 0.5)

            if type == 'First Sessions':
                if s_animalID == animalList[0]:
                    animal_hist_FirstAll = animal_hist
                else:
                    animal_hist_FirstAll = np.vstack((animal_hist_FirstAll,animal_hist))
            else:
                if s_animalID == animalList[0]:
                    animal_hist_LastAll = animal_hist
                else:
                    animal_hist_LastAll = np.vstack((animal_hist_LastAll,animal_hist))

    ax.plot( bins[1:],np.nanmean(animal_hist_LastAll,0 ), 'magenta', linewidth =3, alpha= 1, label = 'Last Session')
    ax.plot( bins[1:],np.nanmean(animal_hist_FirstAll,0), 'k', linewidth =3, alpha= 1, label = 'First Session')
    ax.legend(fontsize='small',loc='upper right', frameon=False)
    ymax = math.ceil(np.nanmax(np.nanmax(animal_hist_LastAll,0 )) / 10.0) * 10

    ax.set_xlim(0, allStimDur*fRate)
    ax.set_xticks (range(0,(allStimDur*fRate)+1,fRate), range(preStimDur*-1,(allStimDur-1),1))
    ax.set_ylim(-0.1,ymax+2)
    ax.set_ylabel('Lick Density (Hz)')
    ax.set_xlabel('Time (sec)')
    ax.axvline(x=preStimDur*fRate, color='magenta', linestyle='-', linewidth=2)
    ax.text((preStimDur-0.6)*fRate, ymax+1, 'Stimuli', fontsize=8, ha='center')
    # add a dashed line for the end of the stimulation
    ax.axvline(x=(preStimDur+visualStimDur) *fRate, color='k', linestyle='--', linewidth=2)
    ax.text((preStimDur+visualStimDur-0.5)*fRate, ymax+1, 'Reward', fontsize=8, ha='center')
    ax.set_title('Trained animals (n = ' + str(len(animalList)) +')')

    ax = axSummary[0]
    # violing plot for each animal average of 1 sec before water delivery and 1 sec before stimuli
    
    OnesecbeforeStimWindow = np.where((bins>(preStimDur-1)*fRate) & (bins<(preStimDur)*fRate))[0]
    mean1secbeforeStim = np.nanmean(animal_hist_LastAll[:,OnesecbeforeStimWindow], axis = 1)
    OnesecbeforeRewardWindow = np.where((bins>(preStimDur-1+visualStimDur)*fRate) & (bins<(preStimDur+visualStimDur)*fRate))[0]
    mean1secbeforeReward = np.nanmean(animal_hist_LastAll[:,OnesecbeforeRewardWindow], axis = 1)
    df = pd.DataFrame({
        'mean1secbeforeStim': mean1secbeforeStim,
        'mean1secbeforeReward': mean1secbeforeReward
        })

    sns.swarmplot(data=[mean1secbeforeStim, mean1secbeforeReward], ax = ax, palette = ['k', 'magenta'])
    for index, row in df.iterrows():
        plt.plot([0, 1], [row['mean1secbeforeStim'], row['mean1secbeforeReward']], color='grey', alpha=0.5)

    ax.set_ylim(-0.5,ymax+2)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel('Lick Density (Hz)')
    ax.set_xticks([0,1])
    ax.set_xticklabels(['before stimuli', 'before reward'], rotation=15)
    return   df['mean1secbeforeReward'] - df['mean1secbeforeStim'] 

def plot_correlationMatrix(stimType,cohort, trainedLevel, responsiveness,  
                           axs=None, savefigname=None, savefigpath=None):
    # Load data
    infoPath = f'C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\analysis\\crossCorrelation_{cohort}_{trainedLevel}_{responsiveness}.pkl'
    linkage_matrices, correlation_matrices, animalID = pd.read_pickle(infoPath)
    # Create or use provided axes
    if axs is None:
        fig, axs = plt.subplots()
    else:
        fig = axs[1].figure
    
    # Select data based on stimType
    Z1, corr_matrix = linkage_matrices[0], correlation_matrices[0]
    Z2, corr_matrix2 = linkage_matrices[1], correlation_matrices[1]
    #Z3, corr_matrix3 = linkage_matrices[2], correlation_matrices[2]

    # Sort correlation matrix
    leaf_order = leaves_list(Z1)
    sorted_corr_matrix = corr_matrix[np.ix_(leaf_order, leaf_order)]
    sorted_corr_matrix2 = corr_matrix2[np.ix_(leaf_order, leaf_order)] # same order with Visual
    cax = axs[0].imshow(sorted_corr_matrix, interpolation='nearest', cmap='coolwarm', aspect='auto')
    cax = axs[1].imshow(sorted_corr_matrix2, interpolation='nearest', cmap='coolwarm', aspect='auto')
    fig.colorbar(cax, ax=axs, label='Cross-Correlation')

    # Dynamic title based on stimType
    axs[0].set_title(f'Visual')
    axs[1].set_title(f'Visual + Opto')
    axs[0].set_ylabel('Cross correlation\nAll responsive cells')
    # Turn off y labels for second plot
    axs[1].set_yticks([])
    axs[0].set_xlabel('All responsive cells')
    axs[1].set_xlabel('All responsive cells')

    # Save figure if requested
    if savefigname and savefigpath:
        plt.savefig(f"{savefigpath}\\{savefigname}_{stimType}_CorrelationMatrix.png")

def plot_correlationMatrix_meanChangeDELETE(stimType,cohort, trainedLevel, responsiveness, params, 
                           axs=None, savefigname=None, savefigpath=None):
    if params == 'All':
        cellresponsivenessList = ['All','None']
        colorpalet = ['grey', 'black'] 
    else:
        print('StimType order is not matched with colors')
        cellresponsivenessList = ['Visual', 'Visual + Opto', 'Opto','All']
        palette1 = sns.color_palette('cividis') # for OPTO - blue
        palette2 = sns.color_palette('viridis') # for BOTH - green
        colorpalet = ['red',palette2[4], palette1[0] , 'grey'] # all, visual, opto, both
        xlabel =   ['Visual','Opto-boosted','Opto','All']


    # Initialize a dictionary to hold average changes for each animal across conditions
    average_changes = {condition: [] for condition in responsiveness}
    for condition in responsiveness:
        # Load data for this condition
        infoPath = f'C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\analysis\\crossCorrelation_{cohort}_{trainedLevel}_{condition}.pkl'

        #Z,Z2,Z3,corr_matrix,corr_matrix2,corr_matrix3 = pd.read_pickle(infoPath)
        linkage_matrices, correlation_matrices, animalID = pd.read_pickle(infoPath)
            # Select data based on stimType
        if stimType[0] == 'Visual':
            Z, corr_matrix = linkage_matrices[0], correlation_matrices[0]
        elif stimType[0] == 'Opto':
            Z, corr_matrix = linkage_matrices[1], correlation_matrices[1]
        elif stimType[0] == 'Visual + Opto':
            Z, corr_matrix = linkage_matrices[2], correlation_matrices[2]
        else:
            raise ValueError(f"Unknown stimType: {stimType[0]}")
        
        if stimType[1] == 'Visual':
            Z, corr_matrix2 = linkage_matrices[0], correlation_matrices[0]
        elif stimType[1] == 'Visual + Opto':
            Z, corr_matrix2 = linkage_matrices[1], correlation_matrices[1]
        else:
            raise ValueError(f"Unknown stimType: {stimType[1]}")
        
        change_matrix = corr_matrix2 - corr_matrix
    # Calculate average change for each matrix (includes all cells from all animal) and store
        animalID = animalID[:int(len(animalID)/len(correlation_matrices))]
        animalList = np.unique(animalID)
       # print(animalList)
        for animal_id in animalList:
            # Extract the correlation matrix for this animal
            animal_indices = animalID == animal_id
            corr_matrix = correlation_matrices[0][animal_indices]
            corr_matrix2 = correlation_matrices[1][animal_indices]
            change_matrix = corr_matrix2 - corr_matrix

            # Calculate the average change in correlation for this animal
            rows, cols = np.tril_indices(change_matrix.shape[0], -1)
            lower_triangle = change_matrix[rows, cols] # Extract the lower triangular part of the matrix
            average_abs_value = np.nanmean(np.abs(lower_triangle)) # Calculate the average of the absolute values
            average_changes[condition].append(average_abs_value) # Store the average change

    # Convert average_changes to a DataFrame for easier plotting
    average_changes_df = pd.DataFrame.from_dict(average_changes, orient='index').transpose()
    # Plotting the results
    
    sns.swarmplot(data=average_changes_df, ax=axs)
    #axs.set_title('Average Change in Correlation for Different Responsiveness Categories')
    axs.set_ylabel('Total average change')
    #axs.set_xlabel('Neural Population Based on Responsiveness')
    axs.set_xticklabels(responsiveness)#, rotation=45)
    axs.set_ylim(0, 0.35)

    # Save figure if requested
    if savefigname and savefigpath:
        plt.savefig(f"{savefigpath}\\{savefigname}_{stimType}_CorrelationMatrixMean.png")
    return average_changes_df

def plot_correlationMatrix_meanChange(compareType,cohort, responsiveness, params, 
                           axs=None, savefigname=None, savefigpath=None):
    if params == 'All':
        cellresponsivenessList = ['All','None']
        colorpalet = [ 'black', 'magenta']
    else:
        print('StimType order is not matched with colors')
        cellresponsivenessList = ['Visual', 'Visual + Opto', 'Opto','All']
        palette1 = sns.color_palette('cividis') # for OPTO - blue
        palette2 = sns.color_palette('viridis') # for BOTH - green
        colorpalet = ['red',palette2[4], palette1[0] , 'grey'] # all, visual, opto, both
        xlabel =   ['Visual','Opto-boosted','Opto','All']


    # Initialize a dictionary to hold average changes for each animal across conditions
    #if in compareType Train 
    if 'Trained' in compareType:
        infoPath = f'C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\analysis\\crossCorrelation_{cohort}_Trained_All.pkl'
    else:
        infoPath = f'C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\analysis\\crossCorrelation_{cohort}_Naive_All.pkl'
    data = pd.read_pickle(infoPath) 
    _, correlation_matrices, animalID = data
    animalID = animalID[:int(len(animalID)/len(correlation_matrices))]
    animalList = np.unique(animalID)
    data_list = []

    # Process both Trained and Naive data
    for condition in responsiveness:
        for group in compareType:
            infoPath = f'C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\analysis\\crossCorrelation_{cohort}_{group}_{condition}.pkl'
            data = pd.read_pickle(infoPath) 
            
            # Extract data based on index - HARD CODE: Visual - Visual+Opto
            linkage_matrices, correlation_matrices, animalID = data
            animalID = animalID[:int(len(animalID)/len(correlation_matrices))]
            change_matrix = correlation_matrices[0] - correlation_matrices[1]

            # Calculate average change for each animal
            for animal_id in animalList:
                animal_indices = animalID == animal_id
                lower_triangle = np.tril_indices_from(change_matrix[animal_indices], -1)
                average_abs_value = np.nanmean(np.abs(change_matrix[lower_triangle]))
                
                data_list.append({
                    'Condition': condition,
                    'Group': group,
                    'Average Change': average_abs_value,
                })

    # Convert list to DataFrame
    df = pd.DataFrame(data_list)

    # Create figure if not provided
    if axs is None:
        fig, axs = plt.subplots(figsize=(12, 6))

    # Plotting the results
    sns.swarmplot(data=df, x='Condition', y='Average Change', hue='Group',
                   palette=colorpalet, ax=axs, dodge=True)
    
    axs.set_ylabel('Average Change in Correlation')
    axs.set_xlabel('Condition')
    axs.legend(loc = 'lower right', fontsize = 10, frameon = False)
    axs.set_ylim(0, 0.45)

    # Optionally save the figure
    if savefigname and savefigpath:
        plt.savefig(f"{savefigpath}/{savefigname}")


## More code here