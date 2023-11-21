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

import os


def set_figure():
    from matplotlib import rcParams
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
            'axes.titlesize':'large',
            'xtick.labelsize':'small',
            'ytick.labelsize':'small',
            'legend.fontsize': 'large'}
    
    rcParams.update(params)


def plot_cell_trace(ax, analysispath,tTypes, s_recDate, s_animalID, s_recID, s_cellID, cell_title, y_min, y_max):
    pathname = analysispath + s_recDate + '_' + str(s_animalID) + '_' + str(s_recID).zfill(3) + '\\'
    cellID = s_cellID

    plot_dff_mean_traces(pathname, cellID, tTypes, ax)
    ax.set_ylim(ymin=y_min, ymax=y_max)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    ax.axis('off')
    ax.legend_.remove()
    current_lines = plt.gca().get_lines()
    for line in current_lines:
        line.set_linewidth(0.5)  # Adjust the value as needed

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
        lineplot_withSEM (data=dffTrace[t][cellID], colorInd = indx, label=t, axis = axis)

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
    palette1 = sns.color_palette('cividis') # for BOTH
    palette2 = sns.color_palette('viridis') # for OPTO
    color = ['red', palette1[0], palette2[4]]
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
    ylim_max = np.ceil (np.nanmax(np.nanmean(data,1)))*1.5
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

def heatmap_comparison(index, dff_traceVis,dff_traceBoth, 
                       dff_meanVis1sec, colormap, axis=None, cbar_ax=None, 
                       savefigname =None ,savefigpath = None, colorbarlimits = None ) :

    if colorbarlimits is None:
        yminValue = -1
        ymaxValue  = 1
    else: 
        yminValue = colorbarlimits[0]
        ymaxValue = colorbarlimits[1]

    step = 30
    set_figure()

    # HEAT PLOT FOR RESPONSIVE SENSORY CELLS
    sortedInd = np.array(dff_meanVis1sec[index]).argsort()[::-1]

    if axis is None:
        ax1 = plt.subplot(1,3,1)
    else:
        ax1 = axis[0]
        
    plot_data = dff_traceVis[index]
    plot_data = plot_data[sortedInd]
    x_labels = np.linspace(-2, 6, plot_data.shape[1], dtype = int)
    xticks = np.arange(0, len(x_labels), step)
    xticklabels = x_labels[::step]
    sns.heatmap(plot_data, vmin = yminValue, vmax = ymaxValue, cbar = False, yticklabels = False,
                     cmap = colormap,  ax = ax1)
    ax1.set_xticks (ticks = xticks, labels= xticklabels)
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=0)
    plot_dataV = plot_data
    ax1.set_xlim(30, 240)
    ax1.set_xlabel('Time (sec)')
    ax1.set_ylabel('Cells n = '  + str(np.sum(index)))
    ax1.set_title('Visual')
    if axis is None:
        ax2 = plt.subplot(1,3,2)
    else:
        ax2 = axis[1]

    plot_data = dff_traceBoth[index]
    plot_data = plot_data[sortedInd]
    sns.heatmap(plot_data, vmin = yminValue, vmax = ymaxValue, cbar = False, yticklabels = False,
                     cmap = colormap,  ax = ax2)
    ax2.set_xticks (ticks = xticks, labels= xticklabels)
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=0)
    plot_dataB = plot_data
    ax2.set_xlim(30, 240)
    ax2.set_xlabel('Time (sec)')
    ax2.set_title('Visual + Opto')

    if axis is None:
        ax3 = plt.subplot(1,3,3)
        cbar_ax = None
    else:
        ax3 = axis[2]
      
    plot_dataDIFF = plot_dataB - plot_dataV

    sns.heatmap(plot_dataDIFF, vmin = yminValue, vmax = ymaxValue, cbar = True, yticklabels = False,
                cmap = colormap,  ax = ax3,  cbar_ax=cbar_ax,   cbar_kws={'label': 'dFF'})
    ax3.set_xticklabels(ax3.get_xticklabels(), rotation=0)
    ax3.set_xticks (ticks = xticks, labels= xticklabels)
    ax3.set_xlim(30, 240)
    ax3.set_xlabel('Time (sec)')
    ax3.set_title('Opto modulation')

    if axis is not None:
        pos = cbar_ax.get_position()
        new_position = [pos.x0 - 0.01, pos.y0, pos.width, pos.height]  # Modify the x0 value to move it closer to ax2
        cbar_ax.set_position(new_position)

    # Create inset axis for colorbar
    #divider = make_axes_locatable(ax4)
    #cax = divider.append_axes("right", size="100%", pad=1)
    #plt.colorbar(ax3.get_children()[0], cax=cax, label='dFF')

    if savefigname != None:
        save_figure(savefigname,savefigpath)
    return plot_dataV, plot_dataB, plot_dataDIFF

def population_plots(index, dff_traceVis,dff_traceBoth, dff_meanVis1sec,
                       analysis_time, colormap,
                       ylimitsforhist,xlimitsforhist,
                       FaceColors,colorbarlimitsForHeatMap,scatterplotlimits,
                       xlimitsforSNR, xlimitsforCV,
                       savefigname =None ,savefigpath = None ) :
    #parameters
    fRate = 1000/30.0
    pre_frames    = 2000.0# in ms
    pre_frames    = int(np.ceil(pre_frames/fRate))
    post_frames   = 6000.0 # in ms
    post_frames   = int(np.ceil(post_frames/fRate))
    #analysis_time = 1500.0 # in ms
    analysis_time = int(np.ceil(analysis_time/fRate))
    simulationDur_ms = 350.0 # in ms
    simulationDur = int(np.ceil(simulationDur_ms/fRate))


    # Define gridspec
    set_figure()
    fig = plt.figure(figsize=(14, 5))
    gs = gridspec.GridSpec(6, 11, height_ratios=[1, 1,1,1,1,1], width_ratios=[0.4,0.4,0.4,0.1,  0.2,  0.3,0.1,  0.3,0.1,  0.4,0.4], wspace = 0.2, hspace = 2)
    # Panel A: Heatmaps
    axes = [fig.add_subplot(gs[:, i]) for i in range(3)]
    cax = fig.add_subplot(gs[:, 3])

    plot_dataV, plot_dataB, plot_dataDIFF = heatmap_comparison(index, dff_traceVis,dff_traceBoth, 
                                                                    dff_meanVis1sec, colormap, axis=axes, cbar_ax=cax,
                                                                    savefigname = savefigname, savefigpath= savefigpath, colorbarlimits = colorbarlimitsForHeatMap )

    # Panel B, C: Histograms
    # Define some common properties
    common_props = {
        'binwidth': 0.02,
    }

    titles = ['Visual', 'Visual + Opto']
    splot = [5, 7,9]
    data_means = [np.nanmean(data[:, pre_frames:(pre_frames + simulationDur + analysis_time)], axis=1) for data in [plot_dataV, plot_dataB, plot_dataDIFF]]
    for idx, (title, color, data_mean) in enumerate(zip(titles, FaceColors, data_means)):
        if idx<2:
            ax = fig.add_subplot(gs[0:3,splot[idx]])
        sns.histplot(data_mean, color=color, ax=ax, **common_props)
        ax.set_title(f'{chr(66 + idx)}', loc='left', fontweight='bold')
        ax.set_ylim(ylimitsforhist)
        ax.set_title(title)
        ax.set_xlabel('DFF')
        ax.set_xlim(xlimitsforhist)

        ax.set_ylabel('Number of cells')

    plot_dataVmean = np.nanmean(plot_dataV [:, pre_frames:(pre_frames + simulationDur + analysis_time)], axis = 1)
    plot_dataBmean = np.nanmean(plot_dataB [:, pre_frames:(pre_frames + simulationDur + analysis_time)], axis = 1)

    # Panel D: Scatter plot
    ax = fig.add_subplot(gs[0:3, 9:10])
    plot_data = pd.DataFrame( {'Visual (dFF)' :plot_dataVmean, 
                            'Opto modulation (dFF)':  plot_dataBmean-plot_dataVmean})
    sns.scatterplot (y = 'Visual (dFF)', x = 'Opto modulation (dFF)', data = plot_data, color='black', ax=ax, linewidth = 0.5, markers='.', s=7)
    plt.axhline(y = 0, color = 'black', linestyle = '--', linewidth = 0.2)
    plt.axvline(color = 'black', linestyle = '--', linewidth = 0.2)
    plt.ylim(scatterplotlimits)
    plt.xlim(scatterplotlimits)
    ax.set_title('D', loc='left', fontweight='bold')
    
    # INSET IN Panel D:  Violin plots
    xsmall = ['V', 'V + O']
    inset_ax = inset_axes(ax, width="10%", height="25%", loc="upper right")
    s_index = np.where((10>plot_dataVmean) &(plot_dataVmean>0))[0] # This 10 does not mean anything right now. 
    plot_data = pd.DataFrame( {'Mean DFF' :np.concatenate((plot_dataVmean[s_index], plot_dataBmean[s_index])), 
                                'Type':  np.concatenate((np.repeat('Visual', len(plot_dataVmean[s_index])), np.repeat('Visual+Opto', len(plot_dataBmean[s_index]))))})

    sns.barplot(x = 'Type', y = 'Mean DFF', data = plot_data, palette=FaceColors, ax=inset_ax, linewidth = 0.1)
    plt.axhline(y = 0, color = 'black', linestyle = '--', linewidth = 0.1)
    plt.annotate('p = {:.5f}'.format(stats.ttest_rel(plot_dataVmean[s_index], plot_dataBmean[s_index], alternative = 'greater')[1], 3), 
                xy=(0.5, 0.7), xytext=(0, 0), textcoords='offset points', 
                ha = 'center', va = 'top', fontsize = 5)
    inset_ax.set_xticklabels(xsmall,rotation=45)
    inset_ax.set_ylabel('')
    inset_ax.set_xlabel('')
    plt.ylim(-0.05,0.7)
    for item in ([inset_ax.title, inset_ax.xaxis.label, inset_ax.yaxis.label] +
             inset_ax.get_xticklabels() + inset_ax.get_yticklabels()):
        item.set_fontsize(5)

    inset_ax = inset_axes(ax, width="10%", height="25%", loc="lower left",
                          bbox_to_anchor=(0.1, 0.1, 1, 1),
                            bbox_transform=ax.transAxes)
    s_index = np.where((0>plot_dataVmean) &(plot_dataVmean>-10))[0] # This 10 does not mean anything right now. 
    plot_data= pd.DataFrame( {'Mean DFF' :np.concatenate((plot_dataVmean[s_index], plot_dataBmean[s_index])), 
                                'Type':  np.concatenate((np.repeat('Visual', len(plot_dataVmean[s_index])), np.repeat('Visual+Opto', len(plot_dataBmean[s_index]))))})

    sns.barplot(x = 'Type', y = 'Mean DFF', data = plot_data, palette=FaceColors, ax=inset_ax,  linewidth = 0.1)
    plt.axhline(y = 0, color = 'black', linestyle = '--',  linewidth = 0.1)
    plt.annotate('p = {:.5f}'.format(stats.wilcoxon(plot_dataVmean[s_index], plot_dataBmean[s_index], alternative = 'less')[1], 3),  
                xy=(0.5, 0.1), xytext=(0, 0), textcoords='offset points', 
                ha = 'center', va = 'top',fontsize = 5)
    plt.ylim(-0.4,0.1)
    inset_ax.set_ylabel('')
    inset_ax.set_xlabel('')
    inset_ax.set_xticklabels(xsmall,rotation=45)
    for item in ([inset_ax.title, inset_ax.xaxis.label, inset_ax.yaxis.label] +
            inset_ax.get_xticklabels() + inset_ax.get_yticklabels()):
        item.set_fontsize(5)


    # Panel E: Cumulative ECDF
    ax7 = fig.add_subplot(gs[3:6, 5:6])
    responses_before = np.nanmean(dff_traceVis[index, pre_frames:(pre_frames + simulationDur + analysis_time)], axis=1)
    x, y = mfun.ecdf(responses_before)
    plt.plot(x, y, marker='.', markersize=1, color= FaceColors[0], linestyle='none',label='Visual', axes=ax7)

    responses_after =  np.nanmean(dff_traceBoth[index, pre_frames:(pre_frames + simulationDur + analysis_time)], axis=1)
    x, y = mfun.ecdf(responses_after)
    plt.plot(x, y, marker='.', markersize=1, linestyle='none', color = FaceColors[1], label='Visual + Opto', axes=ax7)
    plt.annotate('p = {:.5f}'.format(stats.wilcoxon(responses_before, responses_after)[1], 3),  
            xy=(0.1, 0.45), xytext=(0, 0), textcoords='offset points', 
            ha = 'center', va = 'top')

    plt.xlabel('Neural response (dFF)')
    plt.ylabel('ECDF')
    if colormap is None:
        plt.ylim(0.6, 1.1)
    plt.ylim(0.4, 1.1)
    #plt.grid(True)
    plt.xscale('log') 
    plt.legend(fontsize='x-small',loc='upper right', frameon=False)
    ax7.set_title('E', loc='left', fontweight='bold')
    

    # Panel F: SNR
    infoPath = 'C:\\Users\\Huriye\\Documents\\code\\clapfcstimulation\\analysis\\infoForAnalysis-readyForPlotting_moreStats.pkl'
    varianceBoth_pre, varianceVis_pre, varianceOpto_pre, varianceBoth_post, varianceVis_post, varianceOpto_post, snrBoth,snrVis,snrOpto,ccBoth,ccVis,ccOpto  = pd.read_pickle(infoPath) 
    datasets = [dff_traceVis, dff_traceBoth]
    labels = ['Visual', 'Both']

    # Calculate Coefficient of Variation Differences
    indices = np.where(index == True)[0]
    cv_values = []
    snr_values = []
    for data, label in zip(datasets, labels):
        
        if label == 'Visual':  # Get the first element of the tuple, which is the array of indices
            snr = [snrVis[i] for i in indices]
            varianceDiff = [post - pre for post, pre in zip(varianceVis_post, varianceVis_pre)]
            variance = [varianceDiff[i] for i in indices]
        elif label == 'Both':
            snr = [snrBoth[i] for i in indices]
            varianceDiff = [post - pre for post, pre in zip(varianceBoth_post, varianceBoth_pre)]
            variance = [varianceDiff[i] for i in indices]
    
        for value in snr:  # Assuming snr is a numpy array or list
            snr_values.append({'Condition': label, 'SNR': value})
        for value in variance:  # Assuming snr is a numpy array or list
            cv_values.append({'Condition': label, 'CV': value})
    
    ax8 = fig.add_subplot(gs[3:6, 7:8])
    snr_df = pd.DataFrame(snr_values)
    sns.histplot(data=snr_df, x="SNR", hue="Condition", element="step",
                 palette=FaceColors, ax=ax8, **common_props)
    ax8.set_xlabel('SNR')
    ax8.set_xlim(xlimitsforSNR)
    ax8.legend_.remove()

    
    ax8 = fig.add_subplot(gs[3:6, 9:10])
    cv_df = pd.DataFrame(cv_values)
    sns.histplot(data=cv_df, x="CV", hue="Condition", element="step",
                 palette=FaceColors, ax=ax8, **common_props)
    ax8.set_xlabel('CV differences\nlog(postCV)- log(preCV)')
    ax8.set_xlim(xlimitsforCV)
    ax8.legend_.remove()
    
    # save figure
    save_figure(savefigname,savefigpath)