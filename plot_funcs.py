# This code  has plotting functions
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import os


def set_figure():
    from matplotlib import rcParams
        # set the plotting values
    rcParams['figure.figsize'] = [5, 5]
    rcParams['font.size'] = 12
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']

    rcParams['axes.spines.right']  = False
    rcParams['axes.spines.top']    = False
    rcParams['axes.spines.left']   = True
    rcParams['axes.spines.bottom'] = True

    params = {'axes.labelsize': 'large',
            'axes.titlesize':'large',
            'xtick.labelsize':'large',
            'ytick.labelsize':'large',
            'legend.fontsize': 'large'}
    
    rcParams.update(params)

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
    
def lineplot_withSEM (data, colorInd, label):
    #lineplot_matrix(data=pupil_arr[session.outcome=='hit'], x_axis=x_axis, color=COLORS[0], label='hit')
    x_axis = np.linspace(-2, 6, data.shape[0])
    color =  sns.color_palette("Paired")
    color  = [ color[5], (0,0,0),color[2]]
    df = pd.DataFrame(data).melt()
    df['Time (seconds)'] = np.tile(x_axis, data.shape[1])

    ax = sns.lineplot(x='Time (seconds)', y='value', data=df, color=color[colorInd],
                label=label )
    #ylim_min = np.floor(np.min(np.nanmean(data,1)))*1.5
    #ylim_max = np.ceil(np.max(np.nanmean(data,1)))*1.5
    ylim_min = np.floor(np.nanmin(np.nanmean(data,1)))*1.1
    ylim_max = np.ceil (np.nanmax(np.nanmean(data,1)))*1.1
    #if np.isnan(ylim_min): ylim_min = -1
    #if np.isnan(ylim_max): ylim_max = 5
    ylength = np.absolute(ylim_max - ylim_min)
    xlength = 0.25
    #add rectangle to plot
 
    ax.add_patch (Rectangle ((0, ylim_min), xlength, ylength, alpha = 1, facecolor="grey",zorder=10))
    plt.ylabel('DFF')
    #ax.set_ylim( ymin =ylim_min, ymax = ylim_max)
    
def save_figure(name, base_path):
    plt.savefig(os.path.join(base_path, f'{name}.png'), 
                bbox_inches='tight', transparent=False)
   # plt.savefig(os.path.join(base_path, f'{name}.svg'), 
   #             bbox_inches='tight', transparent=True)

def heatmap_comparison(index, dff_traceVis,dff_traceBoth, 
                       dff_meanVis1sec, colormap, savefigname, savefigpath ) :

    yminValue = -1
    ymaxValue = 1
    step = 30


    # HEAT PLOT FOR RESPONSIVE SENSORY CELLS
    sortedInd = np.array(dff_meanVis1sec[index]).argsort()

    plt.subplot(1,3,1)
    plot_data = dff_traceVis[index]
    plot_data = plot_data[sortedInd]
    x_labels = np.linspace(-2, 6, plot_data.shape[1], dtype = int)
    xticks = np.arange(0, len(x_labels), step)
    xticklabels = x_labels[::step]
    ax = sns.heatmap(plot_data, vmin = yminValue, vmax = ymaxValue, cbar = False, yticklabels = False,cmap = colormap)
    plt.xticks (ticks = xticks, labels= xticklabels)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    plot_dataV = plot_data
    plt.xlim(30,240)
    plt.xlabel('Time (sec)')
    plt.ylabel('Cells')
    plt.title('Visual')

    plt.subplot(1,3,2)
    plot_data = dff_traceBoth[index]
    plot_data = plot_data[sortedInd]
    ax = sns.heatmap(plot_data, vmin = yminValue, vmax = ymaxValue, cbar = False, yticklabels = False,cmap = colormap)
    plt.xticks (ticks = xticks, labels= xticklabels)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    plot_dataB = plot_data
    plt.xlim(30,240)
    plt.xlabel('Time (sec)')
    plt.title('Visual + Opto')

    plt.subplot(1,3,3)
    plot_data = plot_dataB - plot_dataV

    ax = sns.heatmap(plot_data, vmin = yminValue, vmax = ymaxValue, cbar = True, yticklabels = False,cmap = colormap,
                    cbar_kws={'label': 'dFF'})
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    plt.xticks (ticks = xticks, labels= xticklabels)
    plot_diff = plot_data
    plt.xlim(30,240)
    plt.xlabel('Time (sec)')
    plt.title('Opto modulations')

    save_figure(savefigpath,savefigpath)