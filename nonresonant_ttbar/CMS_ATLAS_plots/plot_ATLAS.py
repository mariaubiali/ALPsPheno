import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
colors = sns.color_palette("colorblind", 20)
import math
import csv 


font = {'size'   : 14}
import matplotlib
matplotlib.rc('font', **font)

def axion_signal(ct, bins, BR, kfac, kfac_unc):
    '''
    Calculate SM + C^2*sigma + C^4*sigma
    and uncertainties
    '''
    bins_width = bins[1:]-bins[:-1]
    axion_lin = np.divide(BR*np.loadtxt('./Madgraph/pth_nnpdf4p0_alphaS_lin.txt',usecols=(2,),dtype=float),bins_width)
    axion_quad = np.divide(BR*np.loadtxt('./Madgraph/pth_nnpdf4p0_alphaS_quad.txt',usecols=(2,),dtype=float),bins_width)
    sm = np.divide(BR*np.loadtxt('./Madgraph/pth_nnpdf4p0_alphaS_SM.txt',usecols=(2,),dtype=float), bins_width)


    signal = kfac*(sm + axion_lin*ct**2 + axion_quad*ct**4)
    signal_unc = np.zeros(np.size(kfac))
    return signal, signal_unc


if __name__=='__main__':
    
    #ctvals = np.array([0.0,1.0,4.0, 12.0, 20.0])
    ctvals = np.array([0.0, 5.0, 10.0, 13.0, 20.0])
    ctvals = np.array([0.0, 4.8, 7.8, 15.0])
    ctvals = np.array([0.0, 4.0, 12.0, 15.0])
    BR = 0.2877

    #Bin info
    bins_atlas = np.array([355.0, 381.0, 420.0, 478.0, 549.0, 633.0, 720.0, 836.0, 2000.0])
    bins_width = bins_atlas[1:]-bins_atlas[:-1]
    bins_atlas_C = 0.5*(bins_atlas[1:] + bins_atlas[:-1])
    
    #Create digitised kfactor for SM predictions
    sm_lo = np.divide(BR*np.loadtxt('Madgraph/pth_nnpdf4p0_alphaS_SM.txt',usecols=(2,),dtype=float), bins_width)
    sm_nnlo = np.loadtxt('./datathief_sm/nnlo_from_fig19.txt',dtype=float,usecols=(0,))
    kfac_partial = np.divide(sm_nnlo, sm_lo)
    kfac = kfac_partial
    kfac_unc = np.zeros(np.size(kfac))

    #Load preprocessed ATLAS data
    data_atlas = np.loadtxt('./data/ATLAS_data_proc.txt', dtype=float, usecols=(2,))
    data_atlas_unc = np.loadtxt('./data/ATLAS_data_proc.txt', dtype=float, usecols=(3,))
    #create data arrays for plot
    data_atlas_plot = np.concatenate([data_atlas, [data_atlas[-1]]])
    data_atlas_unc_plot = np.concatenate([data_atlas_unc, [data_atlas_unc[-1]]])


    fig = plt.figure(figsize=(7, 7))
    gs = fig.add_gridspec(2, 1, height_ratios=(4, 1),
                      left=0.14, right=0.86, bottom=0.1, top=0.9,
                      wspace=0.55, hspace=0.25)

    # Create the Axes.
    ax = fig.add_subplot(gs[0, 0])
    ax_ratio = fig.add_subplot(gs[1, 0], sharex=ax)

    #Plot ATLAS data
    ax.step(bins_atlas, data_atlas_plot, where='post',label='ATLAS data', color=colors[0])
    ax.fill_between(bins_atlas, data_atlas_plot-data_atlas_unc_plot, data_atlas_plot+data_atlas_unc_plot, color=colors[0],step='post', alpha=0.2)
    ax.set_title('ATLAS data compared to SM+ALP signal')
    ax.set_yscale('log')
    ax.set_xlim([355,2000]);
    ax.set_xlabel(r'$p_{T}(t_{h})$ [GeV]')
    ax.set_ylabel(r'$\frac{d \sigma}{d p_{T}(t_{h})}$ [pb/GeV]')


    ax_ratio.set_xlabel(r'$p_{T}(t_{h})$ [GeV]')
    ax_ratio.set_ylabel('Ratio to data')
    ax_ratio.hlines(xmin=355, xmax=2000, y=1.0,color=colors[0])
    data_atlas_ratio_unc_plot = np.divide(data_atlas_unc_plot,data_atlas_plot)
    ax_ratio.fill_between(bins_atlas, 1-data_atlas_ratio_unc_plot, 1+data_atlas_ratio_unc_plot, color=colors[0],step='post', alpha=0.3)

    
    for ind in np.arange(1,np.size(ctvals)+1):
        axion_data, axion_data_unc = axion_signal(ctvals[ind-1], bins_atlas, BR, kfac, kfac_unc)

        #Create arrays for step plot
        axion_data_plot = np.concatenate([axion_data, [axion_data[-1]]])
        axion_data_unc_plot = np.concatenate([axion_data_unc, [axion_data_unc[-1]]])

        #Create ratio arrays
        axion_data_ratio = np.divide(axion_data_plot,data_atlas_plot)
        ax.step(bins_atlas, axion_data_plot, where='post',label=r'SM + ALP, $c_{t}/f_a = $'+str(int(ctvals[ind-1]))+r' TeV$^{-1}$', color=colors[ind])
        ax_ratio.step(bins_atlas, axion_data_ratio, where='post',color=colors[ind])

    
    
    ax.legend();
    plt.tight_layout()
    plt.savefig('axion_data_ATLAS_ttbar.pdf')
    plt.show()


