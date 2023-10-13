import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
colors = sns.color_palette("colorblind", 8)
import math
import csv 

import sys
mass_alp = sys.argv[1]

font = {'size'   : 13}
import matplotlib
matplotlib.rc('font', **font)


def csv_reader(filename):
    output = []
    with open(filename, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            output.append(row)
        csvfile.close()

    return output


def read_HEPdata_SM():
    data = csv_reader('../data/HEPData-ins1901295-v1-parton_abs_ttm.csv')
    cms_data  = []
    for item in data[9:24]:
        cms_data.append(float(item[3]))

    covdata = csv_reader('../data/HEPData-ins1901295-v1-parton_abs_ttm_covariance.csv')
    covmat = np.zeros(15*15).reshape(15,15)

    covmatlist = []
    count=0
    for item in covdata[9:234]:
        covmatlist.append(float(item[6]))

    for i in range(15):
        for j in range(15):
            covmat[i,j] = covmatlist[count]
            count+=1

    return cms_data, covmat



def axion_signal(ct, bins, BR, kfac, kfac_unc):
    '''
    Calculate SM + C^2*sigma + C^4*sigma
    and uncertainties
    '''
    axion_lin = np.divide(BR*np.loadtxt('../data_from_madgraph/pth_100723_lin_ma'+mass_alp+'.txt',usecols=(2,),dtype=float),bins_width)
    axion_quad = np.divide(BR*np.loadtxt('../data_from_madgraph/pth_100723_quad_ma'+mass_alp+'.txt',usecols=(2,),dtype=float),bins_width)
    sm = np.divide(BR*np.loadtxt('../data_from_madgraph/mtt_SM_ttbar_nnpdf4p0.txt',usecols=(0,),dtype=float), bins_width)
    signal = kfac*(sm + axion_lin*ct**2 + axion_quad*ct**4)
    signal_unc = sm*np.sqrt((0.01)**2 + np.divide(kfac_unc, kfac)**2)

    
    return signal, signal_unc
    




  







if __name__=='__main__':
    
    ctvals = [0.0,4.0,12.0,20.0]
    BR = 0.2877
    kfac = np.loadtxt('../sm/kfac_nnlo_lo_highstats.txt', dtype=float, usecols=(0,))
    kfac_unc = np.loadtxt('../sm/kfac_nnlo_lo_highstats.txt', dtype=float, usecols=(1,))


    #Get CMS data
    bins_cms=np.array([250.0,400.0,480.0,560.0,640.0,720.0,800.0,900.0,1000.0,1150.0,1300.0,1500.0,1700.0,2000.0,2300,3500])
    bins_cms_C = 0.5*(bins_cms[:-1] + bins_cms[1:])
    bins_width = bins_cms[1:]-bins_cms[:-1]
    cms_data, cms_covmat = read_HEPdata_SM();
    cms_unc = list(np.sqrt(np.diag(cms_covmat)))

    #Create arrays for step plot
    cms_data_plot = np.concatenate([cms_data, [cms_data[-1]]])
    cms_unc_plot = np.concatenate([cms_unc, [cms_unc[-1]]])


    fig = plt.figure(figsize=(7, 7))
    gs = fig.add_gridspec(2, 1, height_ratios=(4, 1),
                      left=0.14, right=0.86, bottom=0.1, top=0.9,
                      wspace=0.55, hspace=0.25)
    
    # Create the Axes.
    ax = fig.add_subplot(gs[0, 0])
    ax_ratio = fig.add_subplot(gs[1, 0], sharex=ax)

    ax.step(bins_cms, cms_data_plot, where='post',label='CMS data', color=colors[0])
    ax.fill_between(bins_cms, cms_data_plot-cms_unc_plot, cms_data_plot+cms_unc_plot, color=colors[0],step='post', alpha=0.2)
    ax.set_title('CMS data compared to SM+ALP signal')
    ax.set_yscale('log')
    ax.set_xlim([355,3500]);
    ax.set_xlabel(r'$m_{t \bar{t}}$ [GeV]')
    ax.set_ylabel(r'$\frac{d \sigma}{d m_{t \bar{t}}}$ [pb/GeV]')


    ax_ratio.set_xlabel(r'$m_{t \bar{t}}$ [GeV]')
    ax_ratio.set_ylabel('Ratio to data')
    ax_ratio.hlines(xmin=355, xmax=4000, y=1.0,color=colors[0])


    
    for ind in range(len(ctvals)):
        axion_data, axion_data_unc = axion_signal(ctvals[ind], bins_cms, BR, kfac, kfac_unc)

        #Create arrays for step plot
        axion_data_plot = np.concatenate([axion_data, [axion_data[-1]]])
        axion_data_unc_plot = np.concatenate([axion_data_unc, [axion_data_unc[-1]]])

        #Create ratio arrays
        axion_data_ratio = np.divide(axion_data_plot,cms_data_plot)
        if(ind==0):
            ratio_sm = axion_data_ratio
        data_ratio_unc = np.divide(cms_unc_plot,cms_data_plot)
        axion_data_ratio_unc = axion_data_ratio*np.sqrt(data_ratio_unc**2 + np.divide(axion_data_unc_plot, axion_data_plot)**2)

        ax.step(bins_cms, axion_data_plot, where='post',label=r'SM + ALP, $c_{t}/f_a = $'+str(int(ctvals[ind]))+r' TeV$^{-1}$', color=colors[ind+1])
        ax_ratio.step(bins_cms, axion_data_ratio, where='post',color=colors[ind+1])
    ax_ratio.fill_between(bins_cms, 1-data_ratio_unc, 1+data_ratio_unc, color=colors[0],step='post', alpha=0.4)
    ax.legend();
    plt.tight_layout()
    plt.savefig('axion_data_CMS_ttbar.pdf')
    plt.show()


