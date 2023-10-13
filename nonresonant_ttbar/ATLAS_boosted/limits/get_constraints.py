import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
colors = sns.color_palette("colorblind", 8)
import math
import csv 
from scipy.optimize import minimize, brentq
import sys

font = {'size'   : 12}
import matplotlib
matplotlib.rc('font', **font)

mass_alp = str(sys.argv[1])


def axion_signal(ct, bins, BR, kfac):
    '''
    Calculate SM + C^2*sigma + C^4*sigma
    and uncertainties
    '''
    bins_width = bins[1:]-bins[:-1]
    axion_lin = np.divide(BR*np.loadtxt('../data_from_madgraph/mass_scan_061023/pth_061023_lin_ma'+mass_alp+'.txt',usecols=(2,),dtype=float),bins_width)
    axion_quad = np.divide(BR*np.loadtxt('../data_from_madgraph/mass_scan_061023/pth_061023_quad_ma'+mass_alp+'.txt',usecols=(2,),dtype=float),bins_width)
    sm = np.divide(BR*np.loadtxt('../data_from_madgraph/pth_nnpdf4p0_alphaS_SM.txt',usecols=(2,),dtype=float), bins_width)
    signal = np.c_[kfac*sm, kfac*axion_lin, kfac*axion_quad]
    return signal


def chi2(ct, axion_signal, data, covmat):
    '''
    Compure the chi2 as a function of ct
    '''
    sm = axion_signal[:,0]
    quad = axion_signal[:,1]
    quart = axion_signal[:,2]
    theory = sm + quad*ct**2 + quart*ct**4
    diff = (theory - data)
    Vinv = np.linalg.inv(covmat)
    return ((diff).dot(Vinv)).dot(diff)

  

if __name__=='__main__':
    
    BR = 0.2877
    indexi = 0


    #Get ATLAS data
    bins_atlas = np.array([355.0, 381.0, 420.0, 478.0, 549.0, 633.0, 720.0, 836.0, 2000.0])
    bins_width = bins_atlas[1:]-bins_atlas[:-1]
 
    #Read data
    data_atlas = np.loadtxt('../data/ATLAS_data_proc.txt', dtype=float, usecols=(2,))
    unc_atlas = np.loadtxt('../data/ATLAS_data_proc.txt', dtype=float, usecols=(3,))
    covmat_stat = np.loadtxt('../data/ATLAS_data_covmat.txt', dtype=float)

    #Add systematic uncertainties to covariance matrix (in quadrature)
    covmat = covmat_stat - np.diag(covmat_stat.diagonal()) + np.diag(unc_atlas**2)

    sm_lo = np.divide(BR*np.loadtxt('../data_from_madgraph/pth_nnpdf4p0_alphaS_SM.txt',usecols=(2,),dtype=float), bins_width)
    sm_nnlo = np.loadtxt('../digitised_sm/nnlo_from_fig11.txt',dtype=float,usecols=(0,))
    kfac = np.divide(sm_nnlo, sm_lo)


    #Extract sm, quadratic and quartic pieces of alp signal by setting ct=1
    axion_signal = axion_signal(1, bins_atlas, BR, kfac)
 
    def func_to_solve_deltachi2(ct):
        return chi2(ct, axion_signal[indexi:], data_atlas[indexi:], covmat[indexi:,indexi:])

    
    #Find the chi2 minimum
    min1, min2 = minimize(func_to_solve_deltachi2, x0=-100).x, minimize(func_to_solve_deltachi2, x0=100).x #minima obtained by minimising from either direction
    chi2min = chi2(min1, axion_signal[indexi:], data_atlas[indexi:], covmat[indexi:,indexi:])



    def func_to_solve_68(ct):
        return chi2(ct, axion_signal[indexi:], data_atlas[indexi:], covmat[indexi:,indexi:]) - chi2min - 0.99

    def func_to_solve_95(ct):
        return chi2(ct, axion_signal[indexi:], data_atlas[indexi:], covmat[indexi:,indexi:]) - chi2min - 3.84

    c68a, c68b = brentq(func_to_solve_68, a=-500,b=min1), brentq(func_to_solve_68, a=min2,b=500.0)
    c95a, c95b = brentq(func_to_solve_95, a=-500,b=min1), brentq(func_to_solve_95, a=min2,b=500.0)

    if(abs(min1-min2)/c95a > 1e-1):
        raise ValueError('Possible double minimum structure in chi2, check limits')

    print(c95b)

    #for plotting:
    cvals = np.arange(-100,100,0.01)
    chivals = []
    for ct in cvals:
        chivals.append(func_to_solve_deltachi2(ct) - chi2min)

    plt.plot(cvals, chivals, 'b')

    plt.xlabel(r'$c_{t}$ [TeV$^{-1}$]')
    plt.ylabel(r'$\Delta \chi^{2}$')
    plt.hlines(y=4, xmin=-10,xmax=10,color='grey')
    plt.hlines(y=1, xmin=-10,xmax=10,color='grey')
    plt.vlines(x=c68a, ymin=-10, ymax = 10, color='grey')
    plt.vlines(x=c68b, ymin=-10, ymax = 10, color='grey')
    plt.vlines(x=c95a, ymin=-10, ymax = 10, color='lightgrey')
    plt.vlines(x=c95b, ymin=-10, ymax = 10, color='lightgrey')
    plt.xlim([-10,10])
    plt.ylim([-0.2, 4.5])
    plt.savefig('Chi2.png', dpi=256)


