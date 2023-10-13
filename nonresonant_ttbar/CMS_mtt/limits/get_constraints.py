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


def csv_reader(filename):
    output = []
    with open(filename, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            output.append(row)
        csvfile.close()

    return output


def read_HEPdata_SM():
    '''
    Read data and covmat from hepdata
    '''
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
    '''
    axion_lin = np.divide(BR*np.loadtxt('../data_from_madgraph/pth_100723_lin_ma'+mass_alp+'.txt',usecols=(2,),dtype=float),bins_width)
    axion_quad = np.divide(BR*np.loadtxt('../data_from_madgraph/pth_100723_quad_ma'+mass_alp+'.txt',usecols=(2,),dtype=float),bins_width)
    sm = np.divide(BR*np.loadtxt('../data_from_madgraph/mtt_SM_ttbar_nnpdf4p0.txt',usecols=(0,),dtype=float), bins_width)
    
    axion_signal = np.c_[kfac*sm, kfac*axion_lin, kfac*axion_quad]
    return axion_signal
    

def kfactor(highstats=False):
    '''
    Get SM NNLO kfactor from hightea
    '''
    filename='../sm/kfac_nnlo_lo_highstats.txt'
    kfac = np.loadtxt(filename,dtype=float,usecols=(0,))
    kfac_unc = np.loadtxt(filename,dtype=float,usecols=(1,))

    return kfac, kfac_unc




def chi2(ct, axion_signal, data, covmat):
    sm = axion_signal[:,0]
    quad = axion_signal[:,1]
    quart = axion_signal[:,2]
    theory = sm + quad*ct**2 + quart*ct**4
    diff = (theory - data)[1:]
    Vinv = np.linalg.inv(covmat)[1:,1:]
    return ((diff).dot(Vinv)).dot(diff)

  



if __name__=='__main__':
    
    BR = 0.2877
    index=0
    indexf=-1

    kfac, kfac_unc = kfactor(highstats=True)

    #Get CMS data
    bins_cms=np.array([250.0,400.0,480.0,560.0,640.0,720.0,800.0,900.0,1000.0,1150.0,1300.0,1500.0,1700.0,2000.0,2300,3500])
    bins_width = bins_cms[1:]-bins_cms[:-1]
    
    cms_data, covmat = read_HEPdata_SM();

    axion_signal = axion_signal(1, bins_cms, BR, kfac, kfac_unc)

    #First find minima of the chi profile, such that the delta chi2 can then be calculated
    def func_to_solve_deltachi2(ct):
        return chi2(ct, axion_signal[index:indexf], cms_data[index:indexf], covmat[index:indexf,index:indexf])

    min1, min2 = minimize(func_to_solve_deltachi2, x0=-20).x, minimize(func_to_solve_deltachi2, x0=20).x
    chi2min = chi2(min1, axion_signal[index:indexf], cms_data[index:indexf], covmat[index:indexf,index:indexf])

    plt.figure()
    
    cvals = np.arange(-30,30,1)
    chivals = []
    for ct in cvals:
        chivals.append(chi2(ct, axion_signal[index:indexf], cms_data[index:indexf], covmat[index:indexf,index:indexf]) - chi2min)

    def func_to_solve_68(ct):
        return chi2(ct, axion_signal[index:indexf], cms_data[index:indexf], covmat[index:indexf,index:indexf]) - chi2min - 0.99

    def func_to_solve_95(ct):
        return chi2(ct, axion_signal[index:indexf], cms_data[index:indexf], covmat[index:indexf,index:indexf]) - chi2min - 3.84

    c68a, c68b = brentq(func_to_solve_68, a=-1000,b=min1), brentq(func_to_solve_68, a=min2,b=1000.0)
    c95a, c95b = brentq(func_to_solve_95, a=-1000,b=min1), brentq(func_to_solve_95, a=min2,b=1000.0)

    print(c95b)

