import numpy as np
import math


data = np.loadtxt('../data/ATLAS_data_proc.txt',dtype=float,usecols=(2,))

dtf1 = np.loadtxt('datathief_atlasNR_option1.txt',dtype=float,usecols=(1,))
dtf2 = np.loadtxt('datathief_atlasNR_option2.txt',dtype=float,usecols=(1,))

nnlo2 = dtf2*data
print('nnlo2 predictions',nnlo2)

nnlo1 = np.divide(data,dtf1)
print('nnlo2 predictions',nnlo1)


np.savetxt('nnlo_from_fig19.txt',nnlo1)
np.savetxt('nnlo_from_fig11.txt',nnlo2)
