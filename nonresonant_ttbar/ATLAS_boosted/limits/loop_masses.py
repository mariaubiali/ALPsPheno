import numpy as np
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns

font = {'size'   : 16}
import matplotlib
matplotlib.rc('font', **font)

masses = ['1em3','1em2','1em1','1','10','100','200','300','400','500','600','700','710','800','1000','1075','1225','1400','1600','1850','2150','2900', '3200', '3500','4000','5000','6000','7000']
massvals = np.array([0.001,0.01,0.1,1,10,100,200,300,400,500,600,700,710,800,1000,1075,1225,1400,1600,1850,2150,2900,3200,3500,4000,5000,6000,7000])


#Loop through masses and get limit on c at 95% CL
lim=[]
clim = []
for mass in masses:
    cmd = 'python get_constraints.py '+mass
    out = subprocess.check_output(cmd,shell=True,encoding='UTF-8')
    clim.append(float(out))
    lim.append(1.0/float(out))

clim = np.array(clim,dtype=float)

np.savetxt('ATLAS_ttbar_nonresonant_update.txt',np.c_[massvals, clim])

#Plot limit vs mass
fig,ax = plt.subplots(figsize=(7,7))
plt.plot(massvals,lim,'k+',ms = 10)
ax.set_ylabel(r'$|f_a/c_t|$ [TeV]')
ax.set_xlabel(r'$m_a$ [GeV]')
ax.set_title(r'$|\frac{f_a}{c_t}|$ vs $m_a$')
ax.set_xscale('log')
ax.set_ylim([1e-2,2e3])
ax.set_yscale('log')
plt.tight_layout()
plt.savefig('constraints_vs_ma_ttbar_update.pdf')
plt.show()
