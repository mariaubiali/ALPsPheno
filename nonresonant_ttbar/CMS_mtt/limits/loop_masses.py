import numpy as np
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns

font = {'size'   : 16}
import matplotlib
matplotlib.rc('font', **font)


masses = ['1em3','1em2','1em1','1','10','100','150','200','250','300','350','400','410','420','440','450','500','550','600','650','700','710','750','800','810','820','830','840','850','900','950','1000','1010','1075','1225','1400','1600','1700','1850','2150','2900', '3200', '3500','4000','5000','6000','7000','8000','9000','10000']
massvals = np.array([0.001,0.01,0.1,1,10,100,150,200,'250','300','350','400','410','420','440','450','500','550','600','650','700','710','750','800','810','820','830','840','850','900','950','1000','1010','1075','1225','1400','1600','1700','1850','2150','2900', '3200', '3500','4000','5000','6000','7000','8000','9000','10000'],dtype=float)

start=0
end=-1

lim=[]
clim = []
for mass in masses:
    cmd = 'python get_constraints.py '+mass
    out = subprocess.check_output(cmd,shell=True,encoding='UTF-8')
    clim.append(float(out))
    lim.append(1.0/float(out))

clim = np.array(clim,dtype=float)

np.savetxt('CMS_ttbar_nonresonant.txt',np.c_[massvals*1e9, clim*1e-3])

fig,ax = plt.subplots(figsize=(7,7))
plt.plot(massvals,lim,'ko',ms = 2)
ax.set_ylabel(r'$|f_a/c_t|$ [TeV]')
ax.set_xlabel(r'$m_a$ [GeV]')
ax.set_title(r'$|\frac{f_a}{c_t}|$ vs $m_a$')
ax.set_xscale('log')
ax.set_yscale('log')
plt.tight_layout()
plt.savefig('constraints_vs_ma_ttbar.pdf')
plt.show()
