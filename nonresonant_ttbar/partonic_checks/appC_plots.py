import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

colors = sns.color_palette("colorblind",8)
color1 = colors[0]
color2 = colors[4]

#Load Quad stuff
roots = np.loadtxt('madgraph_cards/ttbar_gg_NR_quad_partonic/mtt_s_partonic_meframe_130323_quad.txt',dtype=float,usecols=(0,))/1000
sigma = np.loadtxt('madgraph_cards/ttbar_gg_NR_quad_partonic/mtt_s_partonic_meframe_130323_quad.txt',dtype=float,usecols=(1,))

#Load Lin stuff
rootsLin = np.loadtxt('madgraph_cards/ttbar_gg_NR_linear_partonic/mtt_s_partonic_meframe_130323_lin.txt',dtype=float,usecols=(0,))/1000
sigmaLin = np.loadtxt('madgraph_cards/ttbar_gg_NR_linear_partonic/mtt_s_partonic_meframe_130323_lin.txt',dtype=float,usecols=(1,))

#Load SM stuff
rootsSM = np.loadtxt('madgraph_cards/ttbar_gg_NR_SM_partonic/mtt_s_partonic_meframe_130323_SM.txt',dtype=float,usecols=(0,))/1000
sigmaSM = np.loadtxt('madgraph_cards/ttbar_gg_NR_SM_partonic/mtt_s_partonic_meframe_130323_SM.txt',dtype=float,usecols=(1,))




#Analytical expressions; constant of proportionality from mathematica notebook ALPSMpartonic.nb
def functionFit(x):
    return 5.36266e-6*(1 - 59512.5/x**2)

func = []
for sv in roots:
    func.append(functionFit(1000*sv))


def functionFitLin(x):
    return (1088.01*np.log(0.0057971*x))/(x**2)

funcLin = []
for sv in rootsLin:
    funcLin.append(functionFitLin(1000*sv))



def functionFitSM(x):
    return (4.32896*10**6*np.log(0.0057971*x))/(x**2)


funcSM = []
for sv in rootsSM:
    funcSM.append(functionFitSM(1000*sv))

#Quad only plot
plt.figure()
plt.plot(roots,func,color=color2,label=r'$\hat{\sigma}_{\mathrm{ALP}} \sim 1-2m_t^2/\hat{s}$')
plt.plot(roots, sigma, color=color1, marker='o', linestyle=' ', ms=3, label=r'Madgraph')
#plt.title(r'ALP contribution to $g g \rightarrow a \rightarrow t \bar{t}$')
plt.xlabel(r'$\sqrt{\hat{s}}$ [TeV]', fontsize=12)
plt.ylabel(r'$\hat{\sigma}_{\mathrm{ALP}}$ [pb]', fontsize=12)
plt.xlim([0.5,10])
plt.legend(fontsize=14, loc='lower right')
plt.tight_layout()
plt.savefig('gg_ax_tt.pdf')


#Lin only plot
plt.figure()
plt.plot(rootsLin,funcLin,color=color2, label=r'$\hat{\sigma}_{\mathrm{SM-ALP}} \sim \frac{1}{\hat{s}} \mathrm{Log}(\frac{\sqrt{\hat{s}}}{m_t})$')
plt.plot(rootsLin, sigmaLin, color=color1,marker='o',linestyle=' ', ms=3, label='Madgraph')
#plt.title('Interfernce term, partonic cross section')
plt.xlabel(r'$\sqrt{\hat{s}}$ [TeV]',fontsize=12)
plt.ylabel(r'$\hat{\sigma}_{\mathrm{SM-ALP}}$ [pb]', fontsize=12)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('gg_ax_tt_lin.pdf')


#SM only plot
plt.figure()
plt.plot(rootsSM,funcSM,color=color2, label=r'$\hat{\sigma}_{\mathrm{SM}} \sim \frac{1}{\hat{s}} \mathrm{Log}(\frac{\sqrt{\hat{s}}}{m_t})$')
plt.plot(rootsSM, sigmaSM,color=color1,marker='o',linestyle=' ', ms=3, label='Madgraph')
#plt.title('SM term, partonic cross section')
plt.xlabel(r'$\sqrt{\hat{s}}$ [TeV]', fontsize=12)
plt.ylabel(r'$\hat{\sigma}_{\mathrm{SM}}$ [pb]', fontsize=12)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('gg_ax_tt_SM.pdf')



plt.show()
