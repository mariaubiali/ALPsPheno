import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#sns.set_style('whitegrid')
colors = sns.color_palette('colorblind',8)

#Bins
#bins_cms =  np.array([250.0,400.0,480.0,560.0,640.0,720.0,800.0,900.0,1000.0,1150.0,1300.0,1500.0,1700.0,2000.0,2300,3500])
bins_cms = np.arange(250,4000,150)
bins_width = bins_cms[1:]-bins_cms[:-1]
bins_cms_C = 0.5*(bins_cms[1:] + bins_cms[:-1])


#Madgraph predictions
BR = 0.2877
axion_lin = np.divide(BR*np.loadtxt('data_from_madgraph/mtt_nnpdf4p0_alphaS_lin.txt',usecols=(2,),dtype=float),bins_width)
axion_quad = np.divide(BR*np.loadtxt('data_from_madgraph/mtt_nnpdf4p0_alphaS_quad.txt',usecols=(2,),dtype=float),bins_width)
sm = np.divide(BR*np.loadtxt('data_from_madgraph/mtt_nnpdf4p0_alphaS_SM.txt',usecols=(2,),dtype=float), bins_width)


sm_plot = np.concatenate([sm, [sm[-1]]])


ax_signal1 = sm + (1**2)*axion_lin + (1**4)*(axion_quad)
ax_signal1_plot = np.concatenate([ax_signal1, [ax_signal1[-1]]])
ax_signal1_ratio_plot = np.divide(ax_signal1_plot,sm_plot)

ax_signal6 = sm + (6**2)*axion_lin + (6**4)*(axion_quad)
ax_signal6_plot = np.concatenate([ax_signal6, [ax_signal6[-1]]])
ax_signal6_ratio_plot = np.divide(ax_signal6_plot,sm_plot)

ax_signal10 = sm + (10**2)*axion_lin + (10**4)*(axion_quad)
ax_signal10_plot = np.concatenate([ax_signal10, [ax_signal10[-1]]])
ax_signal10_ratio_plot = np.divide(ax_signal10_plot,sm_plot)

ax_signal14 = sm + (14**2)*axion_lin + (14**4)*(axion_quad)
ax_signal14_plot = np.concatenate([ax_signal14, [ax_signal14[-1]]])
ax_signal14_ratio_plot = np.divide(ax_signal14_plot,sm_plot)



#First plot: full signal mtt spectrum, and ratio to SM
fig = plt.figure(figsize=(9, 7))
gs = fig.add_gridspec(2, 1, height_ratios=(1, 1),
                  left=0.14, right=0.86, bottom=0.1, top=0.9,
                  wspace=0.55, hspace=0.35)

ax = fig.add_subplot(gs[0, 0])
ax_ratio = fig.add_subplot(gs[1, 0], sharex=ax)

ax.step(bins_cms, sm_plot, where='post', label='SM', color=colors[0])
ax.step(bins_cms, ax_signal1_plot, where='post', label=r'$c_t = 1$', color=colors[1])
ax.step(bins_cms, ax_signal6_plot, where='post', label=r'$c_t = 6$', color=colors[2])
ax.step(bins_cms, ax_signal10_plot, where='post', label=r'$c_t = 10$', color=colors[4])
ax.step(bins_cms, ax_signal14_plot, where='post', label=r'$c_t = 14$', color=colors[3])


ax_ratio.step(bins_cms, ax_signal1_ratio_plot, where='post', label=r'$c_t = 1$', color=colors[1])
ax_ratio.step(bins_cms, ax_signal6_ratio_plot, where='post', label=r'$c_t = 6$', color=colors[2])
ax_ratio.step(bins_cms, ax_signal10_ratio_plot, where='post', label=r'$c_t = 10$', color=colors[4])
ax_ratio.step(bins_cms, ax_signal14_ratio_plot, where='post', label=r'$c_t = 14$', color=colors[3])


ax.ticklabel_format(style='sci')
ax.set_xlabel(r'$m_{t\bar{t}}$ [GeV]')
ax_ratio.set_xlabel(r'$m_{t\bar{t}}$ [GeV]')
ax.set_ylabel(r'$\frac{d \sigma}{d m_{t \bar{t}}}$ [pb/GeV]', fontsize=13)
ax_ratio.set_ylabel('Ratio to SM')
ax.set_xlim([250,3.5e3])
ax_ratio.set_xlim([250,3.5e3])
ax.set_yscale('log')
#ax_ratio.set_yscale('log')
plt.tight_layout()
plt.savefig('appD_mtt.pdf')

ax.legend(loc='upper right')



#Second plot: signal components only

lin_plot = np.concatenate([axion_lin, [axion_lin[-1]]])
lin_ratio = np.divide(axion_lin, sm)
lin_ratio_plot = np.concatenate([lin_ratio, [lin_ratio[-1]]])

quad_plot = np.concatenate([axion_quad, [axion_quad[-1]]])
quad_ratio = np.divide(axion_quad, sm)
quad_ratio_plot = np.concatenate([quad_ratio, [quad_ratio[-1]]])



fig = plt.figure(figsize=(9, 4))
gs = fig.add_gridspec(1, 2,
                  left=0.12, right=0.9, bottom=0.15, top=0.9,
                  wspace=0.3, hspace=0.1)


ax_lin = fig.add_subplot(gs[0, 0])
ax_quad = fig.add_subplot(gs[0, 1])

ax_lin.step(bins_cms, lin_ratio_plot, where='post')
ax_quad.step(bins_cms, quad_ratio_plot, where='post')

ax_lin.ticklabel_format(style='sci')
ax_quad.ticklabel_format(style='sci')

ax_lin.set_xlabel(r'$m_{t\bar{t}}$ [GeV]',fontsize=13)
ax_quad.set_xlabel(r'$m_{t \bar{t}}$ [GeV]',fontsize=13)

ax_lin.set_ylabel(r'$\frac{d \sigma_{\mathrm{SM-ALP}}}{d m_{t \bar{t}}}$/$\frac{d \sigma_{\mathrm{SM}}}{d m_{t \bar{t}}}$',fontsize=13)
ax_quad.set_ylabel(r'$\frac{d \sigma_{\mathrm{ALP}}}{d m_{t \bar{t}}}$/$\frac{d \sigma_{\mathrm{SM}}}{d m_{t \bar{t}}}$',fontsize=13)

ax_lin.set_xlim([250,3.5e3])
ax_quad.set_xlim([250,3.5e3])
ax_lin.set_ylim([1.3e-4,3.5e-4])
ax_quad.set_ylim([1.3e-7,2.3e-6])


ax_lin.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.tight_layout()
plt.savefig('appD_mtt_kinematics.pdf')




plt.show()



