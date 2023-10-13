# Reproducing Section 4 and Appendix D of 2303.17634



The non-resonant ALP analysis with ttbar data for CMS (mtt) and ATLAS (boosted pt) can be found in the following directories:

-  ATLAS_boosted,

-  CMS_mtt.


Each has the same structure:


- limits: the limits subdirectory contains two scripts. "python get_constraints.py 100" will print out the constraint on c/fa in TeV^-1 assuming ma=100 GeV.  "python loop_masses.py" will loop over all masses and compute the constraint at each, for Figure 11.

- plotting:  "python plot_axion_data_comparison.py" will plot the comparison between data and theory (for selected ALP signals) as shown in the paper.  Currently the ALP mass is ma=100 GeV, and ct=[0,4,12,15], but these can all be changed in the code.

- data, data_from_madgraph and sm: contain the different data and theory inputs used. 






Plots and information for appendix D can be found in the directories

- kinematic_plots

- partonic_checks

