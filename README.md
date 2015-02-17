Global N<sub>2</sub>O Ellipses
==============================

Data and code associated with Snider, Venkiteswaran, Schiff, and Spoelstra. 2015. From the Ground Up: Global Nitrous Oxide Sources are Constrained by Stable Isotope Values. PLoS ONE PONE-D-14-15929R2

Manuscript submitted January 2014. Revised manuscript submitted December 2014. Manuscript accepted 30 January 2015.

The raw data and R-script used to perform the analyses and figures are here.

* figures for N2O isotopes.R: *R-script that loads the raw data, performs calculations, and generates the figures*
* mixing model for N2O isotopes ellipses.R: *R-script that creates summaries of the N<sub>2</sub>O isotope data and runs various mixing models combinations via MixSIAR*
* Dataset S1.csv: *Raw N<sub>2</sub>O isotope data with citations*
* Dataset S2.csv: *Emission-weighted average isotopes data with citations, used in Figure 5*
* Figure [1-7].pdf: *Figures as accepted*

The R-code for the analyses, figures, and mixing models requires (at the very least) ggplot2, reshape, RColorBrewer, scales, grid, siar, plyr, R2jags, MASS, reshape, devtools. Thanks to the authors of these packages.
