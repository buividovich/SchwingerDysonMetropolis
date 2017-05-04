# sd_metropolis

This is a Metropolis algorithm for solving Schwinger-Dyson equations, 
used in the paper ArXiv:1705.XXXXX

The core code for solving infinite-dimensional linear systems is implemented 
in the files src/sd_metropolis.c, include/sd_metropolis.h and is described in 
Section IV.B of ArXiv:1705.XXXXX.

The codes in the apps/ folder implement the solution of different 
Schwinger-Dyson (SD) equations for the large-N U(N)xU(N) principal chiral 
model. 

pcm_sc_cspace - SD equations in terms of U(N)-valued variables in coordinate 
space, which generate the strong-coupling expansion 

pcm_sc_mspace - SD equations in terms of U(N)-valued variables in momentum 
space, which generate the strong-coupling expansion 

pcm_wc_mspace - SD equations in terms of Hermitian matrices which 
parameterize perturbative fluctuations, as described in ArXiv:1705.XXXXX. 
This is the production code for ArXiv:1705.XXXXX. 

pcmwcdp - additional small code for processing the data produced by the 
pcm_wc_mspace code, includes the fast Fourier transform for calculating the 
correlators of U(N)-valued variables in coordinate space (see Section V.A of 
ArXiv:1705.XXXXX). 
 
