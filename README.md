# BAcode
# This is the code written during the research project of my Bachelor degree in Physics.  
# I implemented the multicolor blackbody emission model from the accretion disk of an AGN 
# and fitted it to the Optical and UV part of the spectral energy distribution of three high-redshift sources. 
# The devised fitting technique was based on a Markov Chain Monte Carlo method, implemented by myself, 
# that allowed to handle many parameters of the physical system and to explore their mutual correlations.


-> mcLog.cpp reads the SED data of a single source from a text file (data_fileName.txt), 
performs the fit and writes a file (points_filename.txt) with 6 columns: Chi^2, log(M), log(Mr), log(Rin), log(Rout).
Those are the values of the Chi Squared and of the four parameters of the model: 
mass of the black hole, mass accretion rate, inner radius and outer radius of the accretion disk.
Those values represent the points explored during the random walk in the parameters space.

HOW TO RUN THE PROGRAM: ./mcLog fileName.txt
fileName.txt is the file downloaded from SED generator (see the example file S50014813.txt)


-> ss.h defines the functions needed to compute the spectrum of the Shakura-Sunyaev model.


-> normal.h defines the functions needed to generate a random number from a normal standard distribution.


Nearby samples are generally correlated with each other and do not correctly reflect the distribution. 
This means that if we want a set of independent samples, we have to throw away the majority of samples 
and only take every n-th sample, for some value of n.
-> autocorr.cpp finds the value of n and keeps only the uncorrelated data.

//HOW TO RUN THE PROGRAM: ./autocorr fileName.txt
//fileName.txt is the file downloaded from SED generator (see the example file S50014813.txt)


