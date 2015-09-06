# NSwMCMC
Minimal code to implement nested sampling with random walk MCMC moves to explore the likelihood-constrained prior

Aim: To produce a simple (i.e., easily human parse-able), robust implementation of nested sampling with MCMC moves as a benchmark for the second part of the open collaboration at https://astrostatistics.wordpress.com/nested-sampling/breaking-ns-breaking-smc/

The script should follow the basic recipe seen in John Skilling's lighthouse example: http://www.inference.phy.cam.ac.uk/bayesys/r/lighthouse.r
Namely, each replacement 'live particle' should be sought via random walk MCMC.

To facilitate the use of this script across multiple example applications the prior & likelihood should be treated as external functions.  Specifying the input and output requirements of these is the first task on the to-do list here! 
