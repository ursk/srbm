srbm
====

Ising, RBM and sRBM model estimation with MPF and evaluation with AIS

Contains matlab code for Minimum Probability Flow estimation (Sohl-Dickstein et al. 2011) of Ising models, Restricted and Semi-restricted Boltzmann machines. The partition fuction of the model is calculated using Annealed Importance Sampling (Neal et al., 2001).

The package contains an example script call_estim.m that loads a set of artificial test data from demodata.mat and estimates all 3 models. 

Optimization is performed using the minFunc package by Mark Schmidt. It is packaged with the code, but you may need to run mexAll.m in the 3rd_party_code/minFunc directory to compile the code. 