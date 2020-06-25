'''
Global parameters and variables for piecewise linear regression using L1-risk.

Created on 14 Jan 2020

@author: napsu
'''
import numpy as np

#a=0.0

a = np.zeros((1,1))  # Data matrix;
a_scaled = np.zeros((1,1))  # Data matrix;
a_unscaled = np.zeros((1,1))  # Data matrix; # we do not really need both scaled and unscaled but in this first simple version
xprev = np.zeros(1)  # previous solution;
nfea = 0             # Number of features;
nrec = 0             # Number of data points;
ntrain = 0           # Size of the training set;
kmax = 5             # Maximum number of minima under maximum (K);                
nomax = 1            # Current number of minima under maximum (nomax <= kmax);
jk = np.ones((1), dtype=int)  # Number of linear functions under each minimum;

mk = np.zeros((1), dtype=int)          # number of linear functions under each maximum  
min_line = np.zeros((1,1), dtype=int)  # Index of the smallest linear function under the minimum, (k_max * m)-matrix;
max_piece = np.zeros((1), dtype=int)   # Index of the concave piece that gives the maximum, nrec-vector
alpha = np.zeros(1)                    # Index that tells the absolute value of the function rho. No need?
alpha1 = np.zeros((1), dtype=int)      # Index that tells which part of the function rho1 is used. Add 0.5?

clin = np.zeros(1)   # Parameter for scaling;   
dlin = np.zeros(1)   # Parameter for scaling;