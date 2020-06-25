'''
Initialization for the aggregate subgradient method.

Created on 15 Jan 2020.
Modified for PWLRL1 on 29 Feb 2020. 

@author: napsu
'''
# Load library
#import numpy as np

# Load globals
import config as cfg


def init():
    """ Initialization for Aggregate subgradient method.
    
    Input:
        none;
        
    Output:
        eps    - optimality tolerance;
        mit    - maximum number of inner iterations;   
        tau    - proximity parameter, tau > 0;
        sig1   - decrease parameter for tau, sig1 in (0,1);
        delta  - inner iteration tolerance, delta > 0; 
        sig2   - decrease parameter for delta, sig2 in (0,1];
        c1,c2  - line search parameters, 0 < c2 <= c1 < 1; 
        tlimit - time limit for one round of AggSub.
            
    """    
    
    # Parameters, default values are given as comments.
    mit = 1000 # default = 1000.
    
    # Could also depend on nfea (the number of variables at the beginning).
    if cfg.nrec < 100:   
        tau = 2.0        # tau > 0, default = 2, if nrec < 200, 10, otherwise.
        eps = 10.0**(-6)  
        delta = 10.0**(-6)  
        
    elif cfg.nrec < 500:
        tau = 2.0
        eps = 10.0**(-5)  
        delta = 10.0**(-5)  
    
    elif cfg.nrec < 500:
        tau = 2.0
        eps = 10.0**(-4)  
        delta = 10.0**(-4)  
         
    else:
        tau = 10.0
        eps = 10.0**(-4)  
        delta = 10.0**(-4)  
        
    sig1 = 0.5 # sig1 in (0,1), default = 0.8.
    #sig1 = 0.8 # sig1 in (0,1), default = 0.8.
    sig2 = 1.0 # sig2 in (0,1], default = 1.0.
    #c1 = 0.01 # c1 in (0,1), default 0.2.
    c1 = 0.2 # c1 in (0,1), default 0.2.
    c2 = 0.25*c1 # c2 in (0,c1) default = 0.05.  
     
    if cfg.nrec < 50000:
        tlimit=400
    else:        
        tlimit = 1800.00
        
    return eps,mit,tau,delta,sig1,sig2,c1,c2,tlimit
