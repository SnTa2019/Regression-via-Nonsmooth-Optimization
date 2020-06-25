'''
Main program for

L1-based piecewise linear regression using the aggregate 
subgradient method for nonsmooth DC optimization.

Reference:
A. Bagirov, S. Taheri, N. Karmitsa, N. Sultanova, and S. Asadi, "Robust piecewise linear L1-regression via
nonsmooth optimization in data sets with outliers", submitted, 2020.


The work was financially supported by the Academy of Finland 
(Project No. 289500 and 319274).

Created on 13 Jan 2020

@author: Napsu
'''
# Load library
# import time
import time 

# Import initiaization by user and the global variables
import user_init
import config as cfg


# Import Aggregate subgradient method 
import aggsub
#print aggsub.aggsub.__doc__ #This works ok, but if we call aggsub often...
aggsu = aggsub.aggsub
#print aggsu.__doc__ # Prints documentation at the beginning of function

from functions import scaling, rescaling, evaltr, evaltest
from numpy import loadtxt, zeros, ones, ones_like, append, delete #, array, zeros_like




print('PWLR-L1:')

# Data dependent parameters:
# Give your data
#infile = "../../Data/Data for regression/Data with outliers/Data processed/Scottish_new_order.txt"
#cfg.nrec = 35 #
#cfg.nfea = 3 #
#cfg.ntrain = 28 #
#nout = cfg.nfea
#cfg.kmax = 5  # max number of functions (min and/or linear functions) under maximum
#jkmax = 3     # max number of min functions under each maximum

infile     = raw_input('Name of the data file ?  ')
cfg.nrec   = int(input('Number of data points (m)?  '))
cfg.nfea   = int(input('Number of variables (n+1)?  '))
cfg.ntrain = int(input('Size of the training set?  ')) # The last elements in data are supposed to be the test set.
nout       = int(input('Number of the output variable?  '))
cfg.kmax   = int(input('Maximum number of minima under maximum (K)?  '))
jkmax      = int(input('Maximum number of linear functions under a minimum (J)?  \n'))

if nout != cfg.nfea:
# if nout != len(cfg.a)    
    exit('Sorry, in this preliminary version the output variable needs to be the last one.') # Stop the program.

qk=cfg.kmax * jkmax   # this is the maximum 
print('There will be at most %i linear pieces in the model.' %qk)

totalnit  = 0
totalnf   = 0
totalng   = 0   
neval     = 0 

# Load data from infile
datafile = open(infile,"r")
cfg.a = loadtxt(datafile)
#cfg.a = loadtxt(datafile, delimiter=',')
datafile.close()

# Note that we do not really need to ask nfea and nrec, they can be read from confic.a.
# On the other hand, this could be used for error checking.
# Print the number of data points and features
#len(cfg.a)) = number of data points 
tmp=cfg.a.size/len(cfg.a) # number of features
print('There are %i data points and %i features.' % (len(cfg.a),tmp))

cfg.min_line = zeros((cfg.kmax,cfg.nrec), dtype=int) # Index of the smallest linear function under the minimum, (k_max * nrec)-matrix;
cfg.max_piece = zeros((cfg.nrec), dtype=int)         # Index of the largest concave piece, (nrec)-matrix (or smallest convex piece);
cfg.clin = zeros((cfg.nrec))                         # Scaling parameter;
cfg.dlin = zeros((cfg.nrec))                         # Scaling parameter;
cfg.a_scaled = ones_like(cfg.a)                      # Scaled input data;
cfg.alpha1 = zeros((cfg.nrec), dtype=int)            # Index that identifies the part of the function rho1 used;

# Scaling of input matrix: the result will be in cfg.a_scaled.

cfg.a_unscaled=cfg.a
scaling()
cfg.a=cfg.a_scaled

# Just testing
# access to matrix elements: 
# a[i] # i+1:th row
# a[-1] # the last row
# a[:,i] # i+1:th column
# a[0][0] # first element
# a[1][2] # third element of second row


# First compute just linear regression. That is, we have only one linear function under min and only one min under max. 
# The number of variables in the optimization problem is nft.


#print('\n')    
print('\nLinear regression:\n')    
   
cfg.nomax = 1
cfg.jk   = ones(cfg.nomax, dtype=int)
x = ones(cfg.nfea) # Vector of ones. At the beginning the size of x is nfea.
x_reg = zeros(cfg.nfea) # Vector of ones. At the beginning the size of x is nfea.  
x_aux_min = ones(cfg.nfea) # Vector of ones. We always add only one linear function at time.
x_aux_max = ones(cfg.nfea) # Vector of ones. We always add only one linear function at time.


# Starting time
usedtime0 = time.clock()

# Initialization of the method (Data and Step 0).
#print user_init.init.__doc__ # Prints documentation at the beginning of function
eps,mit,tau,delta,sig1,sig2,c1,c2,tlimit = user_init.init()
prob = 0

#taux = tlimit
taux = max(60.0,tlimit/20.0)
#epsaux = eps
epsaux = 0.01
deltaux = 0.01 
mitaux=20
if cfg.nrec < 1000:
    mitaux = 50

 
# Call aggregate subgradient method.
x_reg,f,nit,nf,ng = aggsu(x_reg,prob,eps,mit,tau,sig1,delta,sig2,c1,c2,tlimit)
# This works if we send default parameters. 
# They can be changed in function, but user_init is not working 
# if we use this.
#x,f = aggsu(x,eps= 10.0**(-5),mit=1000,tau=10,sig1= 0.2,
#            delta= 10.0**(-7),sig2=1.0,c1=0.2,c2=0.05,tlimit=18.0)
print "f_reg is", f 
x_unscaled = rescaling(x_reg)
    
# Results:
print('\nResults for linear regression:\n')    
print "nit = ",nit
print "nf  = ",nf
print "ng  = ",ng
print "neval  = ",nf
#print "x   = ",x_reg # this is scaled result
print "x   = ",x_unscaled # this is scaled result

print('f   =  %4.6e' %f)
# CPU time used
usedtime = time.clock() - usedtime0
print('CPU time = %4.2f' %usedtime)
print('\n')

# Indices for training data for linear regression

totalnit +=nit
totalnf  +=nf
totalng  +=ng    
neval    +=nf
    
x = x_reg
    
i_loop = 0    
while True: # loop structure for incremental algorithm
    prob = 1    
    cfg.xprev = x
            
    # Computing aux_min    
    iauxmin = 0
    if cfg.jk[cfg.nomax-1] < jkmax:
        iauxmin=1
        cfg.jk[cfg.nomax-1] += 1
        #print "Maximum number of min-functions. Exit."
        
        # Solve aux-min problem.
        print('\nComputing aux-min:'),cfg.jk
        x_aux_min = x_reg # 
        x_aux_min,f_min,nit,nf,ng = aggsu(x_aux_min,prob,epsaux,mitaux,tau,sig1,deltaux,sig2,c1,c2,taux)
        print "f_aux_min is", f_min #,x_aux_min
        usedtime = time.clock() - usedtime0
        print('CPU time = %4.2f' %usedtime)
        cfg.jk[cfg.nomax-1] -= 1
        totalnit +=nit
        totalnf  +=nf
        totalng  +=ng    
        neval    +=nf

        
        #Stopping with no improvement in min and already full max:
        if (f - f_min)/(f+1) <= eps and cfg.nomax >= cfg.kmax:
            print "Stop. Adding more linear functions under maximum does not improve the model."
            break    
    
    
    # Computing aux_max    
    cfg.nomax += 1
    
    if cfg.nomax > cfg.kmax:
        print "Maximum number of max-functions. Exit."
        cfg.nomax -= 1
        break
    cfg.jk = append(cfg.jk,1) # Adding max as last piece.
        
    # Solve aux-max problem.
    print('\nComputing aux-max:'),cfg.jk
    x_aux_max = x_reg #

    # Call aggregate subgradient method.
    x_aux_max,f_max,nit,nf,ng = aggsu(x_aux_max,prob,epsaux,mitaux,tau,sig1,deltaux,sig2,c1,c2,taux)
    print "f_aux_max is", f_max #,x_aux_max
    usedtime = time.clock() - usedtime0
    print('CPU time = %4.2f' %usedtime)
    totalnit +=nit
    totalnf  +=nf
    totalng  +=ng    
    neval    +=nf

    # Stopping 
    if (f - f_min)/(f+1) <= eps and (f - f_max)/(f+1) <= eps: 
        print "Stop, No improvement in function value by adding more linear pieces."
        cfg.jk = delete(cfg.jk,cfg.nomax-1)    
        cfg.nomax -= 1
        break    
    
    # Choosing min or max:
    if f_max > f_min and iauxmin==1:
        choice_f = 0
        cfg.jk = delete(cfg.jk,cfg.nomax-1)    
        cfg.nomax -= 1
        cfg.jk[cfg.nomax-1] += 1
    
    else:
        choice_f = 1
        
    
    if choice_f == 0:  # Adding min 
        
        # Solve PWLR problem with one additional minimum
        print('\nAdding min:'),cfg.jk
        prob = 0
        x = append(x,x_aux_min) 
    
        # reinitialization
        #eps,mit,tau,delta,sig1,sig2,c1,c2,tlimit = user_init.init()
        
        # Call aggregate subgradient method.
        x,f,nit,nf,ng = aggsu(x,prob,eps,mit,tau,sig1,delta,sig2,c1,c2,tlimit)
        print "f after min is", f #,x

    
    else: # Adding max
            
        print('\nAdding max:'),cfg.jk
        prob = 0

        x = append(x,x_aux_max)
        
        # reinitialization
        #eps,mit,tau,delta,sig1,sig2,c1,c2,tlimit = user_init.init()
 
        # Call aggregate subgradient method.
        x,f,nit,nf,ng = aggsu(x,prob,eps,mit,tau,sig1,delta,sig2,c1,c2,tlimit)
        print "f after max is", f #,x
        # This works if we want send default parameters. 
        # They can be changed in function, but user_init is not working, if we use this.
        #x,f = aggsu(x,eps= 10.0**(-5),mit=1000,tau=10,sig1= 0.2,
    #            delta= 10.0**(-7),sig2=1.0,c1=0.2,c2=0.05)
    
    totalnit +=nit
    totalnf  +=nf
    totalng  +=ng    
    neval    +=nf*sum(cfg.jk)
    x_unscaled = rescaling(x)

    print('\nResults for PWLR:\n')    
    print "nit = ",totalnit 
    print "nf  = ",totalnf
    print "ng  = ",totalng
    print "neval  = ",neval
    #print "x   = ",x # this is scaled
    print "x   = ",x_unscaled
    print('f   =  %4.6e' %f)
# CPU time used
    usedtime = time.clock() - usedtime0
    print('CPU time = %4.2f' %usedtime)


# Final result
print('\n')
print "Results for PWLR with the model",cfg.jk,"."    
print "nit = ",totalnit 
print "nf  = ",totalnf
print "ng  = ",totalng
print "neval  = ",neval
print "x   = ",x_unscaled
print('f   =  %4.6e' %f)
    
# CPU time used
usedtime = time.clock() - usedtime0
print('CPU time = %4.2f' %usedtime)
print('\n')

#x_unscaled = rescaling(x) # Done above

#x=x_unscaled # not efficient but testing


# Print the resulting PWLR-function
# tmp = 0.0
# if cfg.nfea == 2 and cfg.nomax <= 2: # Note! Works for ziczac for others, need to be checked
#     for _ in range(0,100,cfg.nrec):
#     # _ is an invisible variable (i.e. an index which we are not using)
#         #tmp=tmp+320.0/cfg.nrec
#      
#         tmp=tmp+800000.0/cfg.nrec
#      
#         if (cfg.nomax==1 and cfg.jk[0]==1):   # regression
#             y=x[0]*tmp+x[1]
#             print y
#         elif (cfg.nomax==1 and cfg.jk[0]==2): # jk=2 k=1
#             y=min(x[0]*tmp+x[1],x[2]*tmp+x[3])
#             print y
#         elif (cfg.nomax==2 and cfg.jk[0]==1): # jk=[1,1] k=2
#             y=max(x[0]*tmp+x[1],x[2]*tmp+x[3])
#             print y
#         elif (cfg.nomax==2 and cfg.jk[0]==2): # jk=[2,1] k=2
#             y=min(x[0]*tmp+x[1],x[2]*tmp+x[3])
#             z=max(y,x[4]*tmp+x[5])
#             print z
#         else: # jk=2 k=2
#             y=min(x[0]*tmp+x[1],x[2]*tmp+x[3])
#             z=min(x[4]*tmp+x[5],x[6]*tmp+x[7])
#             print(max(y,z))
 
#print "Parameters used:"    
#print "deltaux =",deltaux    
#print "epsaux =",epsaux    
#print "delta =",delta    
#print "eps =",eps    
#print "mitaux =",mitaux
#print "c1 =",c1
#print "sig1 =",sig1
#print "f =", f 

# Indices for training and test data for PWL-regression
ntest = cfg.nrec - cfg.ntrain
if ntest > 0:
    print "Indices (RMSE, MAE, R^2 and r) for training data for PWLR with the model",cfg.jk,"."  
    rmse,mae,ce,r = evaltr(x_unscaled)
    print rmse, mae, ce, r
    print "Indices (RMSE, MAE, R^2 and r) for test set for PWLR with the model",cfg.jk,"."  
    pred = zeros(ntest) # Predicted values for the test set    
    pred,rmse,mae,ce,r = evaltest(x_unscaled,ntest,pred)
    print rmse, mae, ce, r
    #print pred
    