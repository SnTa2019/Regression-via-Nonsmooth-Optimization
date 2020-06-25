'''
Aggregate subgradient method for nonsmooth DC optimization.

Created on 14 Jan 2020.
Modified for PWLRL1 on 29 Feb 2020. 

@author: napsu
'''
# Load library
import numpy as np

# Load global data matrix
import config
#from config import a, tlimit 

# Define aggsub
def aggsub(x,prob,eps,mit,tau,sig1,delta,sig2,c1,c2,tlimit):
    """ Aggregate subgradient method for nonsmooth DC optimization.
    
    Input:
        x       - a starting point;
        prob    - selection of problem:
                0 - PWLRL1 problem
                1 - auxiliary min-problem
        eps     - optimality tolerance;
        mit     - maximum number of inner iterations;   
        tau     - proximity parameter, tau > 0;
        sig1    - decrease parameter for tau, sig1 in (0,1);
        delta   - inner iteration tolerance, delta > 0; 
        sig2    - decrease parameter for delta, sig2 in (0,1];
        c1,c2   - line search parameters, 0 < c2 <= c1 < 1; 
        tlimit  - time limit for aggsub.
        
    Global parameters:    
        a       - an input matrix;
        
    Output:
        x       - a solution obtained;
        f       - the objective value at x;   
        nit     - number of iterations;
        nf      - number of function evaluations;        
        ng      - number of subgradient evaluations;
    """

    # Import DC functions and their subgradients
    import functions
    #print(config.a)   # just testing     
    # import time
    import time     
    
    # Starting time for aggsub
    usedtime0 = time.clock()
    
    
    fdiff = -99.0   
    # Initialization
    if prob == 0:
        f1 = functions.f1(x)
        f2 = functions.f2(x)
    else:
        f1 = functions.auxminf1(x)
        f2 = functions.auxminf2(x)
        
    fold = f1-f2   # Original function value        
    print "f original is", fold
    print "Computing..."
    nf = 1                      # Number of function evaluations.
    ng = 0                      # Number of subgradient evaluations.         
    nit = 0                     # Number of iterations.
    maxii = max(min(x.size,500),50)    # Maximum number of inner iterations.
    stopii = -1                 # reason to stop the inner iteration:
                                # stopii =-1 - not stopped yet;
                                # stopii = 0 - small gradient: canditate solution;
                                # stopii = 1 - decent direction found;
                                # stopii = 2 - too many inner iterations without progress.
    small = 1e-5
    nrmnew = 1e+10
    ii = 0
        
         
    # Outer iteration
    while True:
        # Step 1
        nii = 0                 # Number of inner iterations.
        dd = np.ones_like(x)/np.sqrt(x.size)  # Initialization of the direction dd.
        if prob == 0:
            f1 = functions.f1(x)  # Needed for correct confic tables. 
            g2 = functions.df2(x) # Subgradient of the DC component f2. 
            f1 = functions.f1(x+tau*dd)  # test for bug fix 
            g1 = functions.df1(x+tau*dd) # Subgradient of the DC component f1. 
        else:
            f1 = functions.auxminf1(x)  # Needed for correct confic tables. 
            g2 = functions.dauxminf2(x) # Subgradient of the DC component f2. 
            f1 = functions.auxminf1(x+tau*dd)  # test for bug fix 
            g1 = functions.dauxminf1(x+tau*dd) # Subgradient of the DC component f1. 
        
        # No need to recompute g2 if returning from Step 6.
        ng += 1
        sg = g1-g2            # Approx. subgradient of the function.
        asg = sg              # Aggregate subgradient.
        nrmasg = 0.0          # Norm of the aggregate subgradient.
    
        # Inner iteration
        while True:
            nii += 1
        
            # Step 2
            sgdiff=np.dot(sg-asg,sg-asg)
            if sgdiff > small:
                #lam = (nrmasg*nrmasg - np.dot(sg,asg))/sgdiff
                lam = -np.dot(asg,sg-asg)/sgdiff # this should be the same as above
            else:
                lam = 0.50
                
            if lam > 1.0 or lam < 0.0:
                print "Projecting lambda = ", lam, "to [0,1]."
                if lam < 0:
                    lam = 0.0
                else:
                    lam = 1.0  
                      
            # If lam=0 nothing changes and we end up to inner iteration termination 2.
                    
            # Step 3
            asg = lam*sg + (1.0-lam)*asg  # The new aggregate subgradient.
            nrmasg = np.linalg.norm(asg)  # Norm of the aggregate subgradient.
            
            if nrmasg < delta:  # Inner iteration termination
                stopii = 0
                break           # With this we should go to Step 6.
            
            if nii%5==0:        # Inner iteration termination 2
                nrmold = nrmnew
                nrmnew=nrmasg
                if np.abs(nrmold-nrmnew)<0.0001:
                    #print "Norm is not changing." 
                    stopii = 0
                    break       # With this we should go to Step 6.
                         
            
            # Step 4: Search direction
            dd = - asg / nrmasg
            
            # Step 5
            xtau = x + tau*dd
            if prob == 0:
                f1 = functions.f1(xtau)
                f2 = functions.f2(xtau)
            else:
                f1 = functions.auxminf1(xtau)
                f2 = functions.auxminf2(xtau)
            
            fnew = f1-f2
            nf += 1
            dt = fnew - fold
            
            #if (dt > -c1*tau*nrmasg and -dt/fold > 0.1): # makes results worse (limited testing)
            if (dt > -c1*tau*nrmasg): # Not a descent direction
                #if (dt < 0 and -dt/fold > 0.05): # just for testing purposes 
                #    print "No descent but should it be?",-dt/fold,dt,-c1*tau*nrmasg,fold,fnew
                if prob == 0:    
                    g1 = functions.df1(xtau)
                else:
                    g1 = functions.dauxminf1(xtau)
                sg = g1-g2
                ng += 1
            else:
                stopii = 1
                break # with this we should go to Step 7.    
        
                
            # Additional termination from inner iteration
            if nii > maxii:
                if tau > eps:
                    #print "Too many inner iterations. Adjusting tau."
                    stopii = 2
                else:
                    print "Too many inner iterations with no descent direction found." #,nit,fnew,fold
                    if ii == 0:
                        ii = 1
                        nii = 0               # Number of inner iterations starts again from zero.
                        #dd = -np.ones_like(x)/np.sqrt(x.size)  # Initialization of the direction dd.
                        dd = -dd  # Opposite of the direction dd.
                        if config.nfea < 200: # tau > 0, default = 10, if n<200,  
                            tau = 10.0        #                    50, otherwise.
                        else:
                            tau = 50.0
                        if prob == 0:    
                            f1 = functions.f1(x+tau*dd)  # test for bug fix 
                            g1 = functions.df1(x+tau*dd) # Subgradient of the DC component f1. 
                        else:
                            f1 = functions.auxminf1(x+tau*dd)  # test for bug fix 
                            g1 = functions.dauxminf1(x+tau*dd) # Subgradient of the DC component f1. 
                        
                        ng += 1
                        sg = g1-g2            # Approx. subgradient of the function.
                        asg = sg              # Aggregate subgradient.
                        nrmasg = 0.0          # Norm of the aggregate subgradient.
                        print "Trying opposite direction."                        
                        print "Computing..."
                        continue
                        
                    else:
                        print "Exit."
                        stopii = 3
                break
            
                    
        # Step 6: Stopping criterion
                     
        if stopii == 0: # Small norm 
            if tau <= eps:
                print "Critical point found."
                break # Critical point found   
            else:
                tau *= sig1
                delta *= sig2
                stopii = -1
                nit += 1
                continue
              
            
        # Step 7: Line search
        elif stopii == 1: # Descent direction
            if ii == 1:
                print "Descent direction found." 
                print "Computing..."
                ii = 0
                
            step = tau    # Here fnew = f(xtau), fold = f(x)
            count = 0
            while count < 10: # Maximum of 10 search are made
                count +=1
                step *=2.0
                xnew = x + step*dd
                if prob == 0:
                    f1 = functions.f1(xnew)
                    f2 = functions.f2(xnew)
                else:
                    f1 = functions.auxminf1(xnew)
                    f2 = functions.auxminf2(xnew)
                
                nf += 1
                fdiff = f1-f2 - fold
            
                if fdiff <= -c2 * step * nrmasg:
                    xtau = xnew
                    fnew = f1-f2
                else:
                    break
            
            # Step 8: Update step
            # The last f1 and f2 are computed at xnew =/ xtau, need to compute f1(xtau) to get max_line correctly
        
            x = xtau 
            fold = fnew
        
        if tau > eps:
            tau *= sig1 # tau too small is not good, needs safeguard 
        elif np.abs(fdiff) < 1e-7:
            print "Termination with small change in function values." #,fdiff
            break
            
        delta *= sig2
            
        nit += 1
        
        # Termination with maximum number of iterations.
        if nit >= mit: # Number of iterations > mit
            print "Termination with maximum number of iterations."
            break
        
        # Termination with the time limit for AggSub.
        usedtime = time.clock() - usedtime0
        if usedtime > tlimit:
            print "Termination with the time limit for AggSub."
            break  # with this we should go to Step 7.
            
        if stopii == 3:
            
            #print "No decent direction found in inner iteration."
            break
        
        stopii = -1 # To next inner iteration.
            
    # Outer iteration ends
    
    f=fold
    
    return x,f,nit,nf,ng