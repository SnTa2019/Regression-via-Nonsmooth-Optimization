'''
DC formulations of functions and subgradients for 
piecewise linear regression using L1-risk.

Created on 12 Feb 2020

@author: napsu
'''
# Load library
import numpy as np

# Load globals
import config as cfg

# DC-FUNCTION ###########################################################################################

# Define f1
def f1(x):
    """ First DC component of nonsmooth PWLR-L1 function.
    
    Input:
        x        - current point (including output), x in R^(n*kmax*jk);  
        
    Global parameters:
        ntrain   - number of data points in training set;
        
    Output:
        f        - the first DC component function value at x;   
        
    """
    
    # Sum over data points
    f = 0.0
    for m_ind in range(cfg.ntrain):
        f += f1_part_i(x,m_ind) 
        
    return f


# Define part i of f1
def f1_part_i(x,m_ind):
    """ Part i (i=1,...,m) of the first DC component of nonsmooth PWLR-L1 function.
    
    Input:
        x        - current point (including output), x in R^(n*kmax*jk);  
        m_ind    - current data point;
        
    Global parameters:
        a        - two dimensional data matrix;
        nfea     - number of features (including output);
        alpha
        alpha1   - index that tells which part of the function rho1 is used.
        
    
    Output:
        f        - part i of the first DC component function value at x;   
        
    """

    #f = max(2.0*rho1(x,m_ind)-cfg.a[m_ind,cfg.nfea-1] ,2.0*rho2(x,m_ind)+cfg.a[m_ind,cfg.nfea-1])
    tmp1 = 2.0*rho1(x,m_ind)-cfg.a[m_ind,cfg.nfea-1]
    tmp2 = 2.0*rho2(x,m_ind)+cfg.a[m_ind,cfg.nfea-1]
    
    # checking absolute value of rho-b_i = rho1-rho2-b_i
    #if (tmp1-tmp2 > cfg.a[m_ind,cfg.nfea-1]):
    #    cfg.alpha[m_ind] = 1.0
    #if (tmp1-tmp2 == cfg.a[m_ind,cfg.nfea-1]):
    #    cfg.alpha[m_ind] = 0.5
    #else:
    #    cfg.alpha[m_ind] = 0.0
        
    # checking maximum used in rho1    
    if (tmp1 > tmp2):
        f = tmp1
        cfg.alpha1[m_ind] = 1       
    elif (tmp1 < tmp2):
        f = tmp2
        cfg.alpha1[m_ind] = 0         
    else:
        f = tmp2
        cfg.alpha1[m_ind] = 2         
    
    return f


# Define f2
def f2(x):
    """ Second DC component of nonsmooth PWLR-L1 function.
    
    Input:
        x        - current point (including output), x in R^(n*kmax*jk);  
        
    Global parameters:
        ntrain   - number of data points in the training set;
        
    Output:
        f        - the second DC component function value at x;   
        
    """

    # Sum over data points
    f = 0.0
    for m_ind in range(cfg.ntrain):
        f += rho1(x,m_ind) + rho2(x,m_ind) 
        
    return f


def rho1(x,m_ind):
    """ Subprogram for computing the nonsmooth PWLR-L1 function.
        
    Input:
        x          - current point (including output), x in R^(n*kmax*jk);  
        m_ind      - current data point;
        
    Global parameters:
        nomax      - current number of functions (minima and/or linear functions) under maximum;
        max_piece  - the index of the concave piece that gives the maximum, m vector;
    
    Output:
        f          - auxiliary function value  at x for computing the nonsmooth PWLR-L1 function;   

    """
    cc_sum = rho2(x,m_ind) 
    f = cc_sum + concave_piece(x,0,m_ind)
    cfg.max_piece[m_ind] = 0
                 
    for k_ind in range(1,cfg.nomax):
        f_tmp = cc_sum + concave_piece(x,k_ind,m_ind)   
        if f_tmp > f: 
            f = f_tmp
            cfg.max_piece[m_ind] = k_ind
    
    return f


def rho2(x,m_ind):
    """ Auxiliary function for computing the nonsmooth PWLR-L1 function.
        
    Input:
        x          - current point (including output), x in R^(n*kmax*jk);  
        m_ind      - current data point;
        
    Global parameters:
        nomax      - current number of functions (minima and/or linear functions) under maximum;
        
    Output:
        f          - auxiliary function value  at x for computing the nonsmooth PWLR-L1 function;   

    """
    
    f = 0.0
    for k_ind in range(cfg.nomax):
        f -= concave_piece(x,k_ind,m_ind) 

    return f

# Define concave piece
def concave_piece(x,k_ind,m_ind):
    """ Calculates the concave function f_cc = min_{j=1,...,J_k} { (x^{kj})^T a^{i} + y_{kj} }
    
    Input:
        x        - current point (including output), x in R^(n*kmax*jk);  
        k_ind    - current piece, k_ind in [0,...,kmax-1];
        m_ind    - current data point, m_ind in [0,...,m-1]
        
    Global parameters:
        a        - two dimensional data matrix;
        nfea     - number of features (including output);
        jk       - number of linear functions under each minimum;
        min_line - the index of the smallest value, kmax * m matrix;
    
    Output:
        f_cc     - the concave function value at x;   
        
    """
    line_start=cfg.nfea*sum(cfg.jk[i] for i in range(k_ind))
    f_cc=np.dot(x[line_start:line_start+(cfg.nfea-1)],cfg.a[m_ind,:cfg.nfea-1])+x[line_start+cfg.nfea-1]
    cfg.min_line[k_ind,m_ind] = 0  # a global variable to save the smallest value.
    
    # next lines
    line_start += cfg.nfea
    for j in range(1,cfg.jk[k_ind]): # jk is ok, range does not take the limit itself but jk-1.
        f_tmp = np.dot(x[line_start:line_start+(cfg.nfea-1)],cfg.a[m_ind,:cfg.nfea-1])+x[line_start+cfg.nfea-1]
        
        # Minimum of lines
        if f_tmp < f_cc:
            f_cc = f_tmp
            cfg.min_line[k_ind,m_ind] = j
        line_start += cfg.nfea
    
    return f_cc


# Define subgradient of the concave piece
def grad_cc_piece(x,k_ind,m_ind):
    """ Calculates the subgradient of the concave function 
    
            f_cc = min_{j=1,...,J_k} { (x^{kj})^T a^{i} + y_{kj} }
    
    Input:
        x        - current point (including output), x in R^(n*kmax*jk);  
        k_ind    - current piece, k_ind in [0,...,kmax-1];
        m_ind    - current data point, m_ind in [0,...,m-1];
        
    Global parameters:
        a        - two dimensional data matrix;
        nfea     - number of variables (including output)
        jk       - number of linear functions under each minimum;
        min_line - the index of the smallest value, kmax * m matrix;
        
    Output:
        g_cc     - the subgradient of the concave function at x;   
        
    """
    
    g_cc = np.zeros_like(x)
    i = cfg.min_line[k_ind,m_ind] * cfg.nfea + cfg.nfea*sum(cfg.jk[i] for i in range(k_ind)) # starting of the min line
    g_cc[i:i+cfg.nfea-1] = cfg.a[m_ind,:cfg.nfea-1]
    g_cc[i+cfg.nfea-1] = 1.0
                     
    return g_cc

# Define df1
def df1(x):
    """ Subgradient of f1.
    
    Input:
        x          - current point (including output), x in R^(n*kmax*jk);  
        
    Global parameters:
        ntrain     - number of data points in the training set;
        alpha1     - index that tells which part of the function rho1 is used;
        max_piece  - the index of the concave piece that gives the maximum, 
                     note that f1 need to be computed at x, m vector;
        
    Output:
        g          - the subgradient of the f1 function at x;   
    
    """
    g = np.zeros_like(x)
    for m_ind in range(cfg.ntrain):
        if (cfg.alpha1[m_ind] == 1): # 2*rho1 - b > 2*rho2 + b
            for k_ind in range(cfg.max_piece[m_ind]):
                g -= grad_cc_piece(x,k_ind,m_ind)
                
            for k_ind in range(cfg.max_piece[m_ind]+1,cfg.nomax):
                g -= grad_cc_piece(x,k_ind,m_ind)
                
            
        elif (cfg.alpha1[m_ind] == 0): # 2*rho1 - b < 2*rho2 + b
            for k_ind in range(cfg.nomax):
                g -= grad_cc_piece(x,k_ind,m_ind) #partial(rho2)
            
        else: # 2*rho1 - b == 2*rho2 + b, could be coded, but now just same as above
            for k_ind in range(cfg.nomax):
                g -= grad_cc_piece(x,k_ind,m_ind) # 1/2* partial(rho1) + 1/2* partial(rho2)
                
    g = 2.0 * g
    
    return g


# Define df2
def df2(x):
    """ Subgradient of f2.
    
    Input:
        x        - current point (including output), x in R^(n*kmax*jk);  
        
    Global parameters:
        ntrain     - number of data points in the training set;
        nomax      - current number of functions under maximum (K);
        max_piece  - the index of the concave piece that gives the maximum,
                    note that f1 need to be computed at x, m vector;
         
    Output:
        g          - the subgradient of the f2 function at x;   
        """
    g = np.zeros_like(x)
    for m_ind in range(cfg.ntrain):
        for k_ind in range(cfg.nomax):
            g -= 2.0*grad_cc_piece(x,k_ind,m_ind) 
        g += grad_cc_piece(x,cfg.max_piece[m_ind],m_ind) 
                          
    return g


# AUX-MIN-FUNCTION ###########################################################################################

# Define aux-min-f1
def auxminf1(x):
    """ First DC component of the aux-min-function.
    
    Input:
        x        - current point (including output), x in R^(nfea);  
        
    Global parameters:
        xprev    - previous solution;
        ntrain   - number of data points in training set;
        a        - two dimensional data matrix;
        nfea     - number of features (including output);
        alpha1   - index that tells which part of the function rho1 is used.
        
    Output:
        f        - the first DC component function value at x;   
        
    """
 
# Sum over data points
    f = 0.0
    for m_ind in range(cfg.ntrain):
        f += auxmin_f1_part_i(x,m_ind) 
        
    return f


# Define aux-min-f1
def auxmin_f1_part_i(x,m_ind):
    """ First DC component of the aux-min-function.
    
    Input:
        x        - current point (including output), x in R^(nfea);  
        
    Global parameters:
        xprev    - previous solution;
        ntrain   - number of data points in training set;
        a        - two dimensional data matrix;
        nfea     - number of features (including output);
        alpha1   - index that tells which part of the function rho1 is used.
        
    Output:
        f        - the first DC component function value at x;   
        
    """
     
    tmp1 = 2.0*auxminrho1(x,m_ind)-cfg.a[m_ind,cfg.nfea-1] 
    tmp2 = 2.0*auxminrho2(x,m_ind)+cfg.a[m_ind,cfg.nfea-1]

    # checking maximum used in auxminrho1    
    if (tmp1 > tmp2):
        f = tmp1
        cfg.alpha1[m_ind] = 1  # alpha1 should be ok here. We do not solve aux and real problem at the same time.      
    elif (tmp1 < tmp2):
        f = tmp2
        cfg.alpha1[m_ind] = 0  
    else:
        f = tmp2
        cfg.alpha1[m_ind] = 2  

    return f



# Define aux-min-f2
def auxminf2(x):
    """ Second DC component of the aux-min-function.
    
    Input:
        x        - current point (including output), x in R^(nfea);  
        
    Global parameters:
        ntrain   - number of data points in the training set;
        
    Output:
        f        - the second DC component function value at x;   
        
    """
    # Sum over data points
    f = 0.0
    for m_ind in range(cfg.ntrain):
        f += auxminrho1(x,m_ind) + auxminrho2(x,m_ind) 
        
    return f


def auxminrho1(x,m_ind):
    """ Subprogram for computing the aux-min-function.
        
    Input:
        x          - current point (including output), x in R^(nfea);  
        m_ind      - current data point;
        
    Global parameters:
        xprev      - previous solution;
        nomax      - current number of functions (minima and/or linear functions) under maximum;
        max_piece  - the index of the concave piece that gives the maximum, m vector;
    
    Output:
        f          - auxiliary function value  at x;   

    """
    
    cc_sum = auxminrho2(x,m_ind) 
    f = cc_sum + auxmin_cc_piece(x,0,m_ind) 
    cfg.max_piece[m_ind] = 0 # max_piece should be ok here. We do not solve aux and real problem at the same time.
             
    for k_ind in range(1,cfg.nomax):
    
        f_tmp = cc_sum + auxmin_cc_piece(x,k_ind,m_ind) 
        if f_tmp > f: 
            f = f_tmp
            cfg.max_piece[m_ind] = k_ind
    
    return f


def auxminrho2(x,m_ind):
    """ Subprogram for computing the aux-min-function.
        
    Input:
        x          - current point (including output), x in R^(n*kmax*jk);  
        m_ind      - current data point;
        
    Global parameters:
        nomax      - current number of functions (minima and/or linear functions) under maximum;
        
    Output:
        f          - auxiliary function value  at x for computing the nonsmooth PWLR-L1 function;   

    """
    
    f = 0.0
    for k_ind in range(cfg.nomax):
        f -= auxmin_cc_piece(x,k_ind,m_ind)   

    return f

# Define concave piece for aux-min-function
def auxmin_cc_piece(x,k_ind,m_ind):
    """ Calculates the concave function f_cc = min_{j=1,...,J_k} { (x^{kj})^T a^{i} + y_{kj} }
        only the newest linear function is considered as variable.
    
    Input:
        x        - current point (including output), x in R^(n*kmax*jk);  
        k_ind    - current piece, k_ind in [0,...,kmax-1];
        m_ind    - current data point, m_ind in [0,...,m-1]
        
    Global parameters:
        xprev    - previous solution;
        a        - two dimensional data matrix;
        nfea     - number of features (including output);
        nomax    - number of max-functions;
        jk       - number of linear functions under each minimum;
        min_line - the index of the smallest value, kmax * m matrix;
    
    Output:
        f_cc     - the concave function value at x;   
        
    """
    
    # Adding new linear function as a last function:
    # The first line. If jk = 1 and k_ind = nomax, this is a new line, otherwise an old one.
    line_start=cfg.nfea*sum(cfg.jk[i] for i in range(k_ind))
    #print line_start,cfg.jk,k_ind,cfg.nomax-1,cfg.jk[k_ind], cfg.xprev,x
    if cfg.jk[k_ind]==1 and k_ind==cfg.nomax-1:
        #print "hihu0"
        f_cc=np.dot(x[0:cfg.nfea-1],cfg.a[m_ind,:cfg.nfea-1])+x[cfg.nfea-1]
        cfg.min_line[k_ind,m_ind] = 0  # a global variable to save the smallest value.
        return f_cc
    else:
        #print "hihu1",line_start,k_ind
        f_cc=np.dot(cfg.xprev[line_start:line_start+(cfg.nfea-1)],cfg.a[m_ind,:cfg.nfea-1])+cfg.xprev[line_start+cfg.nfea-1]
        cfg.min_line[k_ind,m_ind] = 0  # a global variable to save the smallest value.
        if cfg.jk[k_ind]==1:
            return f_cc
        
    # Next lines
    line_start += cfg.nfea
    for j in range(1,cfg.jk[k_ind]-1): # Everything but the first and last.
        f_tmp = np.dot(cfg.xprev[line_start:line_start+(cfg.nfea-1)],cfg.a[m_ind,:cfg.nfea-1])+cfg.xprev[line_start+cfg.nfea-1]
        
        # Minimum of lines
        if f_tmp <= f_cc:
            f_cc = f_tmp
            cfg.min_line[k_ind,m_ind] = j
        line_start += cfg.nfea
    
    # The last line.
    if k_ind==cfg.nomax-1:
        #print "hihu3"
        f_tmp = np.dot(x[0:cfg.nfea-1],cfg.a[m_ind,:cfg.nfea-1])+x[cfg.nfea-1]
    else:    
        
        f_tmp = np.dot(cfg.xprev[line_start:line_start+(cfg.nfea-1)],cfg.a[m_ind,:cfg.nfea-1])+cfg.xprev[line_start+cfg.nfea-1]
            
    # Minimum of lines
    if f_tmp <= f_cc:
        f_cc = f_tmp
        cfg.min_line[k_ind,m_ind] = cfg.jk[k_ind]-1 

    return f_cc    


# Define subgradient of the concave piece for aux-min-function
def dauxmin_cc_piece(x,k_ind,m_ind):
    """ Calculates the subgradient of the concave function 
    
            f_cc = min_{j=1,...,J_k} { (x^{kj})^T a^{i} + y_{kj} }
    
    Input:
        x        - current point (including output), x in R^(nrec);  
        k_ind    - current piece, k_ind in [0,...,kmax-1];
        m_ind    - current data point, m_ind in [0,...,m-1];
        
    Global parameters:
        a        - two dimensional data matrix;
        nfea     - number of variables (including output)
        jk       - number of linear functions under each minimum;
        min_line - the index of the smallest value, kmax * m matrix;
        
    Output:
        g_cc     - the subgradient of the concave function at x;   
        
    """
    g_cc = np.zeros_like(x)
    if cfg.min_line[k_ind,m_ind] == cfg.jk[k_ind]-1:
        g_cc[0:cfg.nfea-1]=cfg.a[m_ind,:cfg.nfea-1] 
        g_cc[cfg.nfea-1] = 1.0
                     
    return g_cc


# Define df1
def dauxminf1(x):
    """ Subgradient of aux-min-f1.
    
    Input:
        x          - current point (including output), x in R^(n*kmax*jk);  
        
    Global parameters:
        ntrain     - number of data points in the training set;
        alpha1     - index that tells which part of the function rho1 is used;
        max_piece  - the index of the concave piece that gives the maximum, 
                     note that f1 need to be computed at x, m vector;
        
    Output:
        g          - the subgradient of the f1 function at x;   
    
    """
    
    g = np.zeros_like(x)
    for m_ind in range(cfg.ntrain):
        if (cfg.alpha1[m_ind] == 1): # 2*rho1 - b > 2*rho2 + b
            if cfg.nomax-1 != cfg.max_piece[m_ind]:
                g -= dauxmin_cc_piece(x,cfg.nomax-1,m_ind) 
                
        elif (cfg.alpha1[m_ind] == 0): # 2*rho1 - b < 2*rho2 + b
            g -= dauxmin_cc_piece(x,cfg.nomax-1,m_ind) #partial(rho2)
            
        else: # 2*rho1 - b == 2*rho2 + b
            g -= dauxmin_cc_piece(x,cfg.nomax-1,m_ind) #partial(rho2)
                    
    g = 2.0 * g
    
    return g

# Define df2
def dauxminf2(x):
    """ Subgradient of aux-min-f2.
    
    Input:
        x        - current point (including output), x in R^(n*kmax*jk);  
        
    Global parameters:
        ntrain     - number of data points in the training set;
        nomax      - current number of functions under maximum (K);
        max_piece  - the index of the concave piece that gives the maximum,
                    note that f1 need to be computed at x, m vector;
         
    Output:
        g          - the subgradient of the f2 function at x;   
        """
    g = np.zeros_like(x)
    for m_ind in range(cfg.ntrain):
        if cfg.nomax-1 != cfg.max_piece[m_ind]: 
            g -= 2.0*dauxmin_cc_piece(x,cfg.nomax-1,m_ind) 
        else:  
            g -= dauxmin_cc_piece(x,cfg.nomax-1,m_ind)         
    
    return g

# AUX-MAX-FUNCTION ###########################################################################################

# Define aux-max-f1
def auxmaxf1(x):
    """ First DC component of the aux-max-function.
    
    Input:
        x        - current point (including output), x in R^(nfea);  
        
    Global parameters:
        ntrain   - number of data points in training set;
        
    Output:
        f        - the first DC component function value at x;   
        
    """
 
# Sum over data points
    f = 0.0
    for m_ind in range(cfg.ntrain):
        f += auxmax_f1_part_i(x,m_ind) 
        
    return f


# Define aux-max-f2
def auxmaxf2(x):
    """ Second DC component of the aux-max-function.
    
    Input:
        x        - current point (including output), x in R^(nfea);  
        
    Global parameters:
        ntrain   - number of data points in the training set;
        
    Output:
        f        - the second DC component function value at x;   
        
    """
    # Sum over data points
    f = 0.0
    for m_ind in range(cfg.ntrain):
        f += auxmaxrho1(x,m_ind) + auxmaxrho2(x,m_ind) 
        
    return f


# Define aux-max-f1
def auxmax_f1_part_i(x,m_ind):
    """ Subprogram for computing the aux-max-function.
    
    Input:
        x        - current point (including output), x in R^(nfea);  
        m_ind    - current data point;
                
    Global parameters:
        a        - two dimensional data matrix;
        nfea     - number of features (including output);
        alpha1   - index that tells which part of the function auxmaxrho1 is used.
        
    Output:
        f        - part m_ind of the first DC component function value at x;   
        
    """
     
    tmp1 = 2.0*auxmaxrho1(x,m_ind)-cfg.a[m_ind,cfg.nfea-1] 
    tmp2 = 2.0*auxmaxrho2(x,m_ind)+cfg.a[m_ind,cfg.nfea-1]

    # checking the maximum used in auxmaxrho1    
    if (tmp1 > tmp2):
        f = tmp1
        cfg.alpha1[m_ind] = 1  # alpha1 should be ok here. We do not solve aux and real problems at the same time.      
    elif (tmp1 < tmp2):
        f = tmp2
        cfg.alpha1[m_ind] = 0  
    else:
        f = tmp2
        cfg.alpha1[m_ind] = 2  

    return f


def auxmaxrho1(x,m_ind):
    """ Subprogram for computing the aux-max functions.
        
    Input:
        x          - current point (including output), x in R^(nfea);  
        m_ind      - current data point;
        
    Global parameters:
        nomax      - current number of functions (minima and/or linear functions) under maximum;
        max_piece  - the index of the concave piece that gives the maximum, m vector;
    
    Output:
        f          - auxiliary function value  at x;   

    """
    
    cc_sum = auxmaxrho2(x,m_ind) 
    f = cc_sum + auxmax_cc_piece(x,0,m_ind) 
    cfg.max_piece[m_ind] = 0 # max_piece should be ok here. We do not solve aux and real problem at the same time.
             
    for k_ind in range(1,cfg.nomax):
    
        f_tmp = cc_sum + auxmax_cc_piece(x,k_ind,m_ind) 
        if f_tmp > f: 
            f = f_tmp
            cfg.max_piece[m_ind] = k_ind
    
    return f


def auxmaxrho2(x,m_ind):
    """ Subprogram for computing the aux-max-function.
        
    Input:
        x          - current point (including output), x in R^(n*kmax*jk);  
        m_ind      - current data point;
        
    Global parameters:
        nomax      - current number of functions (minima and/or linear functions) under maximum;
        
    Output:
        f          - auxiliary function value  at x for computing the nonsmooth PWLR-L1 function;   

    """
    
    f = 0.0
    for k_ind in range(cfg.nomax):
        f -= auxmax_cc_piece(x,k_ind,m_ind)   

    return f


# Define concave piece for aux-max-function
def auxmax_cc_piece(x,k_ind,m_ind):
    """ Calculates the concave function f_cc = min_{j=1,...,J_k} { (x^{kj})^T a^{i} + y_{kj} }
        only the newest linear function is considered as variable.
    
    Input:
        x        - current point (including output), x in R^(nfea);  
        k_ind    - current piece, k_ind in [0,...,kmax-1];
        m_ind    - current data point, m_ind in [0,...,m-1]
        
    Global parameters:
        xprev    - previous solution;
        a        - two dimensional data matrix;
        nfea     - number of features (including output);
        nomax    - number of max-functions;
        jk       - number of linear functions under each minimum;
        min_line - the index of the smallest value, kmax * m matrix;
    
    Output:
        f_cc     - the concave function value at x;   
        
    """
    
    # Adding new linear function as a last function:
    # The first line. If k_ind = nomax-1, this is a new line, otherwise an old one.
    line_start=cfg.nfea*sum(cfg.jk[i] for i in range(k_ind))
    if cfg.jk[k_ind]==1 and k_ind==cfg.nomax-1: #
        print "hihu0"
        f_cc=np.dot(x[0:cfg.nfea-1],cfg.a[m_ind,:cfg.nfea-1])+x[cfg.nfea-1]
        return f_cc
    else:
        print "hihu1",line_start
        f_cc=np.dot(cfg.xprev[line_start:line_start+(cfg.nfea-1)],cfg.a[m_ind,:cfg.nfea-1])+cfg.xprev[line_start+cfg.nfea-1]
        
    cfg.min_line[k_ind,m_ind] = 0  # a global variable to save the smallest value.
        
    # Next lines
    line_start += cfg.nfea
    for j in range(1,cfg.jk[k_ind]-1): # Everything but the first and last.
        
        f_tmp = np.dot(cfg.xprev[line_start:line_start+(cfg.nfea-1)],cfg.a[m_ind,:cfg.nfea-1])+cfg.xprev[line_start+cfg.nfea-1]
        
        # Minimum of lines
        if f_tmp <= f_cc:
            f_cc = f_tmp
            cfg.min_line[k_ind,m_ind] = j
        line_start += cfg.nfea
    
    
    # The last line.
    if k_ind==cfg.nomax-1:
         
        f_tmp = np.dot(x[0:cfg.nfea-1],cfg.a[m_ind,:cfg.nfea-1])+x[cfg.nfea-1]
    else:    
        
        f_tmp = np.dot(cfg.xprev[line_start:line_start+(cfg.nfea-1)],cfg.a[m_ind,:cfg.nfea-1])+cfg.xprev[line_start+cfg.nfea-1]
            
    # Minimum of lines
    if f_tmp <= f_cc:
        f_cc = f_tmp
        cfg.min_line[k_ind,m_ind] = cfg.jk[k_ind]-1 
    
    
    return f_cc


# Define subgradient of the concave piece for aux-min-function
def dauxmax_cc_piece(x,k_ind,m_ind):
    """ Calculates the subgradient of the concave function 
    
            f_cc = min_{j=1,...,J_k} { (x^{kj})^T a^{i} + y_{kj} }
    
    Input:
        x        - current point (including output), x in R^(nrec);  
        k_ind    - current piece, k_ind in [0,...,kmax-1];
        m_ind    - current data point, m_ind in [0,...,m-1];
        
    Global parameters:
        a        - two dimensional data matrix;
        nfea     - number of variables (including output)
        jk       - number of linear functions under each minimum;
        min_line - the index of the smallest value, kmax * m matrix;
        
    Output:
        g_cc     - the subgradient of the concave function at x;   
        
    """
    g_cc = np.zeros_like(x)
    
    if cfg.min_line[k_ind,m_ind] == cfg.jk[k_ind]-1:
        g_cc[0:cfg.nfea-1]=cfg.a[m_ind,:cfg.nfea-1] 
        g_cc[cfg.nfea-1] = 1.0
                     
    return g_cc


# Define df1
def dauxmaxf1(x):
    """ Subgradient of aux-min-f1.
    
    Input:
        x          - current point (including output), x in R^(n*kmax*jk);  
        
    Global parameters:
        ntrain     - number of data points in the training set;
        alpha1     - index that tells which part of the function rho1 is used;
        max_piece  - the index of the concave piece that gives the maximum, 
                     note that f1 need to be computed at x, m vector;
        
    Output:
        g          - the subgradient of the f1 function at x;   
    
    """
    
    g = np.zeros_like(x)
    for m_ind in range(cfg.ntrain):
        if (cfg.alpha1[m_ind] == 1): # 2*rho1 - b > 2*rho2 + b
            if cfg.nomax-1 != cfg.max_piece[m_ind]:
                g -= dauxmin_cc_piece(x,cfg.nomax-1,m_ind) 
                
        elif (cfg.alpha1[m_ind] == 0): # 2*rho1 - b < 2*rho2 + b
            g -= dauxmin_cc_piece(x,cfg.nomax-1,m_ind) #partial(rho2) 
            
        else: # 2*rho1 - b == 2*rho2 + b (tata ei viela ole)
            g -= dauxmin_cc_piece(x,cfg.nomax-1,m_ind) #partial(rho2) 
                
    g = 2.0 * g
    
    return g

# Define df2
def dauxmaxf2(x):
    """ Subgradient of aux-min-f2.
    
    Input:
        x        - current point (including output), x in R^(n*kmax*jk);  
        
    Global parameters:
        ntrain     - number of data points in the training set;
        nomax      - current number of functions under maximum (K);
        max_piece  - the index of the concave piece that gives the maximum,
                    note that f1 need to be computed at x, m vector;
         
    Output:
        g          - the subgradient of the f2 function at x;   
        """
    g = np.zeros_like(x)
    for m_ind in range(cfg.ntrain):
        if cfg.nomax-1 != cfg.max_piece[m_ind]: 
            g -= 2.0*dauxmin_cc_piece(x,cfg.nomax-1,m_ind) 
        else:  
            g -= dauxmin_cc_piece(x,cfg.nomax-1,m_ind) 
        
    return g


# SCALING ###########################################################################################

# Define scaling
def scaling():
    """ Scaling of the input data
        
    Global parameters:
        a        - two dimensional data matrix (input of scaling);
        a_scaled - scaled input matrix (output of scaling);   
        nfea     - number of features (including output);
        ntrain   - size of the training set;
            
    """
    
    for i in range(cfg.nfea):
        dm = 0
        var = 0
        for j in range(cfg.ntrain):
            dm += cfg.a[j,i]
        dm = dm/cfg.ntrain
        
        for j in range(cfg.ntrain):
            var += (cfg.a[j,i]-dm)**2

        var = var/cfg.ntrain
        var = np.sqrt(var)
            
        if var >= 10**(-5):
            cfg.clin[i] = 1.0/var   
            cfg.dlin[i] = -dm/var  
            
        else:    
            if np.abs(dm)<=1.0:
                cfg.clin[i] = 1.0
                cfg.dlin[i] = 0.0 
            else:    
                cfg.clin[i] = 1.0/dm
                cfg.dlin[i] = 0.0 
                
        for j in range(cfg.ntrain):
            cfg.a_scaled[j,i] = cfg.clin[i]*cfg.a[j,i] + cfg.dlin[i]
            
    return


# Define rescaling
def rescaling(x_solution):
    """ Rescaling of x.

    Input:
        x_solution  - current point (including output variable);
        
    Global parameters:
        nfea        - number of features (including output);
        nomax       - current number of minima (and/or linear functions) under maximum;
        
    Output:
        x_unscaled  - unscaled solution.
        
    """

    x_unscaled = np.zeros_like(x_solution)
    rabs = np.abs(cfg.clin[cfg.nfea-1]) #
    
    if rabs > 10**(-8):
        i = 0
        for k in range(cfg.nomax):  
            for _ in range (cfg.jk[k]): 
                i += 1 
                d0 = 0.0
                for j1 in range (cfg.nfea-1): 
                    x_unscaled[j1+(i-1)*cfg.nfea] = x_solution[j1+(i-1)*cfg.nfea]*cfg.clin[j1]/rabs      
                    d0 += x_solution[j1+(i-1)*cfg.nfea]*cfg.dlin[j1]
                x_unscaled[i*cfg.nfea-1] = (x_solution[i*cfg.nfea-1]+d0-cfg.dlin[cfg.nfea-1])/rabs
                
    return x_unscaled            


# EVALUATION CRITERIA ###########################################################################################

# Define evaluation criteria for training set
def evaltr(x_solution):
    """ Computing RMSE, MAE, CE and r.

    Input:
        x_solution  - current point (including output variable);
        
    Global parameters:
        a_unscaled  - data matrix;
        nfea        - number of features (including output);
        ntrain      - number of points in training set;
        nomax       - number of maxima;
        jk          - number of linear functions under maximum;
        
    Output:
        rmse,mae,ce,r.
        
    """         
                                     
    large = 10.0**30
    pred = np.zeros(cfg.ntrain)
    e0 = 0.0 # mean of observed values
    y=0.0
    for i in range(cfg.ntrain): # Computation of correct piece
        e0 += cfg.a_unscaled[i][-1]
        pind = 0
        ipbest = 0
        pbest = -large # for max
        
        for j1 in range(cfg.nomax):
            ipmin=pind
            pmin=large # for min
            for _ in range(cfg.jk[j1]):
                piece=x_solution[(pind+1)*cfg.nfea-1] 
                for j3 in range(cfg.nfea-1): #
                    piece += x_solution[pind*cfg.nfea+j3]*cfg.a_unscaled[i][j3]
                if piece < pmin:
                    ipmin = pind
                    pmin  = piece
                pind += 1        
            
            if pmin > pbest:
                ipbest = ipmin
                pbest = pmin
                        
        pred[i] = x_solution[(ipbest+1)*cfg.nfea-1] # Computation of prediction
        for j1 in range(cfg.nfea-1):
            pred[i] += x_solution[ipbest*cfg.nfea+j1]*cfg.a_unscaled[i][j1]
        y += pred[i]
        
    y = y/cfg.ntrain 
    e0 = e0/cfg.ntrain
            
    # Computation of indices
    rmse = 0.0
    mae = 0.0
    e1 = 0.0
    for i in range(cfg.ntrain):
        rmse += (pred[i]-cfg.a_unscaled[i][-1])**2
        mae += np.abs(pred[i]-cfg.a_unscaled[i][-1]) 
        e1 += (cfg.a_unscaled[i][-1] - e0)**2
    ce = 1.0 - rmse/e1 
    rmse = np.sqrt(rmse/cfg.ntrain)
    mae = mae/cfg.ntrain 

    if cfg.ntrain > 1:
        sx=0.0
        sy=0.0
        rcor=0.0
        for i in  range(cfg.ntrain):
            sx += (pred[i]-y)**2
            sy += (cfg.a_unscaled[i][-1]-e0)**2  
            rcor += (pred[i]-y) * (cfg.a_unscaled[i][-1]-e0) 

        r = rcor/np.sqrt(sx*sy)
        
    return rmse,mae,ce,r            

# Define evaluation criteria for training set
def evaltest(x_solution,ntest,pred):
    """ Computing the prediction and the indices RMSE, MAE, CE and r.

    Input:
        x_solution  - current point (including output variable);
        ntest       - number of points in test set;
        
    Global parameters:
        a_unscaled  - data matrix;
        nfea        - number of features (including output);
        ntrain      - number of points in training set;
        nomax       - number of maxima;
        jk          - number of linear functions under maximum;
        
    Output:
        pred,rmse,mae,ce,r.
        
    """
        
    large = 10.0**30
    e0 = 0.0
    y=0.0
    for i in range(ntest): # Computation of correct piece
        e0 += cfg.a_unscaled[cfg.ntrain+i][-1]
        pind = 0
        ipbest = 0
        pbest = -large # for max
        
        for j1 in range(cfg.nomax):
            ipmin=pind
            pmin=large # for min
            for _ in range(cfg.jk[j1]):
                piece=x_solution[(pind+1)*cfg.nfea-1] 
                for j3 in range(cfg.nfea-1): #
                    piece += x_solution[pind*cfg.nfea+j3]*cfg.a_unscaled[cfg.ntrain+i][j3]
                if piece < pmin:
                    ipmin = pind
                    pmin  = piece
                pind += 1        
            
            if pmin > pbest:
                ipbest = ipmin
                pbest = pmin
                        
        pred[i] = x_solution[(ipbest+1)*cfg.nfea-1] # Computation of prediction
        for j1 in range(cfg.nfea-1):
            pred[i] += x_solution[ipbest*cfg.nfea+j1]*cfg.a_unscaled[cfg.ntrain+i][j1]
        y += pred[i]
        
    y = y/ntest 
    e0 = e0/ntest
            
    # Computation of indices
    rmse = 0.0
    mae = 0.0
    e1 = 0.0
    for i in range(ntest):
        rmse += (pred[i]-cfg.a_unscaled[cfg.ntrain+i][-1])**2
        mae += np.abs(pred[i]-cfg.a_unscaled[cfg.ntrain+i][-1]) 
        e1 += (cfg.a_unscaled[cfg.ntrain+i][-1] - e0)**2
    ce = 1.0 - rmse/e1 
    rmse = np.sqrt(rmse/ntest)
    mae = mae/ntest

    if ntest > 1:
        sx=0.0
        sy=0.0
        rcor=0.0
        for i in  range(ntest):
            sx += (pred[i]-y)**2
            sy += (cfg.a_unscaled[cfg.ntrain+i][-1]-e0)**2  
            rcor += (pred[i]-y) * (cfg.a_unscaled[cfg.ntrain+i][-1]-e0) 

        r = rcor/np.sqrt(sx*sy)
        
    return pred,rmse,mae,ce,r            
