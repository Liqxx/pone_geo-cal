import numpy as np

# http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
# Using the above nested radical formula for g=phi_d 
# or you could just hard-code it. 
# phi(1) = 1.6180339887498948482 
# phi(2) = 1.32471795724474602596 

def phi(d): 
    x = 2.0000 
    # numerical approx of root x**(d+1) = x + 1
    # the range gives the precision
    for i in range(10):
        x = pow(1+x, 1./(d+1)) 
    return x
  
def rnew(d, n, seed=0.5):
    '''
    Use quasi-random multi-dimensional sampling
    based on golden-ratio generalization.
    Source: http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
    
    Returns <n> points in <d> dimensions as array of shape (n,d).
    '''
    # get root
    g = phi(d)

    # create coefficients
    alpha = np.zeros(d) 
    for j in range(d): 
        alpha[j] = pow(1/g, j+1) %1 
    z = np.zeros((n, d)) 

    # fill n-d rng table
    for i in range(n): 
        rng = (seed + alpha*(i+1)) %1 
        z[i] = rng 
        
    return z
    

