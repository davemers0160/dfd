import numpy as np
import math

def coc_calc(f_num, f, d_o, limits, step):
   
    CoC_max = (f*f/(f_num*(d_o-f)))

    start = max(limits[0],step)

    r = np.arange(start, limits[1], step)
    CoC = np.zeros(r.shape, dtype=float, order='C')
    tl = (d_o*f*f)/(f_num*(d_o-f))
    
    for d in range(r.shape[0]):
        
        
        if r[d] < d_o:
            CoC[d] = tl*((1.0/r[d]) - (1.0/d_o))
        elif r[d] > d_o:
            CoC[d] = tl*((1.0/d_o)-(1.0/r[d]))        
        # else:
            # CoC[d] = 0
    
    
    
    
    # d_far = d_o:d_step:(limits(2));
    # d_near = limits(1):d_step:d_o-d_step;
    
    # tl = (d_o*f*f)/(f_num*(d_o-f));
    
    # coc_far = tl*((1/d_o)-(1./d_far));
    # coc_near = tl*((1./d_near) - (1/d_o));
    
    # CoC = [coc_near coc_far];  
    # range = [d_near d_far];
    
    
    return r, CoC, CoC_max