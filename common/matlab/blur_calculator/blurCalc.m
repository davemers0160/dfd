function [d_range, CoC, CoC_max] = blurCalc(f_num, f, d_o, limits, d_step)

    % set the range step
    %d_step = 100;
    
    d_far = d_o:d_step:(limits(2));
    d_near = limits(1):d_step:d_o-d_step;
    
    tl = (d_o*f*f)/(f_num*(d_o-f));
    
    coc_far = tl*((1/d_o)-(1./d_far));
    coc_near = tl*((1./d_near) - (1/d_o));
    
    CoC = [coc_near coc_far];  
    d_range = [d_near d_far];
    
    CoC_max = (f*f/(f_num*(d_o-f)));

end


