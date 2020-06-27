function [blur_radius, coc_max, coc_near, coc_far] = calc_blur_radius(d_o, f, f_num, range)

    dn = (range <= d_o);
    
    d_near = range(dn);
    d_far = range(~dn);
    
    tl = (d_o*f*f)/(f_num*(d_o-f));
    
    coc_far = tl*((1/d_o)-(1./d_far));
    coc_near = tl*((1./d_near) - (1/d_o));
    
    blur_radius = [coc_near coc_far]; 
    
    coc_max = (f*f/(f_num*(d_o-f)));
end
