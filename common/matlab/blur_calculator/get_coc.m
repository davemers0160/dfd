function [err] = get_coc(x)
    global pixel px_size range
    
    if (numel(x) > 3)
        f_num = x(1,:);
        f = x(2,:);
        d_o = x(3,:);
    else
        f_num = x(1);
        f = x(2);
        d_o = x(3);
    end
    
    if(numel(d_o) > 1)
        parfor idx=1:numel(d_o)
            dn = (range <= d_o(idx));

            d_near = range(dn);
            d_far = range(~dn);

            tl = (d_o(idx)*f(idx)*f(idx))/(f_num(idx)*(d_o(idx)-f(idx)));

            coc_far = tl*((1/d_o(idx))-(1./d_far));
            coc_near = tl*((1./d_near) - (1/d_o(idx)));

            CoC = [coc_near coc_far];  

            px = ceil(CoC/px_size);
            %px = ceil(px);

            err(idx,1) = sum(abs(pixel - px));
        end

    else
        dn = (range <= d_o);

        d_near = range(dn);
        d_far = range(~dn);

        tl = (d_o*f*f)/(f_num*(d_o-f));

        coc_far = tl*((1/d_o)-(1./d_far));
        coc_near = tl*((1./d_near) - (1/d_o));

        CoC = [coc_near coc_far];  

        px = ceil(CoC/px_size);

        err = sum(abs(pixel - px));    
    end
    
end
