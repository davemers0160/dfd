function print_half_kernel(sigma, kernel_size, dim)

    fprintf('float[,] kernel = new float[,] {\n');

    for idx=1:numel(sigma)
        if(dim == 1)
            kernel = single(create_1D_gauss_kernel(kernel_size, sigma(idx)));            
        else
            kernel = single(create_gauss_kernel(kernel_size, sigma(idx)));
        end
        
        fprintf('{');
        str = '';
        for jdx=floor(kernel_size/2+1):kernel_size
            if(dim == 1)
                % 1-D case
                str = strcat(str, num2str(kernel(jdx), '%11.10ff, '));
            else
                % 2-D case
                str = strcat(str, num2str(kernel(floor(kernel_size/2+1), jdx), '%11.10ff, '));
            end
        end
        str = strcat(str(1:end-1),'},');
        fprintf('%s\n', str);
    
    end
    fprintf('};\n');
    
end