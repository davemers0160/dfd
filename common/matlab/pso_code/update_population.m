function X = update_population(X, V, N, itr, pso_member_limits)
% - Calculate the new population members and check the population limits
%   - x(k+1) = x(k) + v(k+1)

    for idx=1:N
        X(idx, itr+1).iteration = itr+1;
        X(idx, itr+1).number = X(idx, itr).number;
        X(idx, itr+1).con = X(idx,itr).con + V(idx,itr+1).con;
        X(idx, itr+1).bn = X(idx,itr).bn + V(idx,itr+1).bn;
        X(idx, itr+1).act = X(idx,itr).act + V(idx,itr+1).act;
        X(idx, itr+1).cond = X(idx,itr).cond + V(idx,itr+1).cond;
        X(idx, itr+1).input = X(idx,itr).input + V(idx,itr+1).input;
        X(idx, itr+1).input_file = X(idx,itr).input_file + V(idx,itr+1).input_file;
        X(idx, itr+1).crop_size = X(idx,itr).crop_size + V(idx,itr).crop_size;
        X(idx, itr+1).con_map = X(idx, itr).con_map;
        
        X(idx, itr+1) = validate_pso_member(X(idx, itr+1), pso_member_limits);
    end
    
end
