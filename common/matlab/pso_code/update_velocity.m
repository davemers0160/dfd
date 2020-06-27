function V = update_velocity(X, V, P, G, c1, c2, kap, N, itr, pso_velocity_limits)
% - Calculate the updated velocity vector and then check the velocity limits
%   - v(k+1) = kap*(v(k) + c1*r(k) o (P(k) - x(k)) + c2*s(k) o (G - x(k)))

    for idx=1:N
        R = rand(size(X(idx,itr).con));
        S = rand(size(X(idx,itr).con));
        V(idx, itr+1).con = kap * (V(idx,itr).con + (c1*R).*(P(idx,itr).con - X(idx,itr).con) + (c2*S).*(G(itr).con - X(idx,itr).con));
        
        R = rand(size(X(idx,itr).bn));
        S = rand(size(X(idx,itr).bn));
        V(idx, itr+1).bn = kap * (V(idx,itr).bn + (c1*R).*(P(idx,itr).bn - X(idx,itr).bn) + (c2*S).*(G(itr).bn - X(idx,itr).bn));
        
        R = rand(size(X(idx,itr).act));
        S = rand(size(X(idx,itr).act));
        V(idx, itr+1).act = kap * (V(idx,itr).act + (c1*R).*(P(idx,itr).act - X(idx,itr).act) + (c2*S).*(G(itr).act - X(idx,itr).act));
        
        R = rand(size(X(idx,itr).cond));
        S = rand(size(X(idx,itr).cond));
        V(idx, itr+1).cond = kap * (V(idx,itr).cond + (c1*R).*(P(idx,itr).cond - X(idx,itr).cond) + (c2*S).*(G(itr).cond - X(idx,itr).cond));
        
        R = rand(size(X(idx,itr).input));
        S = rand(size(X(idx,itr).input));
        V(idx, itr+1).input = kap * (V(idx,itr).input + (c1*R).*(P(idx,itr).input - X(idx,itr).input) + (c2*S).*(G(itr).input - X(idx,itr).input));
        
        R = rand(size(X(idx,itr).input_file));
        S = rand(size(X(idx,itr).input_file));
        V(idx, itr+1).input_file = kap * (V(idx,itr).input_file + (c1*R).*(P(idx,itr).input_file - X(idx,itr).input_file) + (c2*S).*(G(itr).input_file - X(idx,itr).input_file)); 
        
        R = rand(size(X(idx,itr).crop_size));
        S = rand(size(X(idx,itr).crop_size));
        V(idx, itr+1).crop_size = kap * (V(idx,itr).crop_size + (c1*R).*(P(idx,itr).crop_size - X(idx,itr).crop_size) + (c2*S).*(G(itr).crop_size - X(idx,itr).crop_size));

        V(idx, itr+1) = validate_pso_velocity(V(idx, itr+1), pso_velocity_limits);
    end
       
end
