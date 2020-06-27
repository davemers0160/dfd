function pso_velocity = validate_pso_velocity(pso_velocity, pso_velocity_limits)
    
    pso_velocity.con = max(pso_velocity.con(:,1:3), pso_velocity_limits.con(:,[1,3,5]));
    pso_velocity.con = min(pso_velocity.con(:,1:3), pso_velocity_limits.con(:,[2,4,6]));
    
    pso_velocity.bn = max(pso_velocity.bn(:), pso_velocity_limits.bn(:,1));
    pso_velocity.bn = min(pso_velocity.bn(:), pso_velocity_limits.bn(:,2));

    pso_velocity.act = max(pso_velocity.act(:), pso_velocity_limits.act(:,1));
    pso_velocity.act = min(pso_velocity.act(:), pso_velocity_limits.act(:,2));
    
    pso_velocity.cond = max(pso_velocity.cond, pso_velocity_limits.cond(:,1));
    pso_velocity.cond = min(pso_velocity.cond, pso_velocity_limits.cond(:,2));
    
    pso_velocity.input = max(pso_velocity.input, pso_velocity_limits.input(:,1));
    pso_velocity.input = min(pso_velocity.input, pso_velocity_limits.input(:,2));
    
    pso_velocity.input_file = max(pso_velocity.input_file, pso_velocity_limits.input_file(:,1));
    pso_velocity.input_file = min(pso_velocity.input_file, pso_velocity_limits.input_file(:,2));
    
    pso_velocity.crop_size = max(pso_velocity.crop_size, pso_velocity_limits.crop_size(:,1));
    pso_velocity.crop_size = min(pso_velocity.crop_size, pso_velocity_limits.crop_size(:,2));
    
end