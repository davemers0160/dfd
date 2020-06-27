function pso_member = validate_pso_member(pso_member, pso_member_limits)
    
    pso_member.con = max(floor(pso_member.con(:,1:3)+0.5), pso_member_limits.con(:,[1,3,5]));
    pso_member.con = min(pso_member.con(:,1:3), pso_member_limits.con(:,[2,4,6]));
    
    pso_member.bn = max(floor(pso_member.bn(:)+0.5), pso_member_limits.bn(:,1));
    pso_member.bn = min(pso_member.bn(:), pso_member_limits.bn(:,2));

    %pso_member.act = max(pso_member.act(:), pso_member_limits.act(:,1));
    %pso_member.act = min(pso_member.act(:), pso_member_limits.act(:,2));
    
    pso_member.act = mod(floor(pso_member.act(:)-0.5), pso_member_limits.act(:,2)) + 1;
    
    pso_member.cond = mod(floor(pso_member.cond-0.5), pso_member_limits.cond(:,2)) + 1;
    
    %pso_member.input = max(pso_member.input, pso_member_limits.input(:,1));
    %pso_member.input = min(pso_member.input, pso_member_limits.input(:,2));
    
    pso_member.input = mod(floor(pso_member.input-0.5), pso_member_limits.input(:,2)) + 1;
    
    pso_member.input_file = mod(floor(pso_member.input_file-0.5), pso_member_limits.input_file(:,2)) + 1;
    
    %pso_member.crop_size = max(pso_member.crop_size, pso_member_limits.crop_size(:,1));
    %pso_member.crop_size = min(pso_member.crop_size, pso_member_limits.crop_size(:,2));

    pso_member.crop_size = mod(floor(pso_member.crop_size+0.5-pso_member_limits.crop_size(:,1)), pso_member_limits.crop_size(:,2)-pso_member_limits.crop_size(:,1)+1) + pso_member_limits.crop_size(:,1);
    
end