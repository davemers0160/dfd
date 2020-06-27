function [data_name, num_trials, nmae, nrmse, ssim, nmae_stats, nrmse_stats, ssim_stats] = get_dfd_results(file_name, index)

    % file_name: the full file name of the txt file that contains the
    % results
    % index: this is used based on the type of dfd file that is being read
    % in.  used because files may have silog, depth map variance.  Only
    % looking for nmae, nrmse and ssim for training and testing
    
    data_params = parse_input_parameters(file_name);
    
    % this is a check for a comma in the name only
    if(numel(data_params{1}) > 1)
        data_name = '';
        for idx=1:numel(data_params{1})-1
            data_name = strcat(data_name,data_params{1}{idx},',');
        end
        data_name = strcat(data_name,data_params{1}{end});
    else
        data_name = data_params{1}{1};
    end
    
    data_params(1) = [];

    num_trials = size(data_params,1);

    nmae = zeros(num_trials,2);
    nrmse = zeros(num_trials,2);
    ssim = zeros(num_trials,2);

    for idx=1:num_trials
        nmae(idx,1) = str2double(data_params{idx}{index(1)});
        nrmse(idx,1) = str2double(data_params{idx}{index(2)});
        ssim(idx,1) = str2double(data_params{idx}{index(3)});
        nmae(idx,2) = str2double(data_params{idx}{index(4)});
        nrmse(idx,2) = str2double(data_params{idx}{index(5)});
        ssim(idx,2) = str2double(data_params{idx}{index(6)}); 
    end

    nmae_stats = zeros(4,2);
    nrmse_stats = zeros(4,2);
    ssim_stats = zeros(4,2);

    nmae_stats(1,1) = min(nmae(:,1));
    nmae_stats(2,1) = mean(nmae(:,1));
    nmae_stats(3,1) = max(nmae(:,1));
    nmae_stats(4,1) = std(nmae(:,1));

    nrmse_stats(1,1) = min(nrmse(:,1));
    nrmse_stats(2,1) = mean(nrmse(:,1));
    nrmse_stats(3,1) = max(nrmse(:,1));
    nrmse_stats(4,1) = std(nrmse(:,1));

    ssim_stats(1,1) = min(ssim(:,1));
    ssim_stats(2,1) = mean(ssim(:,1));
    ssim_stats(3,1) = max(ssim(:,1));
    ssim_stats(4,1) = std(ssim(:,1));

    nmae_stats(1,2) = min(nmae(:,2));
    nmae_stats(2,2) = mean(nmae(:,2));
    nmae_stats(3,2) = max(nmae(:,2));
    nmae_stats(4,2) = std(nmae(:,2));

    nrmse_stats(1,2) = min(nrmse(:,2));
    nrmse_stats(2,2) = mean(nrmse(:,2));
    nrmse_stats(3,2) = max(nrmse(:,2));
    nrmse_stats(4,2) = std(nrmse(:,2));

    ssim_stats(1,2) = min(ssim(:,2));
    ssim_stats(2,2) = mean(ssim(:,2));
    ssim_stats(3,2) = max(ssim(:,2));
    ssim_stats(4,2) = std(ssim(:,2));
    
end