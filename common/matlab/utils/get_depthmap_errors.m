function [gt_error_value, combined_hist, cm] = get_depthmap_errors(gt_img, dnn_img, max_depthmap_value)

    dnn_img = double(dnn_img);
    dnn_sz = size(dnn_img);
    gt_img = double(gt_img(1:dnn_sz(1), 1:dnn_sz(2)));

    %% calculate the histogram of the depthmap values
    
    dm_tmp = dnn_img(:);
    gt_tmp = gt_img(:);

    dm_hist = zeros(max_depthmap_value+1,1);
    gt_hist = zeros(max_depthmap_value+1,1);
    
    for idx=1:max_depthmap_value+1
        dm_hist(idx,1) = sum(dm_tmp==(idx-1));   
        gt_hist(idx,1) = sum(gt_tmp==(idx-1));
    end
    combined_hist = cat(2,gt_hist,dm_hist);
    
    %% get the differences between the depthmaps
    error_vals = [1, 2, 3, 4, 5];
    
    gt_img_l = gt_img(:);
    dnn_img_l = dnn_img(:);
    gt_diff_l = dnn_img_l - gt_img_l;
    
    gt_diff = dnn_img - gt_img;
    
    gt_error_value = cell(max_depthmap_value+1,1);
    gt_value_sum = zeros(max_depthmap_value+1,1);
    
    
    dm_error_percent = zeros(max_depthmap_value+1, numel(error_vals));
    dm_value_sum = zeros(max_depthmap_value+1,1);
    dm_error_sum = zeros(numel(error_vals), numel(error_vals));
    
    for idx=1:max_depthmap_value+1
        
        % find the indexes where maps have the depthmap value
        gt_index = gt_img_l == (idx-1);
        dm_index = dnn_img_l == (idx-1);
        
        % get the total number of depthmap values for the gt and dnn
        gt_value_sum(idx) = sum(gt_index);
        dm_value_sum(idx) = sum(dm_index);
        
        gt_error_value{idx, 1} = gt_diff_l(gt_index);
        
        gt_error_hist(idx,:) = histcounts(gt_error_value{idx, 1}, [-max_depthmap_value:1:max_depthmap_value]);
        
        % get the error percentage on a per depthmap value
%         for jdx=1:numel(error_vals)
%             
%             dm_err_index = gt_error_value{idx,1} == error_vals(jdx);
%             dm_error_sum(idx, jdx) = sum(dm_err_index);
%             
%             if(gt_value_sum(idx) == 0)
%                 
%             else
%                 dm_error_percent(idx,jdx) = (dm_error_sum(idx, jdx)/gt_value_sum(idx));
%             end
%             
%         end        
    end
    
    %cm = confusionchart(gt_img_l, dnn_img_l,'RowSummary','row-normalized','ColumnSummary','column-normalized');
    cm = confusionmat(gt_img_l, dnn_img_l, 'Order', [0:1:max_depthmap_value]);
end
