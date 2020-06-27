function [corr_img] = apply_blur_correction(img, kernel, mask_params)

    [img_h,img_w,~] = size(img);
    
    mask = create_ring_mask([img_h,img_w], mask_params.center, mask_params.range(1), 0);
    mask(mask_params.center(2),mask_params.center(1)) = 1;    % do this becasue of the math for removing the center radius

    tmp_blur_img = imfilter(img, kernel{1}, 'symmetric', 'conv');

    %corr_img = [];
    corr_img = tmp_blur_img.*mask;

    for idx=2:numel(mask_params.range)-1

        tmp_blur_img = imfilter(img, kernel{idx}, 'symmetric', 'conv');

        mask = create_ring_mask([img_h,img_w], mask_params.center, mask_params.range(idx), mask_params.range(idx-1));
        corr_img = corr_img + (tmp_blur_img.*mask);

    end
    
end