#ifndef EVAL_DFD_NET_PERFORMANCE
#define EVAL_DFD_NET_PERFORMANCE

#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <array>

// Custom includes
#include "ssim.h"
#include "center_cropper.h"
#include "dlib_matrix_threshold.h"
#include "calc_silog_error.h"

// DLIB includes
#include <dlib/dnn.h>
#include <dlib/matrix.h>
#include <dlib/image_io.h>
#include <dlib/data_io.h>
#include <dlib/image_transforms.h>

//-----------------------------------------------------------------------------
// This function assumes that the ground truth image size and the input image sizes do not have to be the same dimensions
//-----------------------------------------------------------------------------

template <typename net_type>
dlib::matrix<double,1,6> eval_net_performance(
    net_type &net,
    std::array<dlib::matrix<uint16_t>, img_depth> &img_in,
    dlib::matrix<uint16_t> &gt_in,
    dlib::matrix<uint16_t> &map_out, 
    std::pair<uint64_t, uint64_t> crop_size
    //std::pair<uint64_t, uint64_t> scale,
    //float dm_scale = 1.0
)
{
    std::array<dlib::matrix<uint16_t>, img_depth> img_crop;
    dlib::matrix<uint16_t> gt;
    dlib::matrix<uint16_t> gt_crop;
    dlib::matrix<float> gt_cropf, dm_f;

    double ssim_val = 0.0;
    double nrmse_val = 0.0;
    double nmae_val = 0.0;
    double silog_val = 0.0;
    double v_gt = 0.0;
    double v_dm = 0.0;

    uint16_t gt_min = 0;
    uint16_t gt_max = (uint16_t)((dlib::layer<1>(net).layer_details().num_filters() - 1));

    // threshold the ground truth to remove the ignore values from the training
    // this is just incase there are values that are greater than gt_max
    //truncate_threshold(gt_in, gt, gt_max);

    // center crop the ground truth image and use the crop dims to figure out where to crop the input image
    dlib::rectangle rect_gt = get_center_crop_rect(gt_in, crop_size.second, crop_size.first);
    gt_crop = dlib::subm(gt_in, rect_gt);
    gt_cropf = dlib::matrix_cast<float>(gt_crop);

    // get the input image crop info
    //dlib::rectangle rect_td(crop_size.second * scale.second, crop_size.first * scale.first);
    //dlib::rectangle rect_td(rect_gt.left() * scale.second, rect_gt.top() * scale.first, rect_gt.right() * scale.second, rect_gt.bottom() * scale.first);
    dlib::rectangle rect_td(rect_gt.left(), rect_gt.top(), rect_gt.right(), rect_gt.bottom());


    // shift the box around
    //dlib::point offset(rect_gt.left()*scale.second, rect_gt.top()*scale.first);
    //rect_td = dlib::move_rect(rect_td, offset);

    // crop the input image
    for (uint32_t idx = 0; idx<img_depth; ++idx)
    {
        img_crop[idx] = dlib::subm(img_in[idx], rect_td);
    }

    //dlib::find_min_and_max(gt, gt_min, gt_max);

    // get the output results for the network
    map_out = net(img_crop);
    dm_f = dlib::matrix_cast<float>(map_out);

    // calculate the scale invariant log error 
    //silog_val = calc_silog_error(gt_crop, map_out);
    silog_val = calc_silog_error(gt_cropf, dm_f);

    // subtract the two maps
    //dlib::matrix<float> sub_map = dlib::matrix_cast<float>(map_out) - dlib::matrix_cast<float>(gt_crop);
    dlib::matrix<float> sub_map = dm_f - gt_cropf;

    double m1 = dlib::mean(dlib::abs(sub_map));
    double m2 = dlib::mean(dlib::squared(sub_map));

    double rng = (double)std::max(gt_max - gt_min, 1);
    nmae_val = m1 / rng;
    nrmse_val = std::sqrt(m2) / rng;

    dlib::matrix<float> ssim_map;
    //ssim_val = ssim(map_out, gt_crop, ssim_map);
    ssim_val = ssim(dm_f, gt_cropf, ssim_map);

    //v_gt = (double)dlib::variance(dlib::matrix_cast<float>(gt_crop));
    //v_dm = (double)dlib::variance(dlib::matrix_cast<float>(map_out));
    v_gt = (double)dlib::variance(gt_cropf);
    v_dm = (double)dlib::variance(dm_f);

    dlib::matrix<double, 1, 6> res = dlib::zeros_matrix<double>(1, 6);
    res = nmae_val, nrmse_val, ssim_val, silog_val, v_gt, v_dm;
    
    return res;
   
}   // end of eval_net_performance

//-----------------------------------------------------------------------------

template <typename net_type>
dlib::matrix<double, 1, 6> eval_all_net_performance(
    net_type &net,
    std::vector<std::array<dlib::matrix<uint16_t>, img_depth>> &img_in,    
    std::vector<dlib::matrix<uint16_t>> &gt_in,
    std::pair<uint64_t, uint64_t> crop_size
    //std::pair<uint32_t, uint32_t> scale
)
{
    uint32_t idx;
    DLIB_CASSERT(img_in.size() == gt_in.size());

    dlib::matrix<uint16_t> map;
    dlib::matrix<double, 1, 6> results = dlib::zeros_matrix<double>(1,6);

    // cycle through each input image and evaluate
    for (idx = 0; idx < img_in.size(); ++idx)
    {
        results += eval_net_performance(net, img_in[idx], gt_in[idx], map, crop_size);
    }
  
    // return the average values for each metric
    return (results / (double)img_in.size());

}   // end of eval_all_net_performance

//-----------------------------------------------------------------------------
#endif  // EVAL_DFD_NET_PERFORMANCE
