#ifndef _VECT_2_MATRIX_H
#define _VECT_2_MATRIX_H

#include <cstdint>
#include <vector>
#include <algorithm>

// dlib includes
#include <dlib/matrix.h>
#include <dlib/image_transforms.h>


//-----------------------------------------------------------------------------
template <typename img_type, typename gt_type>
inline void vect2matrix(
    uint32_t img_h,
    uint32_t img_w,
    std::vector<uint8_t>& fp1_ptr,
    std::vector<uint8_t>& fp2_ptr,
    std::vector<uint8_t>& dm_ptr,
    std::array<dlib::matrix<img_type>, img_depth>& t,
    dlib::matrix<gt_type>& gt
)
{
    int index = 0, dm_index = 0;
    for (long r = 0; r < img_h; ++r)
    {
        for (long c = 0; c < img_w; ++c)
        {

            dlib::assign_pixel(t[2](r, c), (img_type)fp1_ptr[index]);
            dlib::assign_pixel(t[5](r, c), (img_type)fp2_ptr[index++]);
            dlib::assign_pixel(t[1](r, c), (img_type)fp1_ptr[index]);
            dlib::assign_pixel(t[4](r, c), (img_type)fp2_ptr[index++]);
            dlib::assign_pixel(t[0](r, c), (img_type)fp1_ptr[index]);
            dlib::assign_pixel(t[3](r, c), (img_type)fp2_ptr[index++]);

            dlib::assign_pixel(gt(r, c), (gt_type)dm_ptr[dm_index++]);
        }
    }
}

#endif  // _VECT_2_MATRIX_H
