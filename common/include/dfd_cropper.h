#ifndef DFD_ARRAY_RANDOM_CROPPER_H_
#define DFD_ARRAY_RANDOM_CROPPER_H_

#include <cstdint>
//#include <mutex>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <fstream>


// custom includes
#include "ycrcb_pixel.h"
#include "rot_90.h"


// dlib includes
#include "dlib/threads.h"
#include "dlib/numeric_constants.h"
#include "dlib/image_transforms/interpolation.h"
#include "dlib/rand.h"

// ----------------------------------------------------------------------------------------

struct cropper_stats
{
public:

    uint32_t img_index;     /* index into image dataset */
    uint32_t img_h;         /* image height */
    uint32_t img_w;         /* image width */
    uint32_t x;             /* left corner of the crop */
    uint32_t y;             /* top corner of the crop */
};

// ----------------------------------------------------------------------------------------


class dfd_cropper
{
    
    dlib::chip_dims dims = dlib::chip_dims(32, 32);
    uint32_t scale_x = 1;
    uint32_t scale_y = 1;

    std::vector<cropper_stats> cr_stats;
    std::string recorder_filename = "cropper_stats.bin";
    std::ofstream cropper_recorder;
    bool record_cropper_stats = false;

    dlib::rand rnd;

// ----------------------------------------------------------------------------------------

public:

    void set_seed (time_t seed) { rnd = dlib::rand(seed); }
    
    const dlib::chip_dims& get_chip_dims() const { return dims; }

    void set_chip_dims (const dlib::chip_dims& dims_) { dims = dims_; }

    void set_chip_dims (uint64_t rows, uint64_t cols)
    { set_chip_dims(dlib::chip_dims(rows,cols)); }

    void set_chip_dims (std::pair<uint64_t, uint64_t> p)
    { set_chip_dims(dlib::chip_dims(p.first, p.second)); }
    
    void set_scale_x(uint32_t x) { scale_x = x; }
    void set_scale_y(uint32_t y) { scale_y = y; }

    void set_scale(std::pair<uint32_t, uint64_t> p)
    {
        scale_x = p.second;
        scale_y = p.first;
    }
    
    const uint32_t &get_scale_x() const { return scale_x; }
    const uint32_t &get_scale_y() const { return scale_y; }

    void set_expansion_factor(uint32_t val) { expansion_factor = val; }
    const uint32_t &get_expansion_factor() const { return expansion_factor; }

    std::string get_stats_filename(void) const { return recorder_filename; }
    void set_stats_filename(std::string fn)
    {
        recorder_filename = fn;
        cropper_recorder.open(fn, ios::out | ios::binary);
        record_cropper_stats = cropper_recorder.is_open();
    }

    void close_cropper_stream(void)
    {
        if (record_cropper_stats)
            cropper_recorder.close();
    }
 
// ----------------------------------------------------------------------------------------

        template<typename array_type1, typename array_type2>
        void operator() (
            const uint64_t num_crops,
            const array_type1 &img,
            const array_type2 &gt,
            array_type1& img_crops,
            array_type2& gt_crops
        )
        {

            DLIB_CASSERT(img.size() == gt.size());

            img_crops.clear();
            gt_crops.clear();
            cr_stats.clear();

            append(num_crops, img, gt, img_crops, gt_crops);

        }   // end of operator()    
    
// ----------------------------------------------------------------------------------------

    template<typename array_type1, typename image_type1>
    void append(
        const uint64_t num_crops,
        const array_type1 &img,
        const image_type1 &gt,
        array_type1& img_crops,
        image_type1& gt_crops
    )
    {

        DLIB_CASSERT(img.size() == gt.size());
        DLIB_CASSERT(img_crops.size() == gt_crops.size());

        long original_size = img_crops.size();
        const uint64_t img_count = img.size();

        img_crops.resize(original_size + (num_crops * expansion_factor));
        gt_crops.resize(original_size + (num_crops * expansion_factor));
        cr_stats.resize(original_size + (num_crops));

        dlib::parallel_for(original_size, original_size + num_crops, [&](long idx)
        {
            uint64_t img_index = rnd.get_integer(img_count);

            cr_stats[idx].img_index = img_index;
            cr_stats[idx].img_h = gt[img_index].nr();
            cr_stats[idx].img_w = gt[img_index].nc();

            (*this)(img[img_index], gt[img_index], &img_crops[idx*expansion_factor], &gt_crops[idx*expansion_factor], &cr_stats[idx]);
        });


        save_cropper_stats(cr_stats);

    }   // end of append    
    
// ----------------------------------------------------------------------------------------

    template<typename array_type1, typename image_type1, typename array_type2, typename image_type2, typename cr_struct>
    void operator() (
        const array_type1 &img,
        const image_type1 &gt,
        array_type2 img_crops,
        image_type2 gt_crops,
        cr_struct cr_stats
        )
    {
        uint64_t idx = 0, jdx = 0;
        uint32_t m = 1;
        uint32_t e_f;
        long image_depth = img.size();

        // get a cropping rectangle for a give image
        dlib::rectangle rect_img, rect_gt;        
        make_random_cropping_rect(gt, rect_img, rect_gt);

        cr_stats->x = rect_gt.left();
        cr_stats->y = rect_gt.top();
        
        array_type1 img_t;
        image_type1 gt_t;

        switch (expansion_factor)
        {
            case 4:
                m = 2;
                e_f = 2;
                break;

            case 8:
                m = 1;
                e_f = 4;
                break;

            default:
                m = 2;
                e_f = 1;
                break;
        }
     
        // get the crops
        for (idx = 0; idx < image_depth; ++idx)
        {
            img_t[idx] = dlib::subm(img[idx], rect_img);
        }
        gt_t = dlib::subm(gt, rect_gt);

        // @mem((img_t[0].data).data, UINT16, 1, img_t[0].nc(),img_t[0].nr(),img_t[0].nc()*2)
        // @mem((gt_t.data).data, UINT16, 1, gt_t.nc(),gt_t.nr(),gt_t.nc()*2)

        // create the rotations from the base crop
        array_type1 img_t2;
        for (jdx = 0; jdx < e_f; ++jdx)
        {
            for (idx = 0; idx < image_depth; ++idx)
            {
                img_crops[jdx][idx] = rotate_90(img_t[idx], jdx * m);
            }
            gt_crops[jdx] = rotate_90(gt_t, jdx * m);
        }

        // create the left-right flips from the rotations
        for (jdx = 0; jdx < e_f; ++jdx)
        {
            for (idx = 0; idx < image_depth; ++idx)
            {
                img_crops[jdx + e_f][idx] = dlib::fliplr(img_crops[jdx][idx]);
            }
            gt_crops[jdx + e_f] = dlib::fliplr(gt_crops[jdx]);
        }

    }	// end of operator()    
 
    
// ----------------------------------------------------------------------------------------
private:	

    uint32_t expansion_factor = 8;

    template <typename image_type1>
    void make_random_cropping_rect(const image_type1& img, dlib::rectangle &rect_im, dlib::rectangle &rect_gt)
    {
        uint64_t x = 0, y = 0;
        
        rect_im = dlib::resize_rect(rect_im, dims.cols, dims.rows);
        rect_gt = dlib::resize_rect(rect_gt, (long)(dims.cols/(double)scale_x), (long)(dims.rows/(double)scale_y));
       
        if ((unsigned long)img.nc() <= rect_gt.width())
            x = 0;
        else
            x = (uint64_t)(rnd.get_integer(img.nc() - rect_gt.width()));

        if ((unsigned long)img.nr() <= rect_gt.height())
            y = 0;
        else
            y = (uint64_t)(rnd.get_integer(img.nr() - rect_gt.height()));
            
        // randomly shift the box around
        dlib::point tr_off(x*scale_x, y*scale_y);
        rect_im = dlib::move_rect(rect_im, tr_off);
        
        dlib::point gt_off(x, y);
        rect_gt = dlib::move_rect(rect_gt, gt_off);

        
    }	// end of make_random_cropping_rect    

    // ----------------------------------------------------------------------------------------

    void save_cropper_stats(std::vector<cropper_stats> &cr_stats)
    {

        uint64_t idx;

        if (record_cropper_stats)
        {
            for (idx = 0; idx < cr_stats.size(); ++idx)
            {
                cropper_recorder.write(reinterpret_cast<const char*>(&cr_stats[idx].img_index), sizeof(cr_stats[idx].img_index));
                cropper_recorder.write(reinterpret_cast<const char*>(&cr_stats[idx].img_h), sizeof(cr_stats[idx].img_h));
                cropper_recorder.write(reinterpret_cast<const char*>(&cr_stats[idx].img_w), sizeof(cr_stats[idx].img_w));
                cropper_recorder.write(reinterpret_cast<const char*>(&cr_stats[idx].x), sizeof(cr_stats[idx].x));
                cropper_recorder.write(reinterpret_cast<const char*>(&cr_stats[idx].y), sizeof(cr_stats[idx].y));
            }
        }

    }   // end of save_cropper_stats

    // ----------------------------------------------------------------------------------------
         
};	// end of class

// ----------------------------------------------------------------------------------------
inline std::ostream& operator<< (
    std::ostream& out,
    const dfd_cropper& item
    )
{
    using std::endl;
    out << "dfd_cropper details: " << std::endl;
    out << "  chip_dims.rows:       " << item.get_chip_dims().rows << std::endl;
    out << "  chip_dims.cols:       " << item.get_chip_dims().cols << std::endl;
    out << "  scale_x:              " << item.get_scale_x() << std::endl;
    out << "  scale_y:              " << item.get_scale_y() << std::endl;
    out << "  expansion_factor:     " << item.get_expansion_factor() << std::endl;
    return out;
}

// ----------------------------------------------------------------------------------------	


#endif	// DFD_ARRAY_RANDOM_CROPPER_H_ 
