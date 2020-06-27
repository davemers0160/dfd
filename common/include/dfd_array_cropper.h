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

    uint64_t img_index;     /* index into image dataset */
    uint32_t img_h;         /* image height */
    uint32_t img_w;         /* image width */
    uint32_t x;             /* left corner of the crop */
    uint32_t y;             /* top corner of the crop */
};

// ----------------------------------------------------------------------------------------

class dfd_cropper
{
    
    dlib::chip_dims dims = dlib::chip_dims(32,32);

    //std::mutex rnd_mutex;
    dlib::rand rnd;

    std::vector<cropper_stats> cr_stats;
    std::string recorder_filename = "cropper_stats.bin";
    std::ofstream cropper_recorder;
    bool record_cropper_stats = false;

// ----------------------------------------------------------------------------------------

public:

    void set_seed (time_t seed) { rnd = dlib::rand(seed); }
    
    const dlib::chip_dims& get_chip_dims() const { return dims; }

    void set_chip_dims (const dlib::chip_dims& dims_) { dims = dims_; }

    void set_chip_dims (uint64_t rows, uint64_t cols)
    { set_chip_dims(dlib::chip_dims(rows,cols)); }

    void set_chip_dims (std::pair<uint64_t, uint64_t> p)
    { set_chip_dims(dlib::chip_dims(p.first, p.second)); }
 
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
            const long num_crops,
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

        img_crops.resize(original_size + num_crops*8);
        gt_crops.resize(original_size + num_crops*8);
        cr_stats.resize(original_size + (num_crops));

        const uint64_t img_count = img.size();

        //for (uint64_t idx = original_size; idx < original_size + (num_crops*8); idx+=8)
        //{
        //    //uint32_t img_index = rnd.get_random_32bit_number() % img_count;
        //    uint32_t img_index = rnd.get_integer(img_count);
        //    (*this)(img[img_index], gt[img_index], &img_crops[idx], &gt_crops[idx]);
        //}

        dlib::parallel_for(original_size, original_size + num_crops, [&](long idx) 
        {
            uint64_t img_index = rnd.get_integer(img_count);

            cr_stats[idx].img_index = img_index;
            cr_stats[idx].img_h = gt[img_index].nr();
            cr_stats[idx].img_w = gt[img_index].nc();
            
            (*this)(img[img_index], gt[img_index], &img_crops[idx * 8], &gt_crops[idx * 8], &cr_stats[idx]);
            //(*this)(img[img_index], gt[img_index], &img_crops[idx * 8], &gt_crops[idx * 8]);
        });

        save_cropper_stats(cr_stats);

    }   // end of append    
    
// ----------------------------------------------------------------------------------------

    //template<typename array_type1, typename image_type1, typename array_type2, typename image_type2, typename cr_struct>
    template<typename array_type1, typename image_type1, typename array_type2, typename image_type2>
    void operator() (
        const array_type1 &img,
        const image_type1 &gt,
        array_type2 img_crops,
        image_type2 gt_crops,/*,
        cr_struct cr_stats*/
        cropper_stats *cr_stats
        )
    {
        uint64_t idx = 0;
        long image_depth = img.size();

        // get a cropping rectangle for a give image
        dlib::rectangle crop_rect = make_random_cropping_rect(img[0]);
        array_type1 img_t;
        image_type1 gt_t;

        cr_stats->x = crop_rect.left();
        cr_stats->y = crop_rect.top();

        // get the crops
        for (idx = 0; idx < image_depth; ++idx)
        {
            img_t[idx] = dlib::subm(img[idx], crop_rect);
        }
        gt_t = dlib::subm(gt, crop_rect);

        // create the rotations from the base crop
        array_type1 img_t2;
        for (uint64_t jdx = 0; jdx < 4; ++jdx)
        {
            for (idx = 0; idx < image_depth; ++idx)
            {
                img_crops[jdx][idx] = rotate_90(img_t[idx], jdx);
            }
            gt_crops[jdx] = rotate_90(gt_t, jdx);
        }

        // create the left-right flips from the rotations
        for (uint64_t jdx = 0; jdx < 4; ++jdx)
        {
            for (idx = 0; idx < image_depth; ++idx)
            {
                img_crops[jdx+4][idx] = dlib::fliplr(img_crops[jdx][idx]);
            }
            gt_crops[jdx+4] = dlib::fliplr(gt_crops[jdx]);
        }
            
    }	// end of operator()    
 
    
// ----------------------------------------------------------------------------------------
private:	

    template <typename image_type>
    dlib::rectangle make_random_cropping_rect(const image_type& img)
    {
        
        //dlib::const_image_view<image_type> img(img_);			
        dlib::rectangle rect(dims.cols, dims.rows);

        // randomly shift the box around
        dlib::point offset(rnd.get_random_32bit_number() % (img.nc() - rect.width() - 1), rnd.get_random_32bit_number() % (img.nr() - rect.height() - 1));
        return dlib::move_rect(rect, offset);
        
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
    out << "  chip_dims.rows:            " << item.get_chip_dims().rows << std::endl;
    out << "  chip_dims.cols:            " << item.get_chip_dims().cols << std::endl;
    return out;
}

// ----------------------------------------------------------------------------------------	


#endif	// DFD_ARRAY_RANDOM_CROPPER_H_        
