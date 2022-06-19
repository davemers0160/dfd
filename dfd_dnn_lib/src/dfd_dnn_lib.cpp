#define _CRT_SECURE_NO_WARNINGS

#include <cstdint>
#include <string>
#include <vector>

// Custom includes
//#include "file_ops.h"
#include "dfd_dnn_lib.h"
//#include "prune_detects.h"

// network header file
#include "dfd_net_v16.h"

// dlib includes
#include <dlib/dnn.h>
#include <dlib/data_io.h>
#include <dlib/image_transforms.h>

//----------------------------------------------------------------------------------
// library internal state variables:
adfd_net_type net;
// double pyr_scale;
// unsigned long outer_padding;
// unsigned long padding;
// std::vector<std::string> class_names;
// std::vector<dlib::rgb_pixel> class_color;

//----------------------------------------------------------------------------------
void init_net(const char *net_name, unsigned int *num_classes, struct window_struct* &det_win, unsigned int *num_win)
{
    int gpu = 0;

    dlib::cuda::set_device(gpu);

    dlib::deserialize(net_name) >> net;

}   // end of init_net

//----------------------------------------------------------------------------------
void run_net(unsigned char* image_f1, unsigned char* image_f2, unsigned int nr, unsigned int nc, unsigned unsigned* depth_map)
{

    uint64_t r, c;
    uint64_t idx = 0;
    uint64_t index = 0;
    std::string label;
    dlib::rgb_pixel p;
    
    //dlib::matrix<dlib::rgb_pixel> img(nr, nc);
    //std::array<dlib::matrix<uint8_t>, array_depth> a_img;

    std::array<dlib::matrix<uint8_t>, img_depth> img;
    
    // get the images size and resize the t array
    for (idx = 0; idx < img_depth; ++idx)
    {
        img[idx].set_size(nr, nc);
    }
    
    uint64_t size = (uint64_t)nr * (uint64_t)nc;

    try 
    {
        // copy the pointer into the input image.  The format is assumed to be row-major order
        // and if the array_depth (number of channels) is greater than 1 then channels are not interleaved
        for (idx = 0; idx < img_depth; ++idx)
        {
            img[idx].set_size(nr, nc);
        }
        
        for (r = 0; r < f.nr(); ++r)
        {
            for (c = 0; c < f.nc(); ++c)
            {
                //dlib::assign_pixel(p, f(r, c));
                dlib::assign_pixel(img[0](r, c), *image_f1++);
                dlib::assign_pixel(img[1](r, c), *image_f1++);
                dlib::assign_pixel(img[2](r, c), *image_f1++);
                //dlib::assign_pixel(p, d1(r, c));
                dlib::assign_pixel(img[3](r, c), *image_f2++);
                dlib::assign_pixel(img[4](r, c), *image_f2++);
                dlib::assign_pixel(img[5](r, c), *image_f2++);
            }
        }

        dlib::matrix<uint16_t> dm = net(img);


        det_img = new unsigned char[tmp_img.nr() * tmp_img.nc() * 3L];

        idx = 0;
        for (r = 0; r < nr; ++r)
        {
            for (c = 0; c < nc; ++c)
            {
                det_img[idx++] = tmp_img(r, c).red;
                det_img[idx++] = tmp_img(r, c).green;
                det_img[idx++] = tmp_img(r, c).blue;
            }
        }
    }
    catch (std::exception e)
    {
        std::cout << "Error in run_net function:" << std::endl;
        std::cout << e.what() << std::endl;
    }
}   // end of run_net

/*
//----------------------------------------------------------------------------------
void get_cropped_detections(unsigned char* input_img, 
    unsigned int nr, 
    unsigned int nc, 
    unsigned int x,
    unsigned int y,
    unsigned int w,
    unsigned int h,
    unsigned int* num_dets, 
    struct detection_struct*& dets
)
{
    uint64_t idx = 0;
    std::string label;

    dlib::matrix<dlib::rgb_pixel> img(nr, nc);
    std::array<dlib::matrix<uint8_t>, array_depth> a_img;

    uint64_t size = (uint64_t)nr * (uint64_t)nc;

    dlib::rectangle rect(x, y, w + x - 1, h + y - 1);

    try
    {
        // copy the pointer into the input image.  The format is assumed to be row-major order
        // and if the array_depth (number of channels) is greater than 1 then channels are not interleaved
        for (idx = 0; idx < array_depth; ++idx)
        {
            dlib::matrix<unsigned char> tmp = dlib::mat<unsigned char>((input_img + (idx * (size))), (long)nr, (long)nc);
            a_img[idx] = dlib::subm(tmp, rect);
        }

        std::vector<dlib::mmod_rect> d = net(a_img);

        prune_detects(d, 0.3);

        *num_dets = d.size();
        dets = new detection_struct[d.size()];

        for (idx = 0; idx < d.size(); ++idx)
        {
            // move the rect back to the original image reference frame
            d[idx].rect = dlib::translate_rect(d[idx].rect, x, y);
            // get the center of the rect
            dlib::point c = dlib::center(d[idx].rect);
            label = d[idx].label.substr(0, std::min((size_t)255, label.length()));

            dets[idx] = detection_struct(d[idx].rect.left(), d[idx].rect.top(), d[idx].rect.width(), d[idx].rect.height(), label.c_str());
        }
    }
    catch (std::exception e)
    {
        std::cout << "Error in get_detections function:" << std::endl;
        std::cout << e.what() << std::endl;
    }

}   // end of get_cropped_detections
*/


//----------------------------------------------------------------------------------
void close_lib()
{
    std::cout << "Closing..." << std::endl;
    net.clean();

    //class_names.clear();
    //class_color.clear();

}   // end of close_lib

//----------------------------------------------------------------------------------
void get_layer_01(struct layer_struct *data, const float* &data_params)
{
    auto& lo = dlib::layer<1>(net).get_output();
    data->k = lo.k();
    data->n = lo.num_samples();
    data->nr = lo.nr();
    data->nc = lo.nc();
    data->size = lo.size();
    data_params = lo.host();
}

//----------------------------------------------------------------------------------
void get_input_layer(layer_struct *data, const float* &data_params)
{
    auto& lo = dlib::layer<net_type::num_layers - 2>(net).get_output();
    data->k = lo.k();
    data->n = lo.num_samples();
    data->nr = lo.nr();
    data->nc = lo.nc();
    data->size = lo.size();
    data_params = lo.host();
}

//----------------------------------------------------------------------------------
//void get_layer_05(layer_struct *data, const float **data_params)
//{
//    auto& lo = dlib::layer<5>(net).get_output();
//    data->k = lo.k();
//    data->n = lo.num_samples();
//    data->nr = lo.nr();
//    data->nc = lo.nc();
//    data->size = lo.size();
//    *data_params = lo.host();
//}
//
////----------------------------------------------------------------------------------
//void get_layer_08(layer_struct *data, const float **data_params)
//{
//    auto& lo = dlib::layer<8>(net).get_output();
//    data->k = lo.k();
//    data->n = lo.num_samples();
//    data->nr = lo.nr();
//    data->nc = lo.nc();
//    data->size = lo.size();
//    *data_params = lo.host();
//}
//
////----------------------------------------------------------------------------------
//void get_layer_09(layer_struct *data, const float **data_params)
//{
//    auto& lo = dlib::layer<9>(net).get_output();
//    data->k = lo.k();
//    data->n = lo.num_samples();
//    data->nr = lo.nr();
//    data->nc = lo.nc();
//    data->size = lo.size();
//    *data_params = lo.host();
//}
//
////----------------------------------------------------------------------------------
//void get_layer_12(layer_struct *data, const float **data_params)
//{
//    auto& lo = dlib::layer<12>(net).get_output();
//    data->k = lo.k();
//    data->n = lo.num_samples();
//    data->nr = lo.nr();
//    data->nc = lo.nc();
//    data->size = lo.size();
//    *data_params = lo.host();
//}

//----------------------------------------------------------------------------------
// check to see if we are building the library or a standalone executable
#if !defined(BUILD_LIB)

int main(int argc, char** argv)
{
    uint32_t idx;
    std::string program_root;
    std::string net_directory;
    std::string image_directory;
    std::string test_net_name;

    // timing variables
    typedef std::chrono::duration<double> d_sec;
    auto start_time = std::chrono::system_clock::now();
    auto stop_time = std::chrono::system_clock::now();
    auto elapsed_time = std::chrono::duration_cast<d_sec>(stop_time - start_time);

    //std::vector<std::string> test_images = { "test1.png", "test2.png", "test3.png", "test4.png", "test5.png", "test6.png", "test7.png", "test8.png", "test9.png", "test10.png" };
    std::vector<std::string> test_images = { "mframe_00156.png", "mframe_00163.png", "mframe_00279.png", "mframe_00353.png", "mframe_05042.png" };

    unsigned int num_classes, num_win;

    unsigned char* tiled_img = NULL;
    unsigned char* det_img = NULL;

    unsigned int t_nr = 0, t_nc = 0;
    window_struct* det;
    unsigned int num_dets = 0;
    struct detection_struct* dets;
    detection_center* detects;
    cv::Mat img;
    long nr, nc;

    // setup save variable locations
#if defined(_WIN32) | defined(__WIN32__) | defined(__WIN32) | defined(_WIN64) | defined(__WIN64)
    program_root = path_check(get_path(get_path(get_path(std::string(argv[0]), "\\"), "\\"), "\\"));
#else 
    program_root = get_ubuntu_path();
#endif

    net_directory = program_root + "../common/nets/";
    image_directory = "../images/";

    try
    {
        test_net_name = (net_directory + "tfd_v03_20_20_100_HPC_final_net.dat");

        // initialize the network
        init_net(test_net_name.c_str(), &num_classes, det, &num_win);

        dlib::matrix<uint8_t> ti;

        // run through some images to test the code
        for (idx = 0; idx < test_images.size(); ++idx)
        {
            img = cv::imread(image_directory + test_images[idx], cv::IMREAD_GRAYSCALE);
            nr = img.rows;
            nc = img.cols;

            unsigned char* image = new unsigned char[nr * nc]{ 0 };

            start_time = std::chrono::system_clock::now();

            run_net(img.ptr<unsigned char>(0), nr, nc, det_img, &num_dets, dets);

            get_detections(img.ptr<unsigned char>(0), nr, nc, &num_dets, detects);

            get_cropped_detections(img.ptr<unsigned char>(0), nr, nc, 128, 0, 256, 256, &num_dets, detects);

            layer_struct ls_all, ls_in;
            const float* ld_all, *ld_in;

            // get the input after the pyramid
            get_input_layer(&ls_in, ld_in);

            // get the combined detection map
            get_combined_output(&ls_all, ld_all);

            stop_time = std::chrono::system_clock::now();
            elapsed_time = std::chrono::duration_cast<d_sec>(stop_time - start_time);

            std::cout << "Runtime (s): " << elapsed_time.count() << std::endl;
        }

        close_lib();

    }
    catch (std::exception & e)
    {
        std::cout << std::endl;
        std::cout << e.what() << std::endl;
    }

    std::cout << std::endl << "Program complete.  Press Enter to close." << std::endl;
    std::cin.ignore();
    return 0;

}   // end of main

#endif  // BUILD_LIB
