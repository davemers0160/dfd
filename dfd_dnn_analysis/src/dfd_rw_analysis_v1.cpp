#define _CRT_SECURE_NO_WARNINGS

#if defined(_WIN32) | defined(__WIN32__) | defined(__WIN32) | defined(_WIN64) | defined(__WIN64)
#include <windows.h>
#endif

#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <thread>
#include <sstream>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <string>
#include <utility>

// Custom includes
#include "get_platform.h"
#include "file_ops.h"
#include "file_parser.h"
#include "get_current_time.h"
#include "num2string.h"
#include "dlib_jet_functions.h"
//#include "center_cropper.h"
//#include "gorgon_common.h"
#include "array_image_operations.h"
#include "image_noise_functions.h"

// Net Version
// Things must go in this order since the array size is determined
// by the network header file
#include "dfd_net_v14.h"
#include "dfd_dnn_analysis.h"
#include "load_dfd_data.h"
#include "eval_dfd_net_performance.h"

#include <vs_gen_lib.h>

// dlib includes
#include <dlib/dnn.h>
#include <dlib/image_io.h>
#include <dlib/data_io.h>
#include <dlib/image_transforms.h>

// this is for enabling GUI support, i.e. opening windows
#ifndef DLIB_NO_GUI_SUPPORT
    #include <dlib/gui_widgets.h>
#endif

using namespace std;
using namespace dlib;

// -------------------------------GLOBALS--------------------------------------

extern const uint32_t img_depth;
extern const uint32_t secondary;
std::string platform;

std::string logfileName = "dfd_net_analysis_results_";

// ----------------------------------------------------------------------------
void get_platform_control(void)
{
    get_platform(platform);
}

//-----------------------------------------------------------------------------
void print_usage(void)
{
    std::cout << "Enter the following as arguments into the program:" << std::endl;
    std::cout << "<config file>" << std::endl;
    std::cout << endl;
}

//-----------------------------------------------------------------------------
inline void vect2matrix(
    uint32_t img_h,
    uint32_t img_w,
    std::vector<uint8_t>& fp1_ptr,
    std::vector<uint8_t>& fp2_ptr,
    std::vector<uint8_t>& dm_ptr,
    std::array<dlib::matrix<uint16_t>, img_depth> &t,
    dlib::matrix<uint16_t> &gt
    )
{
    int index = 0, dm_index = 0;
    for (long r = 0; r < img_h; ++r)
    {
        for (long c = 0; c < img_w; ++c)
        {

            dlib::assign_pixel(t[2](r, c), (uint16_t)fp1_ptr[index]);
            dlib::assign_pixel(t[5](r, c), (uint16_t)fp2_ptr[index++]);
            dlib::assign_pixel(t[1](r, c), (uint16_t)fp1_ptr[index]);
            dlib::assign_pixel(t[4](r, c), (uint16_t)fp2_ptr[index++]);
            dlib::assign_pixel(t[0](r, c), (uint16_t)fp1_ptr[index]);
            dlib::assign_pixel(t[3](r, c), (uint16_t)fp2_ptr[index++]);

            dlib::assign_pixel(gt(r, c), (uint16_t)dm_ptr[dm_index++]);
        }
    }
}

//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
    int idx, jdx;

    std::string sdate, stime;

    std::ofstream data_log_stream;
    std::ofstream dm_results_stream;
    std::ofstream cm_results_stream;

    //std::ofstream trainLogStream;
    //std::string train_inputfile;
    std::string test_inputfile;
    std::string net_name;
    std::string parseFilename;
    std::string results_name;
    std::string data_directory;
    //std::string data_home;
    std::string image_num;

    std::vector<std::vector<std::string>> test_file;
    std::vector<std::pair<std::string, std::string>> image_files;
    
    typedef std::chrono::duration<double> d_sec;
    auto start_time = chrono::system_clock::now();
    auto stop_time = chrono::system_clock::now();
    auto elapsed_time = chrono::duration_cast<d_sec>(stop_time - start_time);

    std::vector<std::array<dlib::matrix<uint16_t>, img_depth>> te;
    std::vector<dlib::matrix<uint16_t>> gt_test;

    std::pair<uint64_t, uint64_t > crop_size(32, 32);
    std::pair<uint32_t, uint32_t> scale(1, 1);
    uint16_t gt_max = 0;

   
    //////////////////////////////////////////////////////////////////////////////////

    if (argc == 1)
    {
        print_usage();
        std::cin.ignore();
        return 0;
    }

    get_platform_control();
    uint8_t HPC = 0;
    
    if(platform.compare(0,3,"HPC") == 0)
    {
        std::cout << "HPC Platform Detected." << std::endl;
        HPC = 1;
    }
    
    
    // setup save variable locations
    const std::string os_file_sep = "/";
    std::string program_root;
    std::string output_save_location;
    
#if defined(_WIN32) | defined(__WIN32__) | defined(__WIN32) | defined(_WIN64) | defined(__WIN64)
    program_root = get_path(get_path(get_path(std::string(argv[0]), "\\"), "\\"), "\\") + os_file_sep;

#else    
    if(HPC == 1)
    {
        //HPC version
        program_root = get_path(get_path(get_path(std::string(argv[0]), os_file_sep), os_file_sep), os_file_sep) + os_file_sep;
    }
    else
    {
        // Ubuntu
        program_root = get_ubuntu_path();
    }

#endif

    std::cout << "Reading Inputs... " << std::endl;
    std::cout << "Platform:               " << platform << std::endl;
    std::cout << "program_root:           " << program_root << std::endl;

   try {
        
        ///////////////////////////////////////////////////////////////////////////////
        // Step 1: Read in the training images
        ///////////////////////////////////////////////////////////////////////////////
        //data_home = path_check(get_env_variable("DATA_HOME"));

        parseFilename = argv[1];
        
        // parse through the supplied input file
        parse_dfd_analysis_file(parseFilename, test_inputfile, net_name, results_name, output_save_location, crop_size, scale);
        
        if (test_inputfile == "" | net_name == "")
        {
            std::cout << "------------------------------------------------------------------" << std::endl;
            std::cout << "Error parsing input file.  No test input data file or network file provided." << std::endl;
            std::cout << "test_inputfile: " << test_inputfile << std::endl;
            std::cout << "net_name: " << net_name << std::endl;
            std::cout << "results_name: " << results_name << std::endl;
            std::cout << "Press Enter to continue..." << std::endl;

            std::cin.ignore();
            return 0;
        }

        output_save_location = path_check(output_save_location);
        std::cout << "output_save_location:   " << output_save_location << std::endl;
/*
        // load the test data
#if defined(_WIN32) | defined(__WIN32__) | defined(__WIN32) | defined(_WIN64) | defined(__WIN64)
        parse_csv_file(test_inputfile, test_file);
        //data_directory = data_home + test_file[0][0];
        data_directory = test_file[0][0];
#else
        if (HPC == 1)
        {
            parse_csv_file(test_inputfile, test_file);
            data_directory = test_file[0][2];
        }
        else
        {
            parse_csv_file(test_inputfile, test_file);
            data_directory = test_file[0][1];
        }
#endif

        test_file.erase(test_file.begin());
        std::cout << "data_directory:         " << data_directory << std::endl;
*/

        // make sure that the save location is there and create if not
        mkdir(output_save_location);

        uint32_t num_test_images = 100;
        //-----------------------------------------------------------------------------
        // read in the blur params
        init_from_file(test_inputfile.c_str());
        //-----------------------------------------------------------------------------

        get_current_time(sdate, stime);
        logfileName = logfileName + results_name + "_" + sdate + "_" + stime + ".txt";

        std::cout << "Log File:               " << (output_save_location + logfileName) << std::endl;
        data_log_stream.open((output_save_location + logfileName), ios::out | ios::app);
        dm_results_stream.open((output_save_location + results_name + "_depth_map_result_images.txt"), ios::out);
        cm_results_stream.open((output_save_location + results_name + "_confusion_matrix_results.txt"), ios::out);

        std::cout << "Data Input File:        " << test_inputfile << std::endl << std::endl;

        // Add the date and time to the start of the log file
        data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        data_log_stream << "Version: 2.6    Date: " << sdate << "    Time: " << stime << std::endl << std::endl;

        data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        data_log_stream << "Platform:             " << platform << std::endl;
        data_log_stream << "program_root:         " << program_root << std::endl;
        data_log_stream << "output_save_location: " << output_save_location << std::endl;
        data_log_stream << "data_directory:       " << data_directory << std::endl;
        data_log_stream << "Data Input File:      " << test_inputfile << std::endl << std::endl;

        //std::cout << "Loading test images..." << std::endl;
        //std::cout << "Test image sets to parse: " << num_test_images << std::endl;

        //data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        //data_log_stream << "Test image sets to parse: " << num_test_images << std::endl << std::endl;

        std::cout << "Input Array Depth: " << img_depth << std::endl;
        std::cout << "Secondary data loading value: " << secondary << std::endl << std::endl;

        data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        data_log_stream << "Input Array Depth: " << img_depth << std::endl;
        data_log_stream << "Secondary data loading value: " << secondary << std::endl << std::endl;
        
        // save the eval crop size
        std::cout << "Eval Crop Size: " << crop_size.first << "x" << crop_size.second << std::endl << std::endl;
        data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        data_log_stream << "Eval Crop Size: " << crop_size.first << "x" << crop_size.second << std::endl << std::endl;
        
        ///////////////////////////////////////////////////////////////////////////////
        // Step 2: Load the network
        ///////////////////////////////////////////////////////////////////////////////
        set_dnn_prefer_smallest_algorithms();

        dfd_net_type dfd_net;
        
        std::cout << "Loading " << net_name << std::endl;
        deserialize(net_name) >> dfd_net;

        std::cout << dfd_net << std::endl;

        data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        data_log_stream << "Net Name: " << net_name << std::endl;
        data_log_stream << dfd_net << std::endl;
        //data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        
        ///////////////////////////////////////////////////////////////////////////////
        // Step 3: Analyze the results of the network
        ///////////////////////////////////////////////////////////////////////////////
  
        std::cout << "Ready to analyze the network performance..." << std::endl;
        
#ifndef DLIB_NO_GUI_SUPPORT
        dlib::image_window win0;
        dlib::image_window win1;
#endif
        dlib::matrix<dlib::rgb_pixel> dm_montage, img_montage;
        dlib::matrix<dlib::rgb_pixel> rgb_img1, rgb_img2;

        double nmae_accum = 0.0;
        double nrmse_accum = 0.0;
        double ssim_accum = 0.0;
        double var_gt_accum = 0.0;
        double var_dm_accum = 0.0;
        double silog_accum = 0.0;

        uint64_t count = 0;

        dlib::matrix<uint16_t> map;
        dlib::matrix<double,1,6> results = dlib::zeros_matrix<double>(1,6);

        //-----------------------------------------------------------------------------
        // setup everything for the first run
        // set the size of the first vector
        std::array<dlib::matrix<uint16_t>, img_depth> tmp;
        dlib::matrix<uint16_t> gt_tmp = dlib::zeros_matrix<uint16_t>(crop_size.first, crop_size.second);
        for (int m = 0; m < img_depth; ++m)
        {
            tmp[m].set_size(crop_size.first, crop_size.second);
        }

        // create a temporary container
        std::vector<uint8_t> fp1_ptr(crop_size.first * crop_size.second * 3);
        std::vector<uint8_t> fp2_ptr(crop_size.first * crop_size.second * 3);
        std::vector<uint8_t> dm_ptr(crop_size.first * crop_size.second);

        double vs_scale = 0.1;

        // generate an image 
        generate_scene(vs_scale, crop_size.second, crop_size.first, fp1_ptr.data(), fp2_ptr.data(), dm_ptr.data());

        // convert the vector pointers to dlib::matrix
        vect2matrix(crop_size.first, crop_size.second, fp1_ptr, fp2_ptr, dm_ptr, tmp, gt_tmp);

        // run through the network once.  This primes the GPU and stabilizes the timing
        // don't need the results.
        eval_net_performance(dfd_net, tmp, gt_tmp, map, crop_size, scale);
        dlib::rand rnd(time(NULL));

        float dm_scale = 1.0;
        // assume that the maximum depthmap value is the number of filters of the last layer - 1
        gt_max = (uint16_t)((dlib::layer<1>(dfd_net).layer_details().num_filters() - 1)*dm_scale);

        dlib::matrix<double> cm = dlib::zeros_matrix<double>(gt_max+1, gt_max+1);

        for (idx = 0; idx < num_test_images; ++idx)
        {
            // generate an image 
            generate_scene(vs_scale, crop_size.second, crop_size.first, fp1_ptr.data(), fp2_ptr.data(), dm_ptr.data());

            // convert the vector pointers to dlib::matrix
            vect2matrix(crop_size.first, crop_size.second, fp1_ptr, fp2_ptr, dm_ptr, tmp, gt_tmp);
            
            // add noise
            apply_poisson_noise(tmp, 8.0, rnd, 0.0, 255.0);

            // change lighting intensity
            //for (jdx = 0; jdx < img_depth; ++jdx)
            //{
            //    tmp[jdx] = dlib::matrix_cast<uint16_t>(dlib::matrix_cast<double>(tmp[jdx]) * 0.95);
            //}

            // time and analyze the results
            start_time = chrono::system_clock::now(); 
            results = eval_net_performance(dfd_net, tmp, gt_tmp, map, crop_size, scale, dm_scale);
            stop_time = chrono::system_clock::now();

            elapsed_time = chrono::duration_cast<d_sec>(stop_time - start_time);

            // fill in the confusion matrix for the range of depthmap values
            for (uint32_t r = 0; r < map.nr(); ++r)
            {
                for (uint32_t c = 0; c < map.nc(); ++c)
                {
                    cm(gt_tmp(r,c), map(r,c)) += 1.0;
                }
            }

            // create a depthmap version in RGB 
            dlib::matrix<dlib::rgb_pixel> dm_img = mat_to_rgbjetmat(dlib::matrix_cast<float>(map)*dm_scale, 0.0, (float)gt_max);
            dlib::matrix<dlib::rgb_pixel> gt_img = mat_to_rgbjetmat(dlib::matrix_cast<float>(gt_tmp)*dm_scale, 0.0, (float)gt_max);

            image_num = num2str(idx, "%05d");

#ifndef DLIB_NO_GUI_SUPPORT

            merge_channels(tmp, rgb_img1, 0);
            merge_channels(tmp, rgb_img2, 3);

            img_montage.set_size(rgb_img1.nr(), rgb_img1.nc() * 2);
            dlib::set_subm(img_montage, 0, 0, rgb_img1.nr(), rgb_img1.nc()) = rgb_img1;
            dlib::set_subm(img_montage, 0, rgb_img1.nc(), rgb_img1.nr(), rgb_img1.nc()) = rgb_img2;
            win0.clear_overlay();
            win0.set_image(img_montage);
            win0.set_title("Focus 1 & Focus 2 Images");

            dm_montage.set_size(gt_img.nr(), gt_img.nc() * 2);
            dlib::set_subm(dm_montage, 0, 0, gt_img.nr(), gt_img.nc()) = gt_img;
            dlib::set_subm(dm_montage, 0, gt_img.nc(), gt_img.nr(), gt_img.nc()) = dm_img;
            win1.clear_overlay();
            win1.set_image(dm_montage);
            win1.set_title("Image " + image_num + ": Groundtruth & DFD Depthmaps");

            //std::cin.ignore();
            dlib::sleep(100);
#endif

            std::string dm_filename = output_save_location + "depthmap_image_" + results_name + "_" + image_num + ".png";
            std::string dm_jet_filename = output_save_location + "depthmap_jet_" + results_name + "_" + image_num + ".png";
            std::string gt_filename = output_save_location + "gt_image_" + results_name + "_" + image_num + ".png";

            std::cout << "------------------------------------------------------------------" << std::endl;
            std::cout << "Depthmap generation completed in: " << elapsed_time.count() << " seconds." << std::endl;
            std::cout << "Image Size (h x w): " << map.nr() << " x " << map.nc() << std::endl;
            //std::cout << "Focus File:     " << image_files[idx].first << std::endl;
            //std::cout << "Defocus File:   " << image_files[idx].second << std::endl;
            //std::cout << "Depth Map File: " << dm_filename << std::endl;
            std::cout << "NMAE   " << image_num << ": " << std::fixed << std::setprecision(5) << results(0, 0) << std::endl;
            std::cout << "NRMSE  " << image_num << ": " << std::fixed << std::setprecision(5) << results(0, 1) << std::endl;
            std::cout << "SSIM   " << image_num << ": " << std::fixed << std::setprecision(5) << results(0, 2) << std::endl;
            std::cout << "SILOG  " << image_num << ": " << std::fixed << std::setprecision(5) << results(0, 3) << std::endl;
            std::cout << "Var_GT " << image_num << ": " << std::fixed << std::setprecision(5) << results(0, 4) << std::endl;
            std::cout << "Var_DM " << image_num << ": " << std::fixed << std::setprecision(5) << results(0, 5) << std::endl;

            data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
            data_log_stream << "Depthmap generation completed in: " << elapsed_time.count() << " seconds." << std::endl;
            data_log_stream << "Image Size (h x w): " << map.nr() << " x " << map.nc() << std::endl;
            //data_log_stream << "Focus File:     " << image_files[idx].first << std::endl;
            //data_log_stream << "Defocus File:   " << image_files[idx].second << std::endl;
            //data_log_stream << "Depth Map File: " << dm_filename << std::endl;
            data_log_stream << "NMAE   " << image_num << ": " << std::fixed << std::setprecision(5) << results(0, 0) << std::endl;
            data_log_stream << "NRMSE  " << image_num << ": " << std::fixed << std::setprecision(5) << results(0, 1) << std::endl;
            data_log_stream << "SSIM   " << image_num << ": " << std::fixed << std::setprecision(5) << results(0, 2) << std::endl;
            data_log_stream << "SILOG  " << image_num << ": " << std::fixed << std::setprecision(5) << results(0, 3) << std::endl;
            data_log_stream << "Var_GT " << image_num << ": " << std::fixed << std::setprecision(5) << results(0, 4) << std::endl;
            data_log_stream << "Var_DM " << image_num << ": " << std::fixed << std::setprecision(5) << results(0, 5) << std::endl;

            dm_results_stream << dm_filename << std::endl;

            // add code to save image
            dlib::save_png(dm_img, dm_jet_filename);
            dlib::save_png(dlib::matrix_cast<uint8_t>(map), dm_filename);
            dlib::save_png(gt_img, gt_filename);

            nmae_accum += results(0, 0);
            nrmse_accum += results(0, 1);
            ssim_accum += results(0, 2);
            silog_accum += results(0, 3);
            var_gt_accum += results(0, 4);
            var_dm_accum += results(0, 5);
            ++count;

            //std::cout << "Press Enter to continue..." << std::endl;
            //std::cin.ignore();
            //std::string key;
            //char key;
            //std::getline(cin,key);

            //if(key.compare("q")==0)
            //    break;

        }
        
        // calculate the confusion matrix errors for each of the depthmap values
        dlib::matrix<double> cm_sum = dlib::sum_cols(cm);
        dlib::matrix<double> cm_diag = dlib::diag(cm);
        dlib::matrix<double> cm_error(1,gt_max + 1);

        for (idx = 0; idx < gt_max + 1; ++idx)
        {
            if (cm_sum(idx) == 0)
                cm_error(idx) = 0.0;
            else
                cm_error(idx) = 1.0 - cm_diag(idx) / cm_sum(idx);
        }

        std::cout << "------------------------------------------------------------------" << std::endl;
        std::cout << "Average Image Analysis Results:" << std::endl;
        std::cout << "Average NMAE:   " << nmae_accum / (double)num_test_images << std::endl;
        std::cout << "Average NRMSE:  " << nrmse_accum / (double)num_test_images << std::endl;
        std::cout << "Average SSIM:   " << ssim_accum / (double)num_test_images << std::endl;
        std::cout << "Average SILOG:  " << silog_accum / (double)num_test_images << std::endl;
        std::cout << "Average Var_GT: " << var_gt_accum / (double)num_test_images << std::endl;
        std::cout << "Average Var_DM: " << var_dm_accum / (double)num_test_images << std::endl;

        //std::cout << "Average VIPF Val:  " << vipf_accum / (double)count << std::endl;
        std::cout << "------------------------------------------------------------------" << std::endl;
        std::cout << "Average NMAE, NRMSE, SSIM, SILOG,  Var_GT, Var_DM: " << nmae_accum / (double)num_test_images << ", " << nrmse_accum / (double)num_test_images << ", " << ssim_accum / (double)num_test_images
                  << ", " << silog_accum / (double)num_test_images << ", " << var_gt_accum / (double)num_test_images << ", " << var_dm_accum / (double)num_test_images << std::endl;
        std::cout << std::endl;

        std::cout << "------------------------------------------------------------------" << std::endl;
        std::cout << std::fixed << std::setprecision(0) << dlib::csv << cm << std::endl;

        std::cout << "------------------------------------------------------------------" << std::endl;
        std::cout << std::fixed << std::setprecision(6) << dlib::csv << cm_error << std::endl;

        data_log_stream << std::endl << "#------------------------------------------------------------------------------" << std::endl;
        data_log_stream << "Average Image Analysis Results:" << std::endl;
        data_log_stream << "Average NMAE:   " << nmae_accum / (double)num_test_images << std::endl;
        data_log_stream << "Average NRMSE:  " << nrmse_accum / (double)num_test_images << std::endl;
        data_log_stream << "Average SSIM:   " << ssim_accum / (double)num_test_images << std::endl;
        data_log_stream << "Average SILOG:  " << silog_accum / (double)num_test_images << std::endl;
        data_log_stream << "Average Var_GT: " << var_gt_accum / (double)num_test_images << std::endl;
        data_log_stream << "Average Var_DM: " << var_dm_accum / (double)num_test_images << std::endl;

        // just save everything for easy copying
        data_log_stream << std::endl << "#------------------------------------------------------------------------------" << std::endl;
        data_log_stream << "Average NMAE, NRMSE, SSIM, Var_GT, Var_DM: " << nmae_accum / (double)num_test_images << ", " << nrmse_accum / (double)num_test_images << ", " << ssim_accum / (double)num_test_images
                      << ", " << silog_accum / (double)num_test_images << ", " << var_gt_accum / (double)num_test_images << ", " << var_dm_accum / (double)num_test_images << std::endl;
        data_log_stream << std::endl;

        data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        data_log_stream << "# Confussion Matrix:" << std::endl;
        data_log_stream << std::fixed << std::setprecision(0) << dlib::csv << cm << std::endl;
        cm_results_stream << std::fixed << std::setprecision(0) << dlib::csv << cm;

        data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        data_log_stream << "# Depthmap Error:" << std::endl;
        data_log_stream << std::fixed << std::setprecision(6) << dlib::csv << cm_error << std::endl;

        std::cout << "End of Program.  Press Enter to close!" << std::endl;
        data_log_stream.close();
        dm_results_stream.close();
        cm_results_stream.close();

        std::cin.ignore();

#ifndef DLIB_NO_GUI_SUPPORT
        if(!win0.is_closed())
            win0.close_window();

        if(!win1.is_closed())
            win1.close_window();
#endif

    }
    catch (std::exception& e)
    {
        std::cout << e.what() << std::endl;

        data_log_stream << e.what() << std::endl;
        data_log_stream << "#------------------------------------------------------------------------------" << std::endl;

        data_log_stream.close();
        dm_results_stream.close();

        std::cout << "Press Enter to close..." << std::endl;
        std::cin.ignore();

    }
    return 0;

}    // end of main
       
        
        

