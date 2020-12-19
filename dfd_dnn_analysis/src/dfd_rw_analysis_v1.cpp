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
int main(int argc, char** argv)
{
    int idx, jdx;

    std::string sdate, stime;

    std::ofstream data_log_stream;
    std::ofstream dm_results_stream;
    
    //std::ofstream trainLogStream;
    //std::string train_inputfile;
    std::string test_inputfile;
    std::string net_name;
    std::string parseFilename;
    std::string results_name;
    std::string data_directory;
    //std::string data_home;

    std::vector<std::vector<std::string>> test_file;
    std::vector<std::pair<std::string, std::string>> image_files;
    
    typedef std::chrono::duration<double> d_sec;
    auto start_time = chrono::system_clock::now();
    auto stop_time = chrono::system_clock::now();
    auto elapsed_time = chrono::duration_cast<d_sec>(stop_time - start_time);

    std::vector<std::array<dlib::matrix<uint16_t>, img_depth>> te, te_crop;
    std::vector<dlib::matrix<uint16_t>> gt_train, gt_test, gt_crop;

    std::pair<uint64_t, uint64_t > crop_size(32, 32);
    std::pair<uint32_t, uint32_t> scale(1, 1);
    uint16_t gt_min = 0, gt_max = 0;

    // these are the parameters to load in an image to make sure that it is the correct size
    // for the network.  The first number makes sure that the image is a modulus of the number
    // and the second number is an offest from the modulus.  This is used based on the network
    // structure (downsampling and upsampling tensor sizes).
    std::pair<uint32_t, uint32_t> mod_params(16, 0);  
    
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

        // load the test data
#if defined(_WIN32) | defined(__WIN32__) | defined(__WIN32) | defined(_WIN64) | defined(__WIN64)
        parse_csv_file(test_inputfile, test_file);
        //data_directory = data_home + test_file[0][0];
        data_directory = test_file[0][0];
#else
        if (HPC == 1)
        {
            parse_csv_file(test_inputfile, test_file);
            data_directory = data_home + test_file[0][2];
        }
        else
        {
            parse_csv_file(test_inputfile, test_file);
            data_directory = data_home + test_file[0][1];
        }
#endif

        test_file.erase(test_file.begin());

        // make sure that the save location is there and create if not
        mkdir(output_save_location);
        
        std::cout << "data_directory:         " << data_directory << std::endl;

        get_current_time(sdate, stime);
        logfileName = logfileName + results_name + "_" + sdate + "_" + stime + ".txt";

        std::cout << "Log File:               " << (output_save_location + logfileName) << std::endl;
        data_log_stream.open((output_save_location + logfileName), ios::out | ios::app);
        dm_results_stream.open((output_save_location + "depth_map_result_images.txt"), ios::out);

        std::cout << "Data Input File:        " << test_inputfile << std::endl << std::endl;

        // Add the date and time to the start of the log file
        data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        data_log_stream << "Version: 2.5    Date: " << sdate << "    Time: " << stime << std::endl;
        data_log_stream << "#------------------------------------------------------------------------------" << std::endl;

        data_log_stream << "Platform:             " << platform << std::endl;
        data_log_stream << "program_root:         " << program_root << std::endl;
        data_log_stream << "output_save_location: " << output_save_location << std::endl;
        data_log_stream << "data_directory:       " << data_directory << std::endl;
        data_log_stream << "Data Input File:      " << test_inputfile << std::endl;
        data_log_stream << "#------------------------------------------------------------------------------" << std::endl;

        std::cout << "Test image sets to parse: " << test_file.size() << std::endl;

        data_log_stream << "Test image sets to parse: " << test_file.size() << std::endl;
        data_log_stream << "#------------------------------------------------------------------------------" << std::endl;

        std::cout << "Loading test images..." << std::endl;
        
        start_time = chrono::system_clock::now();
        load_dfd_data(test_file, data_directory, mod_params, te, gt_test, image_files);
        
        stop_time = chrono::system_clock::now();
        elapsed_time = chrono::duration_cast<d_sec>(stop_time - start_time);
        std::cout << "Loaded " << te.size() << " test image sets in " << elapsed_time.count() / 60 << " minutes." << std::endl << std::endl;

        std::cout << "Input Array Depth: " << img_depth << std::endl;
        std::cout << "Secondary data loading value: " << secondary << std::endl << std::endl;
        data_log_stream << "Input Array Depth: " << img_depth << std::endl;
        data_log_stream << "Secondary data loading value: " << secondary << std::endl;
        data_log_stream << "#------------------------------------------------------------------------------" << std::endl << std::endl;
        
        // save the eval crop size
        std::cout << "Eval Crop Size: " << crop_size.first << "x" << crop_size.second << std::endl << std::endl;
        data_log_stream << "Eval Crop Size: " << crop_size.first << "x" << crop_size.second << std::endl;
        data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        
        ///////////////////////////////////////////////////////////////////////////////
        // Step 2: Load the network
        ///////////////////////////////////////////////////////////////////////////////
        set_dnn_prefer_smallest_algorithms();

        dfd_net_type dfd_net;
        
        std::cout << "Loading " << net_name << std::endl;
        deserialize(net_name) >> dfd_net;

        std::cout << dfd_net << std::endl;

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
        //dlib::image_window win2;
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
        
        // run through the network once.  This primes the GPU and stabilizes the timing
        // don't need the results.
        eval_net_performance(dfd_net, te[0], gt_test[0], map, crop_size, scale);
        dlib::rand rnd(time(NULL));

        float dm_scale = 1.0;
        // assume that the maximum depthmap value is the number of filters of the last layer - 1
        gt_max = (uint16_t)((dlib::layer<1>(dfd_net).layer_details().num_filters() - 1)*dm_scale);

        for (idx = 0; idx < te.size(); ++idx)
        {
            // add noise
            //apply_poisson_noise(te[idx], 2.0, rnd, 0.0, 255.0);

            // change lighting intensity
            //for (jdx = 0; jdx < te[idx].size(); ++jdx)
            //{
            //    te[idx][jdx] = dlib::matrix_cast<uint16_t>(dlib::matrix_cast<double>(te[idx][jdx]) * 0.95);
            //}

            // time and analyze the results
            start_time = chrono::system_clock::now(); 
            results = eval_net_performance(dfd_net, te[idx], gt_test[idx], map, crop_size, scale, dm_scale);
            stop_time = chrono::system_clock::now();

            elapsed_time = chrono::duration_cast<d_sec>(stop_time - start_time);

            // create a depthmap version in RGB 
            dlib::matrix<dlib::rgb_pixel> dm_img = mat_to_rgbjetmat(dlib::matrix_cast<float>(map)*dm_scale, 0.0, (float)gt_max);
            dlib::matrix<dlib::rgb_pixel> gt_img = mat_to_rgbjetmat(dlib::matrix_cast<float>(gt_test[idx])*dm_scale, 0.0, (float)gt_max);

#ifndef DLIB_NO_GUI_SUPPORT

            merge_channels(te[idx], rgb_img1, 0);
            merge_channels(te[idx], rgb_img2, 3);

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
            win1.set_title("Groundtruth & DFD DNN Depthmaps");

            //win2.clear_overlay();
            //win2.set_image(dm_img);
            //win2.set_title("DFD DNN Depthmap");

            //std::cin.ignore();
            dlib::sleep(500);
#endif

            std::string dm_filename = output_save_location + "depthmap_image_" + results_name + num2str(idx, "_%05d") + ".png";
            std::string dm_jet_filename = output_save_location + "depthmap_jet_" + results_name + num2str(idx, "_%05d") + ".png";
            std::string gt_filename = output_save_location + "gt_image_" + results_name + num2str(idx, "_%05d") + ".png";

            std::cout << "------------------------------------------------------------------" << std::endl;
            std::cout << "Depthmap generation completed in: " << elapsed_time.count() << " seconds." << std::endl;
            std::cout << "Image Size (h x w): " << map.nr() << " x " << map.nc() << std::endl;
            std::cout << "Focus File:     " << image_files[idx].first << std::endl;
            std::cout << "Defocus File:   " << image_files[idx].second << std::endl;
            std::cout << "Depth Map File: " << dm_filename << std::endl;
            std::cout << "NMAE   " << std::setw(5) << std::setfill('0') << idx << ": " << std::fixed << std::setprecision(5) << results(0, 0) << std::endl;
            std::cout << "NRMSE  " << std::setw(5) << std::setfill('0') << idx << ": " << std::fixed << std::setprecision(5) << results(0, 1) << std::endl;
            std::cout << "SSIM   " << std::setw(5) << std::setfill('0') << idx << ": " << std::fixed << std::setprecision(5) << results(0, 2) << std::endl;
            std::cout << "SILOG  " << std::setw(5) << std::setfill('0') << idx << ": " << std::fixed << std::setprecision(5) << results(0, 3) << std::endl;           
            std::cout << "Var_GT " << std::setw(5) << std::setfill('0') << idx << ": " << std::fixed << std::setprecision(5) << results(0, 4) << std::endl;
            std::cout << "Var_DM " << std::setw(5) << std::setfill('0') << idx << ": " << std::fixed << std::setprecision(5) << results(0, 5) << std::endl;

            data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
            data_log_stream << "Depthmap generation completed in: " << elapsed_time.count() << " seconds." << std::endl;
            data_log_stream << "Image Size (h x w): " << map.nr() << " x " << map.nc() << std::endl;
            data_log_stream << "Focus File:     " << image_files[idx].first << std::endl;
            data_log_stream << "Defocus File:   " << image_files[idx].second << std::endl;
            data_log_stream << "Depth Map File: " << dm_filename << std::endl;
            data_log_stream << "NMAE   " << std::setw(5) << std::setfill('0') << idx << ": " << std::fixed << std::setprecision(5) << results(0, 0) << std::endl;
            data_log_stream << "NRMSE  " << std::setw(5) << std::setfill('0') << idx << ": " << std::fixed << std::setprecision(5) << results(0, 1) << std::endl;
            data_log_stream << "SSIM   " << std::setw(5) << std::setfill('0') << idx << ": " << std::fixed << std::setprecision(5) << results(0, 2) << std::endl;
            data_log_stream << "SILOG  " << std::setw(5) << std::setfill('0') << idx << ": " << std::fixed << std::setprecision(5) << results(0, 3) << std::endl;
            data_log_stream << "Var_GT " << std::setw(5) << std::setfill('0') << idx << ": " << std::fixed << std::setprecision(5) << results(0, 4) << std::endl;
            data_log_stream << "Var_DM " << std::setw(5) << std::setfill('0') << idx << ": " << std::fixed << std::setprecision(5) << results(0, 5) << std::endl;

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
        
        std::cout << "------------------------------------------------------------------" << std::endl;
        std::cout << "Average Image Analysis Results:" << std::endl;
        std::cout << "Average NMAE:   " << nmae_accum / (double)te.size() << std::endl;
        std::cout << "Average NRMSE:  " << nrmse_accum / (double)te.size() << std::endl;
        std::cout << "Average SSIM:   " << ssim_accum / (double)te.size() << std::endl;
        std::cout << "Average SILOG:  " << silog_accum / (double)te.size() << std::endl;
        std::cout << "Average Var_GT: " << var_gt_accum / (double)te.size() << std::endl;
        std::cout << "Average Var_DM: " << var_dm_accum / (double)te.size() << std::endl;

        //std::cout << "Average VIPF Val:  " << vipf_accum / (double)count << std::endl;
        std::cout << "------------------------------------------------------------------" << std::endl;
        std::cout << "Average NMAE, NRMSE, SSIM, SILOG,  Var_GT, Var_DM: " << nmae_accum / (double)te.size() << ", " << nrmse_accum / (double)te.size() << ", " << ssim_accum / (double)te.size()
                  << ", " << silog_accum / (double)te.size() << ", " << var_gt_accum / (double)te.size() << ", " << var_dm_accum / (double)te.size() << std::endl;
        std::cout << std::endl;

        data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        data_log_stream << "Average Image Analysis Results:" << std::endl;
        data_log_stream << "Average NMAE:   " << nmae_accum / (double)te.size() << std::endl;
        data_log_stream << "Average NRMSE:  " << nrmse_accum / (double)te.size() << std::endl;
        data_log_stream << "Average SSIM:   " << ssim_accum / (double)te.size() << std::endl;       
        data_log_stream << "Average SILOG:  " << silog_accum / (double)te.size() << std::endl;
        data_log_stream << "Average Var_GT: " << var_gt_accum / (double)te.size() << std::endl;
        data_log_stream << "Average Var_DM: " << var_dm_accum / (double)te.size() << std::endl;

        // just save everything for easy copying
        data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        data_log_stream << "Average NMAE, NRMSE, SSIM, Var_GT, Var_DM: " << nmae_accum / (double)te.size() << ", " << nrmse_accum / (double)te.size() << ", " << ssim_accum / (double)te.size()
                      << ", " << silog_accum / (double)te.size() << ", " << var_gt_accum / (double)te.size() << ", " << var_dm_accum / (double)te.size() << std::endl;
        data_log_stream << std::endl;


    //#if defined(_WIN32) | defined(__WIN32__) | defined(__WIN32) | defined(_WIN64) | defined(__WIN64)
    //    Beep(500, 1000);
    //#endif

        std::cout << "End of Program.  Press Enter to close!" << std::endl;
        data_log_stream.close();
        dm_results_stream.close();

        std::cin.ignore();

#ifndef DLIB_NO_GUI_SUPPORT

        if(!win0.is_closed())
            win0.close_window();

        if(!win1.is_closed())
            win1.close_window();

        //if(!win2.is_closed())
        //    win2.close_window();
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
       
        
        

