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

// Net Version
// Things must go in this order since the array size is determined
// by the network header file

#include "dfd_dnn_analysis.h"

#include "dfd_dnn_lib.h"

// OpenCV Includes
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>


// -------------------------------GLOBALS--------------------------------------
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
    std::ofstream cm_results_stream;

    std::string test_inputfile;
    std::string net_name;
    std::string parse_filename;
    std::string results_name;
    std::string data_directory;
    std::string image_num;

    std::string input_image_montage = "DfD Input Images";
    std::string depthmap_images = "Depthmap";

    cv::Mat montage;
    cv::Mat img_f1, img_f2, dm;

    std::vector<std::vector<std::string>> test_file;
    std::vector<std::pair<std::string, std::string>> image_files;
    
    typedef std::chrono::duration<double> d_sec;
    auto start_time = chrono::system_clock::now();
    auto stop_time = chrono::system_clock::now();
    auto elapsed_time = chrono::duration_cast<d_sec>(stop_time - start_time);

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

       parse_filename = argv[1];
        
        // parse through the supplied input file
        parse_dfd_analysis_file(parse_filename, test_inputfile, net_name, results_name, output_save_location, crop_size, scale);
        
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

        // make sure that the save location is there and create if not
        mkdir(output_save_location);
        
        std::cout << "data_directory:         " << data_directory << std::endl;

        get_current_time(sdate, stime);
        logfileName = logfileName + results_name + "_" + sdate + "_" + stime + ".txt";

        std::cout << "Log File:               " << (output_save_location + logfileName) << std::endl;
        //data_log_stream.open((output_save_location + logfileName), ios::out | ios::app);
        //dm_results_stream.open((output_save_location + results_name + "_depth_map_result_images.txt"), ios::out);
        //cm_results_stream.open((output_save_location + results_name + "_confusion_matrix_results.txt"), ios::out);

        std::cout << "Data Input File:        " << test_inputfile << std::endl << std::endl;

        // Add the date and time to the start of the log file
        // data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        // data_log_stream << "Version: 2.6    Date: " << sdate << "    Time: " << stime << std::endl << std::endl;

        // data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        // data_log_stream << "Platform:             " << platform << std::endl;
        // data_log_stream << "program_root:         " << program_root << std::endl;
        // data_log_stream << "output_save_location: " << output_save_location << std::endl;
        // data_log_stream << "data_directory:       " << data_directory << std::endl;
        // data_log_stream << "Data Input File:      " << test_inputfile << std::endl << std::endl;

        std::cout << "Loading test images..." << std::endl;
        std::cout << "Test image sets to parse: " << test_file.size() << std::endl;

        // data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        // data_log_stream << "Test image sets to parse: " << test_file.size() << std::endl << std::endl;
      
        start_time = chrono::system_clock::now();
        //load_dfd_data(test_file, data_directory, mod_params, te, gt_test, image_files);
        
        stop_time = chrono::system_clock::now();
        elapsed_time = chrono::duration_cast<d_sec>(stop_time - start_time);
        //std::cout << "Loaded " << te.size() << " test image sets in " << elapsed_time.count() / 60 << " minutes." << std::endl << std::endl;

        //std::cout << "Input Array Depth: " << img_depth << std::endl;
        //std::cout << "Secondary data loading value: " << secondary << std::endl << std::endl;

        // data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        // data_log_stream << "Input Array Depth: " << img_depth << std::endl;
        // data_log_stream << "Secondary data loading value: " << secondary << std::endl << std::endl;
        
        // save the eval crop size
        std::cout << "Eval Crop Size: " << crop_size.first << "x" << crop_size.second << std::endl << std::endl;
        // data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        // data_log_stream << "Eval Crop Size: " << crop_size.first << "x" << crop_size.second << std::endl << std::endl;
        
        ///////////////////////////////////////////////////////////////////////////////
        // Step 2: Load the network
        ///////////////////////////////////////////////////////////////////////////////
        
        std::cout << "Loading " << net_name << std::endl;
        unsigned int num_classes = 0;
        init_net(net_name.c_str(), &num_classes);
        print_net();

        //data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        //data_log_stream << "Net Name: " << net_name << std::endl;
        //data_log_stream << dfd_net << std::endl;
        //data_log_stream << "#------------------------------------------------------------------------------" << std::endl;
        
        ///////////////////////////////////////////////////////////////////////////////
        // Step 3: Analyze the results of the network
        ///////////////////////////////////////////////////////////////////////////////
  
        std::cout << "Ready to analyze the network performance..." << std::endl;

        double nmae_accum = 0.0;
        double nrmse_accum = 0.0;
        double ssim_accum = 0.0;
        double var_gt_accum = 0.0;
        double var_dm_accum = 0.0;
        double silog_accum = 0.0;

        uint64_t count = 0;
        
        // run through the network once.  This primes the GPU and stabilizes the timing
        // don't need the results.


        // create the windows to display the results and the heatmap
        cv::namedWindow(input_image_montage, cv::WINDOW_NORMAL | cv::WINDOW_KEEPRATIO);
        cv::namedWindow(depthmap_images, cv::WINDOW_NORMAL | cv::WINDOW_KEEPRATIO);

        cv::Rect cr = cv::Rect(0, 0, crop_size.second, crop_size.first);
        
        for (idx = 0; idx < test_file.size(); ++idx)
        {
            // add noise
            //apply_poisson_noise(te[idx], 2.0, rnd, 0.0, 255.0);

            // change lighting intensity
            //for (jdx = 0; jdx < te[idx].size(); ++jdx)
            //{
            //    te[idx][jdx] = dlib::matrix_cast<uint16_t>(dlib::matrix_cast<double>(te[idx][jdx]) * 0.95);
            //}
            img_f1 = cv::imread(data_directory + test_file[idx][0], cv::IMREAD_COLOR);
            img_f2 = cv::imread(data_directory + test_file[idx][1], cv::IMREAD_COLOR);
            dm = cv::Mat::zeros(crop_size.first, crop_size.second, CV_16UC1);


            // crop image

            img_f1 = img_f1(cr);
            img_f2 - img_f2(cr);

            cv::hconcat(img_f1, img_f2, montage);

            // time and analyze the results
            start_time = chrono::system_clock::now(); 
            run_net(img_f1.ptr<unsigned char>(0), img_f2.ptr<unsigned char>(0), (unsigned int)crop_size.first, (unsigned int)crop_size.second, dm.ptr<unsigned short>(0));
            stop_time = chrono::system_clock::now();

            elapsed_time = chrono::duration_cast<d_sec>(stop_time - start_time);
            std::cout << "Depthmap generation completed in: " << elapsed_time.count() << " seconds." << std::endl;

            cv::imshow(input_image_montage, montage);
            cv::imshow(depthmap_images, dm*10);
            cv::waitKey(10);

            //std::cout << "Press Enter to continue..." << std::endl;
            //std::cin.ignore();
            //std::string key;
            //char key;
            //std::getline(cin,key);

            //if(key.compare("q")==0)
            //    break;

        }   // end of idx loop

        //std::cin.ignore();

    }
    catch (std::exception& e)
    {
        std::cout << e.what() << std::endl;

        //data_log_stream << e.what() << std::endl;
        //data_log_stream << "#------------------------------------------------------------------------------" << std::endl;

        //data_log_stream.close();
        //dm_results_stream.close();

        std::cout << "Press Enter to close..." << std::endl;
        std::cin.ignore();

    }
    return 0;

}    // end of main
       
        
        

