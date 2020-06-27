#ifndef DFD_DNN_H_
#define DFD_DNN_H_

#include <cstdint>
#include <utility>

#include "file_parser.h"

extern const uint32_t img_depth;
extern const uint32_t secondary;

// ----------------------------------------------------------------------------------------

typedef struct training_params {

    training_params() = default;
    training_params(double ilr, double flr, double lrsf, uint32_t step) : intial_learning_rate(ilr), final_learning_rate(flr), learning_rate_shrink_factor(lrsf), steps_wo_progess(step) {}

    double intial_learning_rate;
    double final_learning_rate;
    double learning_rate_shrink_factor;
    uint32_t steps_wo_progess;

} training_params;

// ----------------------------------------------------------------------------------------

typedef struct crop_info {

    crop_info() = default;
    crop_info(uint64_t n, std::pair<uint64_t, uint64_t> tcs, std::pair<uint64_t, uint64_t> ecs, std::pair<uint32_t, uint32_t> sc) : crop_num(n), train_crop_sizes(tcs), eval_crop_sizes(ecs), scale(sc) {}

    uint64_t crop_num;
    //uint64_t crop_height;
    //uint64_t crop_width;
    std::pair<uint64_t, uint64_t> train_crop_sizes;
    std::pair<uint64_t, uint64_t> eval_crop_sizes;
    std::pair<uint32_t, uint32_t> scale;
} crop_info;

// ----------------------------------------------------------------------------------------

void parse_dnn_data_file(std::string parseFilename, 
    std::string &version, 
    std::vector<int32_t> &gpu,
    std::vector<double> &stop_criteria, 
    training_params& tp,
    std::string &training_file,
    std::string &test_file, 
    crop_info& ci,
    //uint64_t &num_crops, 
    //std::vector<std::pair<uint64_t, uint64_t>> &crop_sizes, 
    //std::pair<uint32_t, uint32_t> &scale, 
    std::array<float, img_depth> &avg_color,
    std::vector<uint32_t> &filter_num
)
{
    /*
    # Version 3.1
    # The file is organized in the following manner:
    # Version (std::string): version name for named svaing of various files
    # Stopping Criteria (double, double) [stop time (hrs), max one step]
    # training_file (std::string): This file contains a list of images and labels used for training
    # test_file (std::string): This file contains a list of images and labels used for testing
    # crop_num (uint64_t): The number of crops to use when using a random cropper
    # crop_size (uint64_t, uint64_t): This is the height and width of the crop size
    # filter_num (uint64_t...): This is the number of filters per layer.  Should be a comma separated list, eg. 10,20,30
    #             if the list does not account for the entire network then the code only uses what is available
    #             and leaves the remaining filter number whatever the default value was.  The order of the filters
    #             goes from outer most to the inner most layer.
    */

    std::vector<std::vector<std::string>> params;
    parse_csv_file(parseFilename, params);

    for (uint64_t idx = 0; idx<params.size(); ++idx)
    {
        switch (idx)
        {

            // get the version
            case 0:
                version = params[idx][0];
                break;

            // get the gpu(s) to target
            case 1:
                try {
                    gpu.clear();
                    for (uint64_t jdx = 0; jdx < params[idx].size(); ++jdx)
                    {
                        gpu.push_back(stoi(params[idx][jdx]));
                    }
                }
                catch (std::exception &e)
                {
                    gpu.clear();
                    gpu.push_back(0);
                    std::cout << e.what() << std::endl;
                    std::cout << "Error getting the GPU(s) to target.  Setting the targeted GPU to: 0" << std::endl;
                }
                break;

            // get the stopping criteria
            case 2:
                try {
                    stop_criteria.clear();
                    for (uint64_t jdx = 0; jdx<params[idx].size(); ++jdx)
                    {
                        stop_criteria.push_back(stod(params[idx][jdx]));
                    }
                }
                catch (std::exception &e) {
                    stop_criteria.clear();
                    stop_criteria.push_back(160.0);
                    stop_criteria.push_back(250000.0);
                    std::cout << e.what() << std::endl;
                    std::cout << "Error getting stopping criteria.  Setting values to default." << std::endl;
                }
                break;

            // get the training parameters
            case 3:
                try {
                    tp = training_params(stod(params[idx][0]), stod(params[idx][1]), stod(params[idx][2]), stol(params[idx][3]));
                }
                catch (std::exception & e) {
                    std::cout << e.what() << std::endl;
                    std::cout << "Using default training parameters..." << std::endl;
                    tp = training_params(0.001, 0.000001, 0.1, 2500);
                }
                break;

            // get the training input file
            case 4:
                training_file = params[idx][0];
                break;

            // get the test input file
            case 5:
                test_file = params[idx][0];
                break;

            // get the crop info
            case 6:
                try {
                    ci = crop_info(stol(params[idx][0]), 
                        std::make_pair(stol(params[idx][1]), stol(params[idx][2])), 
                        std::make_pair(stol(params[idx][3]), stol(params[idx][4])),
                        std::make_pair(stol(params[idx][5]), stol(params[idx][6])));
                }
                catch (std::exception & e) {
                    std::cout << e.what() << std::endl;
                    std::cout << "Setting crop-info to defalut values..." << std::endl;
                    ci = crop_info(40, std::make_pair(32, 32), std::make_pair(352, 352), std::make_pair(1,1));
                }
                break;

            //case 7:
            //    try {
            //        crop_sizes.push_back(std::make_pair(stol(params[idx][0]), stol(params[idx][1])));
            //        crop_sizes.push_back(std::make_pair(stol(params[idx][2]), stol(params[idx][3])));
            //    }
            //    catch (std::exception &e) {
            //        std::cout << e.what() << std::endl;
            //        crop_sizes.push_back(std::make_pair(108, 36));
            //        crop_sizes.push_back(std::make_pair(108, 36));
            //        std::cout << "Setting Training Crop Size to " << crop_sizes[0].first << "x" << crop_sizes[0].second << std::endl;
            //        std::cout << "Setting Evaluation Crop Size to " << crop_sizes[1].first << "x" << crop_sizes[1].second << std::endl;
            //    }
            //    break;

            // 
            //case 8:
            //    try {
            //        scale = std::make_pair(stol(params[idx][0]), stol(params[idx][1]));
            //    }
            //    catch (std::exception &e) {
            //        std::cout << e.what() << std::endl;
            //        scale = std::make_pair(6, 18);
            //        std::cout << "Setting Scale to default: " << scale.first << " x " << scale.second << std::endl;
            //    }
            //    break;

            // get the average colors for the dataset
            case 7:
                try {
                    for (int jdx = 0; jdx < img_depth; ++jdx)
                    {
                        avg_color[jdx] = std::stof(params[idx][jdx]);
                    }
                }
                catch (std::exception & e) {
                    avg_color.fill(128);
                    std::cout << e.what() << std::endl;
                    std::cout << "Error getting average color values.  Setting values to 128." << std::endl;
                }
                break;

            // get the number of conv filters for each layer
            case 8:
                try {
                    filter_num.clear();
                    for (int jdx = 0; jdx<params[idx].size(); ++jdx)
                    {
                        filter_num.push_back(stol(params[idx][jdx]));
                    }
                }
                catch (std::exception &e) {
                    filter_num.clear();
                    std::cout << e.what() << std::endl;
                    std::cout << "Error getting filter numbers.  No values passed on." << std::endl;
                }
                break;

            default:
                break;

        }   // end of switch

    }   // end of for

}   // end of parse_dnn_data_file

#endif  // DFD_DNN_H_
