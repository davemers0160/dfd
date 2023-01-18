#ifndef DFD_DNN_H_
#define DFD_DNN_H_

#include "ryml_all.hpp"

#include <cstdint>
#include <cassert>
#include <utility>
#include <vector>

#include "file_parser.h"

extern const uint32_t img_depth;
extern const uint32_t secondary;

//-----------------------------------------------------------------------------
typedef struct training_params {

    training_params() = default;
    training_params(double ilr, double flr, double lrsf, uint32_t step) : intial_learning_rate(ilr), final_learning_rate(flr), learning_rate_shrink_factor(lrsf), steps_wo_progess(step) {}

    double intial_learning_rate;
    double final_learning_rate;
    double learning_rate_shrink_factor;
    uint32_t steps_wo_progess;

} training_params;

//-----------------------------------------------------------------------------
typedef struct crop_info {

    crop_info() = default;
    crop_info(uint64_t n, std::pair<uint64_t, uint64_t> tcs, std::pair<uint64_t, uint64_t> ecs, std::pair<uint32_t, uint32_t> vs_size_, double std) : crop_num(n), train_crop_sizes(tcs), eval_crop_sizes(ecs), vs_size(vs_size_), noise_std(std) {}

    uint64_t crop_num;
    //uint64_t crop_height;
    //uint64_t crop_width;
    std::pair<uint64_t, uint64_t> train_crop_sizes;
    std::pair<uint64_t, uint64_t> eval_crop_sizes;
    std::pair<uint32_t, uint32_t> vs_size;
    double noise_std;

} crop_info;

//-----------------------------------------------------------------------------
typedef struct input_data {
    uint32_t data_type;
    std::string filename;

    input_data() = default;
    input_data(uint32_t dt, std::string fn) : data_type(dt), filename(fn) {}

} input_data;

//-----------------------------------------------------------------------------
void parse_dnn_data_file(std::string param_filename,
    std::string &version, 
    std::vector<double> &stop_criteria, 
    training_params& tp,
    input_data &training_data,
    input_data &testing_data,
    //std::string &training_file,
    //std::string &test_file, 
    crop_info& ci,
    std::array<float, img_depth> &avg_color,
    std::vector<uint32_t> &filter_num
)
{
    /*
    # file format version: 4.0
    # The file is organized in the following manner:
    # version (std::string): version name for named saving of various files
    # stopping criteria: duration (double), training steps (double)
    # training parameters: initial_learning_rate (double), final_learnig_rate (double), shrink_factor (double), steps_wo_progress (uint32_t)
    # training_file (std::string): This file contains a list of images and labels used for training
    # test_file (std::string): This file contains a list of images and labels used for testing
    # crop info: Number of crops (uint64_t), train_crop_size (h,w) (uint64_t,uint64_t), eval_crop_size (h,w) (uint64_t,uint64_t), 
    #            virtual scene patch size (h,w) (uint64_t,uint64_t), noise std (double)
    # average colors per channel: vector of average color values (float)
    # filter_num (uint32_t...): This is the number of filters per layer.  Should be a comma separated list, eg. 10,20,30
    #             if the list does not account for the entire network then the code only uses what is available
    #             and leaves the remaining filter number whatever the default value was.  The order of the filters
    #             goes from outer most to the inner most layer.
    */

    uint32_t idx;
    double duration, steps;
    uint64_t h, w;
    std::vector <float> average_colors;
    std::string fn;
    uint32_t dt;

    try {
        std::ifstream tmp_stream(param_filename);
        std::stringstream buffer;
        buffer << tmp_stream.rdbuf();
        std::string contents = buffer.str();
        tmp_stream.close();

        //std::cout << contents << std::endl;

        ryml::Tree config = ryml::parse_in_arena(ryml::to_csubstr(contents));

        // version
        config["version"] >> version;

        // stopping criteria: duration, training steps
        ryml::NodeRef stop_criteria_node = config["stop_criteria"];
        stop_criteria_node["hours"] >> duration;
        stop_criteria_node["steps"] >> steps;
        stop_criteria.push_back(duration);
        stop_criteria.push_back(steps);

        // training parameters: initial_learning_rate, final_learnig_rate, shrink_factor, steps_wo_progress
        ryml::NodeRef training_parameters = config["training_parameters"];
        training_parameters["initial_learning_rate"] >> tp.intial_learning_rate;
        training_parameters["final_learning_rate"] >> tp.final_learning_rate;
        training_parameters["shrink_factor"] >> tp.learning_rate_shrink_factor;
        training_parameters["steps_wo_progress"] >> tp.steps_wo_progess;

        // Training/Testing input files
        ryml::NodeRef tng_data = config["training_data"];
        tng_data["data_type"] >> dt;
        tng_data["train_file"] >> fn;
        training_data = input_data(dt, fn);

        ryml::NodeRef tst_data = config["test_data"];
        tst_data["data_type"] >> dt;
        tst_data["test_file"] >> fn;
        testing_data = input_data(dt, fn);

        // Crop Info: Number of crops, train_crop_size (h,w), eval_crop_size (h,w), virtual scene patch size (h,w), noise std
        ryml::NodeRef crop_info_node = config["crop_info"];
        crop_info_node["num_crops"] >> ci.crop_num;
        crop_info_node["crop_height"] >> h;
        crop_info_node["crop_width"] >> w;
        ci.train_crop_sizes = std::make_pair(h, w);
        crop_info_node["eval_height"] >> h;
        crop_info_node["eval_width"] >> w;
        ci.eval_crop_sizes = std::make_pair(h, w);
        crop_info_node["vs_patch_height"] >> h;
        crop_info_node["vs_patch_width"] >> w;
        ci.vs_size = std::make_pair(h, w);
        crop_info_node["noise_std"] >> ci.noise_std;

        // average colors per channel
        config["average_colors"] >> average_colors; //avg_color;
        assert(average_colors.size() == avg_color.size());
        for (idx = 0; idx < avg_color.size(); ++idx)
            avg_color[idx] = average_colors[idx];


        // number of filters per layer
        filter_num.clear();
        config["filter_num"] >> filter_num;

    }
    catch (std::exception& e)
    {
        throw std::runtime_error("Error parsing input file: " + std::string(e.what()));
    }

}   // end of parse_dnn_data_file

#endif  // DFD_DNN_H_
