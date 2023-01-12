#ifndef DFD_ANALYSIS_H_
#define DFD_ANALYSIS_H_

#include "ryml_all.hpp"

#include <vector> 
#include <string>
#include <cstdint>
#include <utility>

#include "file_parser.h"

extern const uint32_t img_depth;
extern const uint32_t secondary;

// ----------------------------------------------------------------------------------------

void parse_dfd_analysis_file(std::string param_filename,
    std::string &data_file, 
    uint32_t &data_type,
    std::string &net_file,
    std::string &results_name, 
    std::string &save_location, 
    std::pair<uint64_t, uint64_t> &crop_size 
)
{
    uint32_t idx;
    uint64_t h, w;

    try {
        std::ifstream tmp_stream(param_filename);
        std::stringstream buffer;
        buffer << tmp_stream.rdbuf();
        std::string contents = buffer.str();
        tmp_stream.close();

        //std::cout << contents << std::endl;

        ryml::Tree config = ryml::parse_in_arena(ryml::to_csubstr(contents));

        // version
        config["results_name"] >> results_name;

        // input file
        config["data_file"] >> data_file;

        // data_type
        config["data_type"] >> data_type;

        // network weight file
        config["net_file"] >> net_file;

        // network weight file
        config["save_location"] >> save_location;


        // Crop Info: Number of crops, train_crop_size (h,w), eval_crop_size (h,w), virtual scene patch size (h,w), noise std
        ryml::NodeRef crop_info_node = config["crop_info"];
        crop_info_node["eval_height"] >> h;
        crop_info_node["eval_width"] >> w;
        crop_size = std::make_pair(h, w);
        //crop_info_node["noise_std"] >> ci.noise_std;

    }
    catch (std::exception& e)
    {
        throw std::runtime_error("Error parsing input file: " + std::string(e.what()));
    }


    //for (uint64_t idx = 0; idx<params.size(); ++idx)
    //{
    //    switch (idx)
    //    {

    //    case 0:
    //        data_file = params[idx][0];
    //        break;

    //    case 1:
    //        net_name = params[idx][0];
    //        break;

    //    case 2:
    //        results_name = params[idx][0];
    //        break;

    //    case 3:
    //        save_location = params[idx][0];
    //        break;

    //    // get the crop size of the input data
    //    case 4:
    //        try {
    //            crop_size = std::make_pair(stol(params[idx][0]), stol(params[idx][1]));
    //        }
    //        catch (std::exception &e) {
    //            std::cout << e.what() << std::endl;
    //            crop_size = std::make_pair(34, 140);
    //            std::cout << "Setting Evaluation Crop Size to " << crop_size.first << "x" << crop_size.second << std::endl;
    //        }
    //        break;

    //    // get the scale for the ground truth depth map size
    //    case 5:
    //        try {
    //            scale = std::make_pair(stol(params[idx][0]), stol(params[idx][1]));
    //        }
    //        catch (std::exception &e) {
    //            std::cout << e.what() << std::endl;
    //            scale = std::make_pair(6, 18);
    //            std::cout << "Setting Scale to default: " << scale.first << " x " << scale.second << std::endl;
    //        }
    //        break;
    //    default:
    //        break;
    //    }   // end of switch

    //}   // end of for    


}   // end of parse_dfd_analysis_file

#endif  // DFD_ANALYSIS_H_

