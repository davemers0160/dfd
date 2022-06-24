#ifndef DFD_ANALYSIS_H_
#define DFD_ANALYSIS_H_

#include <vector> 
#include <string>
#include <cstdint>
#include <utility>

#include "file_parser.h"

//extern const uint32_t img_depth;
//extern const uint32_t secondary;

// ----------------------------------------------------------------------------------------
void parse_dfd_analysis_file(std::string parseFilename, 
    std::string &data_file, 
    std::string &net_name, 
    std::string &results_name, 
    std::string &save_location, 
    std::pair<uint64_t, uint64_t> &crop_size, 
    std::pair<uint32_t, uint32_t> &scale)
{

    std::vector<std::vector<std::string>> params;
    parse_csv_file(parseFilename, params);
    data_file = "";
    net_name = "";

    for (uint64_t idx = 0; idx<params.size(); ++idx)
    {
        switch (idx)
        {

        case 0:
            data_file = params[idx][0];
            break;

        case 1:
            net_name = params[idx][0];
            break;

        case 2:
            results_name = params[idx][0];
            break;

        case 3:
            save_location = params[idx][0];
            break;

        // get the crop size of the input data
        case 4:
            try {
                crop_size = std::make_pair(stol(params[idx][0]), stol(params[idx][1]));
            }
            catch (std::exception &e) {
                std::cout << e.what() << std::endl;
                crop_size = std::make_pair(34, 140);
                std::cout << "Setting Evaluation Crop Size to " << crop_size.first << "x" << crop_size.second << std::endl;
            }
            break;

        // get the scale for the ground truth depth map size
        case 5:
            try {
                scale = std::make_pair(stol(params[idx][0]), stol(params[idx][1]));
            }
            catch (std::exception &e) {
                std::cout << e.what() << std::endl;
                scale = std::make_pair(6, 18);
                std::cout << "Setting Scale to default: " << scale.first << " x " << scale.second << std::endl;
            }
            break;
        default:
            break;
        }   // end of switch

    }   // end of for    


}   // end of parse_dfd_analysis_file

// ----------------------------------------------------------------------------------------
//void load_data(
//    const std::vector<std::vector<std::string>> training_file,
//    const std::string data_directory,
//    std::pair<uint32_t, uint32_t> mod_params,
//
//)
//)


#endif  // DFD_ANALYSIS_H_

