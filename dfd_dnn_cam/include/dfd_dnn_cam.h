#ifndef DFD_DNN_CAM_H_
#define DFD_DNN_CAM_H_

#include <cstdint>
#include <utility>

#include "chameleon_utilities.h"
#include "file_parser.h"

extern const uint32_t img_depth;
extern const uint32_t secondary;

void parse_dnn_cam_file(std::string parse_filename, std::vector<uint8_t> &lens_step, std::vector<uint16_t> &cam_image_params, cam_properties_struct &cam_properties, std::string &net_name)
{

    std::vector<std::vector<std::string>> params;
    parse_csv_file(parse_filename, params);

    for (uint64_t idx = 0; idx<params.size(); ++idx)
    {
        switch (idx)
        {
            // voltage step setting for the lens driver
            case 0:
                try {
                    lens_step.clear();
                    lens_step.push_back(std::stoi(params[idx][0]));
                    lens_step.push_back(std::stoi(params[idx][1]));                    
                }
                catch (std::exception &e) {
                    std::cout << e.what() << std::endl;
                    lens_step.clear();
                    lens_step.push_back(133);
                    lens_step.push_back(130); 
                    std::cout << "Error getting lens step values.  Setting values to default." << std::endl;
                }
                break;
                
            // camera image paramerters (x_offset, y_offset, width, height)
            case 1:
                try {
                    cam_image_params.clear();
                    cam_image_params.push_back(std::stoi(params[idx][0]));
                    cam_image_params.push_back(std::stoi(params[idx][1]));
                    cam_image_params.push_back(std::stoi(params[idx][2]));
                    cam_image_params.push_back(std::stoi(params[idx][3]));

                }
                catch (std::exception &e) {
                    std::cout << e.what() << std::endl;
                    cam_image_params.clear();
                    cam_image_params.push_back(150);   // x offset
                    cam_image_params.push_back(152);   // y offset
                    cam_image_params.push_back(900);   // width
                    cam_image_params.push_back(720);   // height
                    std::cout << "Error getting image capture parameters.  Setting values to default." << std::endl;
                }
                break;
                
            // get the camera capture properties (value, AutoMode, OnOff, absControl)
            case 2:
                try {                   
                    cam_properties.sharpness = cam_prop<uint32_t>(std::stoi(params[idx][0]), false, true, false);
                    cam_properties.fps = cam_prop<float>(std::stof(params[idx][1]), false, true, true);
                    cam_properties.shutter = cam_prop<float>(std::stof(params[idx][2]), false, true, true);
                    cam_properties.gain = cam_prop<float>(std::stof(params[idx][3]), false, true, true);
                }
                catch (std::exception &e) {
                    std::cout << e.what() << std::endl;
                    cam_properties.sharpness = cam_prop<uint32_t>(2500, false, true, false);
                    cam_properties.fps = cam_prop<float>(10.0f, false, true, true);
                    cam_properties.shutter = cam_prop<float>(150.0f, false, true, true);
                    cam_properties.gain = cam_prop<float>(8.0f, false, true, true);
                }
                break;
            // net name
            case 3:
                net_name = params[idx][0];
                break;

            default:
                break;
        }   // end of switch

    }   // end of for

}   // end of parse_dnn_cam_file

#endif  // DFD_DNN_CAM_H_