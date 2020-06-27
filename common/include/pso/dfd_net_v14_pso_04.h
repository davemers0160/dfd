// Auto generated c++ header file for PSO testing
// Iteration Number: 60
// Population Number: 4

#ifndef NET_DEFINITION_H
#define NET_DEFINITION_H

#include <cstdint>
#include <string>

#include "dlib/dnn.h"
#include "dlib/dnn/core.h"

#include "dlib_elu.h"
#include "dlib_srelu.h"

extern const uint32_t img_depth = 6;
extern const uint32_t secondary = 1;
extern const std::string train_inputfile = "dfd_train_data_sm2.txt";
extern const std::string test_inputfile = "dfd_test_data_sm2.txt";
extern const std::string version = "pso_04_60_";
extern const std::vector<std::pair<uint64_t, uint64_t>> crop_sizes = {{32, 32}, {368, 400}};
extern const uint64_t num_crops = 28;
extern const std::pair<uint32_t, uint32_t> crop_scale(1, 1);

//-----------------------------------------------------------------

typedef struct{
    uint32_t iteration;
    uint32_t pop_num;
} pso_struct;
pso_struct pso_info = {60, 4};

//-----------------------------------------------------------------

template <long num_filters, typename SUBNET> using con2d = dlib::con<num_filters, 2, 2, 2, 2, SUBNET>;
template <long num_filters, typename SUBNET> using con33d33 = dlib::con<num_filters, 3, 3, 3, 3, SUBNET>;
template <long num_filters, typename SUBNET> using con32d32 = dlib::con<num_filters, 3, 2, 3, 2, SUBNET>;
template <long num_filters, typename SUBNET> using con21d21 = dlib::con<num_filters, 2, 1, 2, 1, SUBNET>;

template <long num_filters, typename SUBNET> using cont2u = dlib::cont<num_filters, 2, 2, 2, 2, SUBNET>;

template <typename SUBNET> using DTO_0 = dlib::add_tag_layer<200, SUBNET>;
template <typename SUBNET> using DTI_0 = dlib::add_tag_layer<201, SUBNET>;
template <typename SUBNET> using DTO_1 = dlib::add_tag_layer<202, SUBNET>;
template <typename SUBNET> using DTI_1 = dlib::add_tag_layer<203, SUBNET>;
template <typename SUBNET> using DTO_2 = dlib::add_tag_layer<204, SUBNET>;
template <typename SUBNET> using DTI_2 = dlib::add_tag_layer<205, SUBNET>;

//-----------------------------------------------------------------

using dfd_net_type = dlib::loss_multiclass_log_per_pixel<
    dlib::con<256, 3, 3, 1, 1, 
    dlib::prelu<    dlib::bn_con<    dlib::add_prev1<    dlib::con<135, 3, 3, 1, 1, 
    dlib::prelu<    dlib::bn_con<    dlib::con<181, 3, 3, 1, 1, 
    dlib::tag1<    dlib::prelu<    dlib::bn_con<    dlib::con<135, 3, 3, 1, 1, 
    dlib::concat2<DTO_1, DTI_1,
    DTI_1<    dlib::cont<138, 2, 2, 2, 2, 
    dlib::add_prev1<    dlib::con<255, 3, 3, 1, 1, 
    dlib::prelu<    dlib::bn_con<    dlib::con<297, 3, 3, 1, 1, 
    dlib::tag1<    dlib::prelu<    dlib::bn_con<    dlib::con<255, 3, 3, 1, 1, 
    dlib::concat2<DTO_2, DTI_2,
    DTI_2<    dlib::cont<251, 2, 2, 2, 2, 
    dlib::add_prev1<    dlib::con<434, 3, 3, 1, 1, 
    dlib::prelu<    dlib::bn_con<    dlib::con<190, 3, 3, 1, 1, 
    dlib::tag1<    con2d<434,    DTO_2<    dlib::add_prev1<    dlib::con<238, 3, 3, 1, 1, 
    dlib::prelu<    dlib::bn_con<    dlib::con<267, 3, 3, 1, 1, 
    dlib::tag1<    con2d<238,    DTO_1<    dlib::add_prev1<    dlib::con<145, 3, 3, 1, 1, 
    dlib::prelu<    dlib::bn_con<    dlib::con<121, 3, 3, 1, 1, 
    dlib::tag1<    dlib::sig<    dlib::bn_con<    dlib::con<145, 3, 3, 1, 1, 
    dlib::input<std::array<dlib::matrix<uint16_t>, img_depth>>
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>;

//-----------------------------------------------------------------

inline std::ostream& operator<< (std::ostream& out, const pso_struct& item)
{
    out << "------------------------------------------------------------------" << std::endl;
    out << "PSO Info: " << std::endl;
    out << "  Iteration: " << item.iteration << std::endl;
    out << "  Population Number: " << item.pop_num << std::endl;
    out << "------------------------------------------------------------------" << std::endl;
    return out;
}

#endif 
