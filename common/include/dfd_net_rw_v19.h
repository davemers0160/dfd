#ifndef NET_DEFINITION_H
#define NET_DEFINITION_H

#include <cstdint>

// dlib includes
#include "dlib/dnn.h"
#include "dlib/dnn/core.h"

#include "gorgon_capture.h"

extern const uint32_t img_depth = 6;
extern const uint32_t secondary = 1;
 
/* ----------------------------------------------------------------------------------------
dfd_res_NM => con3(N)-> prelu -> bn -> con3(M)
cbpN_blk => prelu -> bn -> con3(N)

input -> con3 -> dfd_res_33 [DTO_0] -------------------------------------------------------------------------------------------------> [DTI_0] con3 -> dfd_res_33 -> prelu -> con3(256) -> output
                       \                                                                                                                        /
                      con2d                                                                                                                  cont2u
                         \                                                                                                                    /
                        con3 -> dfd_res_33 [DTO_1] ---------------------------------------------------------------------> [DTI_1] con3 -> dfd_res_33
                                      \                                                                                            /
                                     con2d                                                                                      cont2u
                                        \                                                                                        /
                                       con3 -> dfd_res_33 [DTO_2] -----------------------------------------> [DTI_2] con3 -> dfd_res_33
                                                     \                                                                /
                                                    con2d                                                          cont2u
                                                       \                                                            /                                 
                                                      con3 -> dfd_res_33 [DTO_3] -------------> [DTI_3] con3 -> dfd_res_33
                                                                    \                                    /
                                                                   con2d                             cont2u
                                                                      \                                /                                 
                                                                       -----> con3 -> dfd_res_33 -----> 
                    
---------------------------------------------------------------------------------------- */

// --------------------------------- Conv Filter Setup ------------------------------------
template <long num_filters, typename SUBNET> using con2d = dlib::con<num_filters, 2, 2, 2, 2, SUBNET>;

template <long num_filters, typename SUBNET> using con1 = dlib::con<num_filters, 1, 1, 1, 1, SUBNET>;
template <long num_filters, typename SUBNET> using con3 = dlib::con<num_filters, 3, 3, 1, 1, SUBNET>;

template <long num_filters, typename SUBNET> using cont2u = dlib::cont<num_filters, 2, 2, 2, 2, SUBNET>;

// --------------------------------- block definitions ------------------------------------

template <int N, typename SUBNET> using cbp3_blk = dlib::prelu<dlib::bn_con<con3<N, SUBNET>>>;

template <int N1, int N2, typename SUBNET> using dfd_block_33 = con3<N1, dlib::prelu<dlib::bn_con<con3<N2, SUBNET>>>>;

template <int N1,int N2,typename SUBNET>
using dfd_res_33 = dlib::add_prev1<dfd_block_33<N1,N2, dlib::tag1<SUBNET>>>;

// --------------------------------- tag definitions --------------------------------------
template <typename SUBNET> using DTO_0 = dlib::add_tag_layer<200, SUBNET>;
template <typename SUBNET> using DTI_0 = dlib::add_tag_layer<201, SUBNET>;
template <typename SUBNET> using DTO_1 = dlib::add_tag_layer<202, SUBNET>;
template <typename SUBNET> using DTI_1 = dlib::add_tag_layer<203, SUBNET>;
template <typename SUBNET> using DTO_2 = dlib::add_tag_layer<204, SUBNET>;
template <typename SUBNET> using DTI_2 = dlib::add_tag_layer<205, SUBNET>;
template <typename SUBNET> using DTO_3 = dlib::add_tag_layer<206, SUBNET>;
template <typename SUBNET> using DTI_3 = dlib::add_tag_layer<207, SUBNET>;

// ----------------------------------------------------------------------------------------

using dfd_net_type = dlib::loss_multiclass_log_per_pixel<
	con3<256, dlib::prelu<dfd_res_33<64,64, con3<64,
    
    dlib::concat2<DTO_0, DTI_0,
    DTI_0<cont2u<64, dfd_res_33<128,128, con3<128,

    dlib::concat2<DTO_1, DTI_1,
    DTI_1<cont2u<128, dfd_res_33<256,256, con3<256,

    dlib::concat2<DTO_2, DTI_2,
    DTI_2<cont2u<256, dfd_res_33<256,256, con3<256,
    
    dlib::concat2<DTO_3, DTI_3,
    DTI_3<cont2u<256, dfd_res_33<512,512, con3<512, con2d<512,
    
    DTO_3<dfd_res_33<256,256, con3<256, con2d<256,
    
    DTO_2<dfd_res_33<256,256, con3<256, con2d<256,
    
    DTO_1<dfd_res_33<128,128, con3<128, con2d<128,
    
    DTO_0<dfd_res_33<64,64, con3<64, dlib::input<std::array<dlib::matrix<uint16_t>, img_depth>>>>
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>;

// ----------------------------------------------------------------------------------------
// Configuration function
// ----------------------------------------------------------------------------------------

template <typename net_type>
void config_net(net_type &net, std::vector<uint32_t> params)
{

    net = net_type(dlib::num_con_outputs(params[0]),
        dlib::num_con_outputs(params[1]),
        dlib::num_con_outputs(params[2]),
        dlib::num_con_outputs(params[3]),
        dlib::num_con_outputs(params[4]),
        dlib::num_con_outputs(params[5]),
        dlib::num_con_outputs(params[6]),
        dlib::num_con_outputs(params[7]),
        dlib::num_con_outputs(params[8]),
        dlib::num_con_outputs(params[9]),
        dlib::num_con_outputs(params[10]),
        dlib::num_con_outputs(params[11]),
        dlib::num_con_outputs(params[12]),
        dlib::num_con_outputs(params[13]),
        dlib::num_con_outputs(params[14]),
        dlib::num_con_outputs(params[15]),
        dlib::num_con_outputs(params[16]),
        dlib::num_con_outputs(params[17]),
        dlib::num_con_outputs(params[18]),
        dlib::num_con_outputs(params[19]),
        dlib::num_con_outputs(params[20]),
        dlib::num_con_outputs(params[21]),
        dlib::num_con_outputs(params[22]),
        dlib::num_con_outputs(params[23]),
        dlib::num_con_outputs(params[24]),
        dlib::num_con_outputs(params[25]),
        dlib::num_con_outputs(params[26]),
        dlib::num_con_outputs(params[27]),
        dlib::num_con_outputs(params[28]),
        dlib::num_con_outputs(params[29]),
        dlib::num_con_outputs(params[30]),
        dlib::num_con_outputs(params[31]),
        dlib::num_con_outputs(params[32]),
        dlib::num_con_outputs(params[33]),
        dlib::num_con_outputs(params[34]),
        dlib::num_con_outputs(params[35]));
        
}   // end of config_net

// ----------------------------------------------------------------------------------------
//  GORGON Functions
// ----------------------------------------------------------------------------------------


//gorgon_capture<1> gc_01(256,128,3,3);
//gorgon_capture<5> gc_02(128,128,3,3);
//gorgon_capture<8> gc_03(128,192,3,3);
//gorgon_capture<14> gc_04(128,128,3,3);
//gorgon_capture<17> gc_05(128,320,3,3);
//gorgon_capture<23> gc_06(256,256,3,3);
//gorgon_capture<26> gc_07(256,576,3,3);
//gorgon_capture<32> gc_08(512,512,3,3);
//gorgon_capture<35> gc_09(512,64,3,3);
//gorgon_capture<40> gc_10(512,512,3,3);
//gorgon_capture<43> gc_11(512,64,3,3);
//gorgon_capture<48> gc_12(256,256,3,3);
//gorgon_capture<51> gc_13(256,64,3,3);
//gorgon_capture<56> gc_14(128,128,3,3);
//gorgon_capture<59> gc_15(128,img_depth,3,3);

void init_gorgon(std::string save_location)
{
    //gc_01.init((save_location + "l01"));
    //gc_02.init((save_location + "l05"));
    //gc_03.init((save_location + "l08"));
    //gc_04.init((save_location + "l14"));
    //gc_05.init((save_location + "l17"));
    //gc_06.init((save_location + "l23"));
    //gc_07.init((save_location + "l26"));
    //gc_08.init((save_location + "l32"));
    //gc_09.init((save_location + "l35"));
    //gc_10.init((save_location + "l40"));
    //gc_11.init((save_location + "l43"));
    //gc_12.init((save_location + "l48"));
    //gc_13.init((save_location + "l51"));
    //gc_14.init((save_location + "l56"));
    //gc_15.init((save_location + "l59"));

}

template<typename net_type>
void save_gorgon(net_type &dfd_net, uint64_t one_step_calls)
{
    //gc_01.save_params(dfd_net, one_step_calls);
    //gc_02.save_params(dfd_net, one_step_calls);
    //gc_03.save_params(dfd_net, one_step_calls);
    //gc_04.save_params(dfd_net, one_step_calls);
    //gc_05.save_params(dfd_net, one_step_calls);
    //gc_06.save_params(dfd_net, one_step_calls);
    //gc_07.save_params(dfd_net, one_step_calls);
    //gc_08.save_params(dfd_net, one_step_calls);
    //gc_09.save_params(dfd_net, one_step_calls);
    //gc_10.save_params(dfd_net, one_step_calls);
    //gc_11.save_params(dfd_net, one_step_calls);
    //gc_12.save_params(dfd_net, one_step_calls);
    //gc_13.save_params(dfd_net, one_step_calls);
    //gc_14.save_params(dfd_net, one_step_calls);
    //gc_15.save_params(dfd_net, one_step_calls);
}


void close_gorgon(void)
{
    //gc_01.close_stream();
    //gc_02.close_stream();
    //gc_03.close_stream();
    //gc_04.close_stream();
    //gc_05.close_stream();
    //gc_06.close_stream();
    //gc_07.close_stream();
    //gc_08.close_stream();
    //gc_09.close_stream();
    //gc_10.close_stream();
    //gc_11.close_stream();
    //gc_12.close_stream();
    //gc_13.close_stream();
    //gc_14.close_stream();
    //gc_15.close_stream();
}


#endif //NET_DEFINITION_H
