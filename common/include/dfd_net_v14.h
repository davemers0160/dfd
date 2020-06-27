#ifndef NET_DEFINITION_H
#define NET_DEFINITION_H

#include <cstdint>
#include <array>

// dlib includes
#include "dlib/dnn.h"
#include "dlib/dnn/core.h"

//#include "ssim_loss.h"
//#include "l1_loss.h"
#include "dfd_dnn_input.h"
//#include "gorgon_capture.h"

extern const uint32_t img_depth = 6;
extern const uint32_t secondary = 1;

//const std::array<float, img_depth> avg_color{ 107.1, 111.9, 132.8, 107.1, 111.9, 132.8 };
//const std::array<float, img_depth> avg_color{ 110.9551, 98.1071, 69.9262 };

// --------------------------------- Conv Filter Setup ------------------------------------
template <long num_filters, typename SUBNET> using con2d = dlib::con<num_filters, 2, 2, 2, 2, SUBNET>;
template <long num_filters, typename SUBNET> using con3d = dlib::con<num_filters, 3, 3, 2, 2, SUBNET>;
template <long num_filters, typename SUBNET> using con5d = dlib::con<num_filters, 5, 5, 2, 2, SUBNET>;
template <long num_filters, typename SUBNET> using con7d = dlib::con<num_filters, 7, 7, 2, 2, SUBNET>;
template <long num_filters, typename SUBNET> using con9d = dlib::con<num_filters, 9, 9, 2, 2, SUBNET>;

template <long num_filters, typename SUBNET> using con1 = dlib::con<num_filters, 1, 1, 1, 1, SUBNET>;
template <long num_filters, typename SUBNET> using con3 = dlib::con<num_filters, 3, 3, 1, 1, SUBNET>;
template <long num_filters, typename SUBNET> using con5 = dlib::con<num_filters, 5, 5, 1, 1, SUBNET>;
template <long num_filters, typename SUBNET> using con7 = dlib::con<num_filters, 7, 7, 1, 1, SUBNET>;
template <long num_filters, typename SUBNET> using con9 = dlib::con<num_filters, 9, 9, 1, 1, SUBNET>;

template <long num_filters, typename SUBNET> using cont2 = dlib::cont<num_filters, 2, 2, 1, 1, SUBNET>;
template <long num_filters, typename SUBNET> using cont3 = dlib::cont<num_filters, 3, 3, 1, 1, SUBNET>;
template <long num_filters, typename SUBNET> using cont5 = dlib::cont<num_filters, 5, 5, 1, 1, SUBNET>;
template <long num_filters, typename SUBNET> using cont7 = dlib::cont<num_filters, 7, 7, 1, 1, SUBNET>;
template <long num_filters, typename SUBNET> using cont9 = dlib::cont<num_filters, 9, 9, 1, 1, SUBNET>;

template <long num_filters, typename SUBNET> using cont2u = dlib::cont<num_filters, 2, 2, 2, 2, SUBNET>;
template <long num_filters, typename SUBNET> using cont3u = dlib::cont<num_filters, 3, 3, 2, 2, SUBNET>;
template <long num_filters, typename SUBNET> using cont4u = dlib::cont<num_filters, 4, 4, 2, 2, SUBNET>;
template <long num_filters, typename SUBNET> using cont5u = dlib::cont<num_filters, 5, 5, 2, 2, SUBNET>;
template <long num_filters, typename SUBNET> using cont7u = dlib::cont<num_filters, 7, 7, 2, 2, SUBNET>;
template <long num_filters, typename SUBNET> using cont9u = dlib::cont<num_filters, 9, 9, 2, 2, SUBNET>;

// ----------------------------------------------------------------------------------------

/*
blockNM => con3(N)->bn->prelu->con3(N)



input -> block55 [dtago0] -------------------------------> [dtagi0] block55 -> con1(256) -> output
            \                                               /
           con3d                                        cont3u
              \                                           /
             block55 [dtago1] -----------------------> [dtagi1] block55
                \                                       /
               con3d                                cont3u
                  \                                   /
                 block53 [dtago2] ---------------> [dtagi2] block35
                    \                               /
                    con3d                       cont3u
                       \                          /
                        --------> block33 -------> 
                    
*/
template<typename SUBNET> using mp2 = dlib::max_pool<2, 2, 2, 2, SUBNET>;
template<typename SUBNET> using up2 = dlib::upsample<2, SUBNET>;

// ----------------------------------------------------------------------------------------
// block definitions with batch norm
// ----------------------------------------------------------------------------------------
template <int N, typename SUBNET> using cbp3_blk = dlib::prelu<dlib::bn_con<con3<N, SUBNET>>>;
template <int N, typename SUBNET> using cbp5_blk = dlib::prelu<dlib::bn_con<con5<N, SUBNET>>>;

template <int N1, int N2, typename SUBNET> using dfd_block_33 = con3<N1, dlib::prelu<dlib::bn_con<con3<N2, SUBNET>>>>;
template <int N1, int N2, typename SUBNET> using dfd_block_35 = con3<N1, dlib::prelu<dlib::bn_con<con5<N2, SUBNET>>>>;
template <int N1, int N2, typename SUBNET> using dfd_block_37 = con3<N1, dlib::prelu<dlib::bn_con<con7<N2, SUBNET>>>>;

template <int N1, int N2, typename SUBNET> using dfd_block_53 = con5<N1, dlib::prelu<dlib::bn_con<con3<N2, SUBNET>>>>;
template <int N1, int N2, typename SUBNET> using dfd_block_55 = con5<N1, dlib::prelu<dlib::bn_con<con5<N2, SUBNET>>>>;
template <int N1, int N2, typename SUBNET> using dfd_block_57 = con5<N1, dlib::prelu<dlib::bn_con<con7<N2, SUBNET>>>>;

template <int N1, int N2, typename SUBNET> using dfd_block_73 = con7<N1, dlib::prelu<dlib::bn_con<con3<N2, SUBNET>>>>;
template <int N1, int N2, typename SUBNET> using dfd_block_75 = con7<N1, dlib::prelu<dlib::bn_con<con5<N2, SUBNET>>>>;
template <int N1, int N2, typename SUBNET> using dfd_block_77 = con7<N1, dlib::prelu<dlib::bn_con<con7<N2, SUBNET>>>>;

template <int N1,int N2,typename SUBNET>
using dfd_res_33 = dlib::add_prev1<dfd_block_33<N1,N2, dlib::tag1<SUBNET>>>;

template <int N1,int N2,typename SUBNET>
using dfd_res_35 = dlib::add_prev1<dfd_block_35<N1,N2, dlib::tag1<SUBNET>>>;

template <int N1,int N2,typename SUBNET>
using dfd_res_37 = dlib::add_prev1<dfd_block_37<N1,N2, dlib::tag1<SUBNET>>>;

template <int N1,int N2,typename SUBNET>
using dfd_res_53 = dlib::add_prev1<dfd_block_53<N1,N2, dlib::tag1<SUBNET>>>;

template <int N1,int N2,typename SUBNET>
using dfd_res_55 = dlib::add_prev1<dfd_block_55<N1,N2, dlib::tag1<SUBNET>>>;

template <int N1,int N2,typename SUBNET>
using dfd_res_57 = dlib::add_prev1<dfd_block_57<N1,N2, dlib::tag1<SUBNET>>>;

// ----------------------------------------------------------------------------------------
// block definitions with affine
// ----------------------------------------------------------------------------------------
template <int N, typename SUBNET> using acbp3_blk = dlib::prelu<dlib::affine<con3<N, SUBNET>>>;

template <int N1, int N2, typename SUBNET> using adfd_block_33 = con3<N1, dlib::prelu<dlib::affine<con3<N2, SUBNET>>>>;

template <int N1,int N2,typename SUBNET>
using adfd_res_33 = dlib::add_prev1<adfd_block_33<N1,N2, dlib::tag1<SUBNET>>>;

// ----------------------------------------------------------------------------------------
// tag definitions
// ----------------------------------------------------------------------------------------
template <typename SUBNET> using dtago0 = dlib::add_tag_layer<200, SUBNET>;
template <typename SUBNET> using dtagi0 = dlib::add_tag_layer<201, SUBNET>;
template <typename SUBNET> using dtago1 = dlib::add_tag_layer<202, SUBNET>;
template <typename SUBNET> using dtagi1 = dlib::add_tag_layer<203, SUBNET>;
template <typename SUBNET> using dtago2 = dlib::add_tag_layer<204, SUBNET>;
template <typename SUBNET> using dtagi2 = dlib::add_tag_layer<205, SUBNET>;

// ----------------------------------------------------------------------------------------
// Training Version with batch norm
// ----------------------------------------------------------------------------------------
using dfd_net_type = dlib::loss_multiclass_log_per_pixel<
    con3<256, dlib::prelu<dlib::bn_con<dfd_res_33<256, 256, cbp3_blk<256,

    dlib::concat2<dtago1, dtagi1,
    dtagi1<cont2u<256, dfd_res_33<512, 512, cbp3_blk<512,

    dlib::concat2<dtago2, dtagi2,
    dtagi2<cont2u<512, dfd_res_33<1024, 1024, con2d<1024,
    
    dtago2<dfd_res_33<512, 512, con2d<512,
    
    dtago1<dfd_res_33<256, 256, cbp3_blk<256, 
    dlib::tag10<dlib::input_dfd_array<uint16_t, img_depth>>
    //dlib::input_dfd_array<uint16_t, img_depth>
    >>> >>> >>>> > >>>> > >>>>> >;
    
// ----------------------------------------------------------------------------------------
// Version with affine layers
// ----------------------------------------------------------------------------------------
using adfd_net_type = dlib::loss_multiclass_log_per_pixel<
    con3<256, dlib::prelu<dlib::affine<adfd_res_33<256, 256, acbp3_blk<256,

    dlib::concat2<dtago1, dtagi1,
    dtagi1<cont2u<256, adfd_res_33<512,512,acbp3_blk<512,

    dlib::concat2<dtago2, dtagi2,
    dtagi2<cont2u<512, adfd_res_33<1024,1024,con2d<1024,
    
    dtago2<adfd_res_33<512,512,con2d<512,
    
    dtago1<adfd_res_33<256, 256, acbp3_blk<256, 
    //dlib::input<std::array<dlib::matrix<uint16_t>, img_depth>
    dlib::tag10<dlib::input_dfd_array<uint16_t, img_depth>>
    >>> >>> >>>> > >>>> > >>>>> >;
 

// ----------------------------------------------------------------------------------------
// Configuration function
// ----------------------------------------------------------------------------------------

template <typename net_type>
//void config_net(net_type& net, std::array<float, img_depth>& avg_color, std::vector<uint32_t>& params)
net_type config_net(std::array<float, img_depth>& avg_color, std::vector<uint32_t>& fn)
{
    net_type net(dlib::num_con_outputs(256),
        dlib::num_con_outputs(fn[1]),
        dlib::num_con_outputs(fn[2]),
        dlib::num_con_outputs(fn[3]),
        dlib::num_con_outputs(fn[4]),
        dlib::num_con_outputs(fn[5]),
        dlib::num_con_outputs(fn[6]),
        dlib::num_con_outputs(fn[7]),
        dlib::num_con_outputs(fn[8]),
        dlib::num_con_outputs(fn[9]),
        dlib::num_con_outputs(fn[10]),
        dlib::num_con_outputs(fn[11]),
        dlib::num_con_outputs(fn[12]),
        dlib::num_con_outputs(fn[13]),
        dlib::num_con_outputs(fn[14]),
        dlib::num_con_outputs(fn[15]),
        dlib::num_con_outputs(fn[16]),
        dlib::num_con_outputs(fn[17])
    );

    net.subnet().layer_details().set_num_filters(fn[0]);
    dlib::layer<net_type::num_layers - 1>(net).set_avg_colors(avg_color);
    dlib::layer<net_type::num_layers - 1>(net).set_scale(1.0/256.0);

    return net;

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
