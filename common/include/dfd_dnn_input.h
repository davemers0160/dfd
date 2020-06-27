#ifndef DLIB_DNN_DFD_INPUT_ARRAY_H_
#define DLIB_DNN_DFD_INPUT_ARRAY_H_

#include <cstdint>
#include <sstream>
#include <array>
#include <algorithm>

#include "dlib/matrix.h"
#include "dlib/pixel.h"
#include "dlib/image_processing.h"
#include "dlib/cuda/tensor_tools.h"

namespace dlib
{

// ----------------------------------------------------------------------------------------

    //template <typename T, long NR, long NC, typename MM, typename L, size_t K>
    template <typename T, size_t K>
    class input_dfd_array
    {
    public:
        typedef std::array<matrix<T>, K> input_type;
        
        input_dfd_array() 
        {
            avg_color.fill(128.0);  // fill in the average color with 128
            set_scale(1.0 / 256.0);
        }

        //input_dfd_array(std::array<float, K> avg_color_)
   //     input_dfd_array(float avg_color_)
   //     {
            ////std::copy(std::begin(avg_color_), std::end(avg_color_), std::begin(avg_color));
   //         //for (uint64_t idx = 0; idx < K; ++idx)
   //         //    avg_color[idx] = avg_color_[idx];
   //         avg_color[0] = avg_color_;

   //     }

        void set_avg_colors(std::array<float, K> avg_color_)
        {
            std::copy(std::begin(avg_color_), std::end(avg_color_), std::begin(avg_color));
        }

        void set_scale(float scale_)
        {
            scale = scale_;
        }

        float get_scale(void) const { return scale; }

        bool image_contained_point ( const tensor& data, const point& p) const { return get_rect(data).contains(p); }
        drectangle tensor_space_to_image_space ( const tensor& /*data*/, drectangle r) const { return r; }
        drectangle image_space_to_tensor_space ( const tensor& /*data*/, double /*scale*/, drectangle r ) const { return r; }

        template <typename forward_iterator>
        void to_tensor (
            forward_iterator ibegin,
            forward_iterator iend,
            resizable_tensor& data
        ) const
        {
            DLIB_CASSERT(std::distance(ibegin,iend) > 0);
            DLIB_CASSERT(ibegin->size() != 0, "When using std::array<matrix> inputs you can't give 0 sized arrays.");
            const auto nr = (*ibegin)[0].nr();
            const auto nc = (*ibegin)[0].nc();
            // make sure all the input matrices have the same dimensions
            for (auto i = ibegin; i != iend; ++i)
            {
                for (size_t k = 0; k < K; ++k)
                {
                    const auto& arr = *i;
                    DLIB_CASSERT(arr[k].nr()==nr && arr[k].nc()==nc,
                        "\t input::to_tensor()"
                        << "\n\t When using std::array<matrix> as input, all matrices in a batch must have the same dimensions."
                        << "\n\t nr: " << nr
                        << "\n\t nc: " << nc
                        << "\n\t k:  " << k 
                        << "\n\t arr[k].nr(): " << arr[k].nr()
                        << "\n\t arr[k].nc(): " << arr[k].nc()
                    );
                }
            }

            
            // initialize data to the right size to contain the stuff in the iterator range.
            data.set_size(std::distance(ibegin,iend), K, nr, nc);

            auto ptr = data.host();
            for (auto i = ibegin; i != iend; ++i)
            {
                for (size_t k = 0; k < K; ++k)
                {
                    for (long r = 0; r < nr; ++r)
                    {
                        for (long c = 0; c < nc; ++c)
                        {
                            //if (is_same_type<T,unsigned char>::value)
                                *ptr++ = ((*i)[k](r,c) - avg_color[k]) * scale;	// this normalizes the input between [-0.5, 0.5)
                            //else
                            //    *ptr++ = (*i)[k](r,c);
                        }
                    }
                }
            }

        }

        friend void serialize(const input_dfd_array& item, std::ostream& out)
        {
            serialize("input_dfd_array", out);
            for (size_t k = 0; k < K; ++k)
                serialize(item.avg_color[k], out);
            serialize(item.scale, out);
        }

        friend void deserialize(input_dfd_array& item, std::istream& in)
        {
            std::string version;
            deserialize(version, in);
            if (version != "input_dfd_array")
                throw serialization_error("Unexpected version found while deserializing dlib::input_dfd_array.");
            for (size_t k = 0; k < K; ++k)
                deserialize(item.avg_color[k], in);	
            deserialize(item.scale, in);
        }

        friend std::ostream& operator<<(std::ostream& out, const input_dfd_array& item)
        {
            out << "input_dfd_array(";
            std::string tmp = "";
            for (size_t k = 0; k < K-1; ++k)
                out << item.avg_color[k] << ",";
            out << item.avg_color[K-1] << ") ";
            out << "array_depth=" << K << ", scale=" << item.scale;
            return out;
        }

        friend void to_xml(const input_dfd_array& item, std::ostream& out)
        {
            out << "<input_dfd_array ";			
            for (size_t k = 0; k < K-1; ++k)
                out << "c" << k << "='" << item.avg_color[k] << "' ";
            out << "c" << K - 1 << "='" << item.avg_color[K - 1] << "' ";
            out << "scale='" << item.scale;
            out << "'/>";
        }

    private:
        std::array<float, K> avg_color;

        float scale;

    };  // end of class

}   // end of namespace

#endif	// DLIB_DNN_DFD_INPUT_ARRAY_H_
