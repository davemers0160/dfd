#ifndef CYCLIC_ANALYSIS_H
#define CYCLIC_ANALYSIS_H

// C++ includes
#include <cstdint>
#include <vector>
#include <cmath>
#include <algorithm>
#include <list>
#include <iterator>

// Custom includes
#include "ycrcb_pixel.h"

// DLIB includes
#include <dlib/image_transforms.h>
#include <dlib/image_transforms/interpolation.h>
#include <dlib/threads.h>

using namespace std;


class cyclic_analysis
{
	dlib::chip_dims dims = dlib::chip_dims(30, 30);
	
	unsigned long x_overlap = 0;    // number of pixels to overlap in the x direction
	unsigned long y_overlap = 0;    // number of pixels to overlap in the y direction
	
// ----------------------------------------------------------------------------------------
public:

	const dlib::chip_dims& get_chip_dims() const { return dims; }

	void set_chip_dims(const dlib::chip_dims& dims_) { 
		dims = dims_; 
		update_step();
	}

	void set_chip_dims(unsigned long rows, unsigned long cols)
	{ set_chip_dims(dlib::chip_dims(rows, cols));  }
	
	void set_chip_dims (std::pair<unsigned long, unsigned long> p)
	{ set_chip_dims(dlib::chip_dims(p.first, p.second)); }

	void set_x_overlap(unsigned long value)
	{ 
		x_overlap = std::max(0UL, value); 
		update_step();
	}

	void set_y_overlap(unsigned long value)
	{  
		y_overlap = std::max(0UL, value);  
		update_step();
	}
	
	void set_overlap(unsigned long x, unsigned long y)
	{ 
		set_x_overlap(x);
		set_y_overlap(y);
	}
	
	unsigned long get_x_overlap() const { return x_overlap; }
	unsigned long get_y_overlap() const { return y_overlap; }
	

	
// ----------------------------------------------------------------------------------------
	
	// this one does a single image and single vector of labels
	template <
		typename net_type,
		typename image_type1,
        typename image_type2
		>
	void operator()
	(
		net_type &net,
		const image_type1 img,
		image_type2 &dm
	)
	{
		
		// get the detects from the net and then prune any overlapping detects
		get_detects(net, img, dm);
		
		// now run the detects through the metrics		
		
		
	}	// end of operator()
	
	
// ----------------------------------------------------------------------------------------
	
private:

	uint64_t row_step;
	uint64_t col_step;
	
// ----------------------------------------------------------------------------------------
	
	void update_step()
	{
		row_step = dims.rows - y_overlap;
        col_step = dims.cols - x_overlap; 
	}	// end of update_step
	
	
// ----------------------------------------------------------------------------------------	
	
	template<
		typename net_type,
		typename image_type1,
        typename image_type2
		>
	void get_detects
	(
		net_type &net,
		const image_type1 img,
		image_type2 &dm
	)        
	{
		uint64_t idx;
		uint64_t row, col;
        uint64_t chip_w, chip_h;

		uint64_t img_width = img[0].nc();
		uint64_t img_height = img[0].nr();
						
		image_type1 tmp_img;
		
		dlib::rectangle rect(dims.cols, dims.rows);

		for (row = 0; row < img_height; row += row_step)
		{
			for (col = 0; col < img_width; col += col_step)
			{
				//rect = dlib::rectangle(dims.cols, dims.rows);

				//rect = dlib::move_rect(rect, dlib::point(col,row));
				//dlib::chip_dims tmp_dims = dlib::chip_dims(dims.rows, dims.cols);

                chip_w = dims.cols;
                chip_h = dims.rows;

				if ((col + chip_w) >= img_width)
				{
					//rect.set_right(img_width - 1);  
					//tmp_dims.cols = rect.width();
                    chip_w = img_width - col;
					//col = img_width;
				}

				if ((row + chip_h) >= img_height)
				{
					//rect.set_bottom(img_height - 1);
					//tmp_dims.rows = rect.height();
					//row = img_height;
                    chip_h = img_height - row;
				}

				// extract the image chip
				//dlib::extract_image_chip(img, dlib::chip_details(rect, tmp_dims, 0.0), tmp_img);
                for (idx = 0; idx < img.size(); ++idx)
                {
                    tmp_img[idx] = dlib::zeros_matrix<uint16_t>(dims.rows, dims.cols);
                    dlib::set_subm(tmp_img[idx], 0, 0, chip_h, chip_w) = dlib::subm(img[idx], row, col, chip_h, chip_w);
                }

				//run the chip through the detector
				image_type2 tmp_dm = net(tmp_img);
                
                
                dlib::set_subm(dm, row, col, chip_h, chip_w) = dlib::subm(tmp_dm, 0, 0, chip_h, chip_w);
                
                // if an overlap is used look at that and then take the average of the overlap and put it back into the depthmap
                // if((row<0) && (col<0))
				// {
                    // image_type2 over_dm = dlib::subm(dm, row, col, dims.rows-y_overlap, dims.cols-x_overlap)*0.5;
                    // dlib::set_subm(dm, row, col, dims.rows-y_overlap, dims.cols-x_overlap) = over_dm;
                // }

				
			}
		}            
		
		
	}   // end of get_detects	
    
// ----------------------------------------------------------------------------------------
	
};  // end of class

#endif	// CYCLIC_ANALYSIS_H
