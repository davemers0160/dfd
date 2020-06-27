#ifndef CALC_SILOG_ERROR_
#define CALC_SILOG_ERROR_

#include <cstdint>

// dlib includes
#include <dlib/dnn.h>
#include <dlib/matrix.h>

template<typename T>
double calc_silog_error(dlib::matrix<T> gt, dlib::matrix<T> map)
{

    // check to make sure that the inputs are the same size
    DLIB_CASSERT((gt.nr() == map.nr()) && (gt.nc() == map.nc()), "gt size : " << gt.nr() << "x" << gt.nc() << "; map size : " << map.nr() << "x" << map.nc());
    
    double sz = (double)(gt.nr()*gt.nc());
    
    dlib::matrix<float> log_map = dlib::log(dlib::matrix_cast<float>(gt)+1) - dlib::log(dlib::matrix_cast<float>(map)+1);

    double p1 = dlib::sum(dlib::squared(log_map)) / sz;
    double p2 = (dlib::sum(log_map)*dlib::sum(log_map)) / (sz*sz);

    double d = p1 - p2;

    return d;

}   // end of calc_silog_error


#endif // CALC_SILOG_ERROR_