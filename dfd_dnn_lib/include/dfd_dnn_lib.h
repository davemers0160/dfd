#ifndef OBJ_DET_DLL_H
#define OBJ_DET_DLL_H

//#define EXTERN_C
//#include <cstdint>
//#include <string>
//#include <vector>

#if defined(_WIN32) | defined(__WIN32__) | defined(__WIN32) | defined(_WIN64) | defined(__WIN64)

#if defined(BUILD_LIB)

#ifdef LIB_EXPORTS
#define LIB_API __declspec(dllexport)
#else
#define LIB_API __declspec(dllimport)
#endif

#else

#define LIB_API

#endif

#else
#define LIB_API

#endif

// ----------------------------------------------------------------------------------------
const int label_size = 256;


// ----------------------------------------------------------------------------------------
struct layer_struct
{
    unsigned int k;
    unsigned int n;
    unsigned int nr;
    unsigned int nc;
    unsigned int size;
};

// ----------------------------------------------------------------------------------------
struct detection_struct
{
    unsigned int x;
    unsigned int y;
    unsigned int w;
    unsigned int h;
    char name[label_size];

    detection_struct()
    {
        x = 0;
        y = 0;
        w = 0;
        h = 0;
        name[0] = 0;
    }

    detection_struct(unsigned int x_, unsigned int y_, unsigned int w_, unsigned int h_, const char name_[])
    {
        x = x_;
        y = y_;
        w = w_;
        h = h_;
        strcpy(name, name_);
    }

};

// ----------------------------------------------------------------------------------------
struct detection_center
{
    unsigned int x;
    unsigned int y;
    char name[label_size];

    detection_center()
    {
        x = 0;
        y = 0;
        name[0] = 0;
    }

    detection_center(unsigned int x_, unsigned int y_, const char name_[])
    {
        x = x_;
        y = y_;
        strcpy(name, name_);
    }

};

// ----------------------------------------------------------------------------------------
struct window_struct
{
    unsigned int w;
    unsigned int h;
    char label[label_size];
};

// ----------------------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif
    // This function will initialize the network and load the required weights
    LIB_API void init_net(const char *net_name, unsigned int *num_classes, struct window_struct* &det_win, unsigned int *num_win);
#ifdef __cplusplus
}
#endif

// ----------------------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif
    // This function will take an grayscale image in unsigned char row major order [r0,c0, r0,c1, r0,c2,...]
    // as an input and will return the bounding boxes of the detections in the image
    LIB_API void run_net(unsigned char* image_f1, unsigned char* image_f2, unsigned int nr, unsigned int nc, unsigned short* depth_map);
#ifdef __cplusplus
}
#endif

// ----------------------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif
    // This function will take an grayscale image in unsigned char row major order [r0,c0, r0,c1, r0,c2,...]
    // as an input and will return the centers of the detections in the image
    LIB_API void get_detections(unsigned char* input_img, unsigned int nr, unsigned int nc, unsigned int* num_dets, struct detection_struct*& dets);
#ifdef __cplusplus
}
#endif

// ----------------------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif
    // This function will take an grayscale image in unsigned char row major order [r0,c0, r0,c1, r0,c2,...]
    // as an input and will return the centers of the detections in the image
    LIB_API void get_cropped_detections(unsigned char* input_img, unsigned int nr, unsigned int nc, unsigned int x, unsigned int y, unsigned int w, unsigned int h, unsigned int* num_dets, struct detection_struct*& dets);
#ifdef __cplusplus
}
#endif

// ----------------------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif
    // This function will output a vector of the output layer for the final classification layer
    LIB_API void close_lib();
#ifdef __cplusplus
}
#endif


// ----------------------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif
    // This function will output a vector of the output layer for the final classification layer
    //OBJ_DLL_API void get_combined_output(unsigned char* input_img, unsigned int nr, unsigned int nc, float* &data_params);
    LIB_API void get_combined_output(struct layer_struct* data, const float*& data_params);
#ifdef __cplusplus
}
#endif
// ----------------------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif
    // This function will output a vector of the output layer for the final classification layer
    LIB_API void get_layer_01(struct layer_struct *data, const float* &data_params);
#ifdef __cplusplus
}
#endif

//// ----------------------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif
    OBJ_DLL_API void get_input_layer(struct layer_struct *data, const float* &data_params);
#ifdef __cplusplus
}
#endif

// ----------------------------------------------------------------------------------------
//#ifdef __cplusplus
//extern "C" {
//#endif
//    //MNIST_DLL_API void get_layer_05(struct layer_struct &data, const float* &data_params);
//    OBJ_DLL_API void get_layer_05(layer_struct *data, const float **data_params);
//#ifdef __cplusplus
//}
//#endif
//
//// ----------------------------------------------------------------------------------------
//#ifdef __cplusplus
//extern "C" {
//#endif
//    //MNIST_DLL_API void get_layer_08(struct layer_struct &data, const float* &data_params);
//    OBJ_DLL_API void get_layer_08(layer_struct *data, const float **data_params);
//#ifdef __cplusplus
//}
//#endif
//
//// ----------------------------------------------------------------------------------------
//#ifdef __cplusplus
//extern "C" {
//#endif
//    //MNIST_DLL_API void get_layer_09(struct layer_struct &data, const float* &data_params);
//    OBJ_DLL_API void get_layer_09(layer_struct *data, const float **data_params);
//#ifdef __cplusplus
//}
//#endif
//
//// ----------------------------------------------------------------------------------------
//#ifdef __cplusplus
//extern "C" {
//#endif
//    //MNIST_DLL_API void get_layer_12(struct layer_struct &data, const float* &data_params);
//    OBJ_DLL_API void get_layer_12(layer_struct *data, const float **data_params);
//#ifdef __cplusplus
//}
//#endif

#endif  // OBJ_DET_DLL_H
