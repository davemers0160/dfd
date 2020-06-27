# Depth from Defocus Repository 
This repository is a combination of several depth from defocus projects based on OpenCV and the dlib library.

## Dependencies

The code in this repository has the following dependecies:

1. [CMake 2.8.12+](https://cmake.org/download/ )
2. [davemers0160 common code repository](https://github.com/davemers0160/Common )
3. [davemers0160 dlib-contrib repository](https://github.com/davemers0160/dlib-contrib )
4. [dlib library](http://dlib.net/ )
5. [OpenCV v4+](https://opencv.org/releases/ )

## Repository Breakdown

### common

This folder contains the project code that is shared between all of the projects.

### dfd_dnn_analysis

This folder contains the project code that runs the performance analysis of a given network against a given dataset.

### dfd_dnn_trainer

This folder contains the project code that runs the deep learning training of the DfD-Net network.

### dfd_graphcuts

This folder contains the project code that applies the graph cuts method to the depth from defocus challenge. 

## References

When using this work please cite the following:

[3-D SCENE RECONSTRUCTION FOR PASSIVE RANGING USING DEPTH FROM DEFOCUS AND DEEP LEARNING](https://hammer.figshare.com/articles/3-D_SCENE_RECONSTRUCTION_FOR_PASSIVE_RANGING_USING_DEPTH_FROM_DEFOCUS_AND_DEEP_LEARNING/8938376/1)

DataCite
```
Emerson, David Ross (2019): 3-D SCENE RECONSTRUCTION FOR PASSIVE RANGING USING DEPTH FROM DEFOCUS AND DEEP LEARNING. figshare. Thesis. https://doi.org/10.25394/PGS.8938376.v1
```

BiBtex
```
@article{Emerson2019,
author = "David Ross Emerson",
title = "{3-D SCENE RECONSTRUCTION FOR PASSIVE RANGING USING DEPTH FROM DEFOCUS AND DEEP LEARNING}",
year = "2019",
month = "10",
url = "https://hammer.figshare.com/articles/3-D_SCENE_RECONSTRUCTION_FOR_PASSIVE_RANGING_USING_DEPTH_FROM_DEFOCUS_AND_DEEP_LEARNING/8938376",
doi = "10.25394/PGS.8938376.v1"
}
```

