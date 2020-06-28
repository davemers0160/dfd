# DFD DNN Analysis Project
Project to analyze the results of the DfD-Net training.

## Dependencies

The code in this repository has the following dependecies:

1. [CMake 2.8.12+](https://cmake.org/download/ )
2. [davemers0160 common code repository](https://github.com/davemers0160/Common )
3. [davemers0160 dlib-contrib repository](https://github.com/davemers0160/dlib-contrib )
4. [dlib library](http://dlib.net/ )

## Build

The project uses CMake as the primary mechanism to build the executables.  There are some modifications that may have to be made to the CMakeLists.txt file in order to get the project to build successfully.

The first thing that must be done is to create an environment variable called "PLATFORM".  The CMakeLists.txt file uses this variable to determine where to look for the other required repositories and/or libraries.  These will be machine specific.

To create an environment variable in Windows (drop the -m if you do not have elevated privileges):
```
setx -m PLATFORM MY_PC
```

In Linux (usually placed in .profile or .bashrc):
```
export PLATFORM=MY_PC
```

In the CMakeLists.txt file make sure to add a check for the platform you've added and point to the right locations for the repositories/libraries.

### Windows

From the directory that contains this file, execute the following commands in a Windows command window:

```
mkdir build
mkdir results
cd build
cmake -G "Visual Studio 15 2017 Win64" -T host=x64 ..
cmake --build . --config Release
```

Or you can use the cmake-gui and set the "source code" location to the location of the CmakeLists.txt file and the set the "build" location to the build folder. Then you can open the project in Visual Studio and compile from there.

### Linux

From the directory that contains this file, execute the following commands in a terminal window:

```
mkdir build
mkdir results
cd build
cmake ..
cmake --build . --config Release -- -j4
```

Or you can use the cmake-gui and set the "source code" location to the location of the CmakeLists.txt file and the set the "build" location to the build folder. Then open a terminal window and navigate to the build folder and execute the follokwing command:

```
cmake --build . --config Release -- -j4
```

The -- -j4 tells the make to use 4 cores to build the code.  This number can be set to as many cores as you have on your PC.

## Running

The code is run by supplying all of the parameters in a single file.  Using this method all of the input parameters must be supplied and they must be in the order outlined in the example file *dfd_dnn_input.txt*

To use the file enter the following:

```
Windows: dfd_dnn_an ../dfd_dnn_input.txt
Linux: ./dfd_dnn_an ../dfd_dnn_input.txt
```

It is important to note that if the output folder specified in the input file does not exist the program will run, but there may not be any indication that the data is not being saved.

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


