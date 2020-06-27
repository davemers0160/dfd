
Common support files for the various DfD projects.

## Repository Breakdown


### include
This is where all of the C++ common header files are located

### inputs
This is where the input files are located for various data sets.  The Middlebury College dataset (mb) is taken from [here](https://vision.middlebury.edu/stereo/data/) and is synthetically blurred using the code in the [blur generation](https://github.com/davemers0160/blur_generation) repository.

The Real World dataset (rw) is a custom dataset taken with a FLIR Chameleon3 camera and a Varioptic microfluidic lens.  The groubnd truth data was collected using an Ouster OS1-64 LIDAR.

### lens_data
This folder contains the blur characterization of the microfluidic lens and camera.

### matlab
This folder contains MATLAB code to perform data analysis, preprocess data and other various tasks.

### nets
This folder conatins the final weights of the DfD-Net trained on the synthetically blurred Middlebury College dataset.

### python
This folder contains general python code that will display an interactive calculator to determine the amount of blur on a knife edge target for the given lens and camera parameters.

### scripts
This folder contains scripts for Photoshop.
