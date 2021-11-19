# ----------------------------------------------------------------------------
#  Root CMake file for DfD_DNN Project
# ----------------------------------------------------------------------------
#
# Example CMake build line for windows & VS2015
# windows w/cuda: cmake -G "Visual Studio 14 2015 Win64" -DUSE_AVX_INSTRUCTIONS=ON ..
# linux w/o cuda: cmake -DUSE_AVX_INSTRUCTIONS=ON -DDLIB_USE_CUDA=OFF ..
#
#
# To compile:
# windows: cmake --build . --config Release
# linux: cmake --build . --config Release -- -j10
#


cmake_minimum_required(VERSION 2.8.12)

PROJECT(dfd_dnn_project)
set(exe_name "dfd_dnn")

# This looks for an environment variable called PLATFORM.  This can be used to point
# cmake to various locations based on the individual machines that the code is complied on
set(PLATFORM "$ENV{PLATFORM}" CACHE INTERNAL "Copied from environment variable")
if("${PLATFORM}" STREQUAL "")
    set(PLATFORM "NOT_SET")
    message("The PLATFORM Environment variable is not set on this system.  Setting it to:")
endif()
message("Platform: " ${PLATFORM})

# find the root of the main project
string(FIND ${PROJECT_SOURCE_DIR} "/" INDEX REVERSE)
string(SUBSTRING ${PROJECT_SOURCE_DIR} 0 ${INDEX} PR_ROOT)

# Build Options
option(USE_AVX_INSTRUCTIONS  "Compile your program with AVX instructions"  ON)

if(${PLATFORM} STREQUAL "Laptop")
    set(COMMON_ROOT "D:/Common")
    set(DLIB_CONTRIB "D:/Projects/dlib-contrib")
    set(VS_GEN_ROOT "D:/Projects/vs_gen")
    set(VS_GEN_LIB "${VS_GEN_ROOT}/build/Release/vs_gen.lib")
    include(D:/dlib-19.22/dlib/cmake)
elseif(${PLATFORM} STREQUAL "MainGear")
    set(_ECLIPSE_VERSION "3.8")
    set(COMMON_ROOT "~/Projects/Common")
    set(DLIB_CONTRIB "~/Projects/dlib-contrib")
    include(~/dlib-19.22/dlib/cmake)
elseif(${PLATFORM} STREQUAL "HPC")
    set(COMMON_ROOT "~/Projects/Common")
    set(DLIB_CONTRIB "~/Projects/dlib-contrib")
    include(~/dlib-19.22/dlib/cmake)
elseif(${PLATFORM} STREQUAL "BEAST")
    set(COMMON_ROOT "E:/Projects/Common")
    set(DLIB_CONTRIB "E:/Projects/dlib-contrib")
    set(VS_GEN_ROOT "E:/Projects/vs_gen")
    set(VS_GEN_LIB "${VS_GEN_ROOT}/build/Release/vs_gen.lib")
    include(E:/dlib-19.22/dlib/cmake)
else()
    set(COMMON_ROOT "D:/Projects/Common")
    set(DLIB_CONTRIB "D:/Projects/dlib-contrib")
    include(D:/dlib-19.22/dlib/cmake)
endif()


# grab all of the *.c and *.cpp files in the src directory
file(GLOB SOURCES "src/*.c*"	
    "${COMMON_ROOT}/src/*.cpp"
    "${DLIB_CONTRIB}/src/*.cpp"
    "${PR_ROOT}/src/*.cpp"
    )            
            
set(project_includes "include" 
    "${PR_ROOT}/common/include"
    "${COMMON_ROOT}/include" 
    "${DLIB_CONTRIB}/include" 
    "${VS_GEN_ROOT}/include"
    )
            
message(" ")
message("--------------------------------------------------------")
message("Source Files: ")
foreach(S ${SOURCES})
    message("  " ${S})
endforeach()
message("--------------------------------------------------------")

message(" ")

message("--------------------------------------------------------")
message("Project Includes: ")
foreach(S ${project_includes})
    message("  " ${S})
endforeach()
message("--------------------------------------------------------")
message(" ")

message("Cuda version: " ${CUDA_VERSION})

# this is for the various cmake versions out there
#if(${CMAKE_VERSION} VERSION_GREATER "3.12.4")
#    add_compile_definitions(CUDA_VERSION=${CUDA_VERSION})
#else()
#    add_definitions(CUDA_VERSION=${CUDA_VERSION})
#endif()

include_directories(${project_includes})
add_executable(${exe_name} ${SOURCES})
target_link_libraries(${exe_name} dlib::dlib ${VS_GEN_LIB})

if(WIN32)
    message("MSVC Version: " ${MSVC_VERSION})
    # uncomment for MSVC and cmake greater than 3.6.0
    if(MSVC AND (${CMAKE_VERSION} VERSION_GREATER "3.6.0"))
        target_compile_options(${exe_name} PRIVATE "/MP")
        set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${exe_name})
    endif()
endif()

message(" ")
message(STATUS "Build Complete!")
message(" ")
