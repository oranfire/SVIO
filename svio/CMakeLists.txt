project(svio)

cmake_minimum_required(VERSION 2.8.3)

# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -std=c++11 -w -g -O3 -fpermissive -mno-avx -mno-avx2")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -w -g -O3 -fpermissive")
SET(CMAKE_BUILD_TYPE "RelWithDebInfo")

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# Set the build type.  Options are:
#  Coverage       : w/ debug symbols, w/o optimization, w/ code-coverage
#  Debug          : w/ debug symbols, w/o optimization
#  Release        : w/o debug symbols, w/ optimization
#  RelWithDebInfo : w/ debug symbols, w/ optimization
#  MinSizeRel     : w/o debug symbols, w/ optimization, stripped binaries
#set(ROS_BUILD_TYPE RelWithDebInfo)

# rosbuild_init()

option(WITH_OPENGV "use opengv" OFF) #ON use opengv or not in hybridPnP 

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

# set(OpenCV_DIR "/home/davidz/work/3rdlibs/OpenCV-2.4.10/build")
# set(OpenCV_DIR "/home/oran/WS/WSa/KU/opencv-3.2.0/opencv/build")
find_package(OpenCV REQUIRED)
# message("OpenCV_INCLUDE_DIR: ${OpenCV_INCLUDE_DIRS}")
# include_directories("/home/hzhang8/work/3rdlib/opencv-3.3.1/install/include")

find_package(catkin REQUIRED COMPONENTS
    roscpp
    cv_bridge
    sensor_msgs
    std_msgs
    pcl_conversions
    pcl_ros
    message_filters
    image_transport
    camera_models_in # copied from VINS-Fusion 
    tf
)

catkin_package(
    # DEPENDS Eigen
    # CATKIN_DEPENDS roscpp std_msgs sensor_msgs cv_bridge
    # INCLUDE_DIRS include
    # LIBRARIES
)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(Eigen3 QUIET)
include_directories(${EIGEN3_INCLUDE_DIR} "/usr/include/eigen3")
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories(${catkin_INCLUDE_DIRS})

find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})


if(WITH_OPENGV)
    add_definitions(-DUSE_OPENGV)
    #### opengv ####    
    set(opengv_DIR "/home/davidz/work/3rdlib/opengv/build")
    MESSAGE(STATUS "USE_OPENGV ")
    find_package(opengv REQUIRED
        PATHS ${opengv_DIR})
    include_directories("${opengv_DIR}/../include")
endif()

# svio library 
add_library(svio state_estimation/svio.cpp state_estimation/svio_init.cpp
            state_estimation/feature_manager.cpp state_estimation/marginalization_factor.cpp state_estimation/depth_factor.cpp
            state_estimation/projection_quat.cpp state_estimation/parameters.cpp state_estimation/plane_factor.cpp
            utility/visualization.cpp utility/utility.cpp utility/CameraPoseVisualization.cpp
            initialization/initial_aligment.cpp initialization/initial_ex_rotation.cpp 
            initialization/initial_sfm.cpp initialization/solve_5pts.cpp
            initialization/solve_opt.cpp initialization/translate_factor.cpp
            initialization/gmm_model.cpp feature_tracking/feature_tracker.cpp
            plane_tools/PlaneExtractor.cpp state_estimation/plane_manager.cpp state_estimation/line_parameterization.cpp
            line_tracking/line_tracker.cpp state_estimation/line_manager.cpp state_estimation/line_projection_factor.cpp
            state_estimation/vanish_point_factor.cpp state_estimation/mns_alg.cpp state_estimation/frontend.cpp)

if(WITH_OPENGV)
    target_link_libraries(svio ${CERES_LIBRARIES} ${catkin_LIBRARIES} ${OpenCV_LIBS} opengv)
else()
    target_link_libraries(svio ${CERES_LIBRARIES} ${catkin_LIBRARIES} ${OpenCV_LIBS})
endif()

add_executable(svio_syn_node state_estimation/svio_syn_node.cpp)
target_link_libraries(svio_syn_node svio)

add_executable(svio_stereo_node state_estimation/svio_stereo_node.cpp)
target_link_libraries(svio_stereo_node svio)


