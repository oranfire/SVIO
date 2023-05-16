# SVIO
## A RGB-D Visual-Inertial Odometry Leveraging Structural Regularity 

**Related Paper**

Pengfei Gu and Ziyang Meng, **"S-VIO: Exploiting Structural Constraints for RGB-D Visual Inertial Odometry"**, *IEEE Robotics and Automation Letters*, 2023. [pdf](https://ieeexplore.ieee.org/document/10107752)

## 1. Prerequisites
### 1.1 **Ubuntu** and **ROS**
Test on Ubuntu 18.04 and ROS Melodic. Follow [ROS Installation](http://wiki.ros.org/ROS/Installation). The following ROS pacakges are needed:
```
    sudo apt-get install ros-YOUR_DISTRO-cv-bridge ros-YOUR_DISTRO-tf ros-YOUR_DISTRO-message-filters ros-YOUR_DISTRO-image-transport
```

### 1.2 **Ceres Solver**
Follow [Ceres Installation](http://ceres-solver.org/installation.html).

## 2. Build SVIO on ROS
Clone the repository and catkin_make:
```
    cd YOUR_PATH_TO_SVIO/SVIO/src
    git clone https://github.com/oranfire/SVIO.git
    cd ../
    catkin_make
```
Before running the SVIO, remember to source the script first:
```
    source YOUR_PATH_TO_SVIO/SVIO/devel/setup.bash
```

## 3. Usage
### 3.1 VCU-RVI handheld dataset 
Download [VCU-RVI Dataset](https://github.com/rising-turtle/VCU_RVI_Benchmark). We use the time-synchronized IMU measurements and RGB-D images for localization. Open two terminals, launch the estimator and play the bag respectively. Take motion_1 for example: 
```
    roslaunch YOUR_PATH_TO_SVIO/SVIO/src/launch/svio_vcu_run.launch
    rosbag play YOUR_PATH_TO_BAG/lab-motion1.bag
```
It will automounsly launch the rviz for visualization.

### 3.2 OpenLORIS-Scene wheeled robot dataset
Download [OpenLORIS-Scene Dataset](https://lifelong-robotic-vision.github.io/dataset/scene). Although it contains multiple sensors, we only use the RGB-D images and IMU measurements from the d400 depth camera. Open two terminals, launch the estimator and play the bag respectively. Take home1-1 for example: 
```
    roslaunch YOUR_PATH_TO_SVIO/SVIO/src/launch/svio_openloris_run.launch
    rosbag play YOUR_PATH_TO_BAG/home1-1.bag
```
Note that the bag files provided by the OpenLORIS-Scene dataset record the acceleration and the angular velocity in two topics. Therefore before playing the bag, we should merge these two topics into one topic by using the script provided by [OpenLORIS-Scene Tools](https://github.com/lifelong-robotic-vision/openloris-scene-tools) first:
```
    python merge_imu_topics.py YOUR_PATH_TO_BAG/home1-1.bag
```  

### 3.3 EuRoC MAV dataset
As a bonus, SVIO also provides interface for a stereo-inertial dataset, [EuRoC Dataset](https://projects.asl.ethz.ch/datasets/doku.php?id=kmavvisualinertialdatasets), where the depth image is computed from the rectified stereo images using the SGM algorithm. Take mh01 for example:
```
    roslaunch YOUR_PATH_TO_SVIO/SVIO/src/launch/svio_euroc_run.launch
    rosbag play YOUR_PATH_TO_BAG/mh01.bag
```

## 4. Credits
Many thanks to the authors of [VINS-Mono](https://github.com/HKUST-Aerial-Robotics/VINS-Mono), [DUI-VIO](https://github.com/rising-turtle/DUI_VIO) and [ManhattanSLAM](https://github.com/razayunus/ManhattanSLAM). Our system is built on the first two projects, and part of codes are borrowed from the third project.

## 5. Citation
```
    @ARTICLE{10107752,
      author={Gu, Pengfei and Meng, Ziyang},
      journal={IEEE Robotics and Automation Letters}, 
      title={S-VIO: Exploiting Structural Constraints for RGB-D Visual Inertial Odometry}, 
      year={2023},
      volume={8},
      number={6},
      pages={3542-3549},
      doi={10.1109/LRA.2023.3270033}}
```

## 6. Licence
The source code is released under [GPLv3](http://www.gnu.org/licenses/) license.