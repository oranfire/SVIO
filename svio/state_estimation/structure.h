#pragma once 

#include <list>
#include <algorithm>
#include <map>
#include <vector>
#include <numeric>
#include <set>
#include "parameters.h"

using namespace std;

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Cholesky>
using namespace Eigen;

typedef enum Strutural{no, horizontal_1, horizontal_2, atlanta_1, atlanta_2, atlanta_3, atlanta_4, vertical} Strutural;

class Axis
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // id
    int id;

    // visualization
    Vector3d color_13, color_24;

    // geo. info
    Vector3d axis_1, axis_2, axis_3, axis_4;
    double theta;

    // obs. info
    int last_observe_cnt, observe_cnt;

    // state
    bool is_mature;

    // used in opt.
    bool is_opting;

    // Axis(int _id, double _theta):last_observe_cnt(0),observe_cnt(1),id(_id),is_mature(false),theta(_theta){}
    Axis(int _id, double _theta):last_observe_cnt(0),observe_cnt(1),id(_id),is_mature(true),theta(_theta){}
    Axis(){}

    void setcolor()
    {
        switch(id%3)
        {
            case 0: 
                color_13 = Eigen::Vector3d(255, 0, 0);
                color_24 = Eigen::Vector3d(0, 255, 0);
                break;
            case 1: 
                color_13 = Eigen::Vector3d(255, 255, 0);
                color_24 = Eigen::Vector3d(255, 0, 255);
                break;
            case 2: 
                color_13 = Eigen::Vector3d(252, 230, 202);
                color_24 = Eigen::Vector3d(160, 32, 240);
                break;
        }
    }

    void setdir()
    {
        axis_1 = Vector3d(cos(theta), sin(theta), 0);
        axis_2 = Vector3d(-sin(theta), cos(theta), 0);
        axis_3 = Vector3d(-cos(theta), -sin(theta), 0);
        axis_4 = Vector3d(sin(theta), -cos(theta), 0);
    }

    void updateObs()
    {
        if (last_observe_cnt != 0)
        {
            last_observe_cnt = 0;
            observe_cnt++;
        }
    }
};

class AxisMap
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    map<int,Axis> axis_l;
    int axis_id;
    Vector3d h_color;

    AxisMap():axis_id(0),h_color(Vector3d(0, 250, 250)){}

    void updateAxisMap()
    {
        for (auto it = axis_l.begin(); it != axis_l.end();)
        {
            auto& axis = it->second;
            if (S_VERIFY > 0 && axis.observe_cnt < 300 && axis.last_observe_cnt > 300)
            {
                axis_l.erase(it++);
                continue;
            }    

            axis.last_observe_cnt++;
            if (axis.observe_cnt > 100)
                axis.is_mature = true;
            it++;
        }
    }

    void printStatus()
    {
        ROS_WARN("AxisMap: total %d axis detected", axis_l.size());
        for (const auto& axis_pair : axis_l)
        {
            const auto& axis = axis_pair.second;
            ROS_WARN("axis: id(%d), theta(%f), mature(%s), observed(%d), last_observe(%d ago)", 
                axis.id, axis.theta*180/M_PI, (axis.is_mature)?"yes":"no", axis.observe_cnt, axis.last_observe_cnt);
        }
    }
};