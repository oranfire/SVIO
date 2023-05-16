#include "visualization.h"

using namespace ros;
using namespace Eigen;
ros::Publisher pub_odometry, pub_latest_odometry, pub_odometry_posestamped;
ros::Publisher pub_path;
ros::Publisher pub_point_cloud, pub_margin_cloud, pub_plane, pub_line;
ros::Publisher pub_key_poses;
ros::Publisher pub_camera_pose;
ros::Publisher pub_camera_pose_visual;
nav_msgs::Path path;

ros::Publisher pub_keyframe_pose;
ros::Publisher pub_keyframe_point;
ros::Publisher pub_extrinsic;

ros::Publisher pub_point_img;
ros::Publisher pub_plane_img;
ros::Publisher pub_line_img;
ros::Publisher pub_line_cluster_img;

CameraPoseVisualization cameraposevisual(0, 1, 0, 1);
CameraPoseVisualization keyframebasevisual(0.0, 0.0, 1.0, 1.0);
static double sum_of_path = 0;
static Vector3d last_path(0.0, 0.0, 0.0);

void registerPub(ros::NodeHandle &n)
{
    pub_latest_odometry = n.advertise<nav_msgs::Odometry>("imu_propagate", 1000);
    pub_path = n.advertise<nav_msgs::Path>("path", 1000);
    pub_odometry = n.advertise<nav_msgs::Odometry>("odometry", 1000);
    pub_odometry_posestamped = n.advertise<geometry_msgs::PoseStamped>("odometry_pose", 1000);
    pub_point_cloud = n.advertise<sensor_msgs::PointCloud>("point_cloud", 1000);
    pub_plane = n.advertise<sensor_msgs::PointCloud2>("plane", 1000);
    pub_line = n.advertise<visualization_msgs::Marker>("line", 1000);
    pub_margin_cloud = n.advertise<sensor_msgs::PointCloud>("margin_cloud", 1000);
    pub_key_poses = n.advertise<visualization_msgs::Marker>("key_poses", 1000);
    pub_camera_pose = n.advertise<nav_msgs::Odometry>("camera_pose", 1000);
    pub_camera_pose_visual = n.advertise<visualization_msgs::MarkerArray>("camera_pose_visual", 1000);
    pub_keyframe_pose = n.advertise<nav_msgs::Odometry>("keyframe_pose", 1000);
    pub_keyframe_point = n.advertise<sensor_msgs::PointCloud>("keyframe_point", 1000);
    // pub_extrinsic = n.advertise<nav_msgs::Odometry>("extrinsic", 1000);
    pub_point_img = n.advertise<sensor_msgs::Image>("point_image", 1000);
    pub_plane_img = n.advertise<sensor_msgs::Image>("plane_image", 1000);
    pub_line_img = n.advertise<sensor_msgs::Image>("line_image", 1000);
    pub_line_cluster_img = n.advertise<sensor_msgs::Image>("line_cluster_image", 1000);
    cameraposevisual.setScale(1);
    cameraposevisual.setLineWidth(0.05);
    keyframebasevisual.setScale(0.1);
    keyframebasevisual.setLineWidth(0.01);
}


void pubLatestOdometry(const Eigen::Vector3d &P, const Eigen::Quaterniond &Q, const Eigen::Vector3d &V, double t)
{
    nav_msgs::Odometry odometry;
    odometry.header.stamp = ros::Time(t);
    odometry.header.frame_id = "world";
    odometry.pose.pose.position.x = P.x();
    odometry.pose.pose.position.y = P.y();
    odometry.pose.pose.position.z = P.z();
    odometry.pose.pose.orientation.x = Q.x();
    odometry.pose.pose.orientation.y = Q.y();
    odometry.pose.pose.orientation.z = Q.z();
    odometry.pose.pose.orientation.w = Q.w();
    odometry.twist.twist.linear.x = V.x();
    odometry.twist.twist.linear.y = V.y();
    odometry.twist.twist.linear.z = V.z();
    pub_latest_odometry.publish(odometry);
}

void pubPointImage(const cv::Mat &imgTrack, const double t)
{
    std_msgs::Header header;
    header.frame_id = "world";
    header.stamp = ros::Time(t);
    sensor_msgs::ImagePtr imgTrackMsg = cv_bridge::CvImage(header, "bgr8", imgTrack).toImageMsg();
    pub_point_img.publish(imgTrackMsg);
}

void pubPlaneSegmentation(const cv::Mat &imgPlane, const double t)
{
    std_msgs::Header header;
    header.frame_id = "world";
    header.stamp = ros::Time(t);
    sensor_msgs::ImagePtr imgPlaneMsg = cv_bridge::CvImage(header, "bgr8", imgPlane).toImageMsg();
    pub_plane_img.publish(imgPlaneMsg);

    // std::string imgName = "/home/oran/WS/WSa/SLAM/VCU-VIOs/SVIO/output/plane_seg/planeSeg_"+to_string(t)+".jpg";
    // cv::imwrite(imgName, imgPlane);
}

void pubLineImg(const cv::Mat &imgLine, const double t)
{
    std_msgs::Header header;
    header.frame_id = "world";
    header.stamp = ros::Time(t);
    sensor_msgs::ImagePtr imgLineMsg = cv_bridge::CvImage(header, "bgr8", imgLine).toImageMsg();
    pub_line_img.publish(imgLineMsg);

    // std::string imgName = "/home/oran/WS/WSa/SLAM/VCU-VIOs/SVIO/output/line/lineTrack_"+to_string(t)+".jpg";
    // cv::imwrite(imgName, imgLine);
}

void pubLineClusterImg(const cv::Mat &imgLineCluster, const double t)
{
    std_msgs::Header header;
    header.frame_id = "world";
    header.stamp = ros::Time(t);
    sensor_msgs::ImagePtr imgLineClusterMsg = cv_bridge::CvImage(header, "bgr8", imgLineCluster).toImageMsg();
    pub_line_cluster_img.publish(imgLineClusterMsg);

    // std::string imgName = "/home/oran/WS/WSa/SLAM/VCU-VIOs/SVIO/output/line_cluster/lineCluster_"+to_string(t)+".jpg";
    // cv::imwrite(imgName, imgLineCluster);
}

void pubOdometry(const SVIO &estimator, const double t)
{
    if (estimator.solver_flag != SolverFlag::INITIAL)
    {
        std_msgs::Header header;
        header.frame_id = "world";
        header.stamp = ros::Time(t);

        nav_msgs::Odometry odometry;
        odometry.header = header;
        odometry.header.frame_id = "world";
        odometry.child_frame_id = "world";
        Quaterniond tmp_Q;
        tmp_Q = Quaterniond(estimator.Rs[WINDOW_SIZE]);
        odometry.pose.pose.position.x = estimator.Ps[WINDOW_SIZE].x();
        odometry.pose.pose.position.y = estimator.Ps[WINDOW_SIZE].y();
        odometry.pose.pose.position.z = estimator.Ps[WINDOW_SIZE].z();
        odometry.pose.pose.orientation.x = tmp_Q.x();
        odometry.pose.pose.orientation.y = tmp_Q.y();
        odometry.pose.pose.orientation.z = tmp_Q.z();
        odometry.pose.pose.orientation.w = tmp_Q.w();
        odometry.twist.twist.linear.x = estimator.Vs[WINDOW_SIZE].x();
        odometry.twist.twist.linear.y = estimator.Vs[WINDOW_SIZE].y();
        odometry.twist.twist.linear.z = estimator.Vs[WINDOW_SIZE].z();
        pub_odometry.publish(odometry);

        geometry_msgs::PoseStamped pose_stamped;
        pose_stamped.header = header;
        pose_stamped.header.frame_id = "world";
        pose_stamped.pose = odometry.pose.pose;
        path.header = header;
        path.header.frame_id = "world";
        path.poses.push_back(pose_stamped);
        pub_path.publish(path);

        pub_odometry_posestamped.publish(pose_stamped);

        // write result to file
        ofstream foutC(VINS_RESULT_PATH, ios::app);
        if (SAVE_TUM == false)
        {
            foutC.setf(ios::fixed, ios::floatfield);
            foutC.precision(0);
            foutC << header.stamp.toSec() * 1e9 << ",";
            foutC.precision(5);
            foutC << estimator.Ps[WINDOW_SIZE].x() << ","
                << estimator.Ps[WINDOW_SIZE].y() << ","
                << estimator.Ps[WINDOW_SIZE].z() << ","
                << tmp_Q.w() << ","
                << tmp_Q.x() << ","
                << tmp_Q.y() << ","
                << tmp_Q.z() << ","
                << estimator.Vs[WINDOW_SIZE].x() << ","
                << estimator.Vs[WINDOW_SIZE].y() << ","
                << estimator.Vs[WINDOW_SIZE].z() << "," << endl;
        }
        else
        {
            foutC.setf(ios::fixed, ios::floatfield);
            foutC.precision(10);
            foutC << header.stamp.toSec() << " ";
            foutC.precision(5);
            foutC << estimator.Ps[WINDOW_SIZE].x() << " "
                << estimator.Ps[WINDOW_SIZE].y() << " "
                << estimator.Ps[WINDOW_SIZE].z() << " "
                << tmp_Q.x() << " "
                << tmp_Q.y() << " "
                << tmp_Q.z() << " " 
                << tmp_Q.w() << endl;
        }
        foutC.close();
        Eigen::Vector3d tmp_T = estimator.Ps[WINDOW_SIZE];
        printf("time: %f, t: %f %f %f q: %f %f %f %f \n", header.stamp.toSec(), tmp_T.x(), tmp_T.y(), tmp_T.z(),
                                                          tmp_Q.w(), tmp_Q.x(), tmp_Q.y(), tmp_Q.z());
    }
}

void pubOldOdometry(const SVIO &estimator, const double t)
{
    if (estimator.solver_flag != SolverFlag::INITIAL)
    {
        std_msgs::Header header;
        header.frame_id = "world";
        header.stamp = ros::Time(t);

        nav_msgs::Odometry odometry;
        odometry.header = header;
        odometry.header.frame_id = "world";
        odometry.child_frame_id = "world";
        Quaterniond tmp_Q;
        tmp_Q = Quaterniond(estimator.Rs[WINDOW_SIZE]);
        odometry.pose.pose.position.x = estimator.Ps[WINDOW_SIZE].x();
        odometry.pose.pose.position.y = estimator.Ps[WINDOW_SIZE].y();
        odometry.pose.pose.position.z = estimator.Ps[WINDOW_SIZE].z();
        odometry.pose.pose.orientation.x = tmp_Q.x();
        odometry.pose.pose.orientation.y = tmp_Q.y();
        odometry.pose.pose.orientation.z = tmp_Q.z();
        odometry.pose.pose.orientation.w = tmp_Q.w();
        odometry.twist.twist.linear.x = estimator.Vs[WINDOW_SIZE].x();
        odometry.twist.twist.linear.y = estimator.Vs[WINDOW_SIZE].y();
        odometry.twist.twist.linear.z = estimator.Vs[WINDOW_SIZE].z();
        pub_odometry.publish(odometry);

        geometry_msgs::PoseStamped pose_stamped;
        pose_stamped.header = header;
        pose_stamped.header.frame_id = "world";
        pose_stamped.pose = odometry.pose.pose;
        path.header = header;
        path.header.frame_id = "world";
        path.poses.push_back(pose_stamped);
        pub_path.publish(path);

        // write result to file
        if (estimator.marginalization_flag == MARGIN_OLD)
        {
            Quaterniond tmp_Q0 = Quaterniond(estimator.Rs[0]);

            ofstream foutC(VINS_RESULT_PATH, ios::app);
            foutC.setf(ios::fixed, ios::floatfield);
            foutC.precision(0);
            foutC << estimator.Headers[0] * 1e9 << ",";
            foutC.precision(5);
            foutC << estimator.Ps[0].x() << ","
                << estimator.Ps[0].y() << ","
                << estimator.Ps[0].z() << ","
                << tmp_Q0.w() << ","
                << tmp_Q0.x() << ","
                << tmp_Q0.y() << ","
                << tmp_Q0.z() << ","
                << estimator.Vs[0].x() << ","
                << estimator.Vs[0].y() << ","
                << estimator.Vs[0].z() << "," << endl;
            foutC.close();
        }
        
        Eigen::Vector3d tmp_T = estimator.Ps[WINDOW_SIZE];
        printf("time: %f, t: %f %f %f q: %f %f %f %f \n", header.stamp.toSec(), tmp_T.x(), tmp_T.y(), tmp_T.z(),
                                                          tmp_Q.w(), tmp_Q.x(), tmp_Q.y(), tmp_Q.z());
    }
}

void pubKeyPoses(const SVIO &estimator, const double t)
{
    if (estimator.key_poses.size() == 0)
        return;
    
    std_msgs::Header header;
    header.frame_id = "world";
    header.stamp = ros::Time(t);

    visualization_msgs::Marker key_poses;
    key_poses.header = header;
    key_poses.header.frame_id = "world";
    key_poses.ns = "key_poses";
    key_poses.type = visualization_msgs::Marker::SPHERE_LIST;
    key_poses.action = visualization_msgs::Marker::ADD;
    key_poses.pose.orientation.w = 1.0;
    key_poses.lifetime = ros::Duration();

    //static int key_poses_id = 0;
    key_poses.id = 0; //key_poses_id++;
    key_poses.scale.x = 0.05;
    key_poses.scale.y = 0.05;
    key_poses.scale.z = 0.05;
    key_poses.color.r = 1.0;
    key_poses.color.a = 1.0;

    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
        geometry_msgs::Point pose_marker;
        Vector3d correct_pose;
        correct_pose = estimator.key_poses[i];
        pose_marker.x = correct_pose.x();
        pose_marker.y = correct_pose.y();
        pose_marker.z = correct_pose.z();
        key_poses.points.push_back(pose_marker);
    }
    pub_key_poses.publish(key_poses);
}

void pubCameraPose(const SVIO &estimator, const double t)
{
    int idx2 = WINDOW_SIZE - 1;

    if (estimator.solver_flag != SolverFlag::INITIAL)
    {
        int i = idx2;
        Vector3d P = estimator.Ps[i] + estimator.Rs[i] * estimator.tic[0];
        Quaterniond R = Quaterniond(estimator.Rs[i] * estimator.ric[0]);

        std_msgs::Header header;
        header.frame_id = "world";
        header.stamp = ros::Time(t);

        nav_msgs::Odometry odometry;
        odometry.header = header;
        odometry.header.frame_id = "world";
        odometry.pose.pose.position.x = P.x();
        odometry.pose.pose.position.y = P.y();
        odometry.pose.pose.position.z = P.z();
        odometry.pose.pose.orientation.x = R.x();
        odometry.pose.pose.orientation.y = R.y();
        odometry.pose.pose.orientation.z = R.z();
        odometry.pose.pose.orientation.w = R.w();

        pub_camera_pose.publish(odometry);

        cameraposevisual.reset();
        cameraposevisual.add_pose(P, R);
        cameraposevisual.publish_by(pub_camera_pose_visual, odometry.header);
    }
}


void pubPointCloud(const SVIO &estimator, const double t)
{
    std_msgs::Header header;
    header.frame_id = "world";
    header.stamp = ros::Time(t);

    sensor_msgs::PointCloud point_cloud;
    point_cloud.header = header;

    for (auto &it_per_id : estimator.f_manager.feature)
    {
        int used_num;
        used_num = it_per_id.feature_per_frame.size();
        if (!(used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2 && it_per_id.estimated_depth > 0))
            continue;
        if (it_per_id.start_frame > WINDOW_SIZE * 3.0 / 4.0 || it_per_id.solve_flag != 1)
            continue;
        int imu_i = it_per_id.start_frame + it_per_id.depth_shift;
        Vector3d pts_i = it_per_id.feature_per_frame[it_per_id.depth_shift].pt * it_per_id.estimated_depth;
        Vector3d w_pts_i = estimator.Rs[imu_i] * (estimator.ric[0] * pts_i + estimator.tic[0]) + estimator.Ps[imu_i];

        geometry_msgs::Point32 p;
        p.x = w_pts_i(0);
        p.y = w_pts_i(1);
        p.z = w_pts_i(2);
        point_cloud.points.push_back(p);
    }
    pub_point_cloud.publish(point_cloud);

    // pub margined potin
    sensor_msgs::PointCloud margin_cloud;
    margin_cloud.header = header;

    for (auto &it_per_id : estimator.f_manager.feature)
    { 
        int used_num;
        used_num = it_per_id.feature_per_frame.size();
        if (!(used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;
        //if (it_per_id->start_frame > WINDOW_SIZE * 3.0 / 4.0 || it_per_id->solve_flag != 1)
        //        continue;

        if (it_per_id.start_frame == 0 && it_per_id.feature_per_frame.size() <= 2 
            && it_per_id.solve_flag == 1 )
        {
            int imu_i = it_per_id.start_frame + it_per_id.depth_shift;
            Vector3d pts_i = it_per_id.feature_per_frame[it_per_id.depth_shift].pt * it_per_id.estimated_depth;
            Vector3d w_pts_i = estimator.Rs[imu_i] * (estimator.ric[0] * pts_i + estimator.tic[0]) + estimator.Ps[imu_i];

            geometry_msgs::Point32 p;
            p.x = w_pts_i(0);
            p.y = w_pts_i(1);
            p.z = w_pts_i(2);
            margin_cloud.points.push_back(p);
        }
    }
    pub_margin_cloud.publish(margin_cloud);
}

void pubPlanePointCloud(const SVIO &estimator, const double t, bool pub_local, bool pub_history)
{
    sensor_msgs::PointCloud2 output_msg;
    pcl::PointCloud<pcl::PointXYZRGB> cloud;

    Vector3d c;
    if (pub_local == true)
    {
        for (const auto& it : estimator.p_manager.plane_l)
        {
            PlanePerId* it_per_id = &(it.second);

            if (it_per_id->plane_ml.size() < 1)
                continue;

            c = it_per_id->color;

            int f_id = it_per_id->start_frame_id;
            int index = 0;
            Matrix3d R0 = estimator.Rs[f_id] * estimator.ric[0];
            Vector3d P0 = estimator.Rs[f_id] * estimator.tic[0] + estimator.Ps[f_id];
            Matrix4d transT0 = Matrix4d::Identity();
            transT0.topLeftCorner<3,3>() = R0.transpose();
            transT0.bottomLeftCorner<1,3>() = P0.transpose();
            Vector4d hesse_in_s = transT0 * it_per_id->hesse;
            for (auto it_per_xyz : it_per_id->plane_ml[0].xyz) // only publish plane points on the oldest frame
            {
                index++;
                if (index%10 != 0)
                    continue;

                Eigen::Vector3d pts_in_world;
                #ifdef VISIUALIZE_RAW_POINT_CLOUD
                {
                    pts_in_world = R0 * it_per_xyz + P0;
                }
                #else
                {
                    Eigen::Vector3d pts_on_plane = it_per_xyz;
                    pts_on_plane(2) = -(hesse_in_s.head<2>().transpose()*it_per_xyz.head<2>()+hesse_in_s(3))/hesse_in_s(2);
                    pts_in_world = R0 * pts_on_plane + P0;    
                }
                #endif

                pcl::PointXYZRGB p;
                p.x = pts_in_world(0);
                p.y = pts_in_world(1);
                p.z = pts_in_world(2);
                p.r = c(0);
                p.g = c(1);
                p.b = c(2);
                cloud.points.push_back(p);
            }

            // ROS_WARN("PUBLISH PLANE %d with %d points", it_per_id->plane_id, cloud.points.size());
        }
    }

    std_msgs::Header header;
    header.frame_id = "world";
    header.stamp = ros::Time(t);

    pcl::toROSMsg(cloud, output_msg);
    output_msg.header = header;
    pub_plane.publish(output_msg);
}

void pubLineCloud(const SVIO &estimator, const double t)
{
    std_msgs::Header header;
    header.frame_id = "world";
    header.stamp = ros::Time(t);

    visualization_msgs::Marker line_cloud;
    line_cloud.header = header;
    line_cloud.ns = "line";
    line_cloud.id = 0;
    line_cloud.scale.x = 0.02;
    line_cloud.type = visualization_msgs::Marker::LINE_LIST;
    line_cloud.action = visualization_msgs::Marker::ADD;
    line_cloud.color.r = 1.0f;
    line_cloud.color.a = 1.0f;
    line_cloud.color.b = line_cloud.color.g = 1.f;

    int cnt = 0;
    for (const auto& it : estimator.l_manager.line_l)
    {
        if (it.second.has_init == false || it.second.line_ml.size() < 2 || it.second.line_ml[0].end_seg_pts.size() == 0)
            continue;
        
        Vector3d nw = it.second.n*it.second.d, vw = it.second.v;
        Matrix4d Lw;
        Lw << Utility::skewSymmetric(nw), vw, -vw.transpose(), 0;  //和普吕克矩阵有微妙差异 

        Vector3d s0 = it.second.line_ml[0].s, e0 = it.second.line_ml[0].e;
        Vector2d delta = s0.cross(e0).head<2>().normalized();
        Vector3d s0_step(s0(0)+delta(0),s0(1)+delta(1),s0(2)), e0_step(e0(0)+delta(0),e0(1)+delta(1),e0(2));
        Vector4d hesse_sc = Vector4d::Zero(), hesse_ec = Vector4d::Zero();
        hesse_sc.head<3>() = s0.cross(s0_step).normalized();
        hesse_ec.head<3>() = e0.cross(e0_step).normalized();

        Vector3d t0 = estimator.Ps[it.second.start_frame_id] + estimator.Rs[it.second.start_frame_id] * estimator.tic[0];
        Matrix3d R0 = estimator.Rs[it.second.start_frame_id] * estimator.ric[0];
        Matrix4d invTransT0 = Matrix4d::Identity();
        invTransT0.block<3,3>(0,0) = R0;
        invTransT0.block<1,3>(3,0) = -t0.transpose()*R0;

        Vector4d hesse_sw = invTransT0*hesse_sc, hesse_ew = invTransT0*hesse_ec;
        Vector4d p0 = Lw*hesse_sw, p1 = Lw*hesse_ew;
        p0 /= p0(3), p1 /= p1(3);

        Vector3d p0_c = R0.transpose() * (p0.head<3>() - t0), p1_c = R0.transpose() * (p1.head<3>() - t0);
        if (p0_c(2) < 0 || p1_c(2) < 0)
            continue;

        geometry_msgs::Point gsw, gew;
        gsw.x = p0[0];
        gsw.y = p0[1];
        gsw.z = p0[2];
        gew.x = p1[0];
        gew.y = p1[1];
        gew.z = p1[2];
        line_cloud.points.push_back(gsw);
        line_cloud.points.push_back(gew);

        cnt++;
    }
    ROS_WARN("publish %d lines", cnt);

    pub_line.publish(line_cloud);
}

void pubTF(const SVIO &estimator, const double t)
{
    if( estimator.solver_flag == SolverFlag::INITIAL)
        return;

    std_msgs::Header header;
    header.frame_id = "world";
    header.stamp = ros::Time(t);

    static tf::TransformBroadcaster br;
    tf::Transform transform;
    tf::Quaternion q;
    // body frame
    Vector3d correct_t;
    Quaterniond correct_q;
    correct_t = estimator.Ps[WINDOW_SIZE];
    correct_q = estimator.Rs[WINDOW_SIZE];

    transform.setOrigin(tf::Vector3(correct_t(0),
                                    correct_t(1),
                                    correct_t(2)));
    q.setW(correct_q.w());
    q.setX(correct_q.x());
    q.setY(correct_q.y());
    q.setZ(correct_q.z());
    transform.setRotation(q);
    br.sendTransform(tf::StampedTransform(transform, header.stamp, "world", "body"));

    // camera frame
    transform.setOrigin(tf::Vector3(estimator.tic[0].x(),
                                    estimator.tic[0].y(),
                                    estimator.tic[0].z()));
    q.setW(Quaterniond(estimator.ric[0]).w());
    q.setX(Quaterniond(estimator.ric[0]).x());
    q.setY(Quaterniond(estimator.ric[0]).y());
    q.setZ(Quaterniond(estimator.ric[0]).z());
    transform.setRotation(q);
    br.sendTransform(tf::StampedTransform(transform, header.stamp, "body", "camera"));

    
    // nav_msgs::Odometry odometry;
    // odometry.header = header;
    // odometry.header.frame_id = "world";
    // odometry.pose.pose.position.x = estimator.tic[0].x();
    // odometry.pose.pose.position.y = estimator.tic[0].y();
    // odometry.pose.pose.position.z = estimator.tic[0].z();
    // Quaterniond tmp_q{estimator.ric[0]};
    // odometry.pose.pose.orientation.x = tmp_q.x();
    // odometry.pose.pose.orientation.y = tmp_q.y();
    // odometry.pose.pose.orientation.z = tmp_q.z();
    // odometry.pose.pose.orientation.w = tmp_q.w();
    // pub_extrinsic.publish(odometry);

}

void pubKeyframe(const SVIO &estimator, const double t)
{
    // pub camera pose, 2D-3D points of keyframe
    if (estimator.solver_flag != SolverFlag::INITIAL && estimator.marginalization_flag == 0)
    {
        int i = WINDOW_SIZE - 2;
        //Vector3d P = estimator.Ps[i] + estimator.Rs[i] * estimator.tic[0];
        Vector3d P = estimator.Ps[i];
        Quaterniond R = Quaterniond(estimator.Rs[i]);

        std_msgs::Header header;
        header.frame_id = "world";
        header.stamp = ros::Time(t);

        nav_msgs::Odometry odometry;
        odometry.header = header; // estimator.Headers[WINDOW_SIZE - 2];
        odometry.header.frame_id = "world";
        odometry.pose.pose.position.x = P.x();
        odometry.pose.pose.position.y = P.y();
        odometry.pose.pose.position.z = P.z();
        odometry.pose.pose.orientation.x = R.x();
        odometry.pose.pose.orientation.y = R.y();
        odometry.pose.pose.orientation.z = R.z();
        odometry.pose.pose.orientation.w = R.w();
        //printf("time: %f t: %f %f %f r: %f %f %f %f\n", odometry.header.stamp.toSec(), P.x(), P.y(), P.z(), R.w(), R.x(), R.y(), R.z());

        pub_keyframe_pose.publish(odometry);

        sensor_msgs::PointCloud point_cloud;
        point_cloud.header = header;
	    point_cloud.header.stamp.fromSec(estimator.Headers[WINDOW_SIZE - 2]);
        for (auto &it_per_id : estimator.f_manager.feature)
        {
            int frame_size = it_per_id.feature_per_frame.size();
            if(it_per_id.start_frame < WINDOW_SIZE - 2 && it_per_id.start_frame + frame_size - 1 >= WINDOW_SIZE - 2 && it_per_id.solve_flag == 1)
            {

                int imu_i = it_per_id.start_frame;
                Vector3d pts_i = it_per_id.feature_per_frame[0].pt * it_per_id.estimated_depth;
                Vector3d w_pts_i = estimator.Rs[imu_i] * (estimator.ric[0] * pts_i + estimator.tic[0])
                                      + estimator.Ps[imu_i];
                geometry_msgs::Point32 p;
                p.x = w_pts_i(0);
                p.y = w_pts_i(1);
                p.z = w_pts_i(2);
                point_cloud.points.push_back(p);

                int imu_j = WINDOW_SIZE - 2 - it_per_id.start_frame;
                sensor_msgs::ChannelFloat32 p_2d;
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].pt.x());
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].pt.y());
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].pt_2d.x());
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].pt_2d.y());
                p_2d.values.push_back(it_per_id.feature_id);
                point_cloud.channels.push_back(p_2d);
            }

        }
        pub_keyframe_point.publish(point_cloud);
    }
}