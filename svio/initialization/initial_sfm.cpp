#include "initial_sfm.h"
#include "../state_estimation/feature_manager.h"

GlobalSFM::GlobalSFM(){}

void GlobalSFM::triangulatePoint(Eigen::Matrix<double, 3, 4> &Pose0, Eigen::Matrix<double, 3, 4> &Pose1,
						Vector2d &point0, Vector2d &point1, Vector3d &point_3d)
{
	Matrix4d design_matrix = Matrix4d::Zero();
	design_matrix.row(0) = point0[0] * Pose0.row(2) - Pose0.row(0);
	design_matrix.row(1) = point0[1] * Pose0.row(2) - Pose0.row(1);
	design_matrix.row(2) = point1[0] * Pose1.row(2) - Pose1.row(0);
	design_matrix.row(3) = point1[1] * Pose1.row(2) - Pose1.row(1);
	Vector4d triangulated_point;
	triangulated_point =
		      design_matrix.jacobiSvd(Eigen::ComputeFullV).matrixV().rightCols<1>();
	point_3d(0) = triangulated_point(0) / triangulated_point(3);
	point_3d(1) = triangulated_point(1) / triangulated_point(3);
	point_3d(2) = triangulated_point(2) / triangulated_point(3);
}


bool GlobalSFM::solveFrameByPnP(Matrix3d &R_initial, Vector3d &P_initial, int i,
								vector<SFMFeature> &sfm_f)
{
	vector<cv::Point2f> pts_2_vector;
	vector<cv::Point3f> pts_3_vector;
	for (int j = 0; j < feature_num; j++)
	{
		if (sfm_f[j].state != true)
			continue;
		Vector2d point2d;
		for (int k = 0; k < (int)sfm_f[j].observation.size(); k++)
		{
			if (sfm_f[j].observation[k].first == i)
			{
				Vector2d img_pts = sfm_f[j].observation[k].second;
				cv::Point2f pts_2(img_pts(0), img_pts(1));
				pts_2_vector.push_back(pts_2);
				cv::Point3f pts_3(sfm_f[j].position[0], sfm_f[j].position[1], sfm_f[j].position[2]);
				pts_3_vector.push_back(pts_3);
				break;
			}
		}
	}
	printf("at frame %d has observations: %d\n", i, pts_2_vector.size());
	if (int(pts_2_vector.size()) < 15)
	{
		printf("initial_sfm.cpp: for frame %d, only %d features are found \
			unstable features tracking, please slowly move you device!\n", i, pts_2_vector.size());
		if (int(pts_2_vector.size()) < 10){
			printf("initial_sfm.cpp: return false since only %d features are found", pts_2_vector.size());
			return false;
		}
	}
	cv::Mat r, rvec, t, D, tmp_r;
	cv::eigen2cv(R_initial, tmp_r);
	cv::Rodrigues(tmp_r, rvec);
	cv::eigen2cv(P_initial, t);
	cv::Mat K = (cv::Mat_<double>(3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);
	bool pnp_succ;
	pnp_succ = cv::solvePnP(pts_3_vector, pts_2_vector, K, D, rvec, t, 1);
	if(!pnp_succ)
	{
		return false;
	}
	cv::Rodrigues(rvec, r);
	//cout << "r " << endl << r << endl;
	MatrixXd R_pnp;
	cv::cv2eigen(r, R_pnp);
	MatrixXd T_pnp;
	cv::cv2eigen(t, T_pnp);
	R_initial = R_pnp;
	P_initial = T_pnp;
	return true;

}

void GlobalSFM::triangulateTwoFramesWithDepthCov(int frame0,  Matrix<double, 3, 4 > & Pose0, int frame1, 
								Eigen::Matrix<double, 3, 4> & Pose1, vector<SFMFeature> & sfm_f)
{
	assert(frame0 != frame1);
	Matrix3d Pose0_R = Pose0.block< 3,3 >(0,0);
	Matrix3d Pose1_R = Pose1.block< 3,3 >(0,0);
	Vector3d Pose0_t = Pose0.block< 3,1 >(0,3);
	Vector3d Pose1_t = Pose1.block< 3,1 >(0,3);
	int cnt_valid = 0; 
	for (int j = 0; j < feature_num; j++)
	{
		if (sfm_f[j].state == true){
			cnt_valid++;
			continue;
		}
		bool has_0 = false, has_1 = false;
		Vector3d point0;
		Vector3d point1;
		double std_z0 = 0; 
		double std_z1 = 0; 
		for (int k = 0; k < (int)sfm_f[j].observation.size(); k++)
		{
			if (sfm_f[j].observation[k].first == frame0)
			{
				point0 = Vector3d(sfm_f[j].observation[k].second.x(),sfm_f[j].observation[k].second.y(), sfm_f[j].observation_depth[k].second);	
				has_0 = (point0.z() > 0.1); // true; 
				std_z0 = sfm_f[j].observation_depth_std[k].second; 
			}
			if (sfm_f[j].observation[k].first == frame1)
			{
				point1 = Vector3d(sfm_f[j].observation[k].second.x(),sfm_f[j].observation[k].second.y(), sfm_f[j].observation_depth[k].second);
				has_1 = (point1.z() > 0.1) ; // true;
				std_z1 = sfm_f[j].observation_depth_std[k].second; 
			}
			if(has_0 && has_1) break;
		}

		if (has_0 && has_1)
		{
			if(std_z0 > 0 && std_z1 > 0){ // Chi2 check 

				Matrix3d C_p0, C_p1, C_p0_p1; 
				sigma_pt3d(C_p0, point0.x(), point0.y(), point0.z(), std_z0); 
				sigma_pt3d(C_p1, point1.x(), point1.y(), point1.z(), std_z1); 

				Matrix3d R10; 
				Vector3d point_3d, point1_reprojected;
				point0.x() = point0.x()*point0.z(); 
				point0.y() = point0.y()*point0.z(); 
				point_3d = Pose0_R.transpose()*point0 - Pose0_R.transpose()*Pose0_t;//shan add:this is point in world;
				point1_reprojected = Pose1_R*point_3d+Pose1_t;

				// approximate pose uncertainty 
				Matrix3d C_pose = 0.01*Matrix3d::Identity(); 
				
				R10 = Pose1_R*Pose0_R.transpose(); 
				C_p0_p1 = R10 * C_p0 * R10.transpose() + C_pose; // covariance forward 

				point1.x() = point1.x()*point1.z(); 
				point1.y() = point1.y()*point1.z();

				Vector3d point_residual_3d = point1_reprojected - point1; 
				Matrix3d C_combined = C_p0_p1 + C_p1; 
				Matrix3d C_inverse = C_combined.inverse(); 
				double chi2_3d = point_residual_3d.transpose()*C_inverse*point_residual_3d; 
			
				if(chi2_3d < 7.815) // 6.25, p-value is 0.9, 7.815, p-value is 0.95 
				{ // pass chi-square test 

					// merge gaussian 
					// Approximated covariance estimation in graphical approaches to SLAM 
					Vector3d point_3d_1; 
					point_3d_1 = Pose1_R.transpose()*point1 - Pose1_R.transpose()*Pose1_t;
					Matrix3d C_p0_w, C_p1_w; 
					C_p0_w = Pose0_R.transpose()*C_p0*Pose0_R; 
					C_p1_w = Pose1_R.transpose()*C_p1*Pose1_R; 

					Matrix3d Omega_p0_w, Omega_p1_w; 
					Omega_p0_w = C_p0_w.inverse(); 
					Omega_p1_w = C_p1_w.inverse(); 

					// cout<<"Omega_p0_w: "<<endl<<Omega_p0_w<<endl; 
					// cout<<"Omega_p0_w.determinant: "<<Omega_p0_w.determinant()<<endl; 
					// cout<<"Omega_p1_w: "<<endl<<Omega_p1_w<<endl;
					// cout<<"Omega_p1_w.determinant: "<<Omega_p1_w.determinant()<<endl; 

					double w0 = 0.5; //  exp(Omega_p0_w.determinant()); 
					double w1 = 0.5; // exp(Omega_p1_w.determinant()); 
					double w01 = w0+w1; 
					w0 = w0/w01; 
					w1 = w1/w01; 

					Matrix3d Omega_p01 = w0*Omega_p0_w + w1*Omega_p1_w; 
					Matrix3d C_p01 = Omega_p01.inverse(); 

					Vector3d p3d_01 = C_p01*(w0*Omega_p0_w*point_3d+w1*Omega_p1_w*point_3d_1); 
				
					sfm_f[j].state = true;
					sfm_f[j].position[0] = p3d_01(0); // point_3d(0);
					sfm_f[j].position[1] = p3d_01(1); // point_3d(1);
					sfm_f[j].position[2] = p3d_01(2); // point_3d(2);
					++cnt_valid;
					// cout<<"good case, chi2_3d: "<<chi2_3d<<endl;
					// cout<<"let'see point_3d: "<<point_3d.transpose()<<" p3d_01: "<<p3d_01.transpose()<<endl; 
					// cout<<"C_p01: "<<endl<<C_p01<<endl;  
					// cout<<"Omega_p0_w: "<<endl<<Omega_p0_w<<endl; 
					// cout<<"Omega_p1_w: "<<endl<<Omega_p1_w<<endl; 
					// cout<<"w0: "<<w0<<" w1: "<<w1<<endl; 
				}else{
				 	// cout<<"bad case, chi2_3d: "<<chi2_3d<<endl; 
					// cout<<"point_residual_3d: "<<point_residual_3d.transpose()<<endl; 
					// cout<<"C_inverse: "<<endl<<C_inverse<<endl; 

				}
				
			}else if (point0.z()<4){ // just use z0, if it is < 4m, accoding to characterization results  
				
				Vector3d point_3d;
				point0.x() = point0.x()*point0.z(); 
				point0.y() = point0.y()*point0.z(); 
				point_3d = Pose0_R.transpose()*point0 - Pose0_R.transpose()*Pose0_t;//shan add:this is point in world;
				sfm_f[j].state = true;
				sfm_f[j].position[0] = point_3d(0);
				sfm_f[j].position[1] = point_3d(1);
				sfm_f[j].position[2] = point_3d(2);
				++cnt_valid;
				// cout<<"for feature: "<< j<<" just use point0 depth measurement: "<<point0.z()<<endl;
			}else if(point1.z()<4){  // just use z1, if it is < 4m, accoding to characterization results  
				
				Vector3d point_3d;
				point1.x() = point1.x()*point1.z(); 
				point1.y() = point1.y()*point1.z();
				point_3d = Pose1_R.transpose()*point1 - Pose1_R.transpose()*Pose1_t;//shan add:this is point in world;
				sfm_f[j].state = true;
				sfm_f[j].position[0] = point_3d(0);
				sfm_f[j].position[1] = point_3d(1);
				sfm_f[j].position[2] = point_3d(2);
				++cnt_valid;
				// cout<<"for feature: "<< j<<" just use point1 depth measurement: "<<point1.z()<<endl;
			}


			}else{

			if(has_0 && point0.z()<4){
				Vector3d point_3d;
				point0.x() = point0.x()*point0.z(); 
				point0.y() = point0.y()*point0.z(); 
				point_3d = Pose0_R.transpose()*point0 - Pose0_R.transpose()*Pose0_t;//shan add:this is point in world;
				sfm_f[j].state = true;
				sfm_f[j].position[0] = point_3d(0);
				sfm_f[j].position[1] = point_3d(1);
				sfm_f[j].position[2] = point_3d(2);
				++cnt_valid;
				// cout<<"for feature: "<< j<<" just use point0 depth measurement: "<<point0.z()<<endl;
			}else if(has_1 && point1.z()<4){
				Vector3d point_3d;
				point1.x() = point1.x()*point1.z(); 
				point1.y() = point1.y()*point1.z(); 
				point_3d = Pose1_R.transpose()*point1 - Pose1_R.transpose()*Pose1_t;//shan add:this is point in world;
				sfm_f[j].state = true;
				sfm_f[j].position[0] = point_3d(0);
				sfm_f[j].position[1] = point_3d(1);
				sfm_f[j].position[2] = point_3d(2);
				++cnt_valid;
				// cout<<"for feature: "<< j<<" just use point1 depth measurement: "<<point1.z()<<endl;
			}
			
		}

			//cout << "trangulated : " << frame1 << "  3d point : "  << j << "  " << point_3d.transpose() << endl;
	}
	// cout<<"cnt_valid sfm features: "<<cnt_valid<<" unused feature number: "<<feature_num - cnt_valid<<endl; 

}

void GlobalSFM::triangulateTwoFramesWithDepth(int frame0, Eigen::Matrix<double, 3, 4> &Pose0,
                                     int frame1, Eigen::Matrix<double, 3, 4> &Pose1,
                                     vector<SFMFeature> &sfm_f)
{
	assert(frame0 != frame1);
	Matrix3d Pose0_R = Pose0.block< 3,3 >(0,0);
	Matrix3d Pose1_R = Pose1.block< 3,3 >(0,0);
	Vector3d Pose0_t = Pose0.block< 3,1 >(0,3);
	Vector3d Pose1_t = Pose1.block< 3,1 >(0,3);
	int cnt_valid = 0; 
	for (int j = 0; j < feature_num; j++)
	{
		if (sfm_f[j].state == true)
			continue;
		bool has_0 = false, has_1 = false;
		Vector3d point0;
		Vector3d point1;
		for (int k = 0; k < (int)sfm_f[j].observation.size(); k++)
		{
			// if (sfm_f[j].observation_depth[k].second < 0.1 || sfm_f[j].observation_depth[k].second >10) //max and min measurement
				// continue;
			if (sfm_f[j].observation[k].first == frame0)
			{
				point0 = Vector3d(sfm_f[j].observation[k].second.x(),sfm_f[j].observation[k].second.y(), sfm_f[j].observation_depth[k].second);	
				has_0 = true; 
			}
			if (sfm_f[j].observation[k].first == frame1)
			{
				// point1 = sfm_f[j].observation[k].second;
				point1 = Vector3d(sfm_f[j].observation[k].second.x(),sfm_f[j].observation[k].second.y(), sfm_f[j].observation_depth[k].second);
				has_1 = true;
			}
			if(has_0 && has_1) break;
		}

		if (has_0 && has_1)
		{
			Vector2d residual;
			// point 0 has depth 
			if(point0.z() > 0.1){
				Vector3d point_3d, point1_reprojected;
				point0.x() = point0.x()*point0.z(); 
				point0.y() = point0.y()*point0.z(); 
				//triangulatePoint(Pose0, Pose1, point0, point1, point_3d);
				point_3d = Pose0_R.transpose()*point0 - Pose0_R.transpose()*Pose0_t;//shan add:this is point in world;
				point1_reprojected = Pose1_R*point_3d+Pose1_t;

				Vector2d pt1(point1.x(), point1.y()); 
				if(point1_reprojected.z() <= 0.1){
					residual = Vector2d(10, 10);
				}else
					residual = pt1 - Vector2d(point1_reprojected.x()/point1_reprojected.z(),point1_reprojected.y()/point1_reprojected.z());

				//std::cout << residual.transpose()<<"norm"<<residual.norm()*460<<endl;
				if (residual.norm() < 1.0/460){
					sfm_f[j].state = true;
					sfm_f[j].position[0] = point_3d(0);
					sfm_f[j].position[1] = point_3d(1);
					sfm_f[j].position[2] = point_3d(2);
					++cnt_valid;
					cout<<"good case, residual: "<<residual.transpose()<<" norm "<<residual.norm()*460<<endl;
				}else{
					cout<<"bad case, residual: "<<residual.transpose()<<" norm "<<residual.norm()*460<<endl;
				}
			}else if(point1.z() > 0.1){
				Vector3d point_3d, point1_reprojected;
				point1.x() = point1.x()*point1.z(); 
				point1.y() = point1.y()*point1.z(); 
				//triangulatePoint(Pose0, Pose1, point0, point1, point_3d);
				point_3d = Pose1_R.transpose()*point1 - Pose1_R.transpose()*Pose1_t;//shan add:this is point in world;
				point1_reprojected = Pose0_R*point_3d+Pose0_t;

				Vector2d pt0(point0.x(), point0.y());
				if(point1_reprojected.z() <= 0.1){
					residual = Vector2d(10,10);
				}else 
					residual = pt0 - Vector2d(point1_reprojected.x()/point1_reprojected.z(),point1_reprojected.y()/point1_reprojected.z());

				//std::cout << residual.transpose()<<"norm"<<residual.norm()*460<<endl;
				if (residual.norm() < 1.0/460){
					sfm_f[j].state = true;
					sfm_f[j].position[0] = point_3d(0);
					sfm_f[j].position[1] = point_3d(1);
					sfm_f[j].position[2] = point_3d(2);
					++cnt_valid;
				}
			}
			//cout << "trangulated : " << frame1 << "  3d point : "  << j << "  " << point_3d.transpose() << endl;
		}
	}
	// cout<<"cnt_valid sfm features: "<<cnt_valid<<" failed feature number: "<<feature_num - cnt_valid<<endl; 
}


void GlobalSFM::triangulateTwoFrames(int frame0, Eigen::Matrix<double, 3, 4> &Pose0, 
									 int frame1, Eigen::Matrix<double, 3, 4> &Pose1,
									 vector<SFMFeature> &sfm_f)
{
	assert(frame0 != frame1);
	for (int j = 0; j < feature_num; j++)
	{
		if (sfm_f[j].state == true)
			continue;
		bool has_0 = false, has_1 = false;
		Vector2d point0;
		Vector2d point1;
		for (int k = 0; k < (int)sfm_f[j].observation.size(); k++)
		{
			if (sfm_f[j].observation[k].first == frame0)
			{
				point0 = sfm_f[j].observation[k].second;
				has_0 = true;
			}
			if (sfm_f[j].observation[k].first == frame1)
			{
				point1 = sfm_f[j].observation[k].second;
				has_1 = true;
			}
		}
		if (has_0 && has_1)
		{
			Vector3d point_3d;
			triangulatePoint(Pose0, Pose1, point0, point1, point_3d);
			sfm_f[j].state = true;
			sfm_f[j].position[0] = point_3d(0);
			sfm_f[j].position[1] = point_3d(1);
			sfm_f[j].position[2] = point_3d(2);
			// cout << "trangulated : " << frame1 << "  3d point : "  << j << "  " << point_3d.transpose() << endl;
		}							  
	}
}

// 	 q w_R_cam t w_R_cam
//  c_rotation cam_R_w 
//  c_translation cam_R_w
// relative_q[i][j]  j_q_i
// relative_t[i][j]  j_t_ji  (j < i)
bool GlobalSFM::construct(int frame_num, Quaterniond* q, Vector3d* T, int l,
			  const Matrix3d relative_R, const Vector3d relative_T,
			  vector<SFMFeature> &sfm_f, map<int, Vector3d> &sfm_tracked_points)
{
	feature_num = sfm_f.size();
	//cout << "set 0 and " << l << " as known " << endl;
	// have relative_r relative_t
	// intial two view
	q[l].w() = 1;
	q[l].x() = 0;
	q[l].y() = 0;
	q[l].z() = 0;
	T[l].setZero();
	q[frame_num - 1] = q[l] * Quaterniond(relative_R);
	T[frame_num - 1] = relative_T;
	//cout << "init q_l " << q[l].w() << " " << q[l].vec().transpose() << endl;
	//cout << "init t_l " << T[l].transpose() << endl;

	//rotate to cam frame
	Matrix3d c_Rotation[frame_num];
	Vector3d c_Translation[frame_num];
	Quaterniond c_Quat[frame_num];
	double c_rotation[frame_num][4];
	double c_translation[frame_num][3];
	Eigen::Matrix<double, 3, 4> Pose[frame_num];

	c_Quat[l] = q[l].inverse();
	c_Rotation[l] = c_Quat[l].toRotationMatrix();
	c_Translation[l] = -1 * (c_Rotation[l] * T[l]);
	Pose[l].block<3, 3>(0, 0) = c_Rotation[l];
	Pose[l].block<3, 1>(0, 3) = c_Translation[l];

	c_Quat[frame_num - 1] = q[frame_num - 1].inverse();
	c_Rotation[frame_num - 1] = c_Quat[frame_num - 1].toRotationMatrix();
	c_Translation[frame_num - 1] = -1 * (c_Rotation[frame_num - 1] * T[frame_num - 1]);
	Pose[frame_num - 1].block<3, 3>(0, 0) = c_Rotation[frame_num - 1];
	Pose[frame_num - 1].block<3, 1>(0, 3) = c_Translation[frame_num - 1];


	//1: trangulate between l ----- frame_num - 1
	//2: solve pnp l + 1; trangulate l + 1 ------- frame_num - 1; 
	for (int i = l; i < frame_num - 1 ; i++)
	{
		// solve pnp
		if (i > l)
		{
			Matrix3d R_initial = c_Rotation[i - 1];
			Vector3d P_initial = c_Translation[i - 1];
			if(!solveFrameByPnP(R_initial, P_initial, i, sfm_f))
				return false;
			c_Rotation[i] = R_initial;
			c_Translation[i] = P_initial;
			c_Quat[i] = c_Rotation[i];
			Pose[i].block<3, 3>(0, 0) = c_Rotation[i];
			Pose[i].block<3, 1>(0, 3) = c_Translation[i];
		}

		// triangulate point based on the solve pnp result
		// triangulateTwoFramesWithDepth(i, Pose[i], frame_num - 1, Pose[frame_num - 1], sfm_f);
		triangulateTwoFramesWithDepthCov(i, Pose[i], frame_num - 1, Pose[frame_num - 1], sfm_f);
		triangulateTwoFrames(i, Pose[i], frame_num - 1, Pose[frame_num - 1], sfm_f)	;
	}
	//3: triangulate l-----l+1 l+2 ... frame_num -2
	for (int i = l + 1; i < frame_num - 1; i++){
// 		triangulateTwoFramesWithDepth(l, Pose[l], i, Pose[i], sfm_f);
		triangulateTwoFramesWithDepthCov(i, Pose[i], frame_num - 1, Pose[frame_num - 1], sfm_f);
		triangulateTwoFrames(l, Pose[l], i, Pose[i], sfm_f);
	}

	//4: solve pnp l-1; triangulate l-1 ----- l
	//             l-2              l-2 ----- l
	for (int i = l - 1; i >= 0; i--)
	{
		//solve pnp
		Matrix3d R_initial = c_Rotation[i + 1];
		Vector3d P_initial = c_Translation[i + 1];
		if(!solveFrameByPnP(R_initial, P_initial, i, sfm_f))
			return false;
		c_Rotation[i] = R_initial;
		c_Translation[i] = P_initial;
		c_Quat[i] = c_Rotation[i];
		Pose[i].block<3, 3>(0, 0) = c_Rotation[i];
		Pose[i].block<3, 1>(0, 3) = c_Translation[i];
		//triangulate
// 		triangulateTwoFramesWithDepth(i, Pose[i],  l, Pose[l], sfm_f);
		triangulateTwoFramesWithDepthCov(i, Pose[i], frame_num - 1, Pose[frame_num - 1], sfm_f);
		triangulateTwoFrames(i, Pose[i], l, Pose[l], sfm_f);
	}
	//5: triangulate all other points
	// for (int j = 0; j < feature_num; j++)
	// {
	// 	if (sfm_f[j].state == true)
	// 		continue;
	// 	if ((int)sfm_f[j].observation.size() >= 2)
	// 	{
	// 		Vector2d point0, point1;
	// 		int frame_0 = sfm_f[j].observation[0].first;
	// 		point0 = sfm_f[j].observation[0].second;
	// 		int frame_1 = sfm_f[j].observation.back().first;
	// 		point1 = sfm_f[j].observation.back().second;
	// 		Vector3d point_3d;
	// 		triangulatePoint(Pose[frame_0], Pose[frame_1], point0, point1, point_3d);
	// 		sfm_f[j].state = true;
	// 		sfm_f[j].position[0] = point_3d(0);
	// 		sfm_f[j].position[1] = point_3d(1);
	// 		sfm_f[j].position[2] = point_3d(2);
	// 		//cout << "trangulated : " << frame_0 << " " << frame_1 << "  3d point : "  << j << "  " << point_3d.transpose() << endl;
	// 	}		
	// }

	//5: triangulate all other points
	for (int j = 0; j < feature_num; j++)
	{
		if (sfm_f[j].state == true)
			continue;
		if ((int)sfm_f[j].observation.size() >= 2)
		{
			Vector3d point0; 
			Vector2d point1;
			int frame_0 = sfm_f[j].observation[0].first;
			double depth_i = sfm_f[j].observation_depth[0].second; 

			if(depth_i <= 0.1)
				continue; // not a valid point 


			point0 = Vector3d(sfm_f[j].observation[0].second.x()*depth_i,sfm_f[j].observation[0].second.y()*depth_i,depth_i);
			int frame_1 = sfm_f[j].observation.back().first;
			point1 = sfm_f[j].observation.back().second;
			Vector3d point_3d;
			//triangulatePoint(Pose[frame_0], Pose[frame_1], point0, point1, point_3d);

			Matrix3d Pose0_R = Pose[frame_0].block< 3,3 >(0,0);
			Matrix3d Pose1_R = Pose[frame_1].block< 3,3 >(0,0);
			Vector3d Pose0_t = Pose[frame_0].block< 3,1 >(0,3);
			Vector3d Pose1_t = Pose[frame_1].block< 3,1 >(0,3);

			Vector2d residual;
			Vector3d point1_reprojected;
			//triangulatePoint(Pose0, Pose1, point0, point1, point_3d);
			point_3d = Pose0_R.transpose()*point0 - Pose0_R.transpose()*Pose0_t;//point in world;
			point1_reprojected = Pose1_R*point_3d+Pose1_t;

			if(point1_reprojected.z()<0.1){
				residual = Vector2d(10,10);
			}else
				residual = point1 - Vector2d(point1_reprojected.x()/point1_reprojected.z(),point1_reprojected.y()/point1_reprojected.z());

			if (residual.norm() < 1.0/460) {//reprojection error
				sfm_f[j].state = true;
				sfm_f[j].position[0] = point_3d(0);
				sfm_f[j].position[1] = point_3d(1);
				sfm_f[j].position[2] = point_3d(2);
			}
			//cout << "trangulated : " << frame_0 << " " << frame_1 << "  3d point : "  << j << "  " << point_3d.transpose() << endl;
		}		
	}

/*
	for (int i = 0; i < frame_num; i++)
	{
		q[i] = c_Rotation[i].transpose(); 
		cout << "solvePnP  q" << " i " << i <<"  " <<q[i].w() << "  " << q[i].vec().transpose() << endl;
	}
	for (int i = 0; i < frame_num; i++)
	{
		Vector3d t_tmp;
		t_tmp = -1 * (q[i] * c_Translation[i]);
		cout << "solvePnP  t" << " i " << i <<"  " << t_tmp.x() <<"  "<< t_tmp.y() <<"  "<< t_tmp.z() << endl;
	}
*/
	//full BA
	ceres::Problem problem;
	ceres::LocalParameterization* local_parameterization = new ceres::QuaternionParameterization();
	//cout << " begin full BA " << endl;
	for (int i = 0; i < frame_num; i++)
	{
		//double array for ceres
		c_translation[i][0] = c_Translation[i].x();
		c_translation[i][1] = c_Translation[i].y();
		c_translation[i][2] = c_Translation[i].z();
		c_rotation[i][0] = c_Quat[i].w();
		c_rotation[i][1] = c_Quat[i].x();
		c_rotation[i][2] = c_Quat[i].y();
		c_rotation[i][3] = c_Quat[i].z();
		problem.AddParameterBlock(c_rotation[i], 4, local_parameterization);
		problem.AddParameterBlock(c_translation[i], 3);
		if (i == l)
		{
			problem.SetParameterBlockConstant(c_rotation[i]);
		}
		if (i == l || i == frame_num - 1)
		{
			problem.SetParameterBlockConstant(c_translation[i]);
		}
	}

	for (int i = 0; i < feature_num; i++)
	{
		if (sfm_f[i].state != true)
			continue;
		for (int j = 0; j < int(sfm_f[i].observation.size()); j++)
		{
			int l = sfm_f[i].observation[j].first;
			ceres::CostFunction* cost_function = ReprojectionError3D::Create(
												sfm_f[i].observation[j].second.x(),
												sfm_f[i].observation[j].second.y());

    		problem.AddResidualBlock(cost_function, NULL, c_rotation[l], c_translation[l], 
    								sfm_f[i].position);	 
		}

	}
	ceres::Solver::Options options;
	options.linear_solver_type = ceres::DENSE_SCHUR;
	// options.minimizer_progress_to_stdout = true;
	options.max_solver_time_in_seconds = 0.2; //*10; // oranfire modify
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	// std::cout << summary.BriefReport() << "\n";
	if (summary.termination_type == ceres::CONVERGENCE || summary.final_cost < 5e-03)
	{
		cout << "initial_sfm.cpp: vision only BA converge" << endl;
	}
	else
	{
		cout << "initial_sfm.cpp: vision only BA not converge " << endl;
		return false;
	}
	for (int i = 0; i < frame_num; i++)
	{
		q[i].w() = c_rotation[i][0]; 
		q[i].x() = c_rotation[i][1]; 
		q[i].y() = c_rotation[i][2]; 
		q[i].z() = c_rotation[i][3]; 
		q[i] = q[i].inverse();
		//cout << "final  q" << " i " << i <<"  " <<q[i].w() << "  " << q[i].vec().transpose() << endl;
	}
	for (int i = 0; i < frame_num; i++)
	{

		T[i] = -1 * (q[i] * Vector3d(c_translation[i][0], c_translation[i][1], c_translation[i][2]));
		//cout << "final  t" << " i " << i <<"  " << T[i](0) <<"  "<< T[i](1) <<"  "<< T[i](2) << endl;
	}
	for (int i = 0; i < (int)sfm_f.size(); i++)
	{
		if(sfm_f[i].state)
			sfm_tracked_points[sfm_f[i].id] = Vector3d(sfm_f[i].position[0], sfm_f[i].position[1], sfm_f[i].position[2]);
	}
	return true;

}

