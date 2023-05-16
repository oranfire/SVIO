#include "solve_5pts.h"
#include "solve_opt.h"

#ifdef USE_OPENGV

#include <opengv/relative_pose/methods.hpp>
#include <opengv/relative_pose/CentralRelativeAdapter.hpp>
#endif 


namespace cv {
    void decomposeEssentialMat( InputArray _E, OutputArray _R1, OutputArray _R2, OutputArray _t )
    {

        Mat E = _E.getMat().reshape(1, 3);
        CV_Assert(E.cols == 3 && E.rows == 3);

        Mat D, U, Vt;
        SVD::compute(E, D, U, Vt);

        if (determinant(U) < 0) U *= -1.;
        if (determinant(Vt) < 0) Vt *= -1.;

        Mat W = (Mat_<double>(3, 3) << 0, 1, 0, -1, 0, 0, 0, 0, 1);
        W.convertTo(W, E.type());

        Mat R1, R2, t;
        R1 = U * W * Vt;
        R2 = U * W.t() * Vt;
        t = U.col(2) * 1.0;

        R1.copyTo(_R1);
        R2.copyTo(_R2);
        t.copyTo(_t);
    }

    int recoverPose( InputArray E, InputArray _points1, InputArray _points2, InputArray _cameraMatrix,
                         OutputArray _R, OutputArray _t, InputOutputArray _mask)
    {

        Mat points1, points2, cameraMatrix;
        _points1.getMat().convertTo(points1, CV_64F);
        _points2.getMat().convertTo(points2, CV_64F);
        _cameraMatrix.getMat().convertTo(cameraMatrix, CV_64F);

        int npoints = points1.checkVector(2);
        CV_Assert( npoints >= 0 && points2.checkVector(2) == npoints &&
                                  points1.type() == points2.type());

        CV_Assert(cameraMatrix.rows == 3 && cameraMatrix.cols == 3 && cameraMatrix.channels() == 1);

        if (points1.channels() > 1)
        {
            points1 = points1.reshape(1, npoints);
            points2 = points2.reshape(1, npoints);
        }

        double fx = cameraMatrix.at<double>(0,0);
        double fy = cameraMatrix.at<double>(1,1);
        double cx = cameraMatrix.at<double>(0,2);
        double cy = cameraMatrix.at<double>(1,2);

        points1.col(0) = (points1.col(0) - cx) / fx;
        points2.col(0) = (points2.col(0) - cx) / fx;
        points1.col(1) = (points1.col(1) - cy) / fy;
        points2.col(1) = (points2.col(1) - cy) / fy;

        points1 = points1.t();
        points2 = points2.t();

        Mat R1, R2, t;
        decomposeEssentialMat(E, R1, R2, t);
        Mat P0 = Mat::eye(3, 4, R1.type());
        Mat P1(3, 4, R1.type()), P2(3, 4, R1.type()), P3(3, 4, R1.type()), P4(3, 4, R1.type());
        P1(Range::all(), Range(0, 3)) = R1 * 1.0; P1.col(3) = t * 1.0;
        P2(Range::all(), Range(0, 3)) = R2 * 1.0; P2.col(3) = t * 1.0;
        P3(Range::all(), Range(0, 3)) = R1 * 1.0; P3.col(3) = -t * 1.0;
        P4(Range::all(), Range(0, 3)) = R2 * 1.0; P4.col(3) = -t * 1.0;

        // Do the cheirality check.
        // Notice here a threshold dist is used to filter
        // out far away points (i.e. infinite points) since
        // there depth may vary between postive and negtive.
        double dist = 50.0;
        Mat Q;
        triangulatePoints(P0, P1, points1, points2, Q);
        Mat mask1 = Q.row(2).mul(Q.row(3)) > 0;
        Q.row(0) /= Q.row(3);
        Q.row(1) /= Q.row(3);
        Q.row(2) /= Q.row(3);
        Q.row(3) /= Q.row(3);
        mask1 = (Q.row(2) < dist) & mask1;
        Q = P1 * Q;
        mask1 = (Q.row(2) > 0) & mask1;
        mask1 = (Q.row(2) < dist) & mask1;

        triangulatePoints(P0, P2, points1, points2, Q);
        Mat mask2 = Q.row(2).mul(Q.row(3)) > 0;
        Q.row(0) /= Q.row(3);
        Q.row(1) /= Q.row(3);
        Q.row(2) /= Q.row(3);
        Q.row(3) /= Q.row(3);
        mask2 = (Q.row(2) < dist) & mask2;
        Q = P2 * Q;
        mask2 = (Q.row(2) > 0) & mask2;
        mask2 = (Q.row(2) < dist) & mask2;

        triangulatePoints(P0, P3, points1, points2, Q);
        Mat mask3 = Q.row(2).mul(Q.row(3)) > 0;
        Q.row(0) /= Q.row(3);
        Q.row(1) /= Q.row(3);
        Q.row(2) /= Q.row(3);
        Q.row(3) /= Q.row(3);
        mask3 = (Q.row(2) < dist) & mask3;
        Q = P3 * Q;
        mask3 = (Q.row(2) > 0) & mask3;
        mask3 = (Q.row(2) < dist) & mask3;

        triangulatePoints(P0, P4, points1, points2, Q);
        Mat mask4 = Q.row(2).mul(Q.row(3)) > 0;
        Q.row(0) /= Q.row(3);
        Q.row(1) /= Q.row(3);
        Q.row(2) /= Q.row(3);
        Q.row(3) /= Q.row(3);
        mask4 = (Q.row(2) < dist) & mask4;
        Q = P4 * Q;
        mask4 = (Q.row(2) > 0) & mask4;
        mask4 = (Q.row(2) < dist) & mask4;

        mask1 = mask1.t();
        mask2 = mask2.t();
        mask3 = mask3.t();
        mask4 = mask4.t();

        // If _mask is given, then use it to filter outliers.
        if (!_mask.empty())
        {
            Mat mask = _mask.getMat();
            CV_Assert(mask.size() == mask1.size());
            bitwise_and(mask, mask1, mask1);
            bitwise_and(mask, mask2, mask2);
            bitwise_and(mask, mask3, mask3);
            bitwise_and(mask, mask4, mask4);
        }
        if (_mask.empty() && _mask.needed())
        {
            _mask.create(mask1.size(), CV_8U);
        }

        CV_Assert(_R.needed() && _t.needed());
        _R.create(3, 3, R1.type());
        _t.create(3, 1, t.type());

        int good1 = countNonZero(mask1);
        int good2 = countNonZero(mask2);
        int good3 = countNonZero(mask3);
        int good4 = countNonZero(mask4);

        if (good1 >= good2 && good1 >= good3 && good1 >= good4)
        {
            R1.copyTo(_R);
            t.copyTo(_t);
            if (_mask.needed()) mask1.copyTo(_mask);
            return good1;
        }
        else if (good2 >= good1 && good2 >= good3 && good2 >= good4)
        {
            R2.copyTo(_R);
            t.copyTo(_t);
            if (_mask.needed()) mask2.copyTo(_mask);
            return good2;
        }
        else if (good3 >= good1 && good3 >= good2 && good3 >= good4)
        {
            t = -t;
            R1.copyTo(_R);
            t.copyTo(_t);
            if (_mask.needed()) mask3.copyTo(_mask);
            return good3;
        }
        else
        {
            t = -t;
            R2.copyTo(_R);
            t.copyTo(_t);
            if (_mask.needed()) mask4.copyTo(_mask);
            return good4;
        }
    }

    int recoverPose( InputArray E, InputArray _points1, InputArray _points2, OutputArray _R,
                         OutputArray _t, double focal, Point2d pp, InputOutputArray _mask)
    {
        Mat cameraMatrix = (Mat_<double>(3,3) << focal, 0, pp.x, 0, focal, pp.y, 0, 0, 1);
        return cv::recoverPose(E, _points1, _points2, cameraMatrix, _R, _t, _mask);
    }
}


namespace{

vector<pair<Vector3d, Vector3d>> getInliers2D(const vector<pair<Vector3d, Vector3d>> &corres, cv::Mat& mask)
{
    assert(mask.rows == corres.size()); 

    vector<pair<Vector3d, Vector3d>> ret;
    for(int i=0; i<corres.size(); i++){

        if(mask.at<unsigned char>(i,0) == 0) continue; 
        
        ret.push_back(make_pair(corres[i].first, corres[i].second)); 
    } 
    return ret; 
}

vector<pair<Vector3d, Vector3d>> getInliersWithValidDepth(const vector<pair<Vector3d, Vector3d>> &corres, cv::Mat& mask)
{
    assert(mask.rows == corres.size()); 

    vector<pair<Vector3d, Vector3d>> ret;
    for(int i=0; i<corres.size(); i++){

        if(mask.at<unsigned char>(i,0) == 0) continue; 

        if(corres[i].first.z() <= 0.1 )
            continue; // invalid depth pair
        
        ret.push_back(make_pair(corres[i].first, corres[i].second)); 
    } 
    return ret; 
}

vector<pair<Vector3d, Vector3d>> getInliersWithValidDepthCov(const vector<pair<Vector3d, Vector3d>> &corres, cv::Mat& mask, vector<Vector3d>& cov_v)
{
    assert(mask.rows == corres.size()); 

    vector<pair<Vector3d, Vector3d>> ret;
    vector<Vector3d> tmp_cv; 
    for(int i=0; i<corres.size(); i++){

        if(mask.at<unsigned char>(i,0) == 0) continue; 

        if(corres[i].first.z() <= 0.1 )
            continue; // invalid depth pair
        
        ret.push_back(make_pair(corres[i].first, corres[i].second)); 
        tmp_cv.push_back(cov_v[i]); 
    } 
    cov_v = tmp_cv; 
    return ret; 
}

}

bool MotionEstimator::solveRelativeRT(const vector<pair<Vector3d, Vector3d>> &corres, Matrix3d &Rotation, Vector3d &Translation)
{
    if (corres.size() >= 15)
    {
        vector<cv::Point2f> ll, rr;
        for (int i = 0; i < int(corres.size()); i++)
        {
            ll.push_back(cv::Point2f(corres[i].first(0), corres[i].first(1)));
            rr.push_back(cv::Point2f(corres[i].second(0), corres[i].second(1)));
        }
        cv::Mat mask;
        cv::Mat E = cv::findFundamentalMat(ll, rr, cv::FM_RANSAC, 0.3 / 460, 0.99, mask);
        cv::Mat cameraMatrix = (cv::Mat_<double>(3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);
        cv::Mat rot, trans;
        int inlier_cnt = cv::recoverPose(E, ll, rr, cameraMatrix, rot, trans, mask);
        cout << "solve_5pts.cpp: in solveRelativeRT(): inlier_cnt " << inlier_cnt << endl;

        Eigen::Matrix3d R;
        Eigen::Vector3d T;
        for (int i = 0; i < 3; i++)
        {   
            T(i) = trans.at<double>(i, 0);
            for (int j = 0; j < 3; j++)
                R(i, j) = rot.at<double>(i, j);
        }

        Rotation = R.transpose();
        Translation = -R.transpose() * T;

        if(inlier_cnt > 12){
            ROS_WARN("---------------5points----------------");
            ROS_WARN("input points %d inliers %d", ll.size(), inlier_cnt); 

            return true;
        }
        else
            return false;
    }
    return false;
}


bool MotionEstimator::solveRelativeHybrid(const vector<pair<Vector3d, Vector3d>> &corres, Matrix3d &Rotation, Vector3d &Translation, vector<Vector3d>* pcov)
{
    OptSolver sopt; 
    if (corres.size() >= 15)
    {
        // solve [R,t] using 2D-2D match 
        vector<cv::Point2f> ll, rr;
        for (int i = 0; i < int(corres.size()); i++)
        {
            ll.push_back(cv::Point2f(corres[i].first(0), corres[i].first(1)));
            rr.push_back(cv::Point2f(corres[i].second(0), corres[i].second(1)));
        }
        cv::Mat mask;
        cv::Mat E = cv::findFundamentalMat(ll, rr, cv::FM_RANSAC, 0.3 / 460, 0.99, mask);
        cv::Mat cameraMatrix = (cv::Mat_<double>(3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);
        cv::Mat rot, trans;
        int inlier_cnt = cv::recoverPose(E, ll, rr, cameraMatrix, rot, trans, mask);
        cout << "solve_5pts.cpp: in solveRelativeHybrid(): inlier_cnt " << inlier_cnt << endl;

	    if(inlier_cnt <= 12){
		  ROS_DEBUG("solve_5pts.cpp: in solveRelativeHybrid(): inlier_cnt %d is too small, return false", inlier_cnt); 
		  return false;
	    }
		
        Eigen::Matrix3d R;
        Eigen::Vector3d T;
        for (int i = 0; i < 3; i++)
        {   
            T(i) = trans.at<double>(i, 0);
            for (int j = 0; j < 3; j++)
                R(i, j) = rot.at<double>(i, j);
        }
         
        // with the speficied number of features 
        vector<pair<Vector3d, Vector3d>> inliers; 

        // get rid of corres with invalid depth  
        if(pcov == NULL){
            inliers = getInliersWithValidDepth(corres, mask); 
        }else{
            inliers = getInliersWithValidDepthCov(corres, mask, *pcov); 
        }

        // if use opengv to refine the rotation 
#ifdef USE_OPENGV

        cout<<"before opengv R: "<<endl<<R<<endl; 
        std::vector<pair<Vector3d, Vector3d>> inliers_2d = getInliers2D(corres, mask);

	 // ofstream ouf("matched_inliers.log"); 
	 

        Eigen::Matrix3d Rij_e = R.transpose();

        opengv::bearingVectors_t bearingVectors1; 
        opengv::bearingVectors_t bearingVectors2;

        for(int j=0; j<inliers_2d.size(); j++){
            Vector3d pi(inliers_2d[j].first(0), inliers_2d[j].first(1), 1.); 
            Vector3d pj(inliers_2d[j].second(0), inliers_2d[j].second(1), 1.);
            bearingVectors1.push_back(pi/pi.norm());
            bearingVectors2.push_back(pj/pj.norm()); 
	     // ouf<<inliers_2d[j].first(0)<<" "<< inliers_2d[j].first(1)<<" "<<inliers_2d[j].first(2)<<" "<<
	// 	inliers_2d[j].second(0)<<" "<< inliers_2d[j].second(1)<<" "<<inliers_2d[j].second(2)<<endl; 
        }

        opengv::relative_pose::CentralRelativeAdapter adapter_rbs(
              bearingVectors1,
              bearingVectors2,
              Rij_e ); // rotation); R12 

        // first use eight pts to compute initial rotation 
        Rij_e = opengv::relative_pose::eigensolver(adapter_rbs);    // bug, return R12, Rij, but R is Rji 
        R = Rij_e.transpose(); 

        cout <<"after opengv R: "<<endl<<R<<endl; 
#endif

        if(inlier_cnt > 12 && inliers.size() >= 0){
            // optimization to solve R, T // R is Rji, T is tji
            // sopt.solveHybrid(inliers, R, T, pcov); 
            sopt.solveTCeres(inliers, R, T, pcov); 
            Rotation = R.transpose(); // Rotation is Rij
            Translation = -R.transpose() * T; // Translation is tij 

            // ROS_WARN("---------------5points----------------");
            ROS_WARN("solve_5pts.cpp: input points %d 2D inliers %d 3D inliers %d", ll.size(), inlier_cnt, inliers.size()); 
            return true;
        }
        else{
            ROS_DEBUG("solve_5pts.cpp: 2D inlier_cnt: %d 3D inliers: %d", inlier_cnt, inliers.size());
            return false;
        }
    }else{
        ROS_DEBUG("solve_5pts.cpp: 2D corresponds cnt: %d", corres.size());
    }
    return false;
}




