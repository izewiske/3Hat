#include "ContourMatcherStructs.h"

#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/nonfree/nonfree.hpp"

#ifndef __OPENCV_CONTOUR_MATCHER__
#define __OPENCV_CONTOUR_MATCHER__
class OpenCVContourMatcher : public ContourMatcher{
	private:
		cv::Mat convertToMatrix(Contour c);
		std::vector<int> convertDMatchesToSetOfPoints(std::vector<cv:DMatch> matches);
		std::vector<int> compareFLANN(Contour contour1, Contour contour2);
		std::vector<int> compareBruteForce(Contour contour1, Contour contour2);
		std::vector<cv::Keypoint> detectFeaturePoints(Contour contour);
		cv::Mat extractDescriptors(Contour contour, std::vector<cv::Keypoint> contourFeatures);
		std::vector<cv::DMatch> matchBruteForce(cv::Mat desc1, cv::Mat desc2);
		std::vector<cv::DMatch> matchFLANN(cv::Mat desc1, cv::Mat desc2);
	public:
		std::vector<int> compare(Contour contour1, Contour contour2, MATCHER_TYPE m);
};

#endif //__OPEN_CV_CONTOUR_MATCHER__