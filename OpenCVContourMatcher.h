#include "ContourMatcher.h"

#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/nonfree/nonfree.hpp"
#include <string>

#ifndef __OPENCV_CONTOUR_MATCHER__
#define __OPENCV_CONTOUR_MATCHER__
class OpenCVContourMatcher : public ContourMatcher {
	private:
		double maxDistance = 0;
		double minDistance = 100;
		int matrixColNum;
		cv::Mat image1;
		cv::Mat image2;

		Plane convertDMatchesToPlane(std::vector<cv:DMatch> matches);
		Plane compareFLANN(Contour contour1, Contour contour2);
		Plane compareBruteForce(Contour contour1, Contour contour2);

		std::vector<cv::Keypoint> detectFeaturePoints(Contour contour);
		cv::Mat extractDescriptors(Contour contour, std::vector<cv::Keypoint> contourFeatures);
		std::vector<cv::DMatch> matchBruteForce(cv::Mat desc1, cv::Mat desc2);
		std::vector<cv::DMatch> matchFLANN(cv::Mat desc1, cv::Mat desc2);
		std::vector<cv::DMatch> refineMatches(std::vector<DMatch> matches);
		void computeMaxAndMinDistances(std::vector<cv::DMatch> matches);
	public:
		Plane compare(Contour contour1, Contour contour2, MATCHER_TYPE m);
		void loadImages(unsigned char* i1, unsigned int height1, unsigned int width1, unsigned char* i2, unsigned int height2, unsigned int width2);
		void loadImages(Image i1, Image i2);
};

#endif //__OPEN_CV_CONTOUR_MATCHER__