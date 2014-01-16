#include "ContourMatcher.h"
#include "ContourMatcherStructs.h"
#include "ContourMatcherExceptions.h"

#include <opencv2/opencv.hpp>
#include "opencv2/highgui/highgui.hpp"
#include <opencv2/features2d/features2d.hpp>
#include "opencv2/imgproc/imgproc.hpp"

#include <opencv2/nonfree/nonfree.hpp>
#include <opencv2/nonfree/features2d.hpp>
#include <opencv2/legacy/legacy.hpp>


#ifndef __OPENCV_CONTOUR_MATCHER__
#define __OPENCV_CONTOUR_MATCHER__
class OpenCVContourMatcher : public ContourMatcher {
	private:
		double maxDistance;
		double minDistance;
		cv::Mat image1;
		cv::Mat image2;

		Plane convertDMatchesToPlane(std::vector< cv::KeyPoint > features, std::vector< cv::KeyPoint > features2, std::vector<cv::DMatch> matches);
		Plane compareFLANN(Contour contour1, Contour contour2);
		Plane compareBruteForce(Contour contour1, Contour contour2);

		//cv::Mat convertToMatrix(Contour contour);

//		cv::Mat sliceContour(std::vector<std::vector<PixelLoc> > contourPixels,cv::Mat image);
		cv::Mat sliceContour(std::vector<PixelLoc> contourPixels,cv::Mat image);



		std::vector<cv::KeyPoint> detectFeaturePoints(cv::Mat contour);
		cv::Mat extractDescriptors(cv::Mat contour, std::vector<cv::KeyPoint> contourFeatures);
		std::vector<cv::DMatch> matchBruteForce(cv::Mat desc1, cv::Mat desc2);
		std::vector<cv::DMatch> matchFLANN(cv::Mat desc1, cv::Mat desc2);
		std::vector<cv::DMatch> refineMatches(std::vector<cv::DMatch> matches);
		void computeMaxAndMinDistances(std::vector<cv::DMatch> matches);
	public:
		Plane compare(Contour contour1, Contour contour2, MATCHER_TYPE m);
		void loadImages(unsigned char* i1, unsigned int height1, unsigned int width1, unsigned char* i2, unsigned int height2, unsigned int width2);
		void loadImages(Image i1, Image i2);
};

#endif //__OPEN_CV_CONTOUR_MATCHER__