#include "ContourMatcher.h"

cv::Mat convertToMatrix(Contour c) {
	cv::Mat contour;
	// convert contour of some format to matrix
	return contour;
}


std::vector<int> convertDMatchesToSetOfPoints(std::vector<cv:DMatch> matches){
	std::vector<int> planePoints;
	// convert DMatches to tbd set of points
	return planePoints;
}


std::vector<cv::Keypoint> detectFeaturePoints(cv::Mat contour) {
	cv::SiftFeatureDetector detector;
	std::vector<cv::Keypoint> contourFeatures;
	detector.detect(contour,contourFeatures);
	return contourFeatures;
}


cv::Mat extractDescriptors(cv::Mat contour, std::vector<cv::Keypoint> contourFeatures){
	cv::SiftDescriptorExtractor extractor;
	cv::Mat contourDescriptors;
	extractor.compute(contour,contourFeatures,contourDescriptors);
	return contourDescriptors;
}

std::vector<int> compare(Contour contour1, Contour contour2, MATCHER_TYPE m = FLANN){
	try {
		if ( m == FLANN) {
			compareFLANN(contour1, contour2);
		} else if ( m == BF ){
			compareBruteForce(contour1, contour2);
		} else {
			throw 2;
		}
	} catch (int e) {
		ERR << "In OpenCVContourMatcher: Unknown matcher defined. Exitcode " << e << "\n";
	}
}



std::vector<int> compareFLANN(Contour contour1, Contour contour2) {
	//Convert from contour class to matrix
	cv::Mat c1 = convertToMatrix(contour1);
	cv::Mat c2 = convertToMatrix(contour2);

	std::vector<cv::Keypoint> feat1 = detectFeaturePoints(c1);
	std::vector<cv::Keypoint> feat2 = detectFeaturePoints(c2);

	cv::Mat desc1 = extractDescriptors(c1,feat1);
	cv::Mat desc2 = extractDescriptors(c2,feat2);

	std::vector<int> planePoints = convertDMatchesToSetOfPoints(matchFLANN(desc1,desc2));


	return planePoints;
}

std::vector<int> compareBruteForce(Contour contour1, Contour contour2){
	cv::Mat c1 = convertToMatrix(contour1);
	cv::Mat c2 = convertToMatrix(contour2);

	std::vector<cv::Keypoint> feat1 = detectFeaturePoints(c1);
	std::vector<cv::Keypoint> feat2 = detectFeaturePoints(c2);

	cv::Mat desc1 = extractDescriptors(c1,feat1);
	cv::Mat desc2 = extractDescriptors(c2,feat2);

	std::vector<int> planePoints = convertDMatchesToSetOfPoints(matchBruteForce(desc1,desc2));

	return planePoints;
}

// Takes two descriptor matrices
std::vector<cv::DMatch> matchBruteForce(cv::Mat desc1, cv::Mat desc2){
	BruteForceMatcher matcher;
	std::vector<cv::DMatch> matches;
	matcher.match(desc1,desc2,matches);
	return matches;
}

// Takes two descriptor matrices
std::vector<cv::DMatch> matchFLANN(cv::Mat desc1, cv::Mat desc2){
	FlannBasedMatcher matcher;
	std::vector<cv::DMatch> matches;
	matcher.match(desc1,desc2,matches)
	return matches;
}