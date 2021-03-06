// This is the file that will be linked into the final application
#include "OpenCVContourMatcher.h"

void OpenCVContourMatcher::loadImages(Image i1, Image i2){
	loadImages((unsigned char*) i1.getData(),i1.getHeight(),i1.getWidth(),(unsigned char*) i2.getData(),i2.getHeight(),i2.getWidth());
	return;
}


void OpenCVContourMatcher::loadImages(unsigned char* i1, unsigned int height1, unsigned int width1, unsigned char* i2, unsigned int height2, unsigned int width2){
	/*	
		!!! This method of constructing cv::Mat objects requires that the 
		image arrays pointed to are managed seperately!!!
	*/
		try{
			if (height1 != height2 | width1 != width2){
				throw OCVMImageSize;
			}
			image1 = cv::Mat(height1,width1,CV_8UC3,(void*) &i1);
			image2 = cv::Mat(height2,width2,CV_8UC3,(void*) &i2);
		} catch (std::exception& e) {
			ERR << e.what() ;
		}
	return;
}


// Takes input from getContour which returns a //vector of vector of PixelLoc
cv::Mat OpenCVContourMatcher::sliceContour(std::vector<PixelLoc> contourPixels,cv::Mat image){
	std::vector<cv::Point> contour;
	//TODO: include boundary margins on region of interest
	for (int j = 0; j < contourPixels.size(); ++j) {
		// Makes a Pixel2f
		cv::Point p(contourPixels[j].x,contourPixels[j].y);
		//adds it to contour
		contour.push_back(p);
	}

	/* 
	From OpenCV on fitEllipse():
	The function calculates the ellipse that fits (in a least-squares sense) a set of 2D points best of all. 
	It returns the rotated rectangle in which the ellipse is inscribed. The algorithm [Fitzgibbon95] is used. 

	!!!
		NOTE: Developer should keep in mind that it is possible that the returned ellipse/rotatedRect data contains 
		negative indices, due to the data points being close to the border of the containing Mat element.
	!!!
	*/
	cv::RotatedRect roi = cv::fitEllipse(contour);
	cv::Point2f vertices[4];
	roi.points(vertices);
	cv::Mat_<uchar> slice(image,vertices);
	return slice;
}

Plane OpenCVContourMatcher::convertDMatchesToPlane(std::vector<cv::KeyPoint> features, std::vector<cv::KeyPoint> features2, std::vector<cv::DMatch> matches){
	Plane plane;
	//	std::vector<PixelLoc> planePoints;
	// convert DMatches to tbd set of points
	// Gets the 
	for (int i = 0; i < matches.size(); i++){
		// queryIdx is in the "left" image 
		PixelLoc pt1(features[matches[i].queryIdx].pt.x,features[matches[i].queryIdx].pt.y);
		plane.leftImage.push_back(pt1);
		// trainIdx is in the "right" image
		PixelLoc pt2(features2[matches[i].trainIdx].pt.x,features2[matches[i].trainIdx].pt.y);
		plane.rightImage.push_back(pt2);
	}
	return plane;
}


std::vector<cv::KeyPoint> OpenCVContourMatcher::detectFeaturePoints(cv::Mat contour){
	cv::SiftFeatureDetector detector;
	std::vector<cv::KeyPoint> contourFeatures;
	detector.detect(contour,contourFeatures);
	return contourFeatures;
}


cv::Mat OpenCVContourMatcher::extractDescriptors(cv::Mat contour, std::vector<cv::KeyPoint> contourFeatures){
	cv::SiftDescriptorExtractor extractor;
	cv::Mat contourDescriptors;
	extractor.compute(contour,contourFeatures,contourDescriptors);
	return contourDescriptors;
}

Plane OpenCVContourMatcher::compare(Contour contour1, Contour contour2, MATCHER_TYPE m = FLANN){
	Plane p;
	try {
		if ( m == FLANN) {
			p = compareFLANN(contour1, contour2);
		} else if ( m == BF ){
			p = compareBruteForce(contour1, contour2);
		} else {
			throw OCVMMatcherNotDef;
		}
	} catch (std::exception& e) {
		ERR << e.what() ;
	}
	return p;
}

Plane OpenCVContourMatcher::compareFLANN(Contour contour1, Contour contour2){
	//Convert from contour class to matrix
	cv::Mat c1 = sliceContour(contour1.getContourBoundary(),image1);
	cv::Mat c2 = sliceContour(contour2.getContourBoundary(),image2);
	// Detect feature points of contour 1 and contour 2
	std::vector<cv::KeyPoint> feat1 = detectFeaturePoints(c1);
	std::vector<cv::KeyPoint> feat2 = detectFeaturePoints(c2);
	// Extract feature descriptors
	cv::Mat desc1 = extractDescriptors(c1,feat1);
	cv::Mat desc2 = extractDescriptors(c2,feat2);
	// Match feature descriptors and refine matches then return as vector of PixelLocations
	std::vector<cv::DMatch> matches = matchFLANN(desc1,desc2);
	Plane planePoints = convertDMatchesToPlane(feat1,feat2,refineMatches(matches));
	return planePoints;
}

Plane OpenCVContourMatcher::compareBruteForce(Contour contour1, Contour contour2){
	//Convert from contour class to matrix
	cv::Mat c1 = sliceContour(contour1.getContourBoundary(),image1);
	cv::Mat c2 = sliceContour(contour2.getContourBoundary(),image2);
	// Detect feature points of contour 1 and contour 2
	std::vector<cv::KeyPoint> feat1 = detectFeaturePoints(c1);
	std::vector<cv::KeyPoint> feat2 = detectFeaturePoints(c2);
	// Extract feature descriptors
	cv::Mat desc1 = extractDescriptors(c1,feat1);
	cv::Mat desc2 = extractDescriptors(c2,feat2);
	// Match feature descriptors and refine matches then return as vector of PixelLocations
	std::vector<cv::DMatch> matches = matchBruteForce(desc1,desc2);
	Plane planePoints = convertDMatchesToPlane(feat1,feat2,refineMatches(matches));
	return planePoints;
}

void OpenCVContourMatcher::computeMaxAndMinDistances(std::vector<cv::DMatch> matches){
	minDistance = 0;
	maxDistance = 100;
	for (int i = 0; i < matches.size(); i++){
		double dist = matches[i].distance;
		if ( dist < minDistance ) minDistance = dist;
		if ( dist > maxDistance ) maxDistance = dist;
	}
	return;
}

std::vector<cv::DMatch> OpenCVContourMatcher::refineMatches(std::vector<cv::DMatch> matches){
	std::vector<cv::DMatch> bestMatches;
	// Determine how close to their expected locations the features 
	computeMaxAndMinDistances(matches);
	for (int i = 0; i < matches.size(); i++ ){
		if (matches[i].distance < minDistance * THRESHOLD){
			bestMatches.push_back(matches[i]);
		}
	}
	return bestMatches;
}

// Takes two descriptor matrices
std::vector<cv::DMatch> OpenCVContourMatcher::matchBruteForce(cv::Mat desc1, cv::Mat desc2){
	cv::BFMatcher matcher;
	std::vector<cv::DMatch> matches;
	matcher.match(desc1,desc2,matches);
	return matches;
}

// Takes two descriptor matrices
std::vector<cv::DMatch> OpenCVContourMatcher::matchFLANN(cv::Mat desc1, cv::Mat desc2){
	cv::FlannBasedMatcher matcher;
	std::vector<cv::DMatch> matches;
	matcher.match(desc1,desc2,matches);
	return matches;
}

