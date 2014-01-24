#ifndef __OPENCV_UTLITY_FUNCTIONS__
#define __OPENCV_UTLITY_FUNCTIONS__

#define _DISTANCE_THRESHHOLD 200
#define _WEIGHT_THRESHOLD 0.50

#include "ContourMatcherStructs.h"
#include "PixelLoc.h"

#include "Contour.h" //eriol header
#include "Image.h" //eriol header`

#include <unordered_map>
#include <math.h>

#include <opencv2/opencv.hpp>
#include <opencv2/features2d/features2d.hpp>
#include "opencv2/imgproc/imgproc.hpp"
#include <opencv2/nonfree/nonfree.hpp>
#include <opencv2/nonfree/features2d.hpp>
#include <opencv2/legacy/legacy.hpp>

cv::RotatedRect defineROI(std::vector<PixelLoc> contourPixels){
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
	return roi;
}

// Takes input from getContour which returns a vector of vector of PixelLoc
cv::Mat sliceContour(std::vector<PixelLoc> contourPixels,cv::Mat image){
	cv::RotatedRect roi = defineROI(contourPixels);
	cv::Point2f vertices[4];
	roi.points(vertices);
	cv::Mat_<uchar> slice(image,vertices);
	return slice;
}

cv::Mat convertImageToMatrix(Image im){
	cv::Mat image(im.getHeight(),im.getWidth(),CV_8UC3,(void *) im.getData());
	return image;
}

cv::Point getCenterOfROI(Contour contour){
	cv::RotatedRect	roi = defineROI(contour.getContourBoundary());
	return roi.center;
}

float distance(int x1,int y1, int x2, int y2){
	return sqrt(pow(x2-x1,2) + pow(y2-y1,2));
}

/*
	Called after initial feature detection run with SURF or SIFT to evaluate 
	the level of texture present in the contour. 
*/
bool contourIsDifficult(Contour contour,Image image){
	cv::Mat c = sliceContour(contour.getContourBoundary(),convertImageToMatrix(image));
	std::vector<cv::KeyPoint> features;
	cv::SiftFeatureDetector detector;
	detector.detect(c,features);

	// Processing features into map according to weight
	std::vector<float> weights;
	std::unordered_map<float,std::vector<Feature> > featureWeightMap;
	for (int i = 0; i < features.size(); i++) {
		cv::Point point = features[i].pt;
		float weight = features[i].response;
		int octave = features[i].octave;
		Feature f;
		f.x = point.x;
		f.y = point.y;
		f.weight = weight;
		f.octave = octave;
		weights.push_back(weight);
		featureWeightMap[weight].push_back(f);
	}
	cv::Point centerPt = getCenterOfROI(contour);

	bool isDifficult = false;
	// Determine the top number of weights we want to use
	// using all weights right now
	for (int i = 0; i < weights.size(); ++i) {
		
		
		// compare some attributs of features
		// If weight is weak then pop from list and don't bother calculating things
		if (weights[i] < _WEIGHT_THRESHOLD){
			//TODO: for speed consider not popping bad weights just skipping them
			featureWeightMap.erase(weights[i]);
			weights.erase(weights.begin()+i);
		}
		for( int j = 0; j<featureWeightMap[i].size(); j++){
			if ( distance(centerPt.x,centerPt.y,featureWeightMap[i][j].x,featureWeightMap[i][j].y) <= _DISTANCE_THRESHHOLD){
				// if features are closer to center then we want to weight them more because they are less likely to be occluded 
			}
		}
	}


	return isDifficult;
}

cv::Mat locsToBool(vector<PixelLoc> contourPixels, cv::Mat img){
	cv::Mat boolMat(img.rows,img.cols,img.type());
	for (int i=0; i<boolMat.cols; i++){
		for (int j=0; j<boolMat.rows; j++){
			//if PixelLoc(i,j) in contourPixels, 1, else 0
		}
	}
}

#endif //__OPENCV_UTLITY_FUNCTIONS__
