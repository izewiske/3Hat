#include "ContourMatcherStructs.h"

#include <unordered_map>
#include <math.h>

#ifndef __UTILITY_FUNCTIONS__
#define __UTILITY_FUNCTIONS__


// Commonly mentioned in literature is a JND of 1.0
// Mahy et al. (1994) assessed a JND of 2.3 
#define _COLOUR_DIFFERENCE_THRESHOLD 2.3
#define _ACCEPTABLE_PERCENTAGE_DEVIANT 0.5

// approximation of pi for use in CIEDE 2000 calculations for color difference
#define pi 3.141592653589793238462643383279

#ifndef __OPENCV_UTLITY_FUNCTIONS__
#define __OPENCV_UTLITY_FUNCTIONS__

#define _DISTANCE_THRESHHOLD 200
#define _WEIGHT_THRESHOLD 0.50

#include <opencv2/opencv.hpp>
#include <opencv2/features2d/features2d.hpp>
#include "opencv2/imgproc/imgproc.hpp"
#include <opencv2/legacy/legacy.hpp>

#ifdef NONFREE_ENABLED
#include <opencv2/nonfree/nonfree.hpp>
#include <opencv2/nonfree/features2d.hpp>
#endif

float distance(int x1,int y1, int x2, int y2){
	return sqrt(pow(x2-x1,2) + pow(y2-y1,2));
}

bool pixelNeighbors(int x, int y, cv::Mat contour){
	int numNeighbors = 0;
	if (x == 828 && y == 540 ){
			OUT << "Contour - xrange: " << 0 << "-" << contour.cols << " yrange: "<<0<<"-"<<contour.rows << "\n";
		}
	if (x >= 1 && x <= contour.cols -1 ){
		// to the left 
		if (contour.at<cv::Vec3b>(x-1,y)[0]==1 && contour.at<cv::Vec3b>(x-1,y)[1]==1 && contour.at<cv::Vec3b>(x-1,y)[2]==1 ){
			numNeighbors++;
		}
		// to the right
		if (contour.at<cv::Vec3b>(x+1,y)[0]==1 && contour.at<cv::Vec3b>(x+1,y)[1]==1 && contour.at<cv::Vec3b>(x-1,y)[2]==1 ){
			numNeighbors++;
		}
	}
	if (x == 828 && y == 540 ){
		OUT << "Left and right.\n";
	}
	if (y >= 1 && y <= contour.cols -1 ){
		// above
		if (contour.at<cv::Vec3b>(x,y+1)[0]==1 && contour.at<cv::Vec3b>(x,y+1)[1]==1 && contour.at<cv::Vec3b>(x,y+1)[2]==1 ){
			numNeighbors++;
		}
		if (x == 828 && y == 540 ){
			OUT << "Above.\n";
		}
		// below 
		if (contour.at<cv::Vec3b>(x,y-1)[0]==1 && contour.at<cv::Vec3b>(x,y-1)[1]==1 && contour.at<cv::Vec3b>(x,y-1)[2]==1 ){
			numNeighbors++;
		}
		if (x == 828 && y == 540 ){
			OUT << "Below.\n";
		}
	}
	if (numNeighbors >= 2){
		return true;
	} else {
		return false;
	}
}

cv::Mat locsToBool(vector<PixelLoc> contourPixels, cv::Mat img, int pixelBuffer = 3){
	cv::Mat boolMat(img.rows,img.cols,CV_8UC3,cvScalar(0,0,0));


	for (unsigned int i=0; i<contourPixels.size(); i++){
		for(int k=-pixelBuffer; k<pixelBuffer; k++) {
			for(int h=-pixelBuffer; h<pixelBuffer;h++){
				for (int j=0; j<3; j++){
					if (contourPixels[i].y+h >= 0 && contourPixels[i].y+h < img.rows && contourPixels[i].x+k >= 0 && contourPixels[i].x+k < img.rows) {
						boolMat.at<cv::Vec3b>(contourPixels[i].y+h,contourPixels[i].x+k)[j]=1;
					}
				}
			}
		}
	}
	/*
	 *	More preparation of the boolean contour array is possible but not implemented at this time.
	 */ 

	/*
	// iterate over boolMat and find pits
	for(int i=x1-pixelBuffer; i<= x2+pixelBuffer; i++){
		for(int j = y1-pixelBuffer; j <= y2+pixelBuffer; j++){
			// if it's on the edge of matrix check for invalid memory
			if (boolMat.cols > i+1 && boolMat.rows > j+1 && i >= 1 && j >= 1) {
				// check if pixel is already part of contour
				if ( boolMat.at<cv::Vec3b>(i,j)[0]!=1 && boolMat.at<cv::Vec3b>(i,j)[1]!=1 && boolMat.at<cv::Vec3b>(i,j)[2]!=1) {
					// if it isn't then check that it has 2 or more neighbors
					if (pixelNeighbors(i,j,boolMat)) {
						// if it has 2 or more neighbors then it belongs in the contour
						boolMat.at<cv::Vec3b>(contourPixels[i].y,contourPixels[i].x)[j]=1;
					}
				}
			}
		}
	}
	*/

	return boolMat;
}

cv::Point getCenterOfROI(cv::Rect roi){
	cv::Point center(roi.x + roi.width/2, roi.y + roi.height/2);
	return center;
}


// function declaration to use non-openCV function in the openCV utility function section
double calculateColorDifference(double * lab1, double * lab2);

// Determines if a contour has uniform color by calculating the average LAB value for the contour and then 
// determining the percent of pixels that deviate strongly from that average
bool contourHasUniformColor(std::string tileID,std::string image){
	Image imgPrime(image.c_str());
	cv::Mat img(imgPrime.getHeight(),imgPrime.getWidth(),CV_8UC3,(void *) imgPrime.getData());
	if(! img.data ) {
				ERR <<	"Could not open or find the image (1)\n";
				return -1;
	}

	std::vector<PixelLoc> pixels = getContour(tileID,image);
	std::vector<cv::Point2f> contour;
    for (int j = 0; j < pixels.size(); ++j) {
        cv::Point2f p(pixels[j].x,pixels[j].y);
        contour.push_back(p);
    }
 	cv::Rect roi = cv::boundingRect(contour);
	cv::Mat slice(img,roi);
	cv::Mat contourMatrix = slice.clone();
	cv::Mat bools = locsToBool(pixels,img);
	cv::Mat sliceB(bools,roi);
	cv::Mat contourBools = sliceB.clone();
	cv::Mat contourOnly = contourMatrix.mul(contourBools);
	cv::Mat c;
	cv::cvtColor(contourOnly, c, CV_BGR2Lab);
	//cv::cvtColor(src, c, CV_RGB2HLS);
	int step = c.step;
	int channels = c.channels();
	double contourAvgLAB[3] = {0.0,0.0,0.0};


	for (int i = 0; i < c.rows; i++) {
	    for (int j = 0; j < c.cols; j++) {
	    	double pixelLAB[3] = {c.data[step*i + channels*j + 0], c.data[step*i + channels*j + 1], c.data[step*i + channels*j + 2]};
	    	contourAvgLAB[0] = pixelLAB[0] + contourAvgLAB[0];
	  		contourAvgLAB[1] = pixelLAB[1] + contourAvgLAB[1];
	  		contourAvgLAB[2] = pixelLAB[2] + contourAvgLAB[2];
    	}
	}
	contourAvgLAB[0] = contourAvgLAB[0]/(c.rows * c.cols);
	contourAvgLAB[1] = contourAvgLAB[1]/(c.rows * c.cols);
	contourAvgLAB[2] = contourAvgLAB[2]/(c.rows * c.cols);

	int totalDeviantPixels = 0;
	#pragma omp parallel for reduction (+:totalDeviantPixels){
	for (int i = 0; i < c.rows; i++) {
	    for (int j = 0; j < c.cols; j++) {
	    	double pixelLAB[3] = {c.data[step*i + channels*j + 0], c.data[step*i + channels*j + 1], c.data[step*i + channels*j + 2]};
	    	double contourDifference = calculateColorDifference(pixelLAB,contourAvgLAB);
	    	if (contourDifference >= _COLOUR_DIFFERENCE_THRESHOLD){
	    		totalDeviantPixels++;
	    	}
    	}
	}
	if (totalDeviantPixels/(c.cols*c.rows) >= _ACCEPTABLE_PERCENTAGE_DEVIANT) {
		return false;
	} else {
		return true;
	}
}


/*
	Called after initial feature detection run with SURF or SIFT to evaluate 
	the level of texture present in the contour. 
*/
bool contourIsDifficult(std::string tileID,std::string image){
	Image imgPrime(image.c_str());
	cv::Mat img(imgPrime.getHeight(),imgPrime.getWidth(),CV_8UC3,(void *) imgPrime.getData());
	if(! img.data ) {
		ERR <<	"Could not open or find the image (1)\n";
		return -1;
	}

	std::vector<PixelLoc> pixels = getContour(tileID,image);
	std::vector<cv::Point2f> contour;
    for (int j = 0; j < pixels.size(); ++j) {
        cv::Point2f p(pixels[j].x,pixels[j].y);
        contour.push_back(p);
    }
 	cv::Rect roi = cv::boundingRect(contour);
	cv::Mat slice(img,roi);
	cv::Mat contourMatrix = slice.clone();
	cv::Mat bools = locsToBool(pixels,img);
	cv::Mat sliceB(bools,roi);
	cv::Mat contourBools = sliceB.clone();
	cv::Mat c = contourMatrix.mul(contourBools);
	std::vector<cv::KeyPoint> features;
#ifdef NONFREE_ENABLED
	ERR << "Nonfree is enabled, using the Sift feature detector.\n";
	cv::SiftFeatureDetector detector;
#else 
	ERR << "Nonfree is disabled, using the STAR feature detector.\n";
	cv::StarFeatureDetector detector;
#endif 
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
	cv::Point centerPt = getCenterOfROI(roi);

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

#endif //__OPENCV_UTLITY_FUNCTIONS__

double calculateColorDifference(double * lab1, double * lab2) {
	double L1 = lab1[0];
	double a1 = lab1[1];
	double b1 = lab1[2];

	double L2 = lab2[0];
	double a2 = lab2[1];
	double b2 = lab2[2];

	// Calculate the magnitude of a*b
	double magab1 = sqrt(pow(a1,2)+pow(b1,2));
	double magab2 = sqrt(pow(a2,2)+pow(b2,2));

	// Calculate arithmetic mean 
	double mag_mean = (magab1 + magab2)/2.0;

	// Intermediary value used to calculate a_i'
	double G = 0.5*( 1.0 - sqrt( pow(mag_mean,7.0)/(pow(mag_mean,7.0) + pow(25.0,7.0))));

	double aPrime1 = (1.0+G)*a1;
	double aPrime2 = (1.0+G)*a2;

	// Calculate magnitude of a_i'
	double mag_aPrime1 = sqrt(pow(aPrime1,2)+pow(b1,2));
	double mag_aPrime2 = sqrt(pow(aPrime2,2)+pow(b2,2));

	// Combined color
	double colorProd = mag_aPrime1 * mag_aPrime2;
	// Hue for colors mod 2pi
	double huePrime1 = atan2(b1,aPrime1);
	if (huePrime1 < 0) {
		huePrime1 += 2.0 * pi;
	}
	double huePrime2 = atan2(b2,aPrime2);
	if (huePrime2<0) { 
		huePrime2+= 2.0*pi;
	}
	// Absolute value of magnitude of aPrime2 and b2
	if ((fabs(mag_aPrime2) + fabs(b2)) == 0.0) {
		huePrime2= 0.0;	
	} 

	// calculate intermediary values
	double luminosityDifference = L2-L1;
	double primeDifference = aPrime2-aPrime1;
	double huePrimeDifference = huePrime2-huePrime1;

	// Make sure hueDifference and colorDifference conform to algebraic contraints from paper
	if (huePrimeDifference > pi)  {
		huePrimeDifference -= 2.0 * pi;
	} else if (huePrimeDifference < -pi) {
		huePrimeDifference += 2.0 * pi;
	}
	if (colorProd == 0.0) { 
		huePrimeDifference = 0.0;
	}

	double hueDifference = 2.0 * sqrt(colorProd) * sin(huePrimeDifference/2.0);
	double averageLuminosity= (L2 + L1)/2.0;
	double averageAPrime= (aPrime1 + aPrime2)/2.0;

	// Calculate average hue - huePrime and conform to algebraic constraints
	double huePrime = (huePrime2 + huePrime1)/2.0;
	if ( fabs(huePrime1-huePrime2)  > pi ) {
		huePrime -= pi;
	}
	if (huePrime <0) { 
		huePrime += 2.0*pi;
	}
	if (colorProd==0.0) {
		huePrime = huePrime2 + huePrime1;
	}

	// Primary calculations that are composed to calculate difference - variable names as seen in paper
	double S_l= 1.0+0.015 * (averageLuminosity-50.0)*(averageLuminosity-50.0)/sqrt(20.0+(averageLuminosity-50.0)*(averageLuminosity-50.0));
	double S_c= 1.0+0.045*averageAPrime;
	double T= 1.0 - 0.17 * cos(huePrime - pi/6.0) + 0.24*cos(2.0*huePrime) + 0.32*cos(3.0*huePrime+pi/30.0) - 0.20*cos(4.0*huePrime-63.0*pi/180.0);
	double S_h= 1.0 + 0.015 * averageAPrime * T;
	double R_c=  2.0*sqrt(pow(averageAPrime,7.0)/(pow(averageAPrime,7.0) + pow(25.0,7.0)));
	double R_T= -sin(2.0*(30.0*pi/180.0)*exp(- pow(( (180.0/pi*huePrime-275.0)/25.0),2.0)))*R_c;

	// The CIEDE 2000 color difference (double)
	return sqrt( pow((luminosityDifference/S_l),2.0) + pow((primeDifference/S_c),2.0) + pow((hueDifference/S_h),2.0) + R_T*(primeDifference/S_c)*(hueDifference/S_h) );
}


#endif //__UTILITY_FUNCTIONS__
