/*
*
*	This program will load multiple images and process them then determine how well they work
*
*/
#include <vector>
#include <string>
#include <iostream>

#include <opencv2/opencv.hpp>
#include "opencv2/highgui/highgui.hpp"
#include <opencv2/features2d/features2d.hpp>
#include "opencv2/imgproc/imgproc.hpp"

#ifdef USE_GL
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#endif

#include "surflib.h"
#include "eriolHeader.h"
#include "UtilityFunctions.h"
#include "ContourMatcherStructs.h"
#include "GetMatchesSURF.h"


/*
 * Approx threshold is the threshold for "close-enough" that is used to allow a perfect left image match and a slightly-off
 * right image match to still supercede the "best-fit" match approximation.
 */
#define APPROX_THRESHOLD 10
#define _HESSIAN_THRESH 15
// function finds the best possible plane we could have computed with a series of points 
Plane bestPossibleComputedPlane(Plane computedPoints, Plane actualPoints){ 
	// determines the best fitting plane computed - best possible scenario
	Plane bestComputedPlane;
	for (int j=0; j< actualPoints.leftImage.size(); j++){
		// for each computed point search for a closest match to actuality
		int closestMatchIndex = 0;
		float closestDist = 10000.0;
		for(int i=0;i< computedPoints.leftImage.size(); i++){
			float distL = distance(computedPoints.leftImage[i].x,computedPoints.leftImage[i].y,actualPoints.leftImage[j].x,actualPoints.leftImage[j].y);
			float distR = distance(computedPoints.rightImage[i].x,computedPoints.rightImage[i].y,actualPoints.rightImage[j].x,actualPoints.rightImage[j].y);
			// if the left image point is darn close and the right image point is darn close
			//if ( distL <= closestDistL + APPROX_THRESHOLD && distR <= closestDistR + APPROX_THRESHOLD ) {
			// replaced threshold with a notion of average distance from desired
			if ( ((distL + distR) / 2) < closestDist){
				closestDist = ((distL + distR) / 2);
				closestMatchIndex = i;
			}
		}
		bestComputedPlane.leftImage.push_back(PixelLoc(computedPoints.leftImage[closestMatchIndex].x,computedPoints.leftImage[closestMatchIndex].y));
		bestComputedPlane.rightImage.push_back(PixelLoc(computedPoints.rightImage[closestMatchIndex].x,computedPoints.rightImage[closestMatchIndex].y));
	}
	return bestComputedPlane;
}


// function returns a summative overall quality of all matches 
double compareFeaturePoints(Plane computedPoints, Plane actualPoints){
	// This determines how well the features we have computed match those chosen by humans
	float sumOfDistances = 0;
	// TODO: consider using RANSAC here
	for(int i=0;i< computedPoints.leftImage.size(); i++){
		// for each computed point search for a closest match to actuality
		int closestMatchIndex = 0;
		float closestDist = 10000.0;
		for(int j=0; j< actualPoints.leftImage.size(); j++){
			float distL = distance(computedPoints.leftImage[i].x,computedPoints.leftImage[i].y,actualPoints.leftImage[j].x,actualPoints.leftImage[j].y);
			float distR = distance(computedPoints.rightImage[i].x,computedPoints.rightImage[i].y,actualPoints.rightImage[j].x,actualPoints.rightImage[j].y);
			// if the left image point is darn close and the right image point is darn close
			//if ( distL <= closestDistL + APPROX_THRESHOLD && distR <= closestDistR + APPROX_THRESHOLD ) {
			// replaced threshold with a notion of average distance from desired
			if ( ((distL + distR) / 2) < closestDist){
				closestDist = ((distL + distR) / 2);
				closestMatchIndex = i;
			}
		}
		sumOfDistances = sumOfDistances + closestDist;
	}

	// for now quality is determined by the sum of differences.
	// TODO: consider counting number of matches with distance of under threshold
	double quality = sumOfDistances;
	return quality;
}

Plane getUserDefinedPlane(std::string tileID,std::string imageID){
	// TODO: I don't know if these are ordered in a fashion such that coresponding points have same index
	std::string imageIDL = imageID + "L";
	std::string imageIDR = imageID + "R";
	std::vector<Coord> leftStandard = getFeaturePoints(tileID,imageIDL);
	std::vector<Coord> rightStandard = getFeaturePoints(tileID,imageIDR);
	if (leftStandard.size() != rightStandard.size()){
		ERR << "That's odd, the number of points in the right and left images of the gold standard differs. You broke science.\n";
		exit;
	}
	Plane actual;
	for (int i = 0; i < leftStandard.size() ; ++i) {
		actual.leftImage.push_back(PixelLoc(leftStandard[i].x,leftStandard[i].y));
		actual.rightImage.push_back(PixelLoc(rightStandard[i].x,rightStandard[i].y));
	}
	return actual;
}


int main(int argc, char** argv){
	if( argc == 1) {
			// scottt100 1149L 1149R
		 	ERR <<" Usage: imageID1 imageID2 ... imageIDn" << std::endl;
		 	return -1;
		}
	//loop through images
	for (int i = 1; i < argc; i++){
		try {
		std::string imageID = argv[i];
		std::string imageIDL = imageID + "L";
		std::string imageIDR = imageID + "R";
		Image im1( imageIDL.c_str());
		Image im2( imageIDR.c_str());
		cv::Mat image1(im1.getHeight(),im1.getWidth(),CV_8UC3,(void *) im1.getData());
		cv::Mat image2(im2.getHeight(),im2.getWidth(),CV_8UC3,(void *) im2.getData());
		if(! image1.data | !image2.data) {
				OUT  <<	"Could not open one or more images.\n";
				return -1;
		}
		//get list of tiles in image pair
		std::vector<std::string> listOfTiles;
		std::vector<std::string> listOfTilesL = getTileIDsForImage(imageIDL);
		std::vector<std::string> listOfTilesR = getTileIDsForImage(imageIDR);
		for (int w=0; w<listOfTilesR.size();w++){
			if(std::find(listOfTilesL.begin(), listOfTilesL.end(), listOfTilesR[w]) != listOfTilesL.end()) {
	    		listOfTiles.push_back(listOfTilesR[w]);
			} else {
			    continue;
			}
		}
		// iterate over set of contours that appear in both images
		for(int k = 0; k < listOfTiles.size(); k++) {
			std::string tileID = listOfTiles[k];
			std::vector<PixelLoc> pixels1 = getContour(tileID,imageIDL);
			std::vector<PixelLoc> pixels2 = getContour(tileID,imageIDR);
			if(pixels1.empty()) continue;
			if(pixels2.empty()) continue;
			OUT << "Image: " << imageID<< "\n";
			OUT <<"Tile ID: " << tileID << " #: " << k+1 << "\n";

			std::vector<cv::Point2f> contour1;
	        for (int j = 0; j < pixels1.size(); ++j) {
				if ( pixels1[j].x > 0 && pixels1[j].y > 0 && pixels1[j].x < image1.cols && pixels1[j].y < image1.rows) {
	        		cv::Point2f p1(pixels1[j].x,pixels1[j].y);
	        		contour1.push_back(p1);
				}
	        }
			std::vector<cv::Point2f> contour2;
	        for (int j = 0; j < pixels2.size(); ++j) {
				if ( pixels2[j].x > 0 && pixels2[j].y > 0 && pixels2[j].x < image2.cols && pixels2[j].y < image2.rows) {
	       			cv::Point2f p2(pixels2[j].x,pixels2[j].y);
	        		contour2.push_back(p2);
	        	}
			}

			if(contour1.empty()){
				continue;
			}
			if(contour2.empty()){
				continue;
			}
	
	        cv::Rect roi1 = cv::boundingRect(contour1);
			cv::Mat slice1(image1,roi1);
			//OUT <<"Slice 1: "<<slice1.rows << " "<<slice1.cols << "\n";

			cv::Rect roi2 = cv::boundingRect(contour2);
			cv::Mat slice2(image2,roi2);
			//OUT <<"Slice 2: "<<slice2.rows << " "<<slice2.cols << "\n";
			if(slice1.rows <= _HESSIAN_THRESH || slice2.rows<=_HESSIAN_THRESH || slice1.cols <= _HESSIAN_THRESH || slice2.cols <= _HESSIAN_THRESH){
                                //ERR << "Tile: " << tileID << " is too small for a Hessian matrix.\n";
                                continue;
                            }		
			 			cv::Mat contourMatrix1 = slice1.clone();
                        //OUT << "Got to locsToBool\n";
                        cv::Mat bools1 = locsToBool(pixels1,image1);
                        //OUT << "Got passed locsToBool\n";
                        cv::Mat sliceB1(bools1,roi1);
                        cv::Mat contourBools1 = sliceB1.clone();
                        cv::Mat contourOnly1 = contourMatrix1.mul(contourBools1);


			cv::Mat contourMatrix2 = slice2.clone();
			cv::Mat bools2 = locsToBool(pixels2,image2);
			cv::Mat sliceB2(bools2,roi2);
			cv::Mat contourBools2 = sliceB2.clone();
			cv::Mat contourOnly2 = contourMatrix2.mul(contourBools2);
			// find feature points

			// get actual features from human input 
			Plane goldStandard = getUserDefinedPlane(tileID,imageID);

			Plane surfMatches = matchStrengthsContour(true,contourOnly1, contourOnly2, contourBools1, contourBools2,true);
			// compare with stats
			OUT  << "Overall (all enabled) match quality: " << compareFeaturePoints(surfMatches,goldStandard) << ".\n";
		
			// surfMatches = matchStrengthsSimpleBoundsInContour(true,contourOnly1, contourOnly2, contourBools1, contourBools2);
			// // compare with stats
			// OUT  << "Overall (no partial) match quality: " << compareFeaturePoints(surfMatches,goldStandard) << ".\n";
	
			// surfMatches = matchStrengthsSimpleBoundsInContour(false,contourOnly1, contourOnly2, contourBools1, contourBools2);
			// // compare with stats
			// OUT  << "Overall (local) match quality: " << compareFeaturePoints(surfMatches,goldStandard) << ".\n";

			// surfMatches = matchStrengthsSimpleBoundsInContour(false,contourMatrix1, contourMatrix2, contourBools1, contourBools2);
			// // compare with stats
			// OUT  << "Overall (local, full image) match quality: " << compareFeaturePoints(surfMatches,goldStandard) << ".\n";

			// matching fewer than three points is useless.
			if (surfMatches.leftImage.size() < 3 || surfMatches.rightImage.size() < 3) {
				//OUT << "Tile: " << tileID << " does not have any strong SURF features. Consider alternative methods.\n";
				continue;
			}
			Plane best = bestPossibleComputedPlane(surfMatches,goldStandard);
			//OUT <<"Best possible computed plane:";
			for (int j=0; j< best.leftImage.size(); j++){
				//OUT  << "\t (" << best.leftImage[j].x << ", " << best.leftImage[j].y << ")->(" << best.rightImage[j].x << ", " << best.rightImage[j].y << ")";			
			}
			//OUT <<"\nGold standard plane:";
			for (int j=0; j< goldStandard.leftImage.size(); j++){
				//OUT  << "\t (" << goldStandard.leftImage[j].x << ", " << goldStandard.leftImage[j].y << ")->(" << goldStandard.rightImage[j].x << ", " << goldStandard.rightImage[j].y << ")";		
			}
		}
		} catch (int e) {
			ERR << "An unknown error occured. Probably with some negative in index in OpenCV because they haunt my dreams.\n";
			continue;
		}
	}
	return 0;
}
