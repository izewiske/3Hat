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


/*
 * Approx threshold is the threshold for "close-enough" that is used to allow a perfect left image match and a slightly-off
 * right image match to still supercede the "best-fit" match approximation.
 */
#define APPROX_THRESHOLD 10


Plane matchStrengths(cv::Mat mimg1, cv::Mat mimg2, cv::Mat bools1, cv::Mat bools2) {
	bool matchGlobalOrientations = true;
	//OUT <<"Running with matchGlobalOrientations = "<<matchGlobalOrientations<<" first."<<std::endl;

	// Make images as Mats; convert to IplImage for OpenSURF library actions
	cv::Mat mc1 = mimg1.clone();
	cv::Mat mc2 = mimg2.clone();


	IplImage iimg1, iimg2, bi1, bi2;
	iimg1=mc1;
	iimg2=mc2;
	bi1 = bools1;
	bi2 = bools2;

	IplImage *img1, *img2, *b1, *b2;
	img1 = &iimg1;
	img2 = &iimg2;
	b1 = &bi1;
	b2 = &bi2;

	IpVec ipts1, ipts2;
	surfDetDes(img1,ipts1,false,4,4,2,0.0001f,matchGlobalOrientations,b1);
	surfDetDes(img2,ipts2,false,4,4,2,0.0001f,matchGlobalOrientations,b2);

	MatchVec matches;
	getMatchesSymmetric(ipts1,ipts2,matches,true);

	IpVec mpts1, mpts2;
	Plane matchesVector;

	for (unsigned int i = 0; i < matches.size(); ++i)	{
		float strengthOverThreshold = 1 - matches[i].second; // /MATCH_THRESHOLD;
		strengthOverThreshold*=255;
		CvScalar clr = cvScalar(strengthOverThreshold,strengthOverThreshold,strengthOverThreshold);
		clr = cvScalar(255,255,255);
		
		//mpts1.push_back(matches[i].first.first);
		//mpts2.push_back(matches[i].first.second);
	
		//cvLine(img1,cvPoint(matches[i].first.first.x,matches[i].first.first.y),cvPoint(matches[i].first.second.x+w,matches[i].first.second.y), clr,1);
		//cvLine(img2,cvPoint(matches[i].first.first.x-w,matches[i].first.first.y),cvPoint(matches[i].first.second.x,matches[i].first.second.y), clr,1);

		matchesVector.leftImage.push_back(PixelLoc(matches[i].first.first.x,matches[i].first.first.y));
		matchesVector.rightImage.push_back(PixelLoc(matches[i].first.second.x,matches[i].first.second.y));
	}

	OUT << "Number of OpenSURF Matches: " << matches.size() << std::endl;
	return matchesVector;
}

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

		/*for (int w=0; w<listOfTilesR.size();w++){
			listOfTiles.push_back(listOfTilesR[w]);
		}
		// make unique list of tiles
		std::sort(listOfTiles.begin(), listOfTiles.end());
		listOfTiles.erase(std::unique(listOfTiles.begin(), listOfTiles.end()), listOfTiles.end());
		*/
		//start that iteration
		for(int k = 0; k < listOfTiles.size(); k++) {
			std::string tileID = listOfTiles[k];
			std::vector<PixelLoc> pixels1 = getContour(tileID,imageIDL);
			std::vector<PixelLoc> pixels2 = getContour(tileID,imageIDR);
			std::vector<cv::Point2f> contour1;
	        for (int j = 0; j < pixels1.size(); ++j) {
	        	cv::Point2f p1(pixels1[j].x,pixels1[j].y);
	        	contour1.push_back(p1);
	        }
			std::vector<cv::Point2f> contour2;
	        for (int j = 0; j < pixels2.size(); ++j) {
	       		cv::Point2f p2(pixels2[j].x,pixels2[j].y);
	        	contour2.push_back(p2);
	        }

	        cv::Rect roi1 = cv::boundingRect(contour1);
			cv::Mat slice1(image1,roi1);
			cv::Mat contourMatrix1 = slice1.clone();
			cv::Mat bools1 = locsToBool(pixels1,image1);
			cv::Mat sliceB1(bools1,roi1);
			cv::Mat contourBools1 = sliceB1.clone();
			cv::Mat contourOnly1 = contourMatrix1.mul(contourBools1);

			cv::Rect roi2 = cv::boundingRect(contour2);
			cv::Mat slice2(image2,roi2);
			cv::Mat contourMatrix2 = slice2.clone();
			cv::Mat bools2 = locsToBool(pixels2,image2);
			cv::Mat sliceB2(bools2,roi2);
			cv::Mat contourBools2 = sliceB2.clone();
			cv::Mat contourOnly2 = contourMatrix2.mul(contourBools2);
			// find feature points
			Plane surfMatches = matchStrengths(contourOnly1, contourOnly2, contourBools1, contourBools2);

			// get actual features from human input 
			Plane goldStandard = getUserDefinedPlane(tileID,imageID);

			// compare with stats
			OUT  << "Overall match quality: " << compareFeaturePoints(surfMatches,goldStandard) << ".\n";
			// matching fewer than three points is useless.
			if (surfMatches.leftImage.size() < 3 || surfMatches.rightImage.size() < 3) {
				OUT << "Tile: " << tileID << " does not have any strong SURF features. Consider alternative methods.\n";
				continue;
			}
			Plane best = bestPossibleComputedPlane(surfMatches,goldStandard);
			OUT <<"Best possible computed plane:";
			for (int j=0; j< best.leftImage.size(); j++){
				OUT  << "\t (" << best.leftImage[j].x << ", " << best.leftImage[j].y << ")->(" << best.rightImage[j].x << ", " << best.rightImage[j].y << ")";			
			}
			OUT <<"\nGold standard plane:";
			for (int j=0; j< goldStandard.leftImage.size(); j++){
				OUT  << "\t (" << goldStandard.leftImage[j].x << ", " << goldStandard.leftImage[j].y << ")->(" << goldStandard.rightImage[j].x << ", " << goldStandard.rightImage[j].y << ")";		
			}
		}
	}
	return 0;
}