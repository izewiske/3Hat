/*
*
*	This program will load multiple images and process them then determine how well the 
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



std::vector<string> getTileIDS(Image im);
Plane getFeaturesFromTile(std::string tildID);




Plane matchStrengths(cv::Mat mimg1, cv::Mat mimg2, cv::Mat bools1, cv::Mat bools2) {
	bool matchGlobalOrientations = true;
	std::cout<<"Running with matchGlobalOrientations = "<<matchGlobalOrientations<<" first."<<std::endl;

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

	drawIpoints(img1, ipts1);
	drawIpoints(img2, ipts2);

	MatchVec matches;
	getMatchesSymmetric(ipts1,ipts2,matches,true);

	IpVec mpts1, mpts2;

	const int & w = img1->width;

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

	std::cout<< "OpenSURF Matches: " << matches.size() << std::endl;
	return matchesVector;
}

// function finds the best possible plane we could have computed with a series of points 
Plane bestPossibleComputedPlane(Plane computedPoints, Plane actualPoints){ 
	// determines the best fitting plane computed - best possible scenario
	Plane bestComputedPlane;
	for (int j=0; j< actualPoints.size(); j++){
		// for each computed point search for a closest match to actuality
		int closestMatchIndex = 0;
		float closestDist = 10000.0;
		for(int i=0;i< computedPoints.size(); i++){
			float dist = distance(computedPoints[i].x,computedPoints[i].y,actualPoints[j].x,actualPoints[j].y);
			if ( dist < closestDist ) {
				closestDist = dist;
				closestMatchIndex = i;
			}
		}
		bestComputedPlane.push_back(PixelLoc(computedPoints[closestMatchIndex].x,computedPoints[closestMatchIndex].y));
	}
	return bestComputedPlane;
}


// function returns a summative overall quality of all matches 
double compareFeaturePoints(Plane computedPoints, Plane actualPoints){
	// This determines how well the features we have computed match those chosen by humans
	float sumOfDistances = 0;
	// TODO: consider using RANSAC here
	for(int i=0;i< computedPoints.size(); i++){
		// for each computed point search for a closest match to actuality
		int closestMatchIndex = 0;
		float closestDist = 10000.0;
		for (int j=0; j< actualPoints.size(); j++){
			float dist = distance(computedPoints[i].x,computedPoints[i].y,actualPoints[j].x,actualPoints[j].y);
			if (dist < closestDist){
				closestDist = dist;
				closestMatchIndex = j;
			}
		}
		sumOfDistances = sumOfDistances + closestDist;
	}

	// for now quality is determined by the sum of differences.
	// TODO: consider counting number of matches with distance of under threshold
	double quality = sumOfDistances;
	return quality;
}


int main(int argc, char** argv){
	if( argc == 1) {
			// scottt100 1149L 1149R
		 	std::err <<" Usage: image1ID image2ID ... " << std::endl;
		 	return -1;
		}
	//loop through images
	for (int i = 0; i < argc; i++){
		std::string imageID = argv[i];
		Image im1(imageID + "L");
		Image im2(imageID + "R");
		cv::Mat image1(im1.getHeight(),im1.getWidth(),CV_8UC3,(void *) im1.getData());
		cv::Mat image2(im2.getHeight(),im2.getWidth(),CV_8UC3,(void *) im2.getData());
		if(! image1.data | !image2.data) {
				std::cout <<	"Could not open one or more images.\n";
				return -1;
		}
		//get list of tiles in image pair
		std::vector<std::string> listOfTiles = getTileIDS(imageID);
		for(int k = 0; k < listOfTiles.size(); k++) {
			std::string tileID = listOfTiles[k];
			std::vector<PixelLoc> pixels1 = getContour(tileID,imageID+"L");
			std::vector<PixelLoc> pixels2 = getContour(tileID,imageID+"R");
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
			Plane goldStandard = getFeaturesFromTile(tileID);
			// compare with stats
			std::cout << "Overall match quality: " << compareFeaturePoints(surfMatches,goldStandard) << ".\n";
			Plane best = bestPossibleComputedPlane(surfMatches,goldStandard);
			std::cout<<"Best possible computed plane:" 
			for (int j=0; j< best.size(); j++){
				std::cout << "\t (" << best[j].x << ", " << best[j].y << ")";			
			}
			std::cout<<"\nGold standard plane:" 
			for (int j=0; j< goldStandard.size(); j++){
				std::cout << "\t (" << goldStandard[j].x << ", " << goldStandard[j].y << ")";			
			}
		}
	}
	return 0;
}