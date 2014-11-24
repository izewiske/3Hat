#include <vector>
#include <string>
#include <stdlib.h>
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
#include "3Hat.h"

using namespace std;

Plane singleContourDriver(std::string tileID, std::string imageIDL, std::string imageIDR){
	try {
		Image im1( imageIDL.c_str());
		Image im2( imageIDR.c_str());
		cv::Mat image1(im1.getHeight(),im1.getWidth(),CV_8UC3,(void *) im1.getData());
		cv::Mat image2(im2.getHeight(),im2.getWidth(),CV_8UC3,(void *) im2.getData());
		if(! image1.data | !image2.data) {
			OUT  <<	"Could not open one or more images.\n";
			exit(EXIT_FAILURE);
		}
		std::vector<PixelLoc> pixels1 = getContour(tileID,imageIDL);
		std::vector<PixelLoc> pixels2 = getContour(tileID,imageIDR);
		if(pixels1.empty()) throw 4;
		if(pixels2.empty()) throw 4;
		//OUT << "Image: " << imageID<< "\n";
		//OUT <<"Tile ID: " << tileID << " #: " << k+1 << "\n";
			
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
			throw 2;
		}
		if(contour2.empty()){
			throw 2;
		}

	    cv::Rect roi1 = cv::boundingRect(contour1);
		cv::Mat slice1(image1,roi1);
		//OUT <<"Slice 1: "<<slice1.rows << " "<<slice1.cols << "\n";

		cv::Rect roi2 = cv::boundingRect(contour2);
		cv::Mat slice2(image2,roi2);
		//OUT <<"Slice 2: "<<slice2.rows << " "<<slice2.cols << "\n";
		if(slice1.rows <= _HESSIAN_THRESH || slice2.rows<=_HESSIAN_THRESH || slice1.cols <= _HESSIAN_THRESH || slice2.cols <= _HESSIAN_THRESH){
			//ERR << "Tile: " << tileID << " is too small for a Hessian matrix.\n";
	    	throw 3;
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
		//Plane goldStandard = getUserDefinedPlane(tileID,imageID);

		Plane surfMatches = matchStrengthsContour(true,contourOnly1, contourOnly2, contourBools1, contourBools2,true);
		// compare with stats
		//OUT  << "Overall (all enabled) match quality: " << compareFeaturePoints(surfMatches,goldStandard) << ".\n";
		
		// matching fewer than three points is useless.
		if (surfMatches.leftImage.size() < 3 || surfMatches.rightImage.size() < 3) {
			//OUT << "Tile: " << tileID << " does not have any strong SURF features. Consider alternative methods.\n";
			throw 5;
		}
		return surfMatches;
	} catch (int e){
		if(e==2) {
			ERR << "Tile: " << tileID << " is empty.\n";
		} else 
		if (e==3){
			ERR << "Tile: " << tileID << " is too small for a Hessian matrix.\n";
		} else 
		if (e==5){
			ERR << "Tile: " << tileID << " does not have enough strong SURF features.\n";
		}
	}	
}

int singleImageDriver(std::string imageID){
	try{
		OUT << "Image: " << imageID << "\n";
		std::string imageIDL = imageID + "L";
		std::string imageIDR = imageID + "R";
		//get list of tiles in image pair
		std::vector<std::string> listOfTiles;
		std::vector<std::string> listOfTilesL = getTileIDsForImage(imageIDL);
		std::vector<std::string> listOfTilesR = getTileIDsForImage(imageIDR);
		for (int w=0; w<listOfTilesR.size();w++){
			if(w%1000==0) ERR << w << "\n";
			if(std::find(listOfTilesL.begin(), listOfTilesL.end(), listOfTilesR[w]) != listOfTilesL.end()) {
		   		listOfTiles.push_back(listOfTilesR[w]);
			} else {
				//ERR << "No tiles found.\n";
			    continue;
			}
		}
		// iterate over set of contours that appear in both images
		for(int k = 0; k < listOfTiles.size(); k++) {
			std::string tileID = listOfTiles[k];
			OUT << "Tile: " << tileID << " ";
			singleContourDriver(tileID,imageIDL,imageIDR);
			if (k%20==19) OUT <<"\n";
		}
	} catch (int e) {
		ERR << "An unknown error occured. Probably with some negative in index in OpenCV because they haunt my dreams.\n";
		exit(EXIT_FAILURE);
	}
	return 0;
}