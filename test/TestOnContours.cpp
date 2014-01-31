//#include "ContourMatcher.h"
//#include "ContourMatcherStructs.h"
//#include "ContourMatcherExceptions.h"

#include <string>
#include <iostream>

#include "surflib.h"
#include "eriolHeader.h"
#include "GetMatchesSURF.h"

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

cv::Rect defineROI(std::vector<PixelLoc> contourPixels){
	std::vector<cv::Point2f> contour;
	//TODO: include boundary margins on region of interest
	for (int j = 0; j < contourPixels.size(); ++j) {
		// Makes a Pixel2f
		cv::Point2f p(contourPixels[j].x,contourPixels[j].y);
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
	cv::Rect roi = cv::boundingRect(contour);
	return roi;
}

// Takes input from getContour which returns a vector of vector of PixelLoc
void sliceContour(std::vector<PixelLoc> contourPixels, cv::Mat& image,cv::Mat& contour){
	cv::Rect roi = defineROI(contourPixels);
	std::cout<<"roi: " << roi << std::endl;
	cv::Mat slice(image,roi);
	std::cout << slice.data << std::endl;
	std::cout << "slice:" << slice << std::endl;
	contour = slice;
	std::cout << "slice:" << slice << std::endl;
	return ;
}

void convertImageToMatrix(Image im,cv::Mat& image){
	cv::Mat img(im.getHeight(),im.getWidth(),CV_8UC3,(void *) im.getData());
	image = img;
	return ;
}

//! Generate a boolean mat for an image 
// TODO: put this in a better place
cv::Mat locsToBool(vector<PixelLoc> contourPixels, cv::Mat img){
	cv::Mat boolMat(img.rows,img.cols,CV_8UC3,cvScalar(0,0,0));

        for (unsigned int i=0; i<contourPixels.size(); i++){
		for (int j=0; j<3; j++){
			boolMat.at<cv::Vec3b>(contourPixels[i].y,contourPixels[i].x)[j]=1;
		}
	}

	return boolMat;
}


//-------------------------------------------------------

int matchStrengths(cv::Mat mimg1, cv::Mat mimg2, cv::Mat cimg1, cv::Mat cimg2, cv::Mat bools1, cv::Mat bools2)
{
  bool matchGlobalOrientations = true;

  // Make images as Mats; convert to IplImage for OpenSURF library actions
  cv::Mat mc1 = cimg1.clone();
  cv::Mat mc2 = cimg2.clone();
  cv::Mat mc3 = cimg1.clone();
  cv::Mat mc4 = cimg2.clone();

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
  surfDetDes(img1,ipts1,false,4,4,2,0.0001f,!matchGlobalOrientations,b1);
  surfDetDes(img2,ipts2,false,4,4,2,0.0001f,!matchGlobalOrientations,b2);

  drawIpoints(img1, ipts1);
  drawIpoints(img2, ipts2);

  MatchVec matches;
  getMatchesSymmetric(ipts1,ipts2,matches,false);

  IpVec mpts1, mpts2;

  const int & w = img1->width;

  for (unsigned int i = 0; i < matches.size(); ++i)
  {
    float strengthOverThreshold = 1 - matches[i].second; // /MATCH_THRESHOLD;
    strengthOverThreshold*=255;
    CvScalar clr = cvScalar(strengthOverThreshold,strengthOverThreshold,strengthOverThreshold);
    clr = cvScalar(255,255,255);
    
    mpts1.push_back(matches[i].first.first);
    mpts2.push_back(matches[i].first.second);
  
    cvLine(img1,cvPoint(matches[i].first.first.x,matches[i].first.first.y),cvPoint(matches[i].first.second.x+w,matches[i].first.second.y), clr,1);
    cvLine(img2,cvPoint(matches[i].first.first.x-w,matches[i].first.first.y),cvPoint(matches[i].first.second.x,matches[i].first.second.y), clr,1);
  }

  drawIpoints(img1,mpts1);
  drawIpoints(img2,mpts2);

  std::cout<< "Matches: " << matches.size() << std::endl;

  cvNamedWindow("1", CV_WINDOW_AUTOSIZE );
  cvNamedWindow("2", CV_WINDOW_AUTOSIZE );
  cvShowImage("1", img1);
  cvShowImage("2",img2);
  cvWaitKey(0);


  // NOW DO IT AGAIN!
  IplImage iimg3, iimg4;
  iimg3=mc3;
  iimg4=mc4;

  IplImage *img3, *img4;
  img3 = &iimg3;
  img4 = &iimg4;

  IpVec ipts3, ipts4;
  surfDetDes(img3,ipts3,false,4,4,2,0.0001f,matchGlobalOrientations,b1);
  surfDetDes(img4,ipts4,false,4,4,2,0.0001f,matchGlobalOrientations,b2);

  drawIpoints(img3,ipts3);
  drawIpoints(img4,ipts4);

  matches.clear();
  getMatchesSymmetric(ipts3,ipts4,matches,false);

  IpVec mpts3, mpts4;

  for (unsigned int i = 0; i < matches.size(); ++i)
  {
    float strengthOverThreshold = 1 - matches[i].second; // /MATCH_THRESHOLD;
    strengthOverThreshold*=255;
    CvScalar clr = cvScalar(strengthOverThreshold,strengthOverThreshold,strengthOverThreshold);
    clr = cvScalar(255,255,255);
    
    mpts3.push_back(matches[i].first.first);
    mpts4.push_back(matches[i].first.second);
  
    cvLine(img3,cvPoint(matches[i].first.first.x,matches[i].first.first.y),cvPoint(matches[i].first.second.x+w,matches[i].first.second.y), clr,1);
    cvLine(img4,cvPoint(matches[i].first.first.x-w,matches[i].first.first.y),cvPoint(matches[i].first.second.x,matches[i].first.second.y), clr,1);
  }

  drawIpoints(img3,mpts3);
  drawIpoints(img4,mpts4);

  std::cout<< "Matches: " << matches.size() << std::endl;

  cvNamedWindow("3", CV_WINDOW_AUTOSIZE );
  cvNamedWindow("4", CV_WINDOW_AUTOSIZE );
  cvShowImage("3", img3);
  cvShowImage("4",img4);
  cvWaitKey(0);

  return 0;
}

//-------------------------------------------------------

// couldn't get main to work with functions so I flattened it.
int main( int argc, char** argv ) {

		if( argc != 4) {
			// scottt100 1149L 1149R
		 	std::cout <<" Usage: tileID image1 image2" << std::endl;
		 	return -1;
		}

		string tileID = argv[1];
		string imageID1 = argv[2];
		string imageID2 = argv[3];

		Image im1(imageID1.c_str());
		cv::Mat image1(im1.getHeight(),im1.getWidth(),CV_8UC3,(void *) im1.getData());
		if(! image1.data ) {
				std::cout <<	"Could not open or find the image (1)\n";
				return -1;
		}
		Image im2(imageID2.c_str());
		cv::Mat image2(im2.getHeight(),im2.getWidth(),CV_8UC3,(void *) im2.getData());
		if(! image2.data ) {
				std::cout <<	"Could not open or find the image (2)\n";
				return -1;
		}
		
		std::vector<PixelLoc> pixels1 = getContour(tileID,imageID1,true);
		std::cerr<< "Number of pixels (1): "  << pixels1.size() << std::endl;
	        std::vector<cv::Point2f> contour1;
        	for (int j = 0; j < pixels1.size(); ++j) {
        	        cv::Point2f p1(pixels1[j].x,pixels1[j].y);
        	        contour1.push_back(p1);
        	}
		std::vector<PixelLoc> pixels2 = getContour(tileID,imageID2,true);
		std::cerr<< "Number of pixels (2): "  << pixels2.size() << std::endl;
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

		//detect and match features using (modified) OpenSURF
		matchStrengths(contourMatrix1, contourMatrix2, contourOnly1, contourOnly2, contourBools1, contourBools2);

		return 0;
}

