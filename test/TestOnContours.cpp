//#include "ContourMatcher.h"
//#include "ContourMatcherStructs.h"
//#include "ContourMatcherExceptions.h"

#include <string>
#include <iostream>

#include "surflib.h"
#include "eriolHeader.h"

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

//-------------------------------------------------------

int matchStrengths(cv::Mat mimg1, cv::Mat mimg2)
{
  bool matchGlobalOrientations = true;

  // Make images as Mats; convert to IplImage for OpenSURF library actions

  IplImage iimg1, iimg2;
  iimg1=mimg1;
  iimg2=mimg2;

  IplImage *img1, *img2;
  img1 = &iimg1;
  img2 = &iimg2;

  IpVec ipts1, ipts2;
  surfDetDes(img1,ipts1,false,4,4,2,0.0001f,matchGlobalOrientations);
  surfDetDes(img2,ipts2,false,4,4,2,0.0001f,matchGlobalOrientations);

  MatchVec matches;
  getMatchesSymmetric(ipts1,ipts2,matches);

  IpVec mpts1, mpts2;

  const int & w = img1->width;

  for (unsigned int i = 0; i < matches.size(); ++i)
  {
    float strengthOverThreshold = 1 - matches[i].second; // /MATCH_THRESHOLD;
    strengthOverThreshold*=255;
    CvScalar clr = cvScalar(strengthOverThreshold,strengthOverThreshold,strengthOverThreshold);
    clr = cvScalar(255,255,255);
    
    //drawPoint(img1,matches[i].first.first,clr);
    //drawPoint(img2,matches[i].first.second,clr),
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
  cv::Mat mimg3, mimg4;
  mimg3=cv::imread("OpenSURF/imgs/img1.jpg", CV_LOAD_IMAGE_COLOR);
  mimg4=cv::imread("OpenSURF/imgs/img2.jpg", CV_LOAD_IMAGE_COLOR);

  IplImage iimg3, iimg4;
  iimg3=mimg3;
  iimg4=mimg4;

  IplImage *img3, *img4;
  img3 = &iimg3;
  img4 = &iimg4;

  IpVec ipts3, ipts4;
  surfDetDes(img3,ipts3,false,4,4,2,0.0001f,!matchGlobalOrientations);
  surfDetDes(img4,ipts4,false,4,4,2,0.0001f,!matchGlobalOrientations);

  matches.clear();
  getMatchesSymmetric(ipts3,ipts4,matches);

  IpVec mpts3, mpts4;

  for (unsigned int i = 0; i < matches.size(); ++i)
  {
    float strengthOverThreshold = 1 - matches[i].second; // /MATCH_THRESHOLD;
    strengthOverThreshold*=255;
    CvScalar clr = cvScalar(strengthOverThreshold,strengthOverThreshold,strengthOverThreshold);
    clr = cvScalar(255,255,255);
    
    //drawPoint(img1,matches[i].first.first,clr);
    //drawPoint(img2,matches[i].first.second,clr),
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

		if( argc != 3) {
			// scottt100 1149L
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
		
		std::vector<PixelLoc> pixels1 = getContour(tileID,imageID1);
		std::cerr<< "Number of pixels (1): "  << pixels1.size() << std::endl;
	        std::vector<cv::Point2f> contour1;
        	for (int j = 0; j < pixels1.size(); ++j) {
        	        cv::Point2f p1(pixels1[j].x,pixels1[j].y);
        	        contour1.push_back(p1);
        	}
		std::vector<PixelLoc> pixels2 = getContour(tileID,imageID2);
		std::cerr<< "Number of pixels (2): "  << pixels2.size() << std::endl;
	        std::vector<cv::Point2f> contour2;
        	for (int j = 0; j < pixels2.size(); ++j) {
        	        cv::Point2f p2(pixels2[j].x,pixels2[j].y);
        	        contour2.push_back(p2);
        	}

	        cv::Rect roi1 = cv::boundingRect(contour1);
		cv::Mat slice1(image1,roi1);
		cv::Mat contourMatrix1 = slice1.clone();

	        cv::Rect roi2 = cv::boundingRect(contour2);
		cv::Mat slice2(image2,roi2);
		cv::Mat contourMatrix2 = slice2.clone();

		//detect and match features using (modified) OpenSURF
		matchStrengths(contourMatrix1, contourMatrix2);

		return 0;
}

