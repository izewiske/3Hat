//#include "ContourMatcher.h"
//#include "ContourMatcherStructs.h"
//#include "ContourMatcherExceptions.h"

#include <string>
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


//#include "../src/UtilityFunctions.h"


using namespace std;


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
cv::Mat sliceContour(std::vector<PixelLoc> contourPixels, cv::Mat image){
	cv::RotatedRect roi = defineROI(contourPixels);
	cv::Point2f vertices[4];
	roi.points(vertices);
	cv::Mat_<uchar> slice(image,vertices);
	return slice;
}


int main( int argc, char** argv ) {

		if( argc != 3) {
			// scottt100 1149L
		 cout <<" Usage: tileID image" << endl;
		 return -1;
		}

		string tileID = argv[1];
		string imageID = argv[2];

		cv::Mat image;
		image = cv::imread(argv[2], CV_LOAD_IMAGE_COLOR); 

		if(! image.data ) {
				cout <<	"Could not open or find the image" << std::endl ;
				return -1;
		}
		std::vector<PixelLoc> pixels = getContour(tileID,imageID);
		std::cerr<< "Number of pixels: "  << pixels.size() << endl;



		cv::Mat contour = sliceContour(pixels,image);
		cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE );
		cv::imshow( "Display window", contour );
		cv::waitKey(0);
		return 0;
}

