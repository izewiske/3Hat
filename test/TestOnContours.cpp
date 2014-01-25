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


// couldn't get main to work with functions so I flattened it.
int main( int argc, char** argv ) {

		if( argc != 3) {
			// scottt100 1149L
		 	std::cout <<" Usage: tileID image" << std::endl;
		 	return -1;
		}

		string tileID = argv[1];
		string imageID = argv[2];

		Image im(imageID.c_str());
		cv::Mat image(im.getHeight(),im.getWidth(),CV_8UC3,(void *) im.getData());
		if(! image.data ) {
				std::cout <<	"Could not open or find the image\n";
				return -1;
		}
		std::vector<PixelLoc> pixels = getContour(tileID,imageID);
		std::cerr<< "Number of pixels: "  << pixels.size() << std::endl;
	        std::vector<cv::Point2f> contour;
        	for (int j = 0; j < pixels.size(); ++j) {
        	        cv::Point2f p(pixels[j].x,pixels[j].y);
        	        contour.push_back(p);
        	}

	        cv::Rect roi = cv::boundingRect(contour);
		cv::Mat slice(image,roi);
		cv::Mat contourMatrix = slice.clone();



//		std::cout << "Contour = "<< std::endl << " "  << contourMatrix << std::endl << std::endl;
		cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE );
		cv::imshow( "Display window", contourMatrix );
		cv::waitKey(0);
		return 0;
}

