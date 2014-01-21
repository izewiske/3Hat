//#include "ContourMatcher.h"
//#include "ContourMatcherStructs.h"
//#include "ContourMatcherExceptions.h"

#include <string.h>
#include "eriolHeader.h"
#include <opencv2/opencv.hpp>
#include "opencv2/highgui/highgui.hpp"
#include <opencv2/features2d/features2d.hpp>
#include "opencv2/imgproc/imgproc.hpp"

using namespace std;
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


// int main( int argc, char** argv ) {

// 		if( argc != 3) {
// 			// scottt100 1149L
// 		 cout <<" Usage: tileID image" << endl;
// 		 return -1;
// 		}

// 		cv::Mat image;
// 		image = cv::imread(argv[2], CV_LOAD_IMAGE_COLOR); 

// 		if(! image.data ) {
// 				cout <<	"Could not open or find the image" << std::endl ;
// 				return -1;
// 		}

// 		cv::Mat contour = sliceContour(getContour(argv[1],argv[2]),image);
// 		cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE );
// 		cv::imshow( "Display window", contour );
// 		cv::waitKey(0);
// 		return 0;
// }

#ifdef USE_GL
void commonKeyboard( unsigned char c, int x, int y, int windowIndex)
{
	switch( c ) {
		case 'w':
			resetZoom(windowIndex); break;
		case 3: case 27:
			if ( -1 != myWindowID[1] ) glutDestroyWindow(myWindowID[1]);
			glutDestroyWindow(myWindowID[0]);
			exit(0);
	}
}

// the mouse function is called when a mouse button is pressed down or released
void mouse(int button, int state, int x, int y, int windowIndex) {
	const bool v = false;
	if ( v ) cerr << "mouse " << button << " " << state << endl;
	shiftPressed = shiftIsPressed();
	ctrlPressed = ctrlIsPressed();
	if ( GLUT_LEFT_BUTTON == button ) {
		if ( GLUT_DOWN == state ) {
			leftMouseButtonIsDown = true;
			currXY = PixelLoc(x,y);
			startImagePt = ImageUI::asImageCoord(x,y, windowIndex);
			drawMouseInfo(x,y, shiftPressed, ctrlPressed, windowIndex);
		} else if ( GLUT_UP == state ) {
			leftMouseButtonIsDown = false;
			undrawMouseInfo();
		}
	} else if ( GLUT_MIDDLE_BUTTON == button ) {
		if ( GLUT_DOWN == state ) { // start to zoom in
			middleMouseButtonIsDown = true;
			zoomInBy(sqrt(2.0), x, y, windowIndex);
		} else {
			middleMouseButtonIsDown = false;
		}
	} else if ( GLUT_RIGHT_BUTTON == button ) {
		if ( GLUT_DOWN == state ) { // start to zoom out
			rightMouseButtonIsDown = true;
			zoomInBy(sqrt(0.5), x, y, windowIndex);
		} else {
			rightMouseButtonIsDown = false;
		}
	} else if ( 3 == button ) { // mouse wheel up (doesn't work on Mac)
		zoomInBy(sqrt(2.0), x, y, windowIndex);
	} else if ( 4 == button ) { // mouse wheel down (doesn't work on Mac)
		zoomInBy(sqrt(0.5), x, y, windowIndex);
	}
	glutPostRedisplay();
}
#endif

// simple usage message for how to use this program
void usage(const char *progname) {
	cerr << "Usage: " << progname << " tileID pair#" << endl;
	exit(-1);
}

int main(int argc, char **argv) {
	// this is all code made necessary by eriol to load the contour
	const char *progName = argv[0];
	initAllSides();	// helper structure!
	if ( argc > 1 && !strcmp(argv[1], "-no") ) {
		NO_DISPLAY = true;
		--argc; ++argv;
	}
	if ( argc < 3 ) usage(progName);
	Wright::selectedTileID = argv[1];
	--argc; ++argv;
	Wright::oneTileMode = true;
	vector<string> toAdd;
	string imgName(argv[1]);
	if ( imgName.size() == 4 && 
			 isDigit(imgName[0]) && isDigit(imgName[1]) && 
			 isDigit(imgName[2]) && isDigit(imgName[3]) ) { // assume image ID
		toAdd.push_back(imgName + 'L');
		toAdd.push_back(imgName + 'R');
	} else toAdd.push_back(imgName);
	for ( unsigned i=0,i_end=toAdd.size(); i<i_end; ++i )
		ImageUI::addImage(toAdd[i].c_str(), false);	// don't start active
	if ( !ImageUI::allImages.empty() ) {
		if ( !ImageUI::currImage()->active ) ImageUI::currImage()->makeActive();
		ImageUI::bbox.w = ImageUI::currImageWidth(false);
		ImageUI::bbox.h = ImageUI::currImageHeight(false);
		for ( unsigned i=0,i_end=ImageUI::allImages.size(); i<i_end; ++i ) {
			Image *img = ImageUI::allImages[i];
			if ( img->active )
				estimateImageScale(img->getWidth(), img->getHeight(), img->imageScale);
		}
		ImageUI::windowPixelWidth[0] = ImageUI::imageScale(false) * ImageUI::bbox.xLim();
		ImageUI::windowPixelHeight[0] = ImageUI::imageScale(false) * ImageUI::bbox.yLim();
		ImageUI::sourceColorImage = ImageUI::currImage(false);
	}
#ifdef USE_GL
	if ( ! NO_DISPLAY ) init_gl_window();	// does not return
#endif
	// examples of loading and saving information about contours
	vector<PixelLoc> in = getContour("100scottt", "1149L");
	if ( in.empty() ) cerr << " in.empty()" << endl;
	else {
		cerr << "# of Pixels " << in.size() << " First pixel: " << in[0] << " Last pixel: " << in.back() << endl;

		cv::Mat image;
		image = cv::imread(imgName, CV_LOAD_IMAGE_COLOR); 

		if(! image.data ) {
				cout <<	"Could not open or find the image" << std::endl ;
				return -1;
		}
		std::vector<PixelLoc> pixels = getContour("100scottt", "1149L");
		std::cerr<< "Number of pixels: "  << pixels.size() << endl;



		cv::Mat contour = sliceContour(pixels,image);
		cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE );
		cv::imshow( "Display window", contour );
		cv::waitKey(0);
		return 0;
	}
}
