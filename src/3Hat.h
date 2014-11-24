#ifndef __3HAT__
#define __3HAT__

#include <iostream>
#include <vector>
#include <string>
#include "eriolHeader.h"
#include <stdlib.h>


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
#include "3HatDefinitions.h"

// turns contour pixel locations into a boolean matrix
cv::Mat locsToBool(vector<PixelLoc> contourPixels, cv::Mat img, int pixelBuffer=3);

// function finds the best possible plane we could have computed with a series of points 
Plane bestPossibleComputedPlane(Plane computedPoints, Plane actualPoints);

// function returns a summative overall quality of all matches 
double compareFeaturePoints(Plane computedPoints, Plane actualPoints);

// Interfaces with Eriol to get the Gold standard
Plane getUserDefinedPlane(std::string tileID,std::string imageID);

Plane matchStrengthsSimpleBounds(bool globalOri, cv::Mat mimg1, cv::Mat mimg2);

Plane matchStrengthsSimpleBoundsInContour(bool globalOri, cv::Mat mimg1, cv::Mat mimg2, cv::Mat bools1, cv::Mat bools2);

Plane matchStrengthsContour(bool globalOri, cv::Mat cimg1, cv::Mat cimg2, cv::Mat bools1, cv::Mat bools2, bool partialFeatureAdjust);

Plane singleContourDriver(std::string tileID, std::string imageIDL, std::string imageIDR);

int singleImageDriver(std::string imageID);

#endif //__3HAT__