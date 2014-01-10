#include "Image.h"
#include "PixelLoc.h"
#include "Contour.h"
#include "ContourMatcherStructs.h"
#include <iostream>
#include <vector>

#ifndef __CONTOUR_MATCHER__
#define __CONTOUR_MATCHER__
class ContourMatcher{
	public:
		Plane compare(Contour contour1, Contour contour2, MATCHER_TYPE m);
		void loadImages(unsigned char* i1, unsigned int height1, unsigned int width1, unsigned char* i2, unsigned int height2, unsigned int width2);
		void loadImages(Image i1, Image i2);
};

#endif //__CONTOUR_MATCHER__