#include "PixelLoc.h"
#include "Contour.h"
#include "ContourMatcherStructs.h"


#ifndef __CONTOUR_MATCHER__
#define __CONTOUR_MATCHER__
class ContourMatcher{
	public:
		std::vector<PixelLoc> compare(Contour contour1, Contour contour2, MATCHER_TYPE m);
		void loadImages(unsigned char* image, unsigned char* image);
};

#endif //__CONTOUR_MATCHER__