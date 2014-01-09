#include <iostream>

#ifndef __CONTOUR_MATCHER_STRUCTS__
#define __CONTOUR_MATCHER_STRUCTS__

// IO macros for easy redirecting when putting in full application
#define ERR std::cerr
#define OUT std::cout

#define THRESHOLD 3

enum MATCHER_TYPE { FLANN, BF };

struct Plane {
	std::vector<PixelLoc> leftImage;
	std::vector<PixelLoc> rightImage;
};
#endif //__CONTOUR_MATCHER_STRUCTS__