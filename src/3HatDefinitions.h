#ifndef __3HAT_DEFINITIONS__
#define __3HAT_DEFINITIONS__

#include <vector>
#include "eriolHeader.h"

/*
 * Approx threshold is the threshold for "close-enough" that is used to allow a perfect left image match and a slightly-off
 * right image match to still supercede the "best-fit" match approximation.
 */
#define APPROX_THRESHOLD 10
#define _HESSIAN_THRESH 15

// IO macros for easy redirecting when putting in full application
#define ERR std::cerr
#define OUT std::cout

#define THRESHOLD 3


// Structs
 struct Plane {
	std::vector<PixelLoc> leftImage;
	std::vector<PixelLoc> rightImage;
};


struct Feature {
	int x,y;
	float weight;
	int octave;
};

#endif __3HAT_DEFINITIONS__