#include "ContourMatcherStructs.h"

#ifndef __CONTOUR_MATCHER__
#define __CONTOUR_MATCHER__
class ContourMatcher{
	public:
		std::vector<int> compare(Contour contour1, Contour contour2, MATCHER_TYPE m = NULL);
};

#endif //__CONTOUR_MATCHER__