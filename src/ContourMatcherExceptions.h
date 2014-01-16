#include <exception>
#include <iostream>

#ifndef __CONTOUR_MATCHER_EXCEPTIONS__
#define __CONTOUR_MATCHER_EXCEPTIONS__

class OCVMMatcherNotDefined: public std::exception {
 	virtual const char* what() const throw() {
    	return "ERROR: In OpenCVContourMatcher::compare: Matcher not defined for OpenCVContourMatcher.\n";
  	}
} OCVMMatcherNotDef;

class OCVMImageSizeError: public std::exception {
 	virtual const char* what() const throw() {
    	return "ERROR: In OpenCVContourMatcher::loadImages: Cannot compare images of different sizes.\n";
  	}
} OCVMImageSize;


class OCVMImageOpenError: public std::exception {
 	virtual const char* what() const throw() {
    	return	"ERROR: Couldn't open/locate test images. \n";
  	}
} OCVMImageOpen;

class OCVMUnknownError: public std::exception {
	virtual const char* what() const throw() {
		return "ERROR: In OpenCVContourMatcher: An unknown error has occured. \n";
	}
} OCVMUnknown;


#endif //__CONTOUR_MATCHER_EXCEPTIONS__