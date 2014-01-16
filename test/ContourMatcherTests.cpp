#include "../src/OpenCVContourMatcher.h"
#include "../src/UtilityFunctions.h"
#include <iostream>
#include <fstream>
#include <string>
#include <exception>

#ifndef __OCVM_Tests__
#define __OCVM_Tests__

static Plane expected;

void printPlane(Plane p){
	OUT << "Plane:";
	for (int i = 0; i<p.leftImage.size();i++){
		OUT << "Pixel-Locations: ( " <<p.leftImage[i] <<" , "<<p.rightImage[i] << " )\n";
	}
	return ;
}

bool planeIsAsExpected(Plane p){
	bool valid = true;
	//check that plane points coincide with what we expect
	return valid;
}


int main(){



	try {
		OUT << "Executing OpenCVContourMatcher test suite..." << "\n";
		OUT << "Opening test images..." << "\n";
		std::fstream imageFile;
		//Image img1(imageFile.open("openCVTestImage1.ppm"));
		//Image img2(imageFile.open("openCVTestImage2.ppm"));
		Image img1("openCVTestImage1.ppm");
		Image img2("openCVTestImage2.ppm");
		imageFile.close();
		OUT << "Images opened.\n";
		OUT << "Initializing matcher...\n";
		OpenCVContourMatcher matcher;
		matcher.loadImages(img1, img2);
		OUT << "Matcher initilized...\n";

		OUT << "Initializing contours...\n";
		Contour c1("L");
		Contour c2("R");
		OUT <<  "Contours initilized.\n";

		OUT << "Executing FLANN matcher test...\n";
		Plane plane = matcher.compare(c1,c2,FLANN);
		if (!planeIsAsExpected(plane)) {
			printPlane(expected);
			printPlane(plane);
			OUT << "FLANN matcher failed. \n";
			throw OCVMUnknown;
		}
		OUT << "Executing Brute force matcher test...\n";
		plane = matcher.compare(c1,c2,BF);
		if (!planeIsAsExpected(plane)) {
			printPlane(expected);
			printPlane(plane);
			OUT << "Brute force matcher failed. \n";
			throw OCVMUnknown;
		}
	} catch (std::exception& e) {
		ERR << e.what() ;
	}

	return 0;
}

#endif //__OCVM_Tests__