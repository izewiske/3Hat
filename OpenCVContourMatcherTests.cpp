#include "Image.h"
#include "OpenCVContourMatcher.h"
#include <iostream>

using namespace std;

int main(){
	try {
		OUT << "Executing OpenCVContourMatcher test suite..." << "\n";
		OUT << "Loading test images..." << "\n";
		OpenCVContourMatcher matcher;

	} catch (int e) {
		ERR << "One or more tests failed. Exit code " << e << "\n";
	}

}