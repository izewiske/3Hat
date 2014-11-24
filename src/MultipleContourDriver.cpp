/*
*
*	This program will load multiple images and process them then determine how well they work
*
*/
#include "3Hat.h"
using namespace std;

int main(int argc, char** argv){
	if( argc == 1) {
		// scottt100 1149L 1149R
	 	ERR <<" Usage: imageID1 imageID2 ... imageIDn" << std::endl;
	 	return -1;
	}
	//loop through images
	for (int i = 1; i < argc; i++){
		singleImageDriver(argv[i]);
	}
	return 0;
}
