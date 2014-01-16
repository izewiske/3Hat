// Contour.h by Olaf Hall-Holt, 2007-2011
#include "Color.h"
#include "Coord.h"
#include <string>
#include <vector>
#include <iostream>

#ifndef __CONTOUR__
#define __CONTOUR__

class Contour {
	std::string camName;
	std::vector<int> superPixelID;
	public:
		std::vector<PixelLoc> bdy;  // cached set of locations of the contour boundary
		Color col;
		Contour(const std::string &cam) : camName(cam) {}
		std::string getCamName() const { return camName; }
		std::vector<int> getSuperPixels() const { return superPixelID; }
		void addSuperPixel(const std::string &tID, int i);
		void removeSuperPixel(const std::string &tID, int i);
		void drawOnOverlayImg(bool onRight);
		void computeAvgCol(bool onRight);
		void invertCol();
		std::vector< std::vector<Coord> > getSimplifiedBdy(bool onRight);
		std::vector<PixelLoc> getContourBoundary(){ return bdy;}
};
extern Contour noContour;
std::ostream &operator<<(std::ostream &os, const Contour &m);

#endif //__CONTOUR__