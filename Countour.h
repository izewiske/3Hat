#ifndef __CONTOUR__
#define __CONTOUR__

class Contour {
	string camName;
	vector<int> superPixelID;
	public:
		vector<PixelLoc> bdy;  // cached set of locations of the contour boundary
		Color col;
		Contour(const string &cam) : camName(cam) {}
		string getCamName() const { return camName; }
		vector<int> getSuperPixels() const { return superPixelID; }
		void addSuperPixel(const string &tID, int i);
		void removeSuperPixel(const string &tID, int i);
		void drawOnOverlayImg(bool onRight);
		void computeAvgCol(bool onRight);
		void invertCol();
		vector< vector<Coord> > getSimplifiedBdy(bool onRight);
};
extern Contour noContour;
ostream &operator<<(ostream &os, const Contour &m);

#endif //__CONTOUR__