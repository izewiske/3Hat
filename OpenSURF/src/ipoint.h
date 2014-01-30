/*********************************************************** 
*  --- OpenSURF ---                                       *
*  This library is distributed under the GNU GPL. Please   *
*  use the contact form at http://www.chrisevansdev.com    *
*  for more information.                                   *
*                                                          *
*  C. Evans, Research Into Robust Visual Features,         *
*  MSc University of Bristol, 2008.                        *
*                                                          *
************************************************************/

#ifndef IPOINT_H
#define IPOINT_H

#include <vector>
#include <math.h>

#define MATCH_THRESHOLD 0.65

//-------------------------------------------------------

class Ipoint; // Pre-declaration
typedef std::vector<Ipoint> IpVec;
typedef std::vector<std::pair<Ipoint, Ipoint> > IpPairVec;
typedef std::vector<std::pair<std::pair<Ipoint, Ipoint>, float> > MatchVec;

//-------------------------------------------------------

//! Ipoint operations
void getMatches(IpVec &ipts1, IpVec &ipts2, IpPairVec &matches);
void getMatches(IpVec &ipts1, IpVec &ipts2, MatchVec &matches);
void getMatchesSymmetric(IpVec &ipts1, IpVec &ipts2, MatchVec &matches, bool partial = false);
int translateCorners(IpPairVec &matches, const CvPoint src_corners[4], CvPoint dst_corners[4]);

//-------------------------------------------------------

class Ipoint {

public:

  //! Destructor
  ~Ipoint() {};

  //! Constructor
  Ipoint() : orientation(0) {};

  //! Gets the distance in descriptor space between Ipoints
  float operator-(const Ipoint &rhs){
    float sum=0.f;
    for(int i=0; i < 64; ++i)
      sum += (this->descriptor[i] - rhs.descriptor[i])*(this->descriptor[i] - rhs.descriptor[i]);
    return sqrt(sum);
  };

  //! Gets the distance in descriptor space between partial Ipoint descriptors
  float partialDistance(const Ipoint &rhs){
    float sum=0.f;
    int count = 0;
    for(int i=0; i < 64; ++i){
      if(std::isfinite(this->descriptor[i]) && std::isfinite(rhs.descriptor[i])){
        sum += (this->descriptor[i] - rhs.descriptor[i])*(this->descriptor[i] - rhs.descriptor[i]);
        count++;
      }
    }
    
    if (count==0)
      return FLT_MAX;
    else
      return sqrt(sum/((float)count));
  };

  //! Compares two Ipoints
  bool operator==(const Ipoint &rhs) const {
    for (int i=0; i<64; i++)
      if (this->descriptor[i]!=rhs.descriptor[i])
        return false;
    return true;
  }

  //! Coordinates of the detected interest point
  float x, y;

  //! Detected scale
  float scale;

  //! Orientation measured anti-clockwise from +ve x-axis
  float orientation;

  //! Sign of laplacian for fast matching purposes
  int laplacian;

  //! Vector of descriptor components
  float descriptor[64];

  //! Placeholds for point motion (can be used for frame to frame motion analysis)
  float dx, dy;

  //! Used to store cluster index
  int clusterIndex;
};

//-------------------------------------------------------


#endif
