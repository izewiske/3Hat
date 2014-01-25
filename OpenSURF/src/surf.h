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

#ifndef SURF_H
#define SURF_H

#include <cv.h>
#include "ipoint.h"
#include "integral.h"
#include "fasthessian.h"

#include <vector>
#include <unordered_map>

class Surf {
  
  public:
    
    //! Standard Constructor (img is an integral image)
    Surf(IplImage *img, std::vector<Ipoint> &ipts);

    //! Describe all features in the supplied vector
    void getDescriptors(bool bUpright = false, IplImage* contour=NULL);
    void getDescriptorsGlobal(bool bUpright = false, const int init_sample=INIT_SAMPLE, IplImage* contour=NULL);

  private:
    
    //---------------- Private Functions -----------------//

    //! Assign the current Ipoint an orientation
    void getOrientation(IplImage* contour=NULL);
    
    //! Determine global orientation of the image at the given scale
    void getOrientationGlobal(const int init_sample=INIT_SAMPLE, IplImage* contour=NULL);

    //! Get the descriptor. See Agrawal ECCV 08
    void getDescriptor(bool bUpright = false, IplImage* contour=NULL);
    void getDescriptorGlobal(bool bUpright = false, IplImage* contour=NULL);

    //! Calculate the value of the 2d gaussian at x,y
    inline float gaussian(int x, int y, float sig);
    inline float gaussian(float x, float y, float sig);

    //! Calculate Haar wavelet responses in x and y directions
    inline float haarX(int row, int column, int size);
    inline float haarY(int row, int column, int size);
    inline float haarXContour(int row, int column, int size, IplImage* contour);
    inline float haarYContour(int row, int column, int size, IplImage* contour);

    //! Get the angle from the +ve x-axis of the vector given by [X Y]
    float getAngle(float X, float Y);

    //! Calculate the weight mask for the given offsets and orientation
    float weightMask(int x, int y, float ori);


    //---------------- Private Variables -----------------//

    //! Integral image where Ipoints have been detected
    IplImage *img;

    //! Ipoints vector
    IpVec &ipts;

    //! Index of current Ipoint in the vector
    int index;

    //! Map of scales to global orientations
    std::unordered_map<int,float> oris;
};


#endif
