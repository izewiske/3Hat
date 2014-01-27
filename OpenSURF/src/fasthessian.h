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

#ifndef FASTHESSIAN_H
#define FASTHESSIAN_H

#include <cv.h>
#include "ipoint.h"
#include "ContourMat.h"

#include <vector>

class ResponseLayer;
static const int OCTAVES = 10;
static const int INTERVALS = 4;
static const float THRES = 0.0004f;
static const int INIT_SAMPLE = 2;


class FastHessian {
  
  public:
   
    //! Constructor without image
    FastHessian(std::vector<Ipoint> &ipts, 
                const int octaves = OCTAVES, 
                const int intervals = INTERVALS, 
                const int init_sample = INIT_SAMPLE, 
                const float thres = THRES,
		IplImage* contour = NULL,
		IplImage* int_con = NULL);

    //! Constructor with image
    FastHessian(IplImage *img, 
                std::vector<Ipoint> &ipts, 
                const int octaves = OCTAVES, 
                const int intervals = INTERVALS, 
                const int init_sample = INIT_SAMPLE, 
                const float thres = THRES,
		IplImage* contour = NULL,
		IplImage* int_con = NULL);

    //! Destructor
    ~FastHessian();

    //! Save the parameters
    void saveParameters(const int octaves, 
                        const int intervals,
                        const int init_sample, 
                        const float thres);

    //! Set or re-set the integral image source
    void setIntImage(IplImage *img);

    //! Set or re-set the contour map source
    void setConMap(IplImage *contour);

    //! Set or re-set the contour integral image source
    void setConImage(IplImage *int_con);

    //! Find the image features and write into vector of features
    void getIpoints();

  private:

    //---------------- Private Functions -----------------//

    //! Build map of DoH responses
    void buildResponseMap();

    //! Calculate DoH responses for supplied layer
    void buildResponseLayer(ResponseLayer *r);

    //! 3x3x3 Extrema test
    int isExtremum(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b);    
    int isExtremumInContour(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b, ContourMat* con); 
    
    //! Interpolation functions - adapted from Lowe's SIFT implementation
    void interpolateExtremum(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b);
    void interpolateExtremumInContour(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b, ContourMat* con);
    void interpolateStep(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b,
                          double* xi, double* xr, double* xc );
    void interpolateStepInContour(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b,
                          double* xi, double* xr, double* xc, ContourMat* con );
    CvMat* deriv3D(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b);
    CvMat* deriv3DInContour(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b, ContourMat* con);
    CvMat* hessian3D(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b);
    CvMat* hessian3DInContour(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b, ContourMat* con);

    //---------------- Private Variables -----------------//

    //! Pointer to the integral Image, and its attributes 
    IplImage *img, *contour, *int_con;
    int i_width, i_height;

    //! Reference to vector of features passed from outside 
    std::vector<Ipoint> &ipts;

    //! Response stack of determinant of hessian values
    std::vector<ResponseLayer *> responseMap;

    //! Number of Octaves
    int octaves;

    //! Number of Intervals per octave
    int intervals;

    //! Initial sampling step for Ipoint detection
    int init_sample;

    //! Threshold value for blob resonses
    float thresh;
};


#endif
