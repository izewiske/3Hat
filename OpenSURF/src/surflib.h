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

#ifndef SURFLIB_H
#define SURFLIB_H

#include <cv.h>
#include <highgui.h>

#include "integral.h"
#include "fasthessian.h"
#include "surf.h"
#include "ipoint.h"
#include "utils.h"


//! Library function builds vector of described interest points
inline void surfDetDesContour(IplImage *con_img,  /* image to find Ipoints in */
                       std::vector<Ipoint> &ipts, /* reference to vector of Ipoints */
                       bool upright = false, /* run in rotation invariant mode? */
                       int octaves = OCTAVES, /* number of octaves to calculate */
                       int intervals = INTERVALS, /* number of intervals per octave */
                       int init_sample = INIT_SAMPLE, /* initial sampling step */
                       float thres = THRES, /* blob response threshold */
		       bool global = false, /* run in global orientation mode? */
		       IplImage* contour = NULL)
{
  // Create integral-image representation of the image
  IplImage *int_img = Integral(con_img);
  
  // Creat integral-image representation of the contour map
  IplImage *int_con;
  if (contour!=NULL)
    int_con = Integral(contour);
  else
    int_con = NULL;
  
  // Create Fast Hessian Object
  FastHessian fh(int_img, ipts, octaves, intervals, init_sample, thres, contour, int_con);
 
  // Extract interest points and store in vector ipts
  fh.getIpoints();
  
  // Create Surf Descriptor Object
  Surf des(int_img, ipts);

  // Extract the descriptors for the ipts
  if (global)
    des.getDescriptorsGlobal(upright, int_con);
  else
    des.getDescriptors(upright, int_con);

  // Deallocate the integral images
  cvReleaseImage(&int_img);
  cvReleaseImage(&int_con);
}

//! Library function builds vector of described interest points
//  Only use contour for throwing out interest points (in fasthessian)
inline void surfDetDes(IplImage *img,  /* image to find Ipoints in */
                       std::vector<Ipoint> &ipts, /* reference to vector of Ipoints */
                       bool upright = false, /* run in rotation invariant mode? */
                       int octaves = OCTAVES, /* number of octaves to calculate */
                       int intervals = INTERVALS, /* number of intervals per octave */
                       int init_sample = INIT_SAMPLE, /* initial sampling step */
                       float thres = THRES, /* blob response threshold */
		       bool global = false, /* run in global orientation mode? */
		       IplImage* contour = NULL)
{
  // Create integral-image representation of the image
  IplImage *int_img = Integral(img);
  
  // Creat integral-image representation of the contour map
  IplImage *int_con;
  if (contour!=NULL)
    int_con = Integral(contour);
  else
    int_con = NULL;
  
  // Create Fast Hessian Object
  FastHessian fh(int_img, ipts, octaves, intervals, init_sample, thres, contour, int_con);
 
  // Extract interest points and store in vector ipts
  fh.getIpoints();
  
  // Create Surf Descriptor Object
  Surf des(int_img, ipts);

  // Extract the descriptors for the ipts
  if (global)
    des.getDescriptorsGlobal(upright);
  else
    des.getDescriptors(upright);

  // Deallocate the integral images
  cvReleaseImage(&int_img);
  cvReleaseImage(&int_con);
}

//! Library function builds vector of interest points
inline void surfDet(IplImage *img,  /* image to find Ipoints in */
                    std::vector<Ipoint> &ipts, /* reference to vector of Ipoints */
                    int octaves = OCTAVES, /* number of octaves to calculate */
                    int intervals = INTERVALS, /* number of intervals per octave */
                    int init_sample = INIT_SAMPLE, /* initial sampling step */
                    float thres = THRES, /* blob response threshold */
		    IplImage* contour = NULL)
{
  // Create integral image representation of the image
  IplImage *int_img = Integral(img);

 // Creat integral-image representation of the contour map
  IplImage *int_con;
  if (contour!=NULL)
    int_con = Integral(contour);
  else
    int_con = NULL;

  // Create Fast Hessian Object
  FastHessian fh(int_img, ipts, octaves, intervals, init_sample, thres, contour, int_con);

  // Extract interest points and store in vector ipts
  fh.getIpoints();

  // Deallocate the integral images
  cvReleaseImage(&int_img);
  cvReleaseImage(&int_con);
}




//! Library function describes interest points in vector
inline void surfDes(IplImage *img,  /* image to find Ipoints in */
                    std::vector<Ipoint> &ipts, /* reference to vector of Ipoints */
                    bool upright = false, /* run in rotation invariant mode? */
		    bool global = false, /* use global orientations? */
		    IplImage* contour = NULL)
{ 
  // Create integral image representation of the image
  IplImage *int_img = Integral(img);

  // Creat integral-image representation of the contour map
  IplImage *int_con;
  if (contour!=NULL)
    int_con = Integral(contour);
  else
    int_con = NULL;

  // Create Surf Descriptor Object
  Surf des(int_img, ipts);

  // Extract the descriptors for the ipts
  if (global)
    des.getDescriptorsGlobal(upright, int_con);
  else
    des.getDescriptors(upright, int_con);
  
  // Deallocate the integral images
  cvReleaseImage(&int_img);
  cvReleaseImage(&int_con);
}


#endif
