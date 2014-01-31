#ifndef GET_MATCHES_H
#define GET_MATCHES_H
#include "surflib.h"
#include "ContourMatcherStructs.h"

//! We have the following options for SURF feature matching:
//     1. global orientation | local orientation
//     2. black outside contour | simple bounding box | ignore outside contour
//     3. adjust difference for partial features | no adjustment
//     4. only choose features inside contour | choose any features

/* mimg1 and mimg2 are slices of images within a bounding box
     around the contour; further processing optional
   globalOrientation indicates whether we will assign an orientation
     for the contour globally
  
   there will be no partial features.  */
Plane matchStrengthsSimpleBounds(bool globalOri, cv::Mat mimg1, cv::Mat mimg2){
  // Copy the images to avoid screwing up other stuff
  cv::Mat mc1 = mimg1.clone();
  cv::Mat mc2 = mimg1.clone();
  
  // Convert to IplImages for OpenSURF library use
  IplImage iimg1, iimg2;
  iimg1=mc1;
  iimg2=mc2;

  IplImage *img1, *img2;
  img1 = &iimg1;
  img2 = &iimg2;
  
  // Vectors to fill with interest points from each image
  IpVec ipts1, ipts2;
 
  surfDetDes(img1,ipts1,false,4,4,2,0.0001f,globalOri);
  surfDetDes(img2,ipts2,false,4,4,2,0.0001f,globalOri);

  MatchVec matches;
  getMatchesSymmetric(ipts1,ipts2,matches,false);

  Plane matchesVector;
  for (unsigned int i = 0; i < matches.size(); ++i){
    matchesVector.leftImage.push_back(PixelLoc(matches[i].first.first.x,matches[i].first.first.y));
    matchesVector.rightImage.push_back(PixelLoc(matches[i].first.second.x,matches[i].first.second.y));
  }
  
  if (globalOri)
    cout<<"global ";
  std::cout<<"SURF matches: " << matches.size() << std::endl;
  
  return matchesVector;
}

/* mimg1 and mimg2 are slices of images within a bounding box
     around the contour; further processing optional
   bools1 and bools2 are boolean maps of the contour's location;
     we will only choose features within the contour
   globalOrientation indicates whether we will assign an orientation
     for the contour globally
  
   there will be no partial features.  */
Plane matchStrengthsSimpleBoundsInContour(bool globalOri, cv::Mat mimg1, cv::Mat mimg2, cv::Mat bools1, cv::Mat bools2){
  // Copy the images to avoid screwing up other stuff
  cv::Mat mc1 = mimg1.clone();
  cv::Mat mc2 = mimg1.clone();
  
  // Convert to IplImages for OpenSURF library use
  IplImage iimg1, iimg2;
  iimg1=mc1;
  iimg2=mc2;

  IplImage *img1, *img2;
  img1 = &iimg1;
  img2 = &iimg2;
  
  // Vectors to fill with interest points from each image
  IpVec ipts1, ipts2;
 
  // Convert/initialize boolean contour maps
  // Ignore features outside the contour, but process the pixels 
  std::cout<<"Contour-only ";
  cv::Mat bc1 = bools1.clone();
  cv::Mat bc2 = bools2.clone();

  IplImage bi1, bi2;
  bi1 = bc1;
  bi2 = bc2;

  IplImage *b1, *b2;
  b1 = &bi1;
  b2 = &bi2;

  surfDetDes(img1,ipts1,false,4,4,2,0.0001f,globalOri,b1);
  surfDetDes(img2,ipts2,false,4,4,2,0.0001f,globalOri,b2);

  MatchVec matches;
  getMatchesSymmetric(ipts1,ipts2,matches,false);

  Plane matchesVector;
  for (unsigned int i = 0; i < matches.size(); ++i){
    matchesVector.leftImage.push_back(PixelLoc(matches[i].first.first.x,matches[i].first.first.y));
    matchesVector.rightImage.push_back(PixelLoc(matches[i].first.second.x,matches[i].first.second.y));
  }
  
  if (globalOri)
    cout<<"global ";
  std::cout<<"SURF matches: " << matches.size() << std::endl;
  
  return matchesVector;
}


/* cimg1 and cimg2 are slices of images within a bounding box
     around the contour; pixels outside the contour are (0,0,0)
   bools1 and bools2 are boolean maps of the contour's location
     if assigned, we will only choose features within the contour
   globalOrientation indicates whether we will assign an orientation
     for the contour globally
  
   there will be no partial features.  */
Plane matchStrengthsContour(bool globalOri, cv::Mat cimg1, cv::Mat cimg2, cv::Mat bools1, cv::Mat bools2, bool partialFeatureAdjust=true){
  // Copy the images to avoid screwing up other stuff
  cv::Mat mc1 = cimg1.clone();
  cv::Mat mc2 = cimg1.clone();
  
  // Convert to IplImages for OpenSURF library use
  IplImage iimg1, iimg2;
  iimg1=mc1;
  iimg2=mc2;

  IplImage *img1, *img2;
  img1 = &iimg1;
  img2 = &iimg2;
  
  // Vectors to fill with interest points from each image
  IpVec ipts1, ipts2;
 
  // Convert/initialize boolean contour maps
  std::cout<<"Contour-only ";
  cv::Mat bc1 = bools1.clone();
  cv::Mat bc2 = bools2.clone();

  IplImage bi1, bi2;
  bi1 = bc1;
  bi2 = bc2;

  IplImage *b1, *b2;
  b1 = &bi1;
  b2 = &bi2;

  surfDetDesContour(img1,ipts1,false,4,4,2,0.0001f,globalOri,b1);
  surfDetDesContour(img2,ipts2,false,4,4,2,0.0001f,globalOri,b2);

  MatchVec matches;
  getMatchesSymmetric(ipts1,ipts2,matches,partialFeatureAdjust);

  Plane matchesVector;
  for (unsigned int i = 0; i < matches.size(); ++i){
    matchesVector.leftImage.push_back(PixelLoc(matches[i].first.first.x,matches[i].first.first.y));
    matchesVector.rightImage.push_back(PixelLoc(matches[i].first.second.x,matches[i].first.second.y));
  }
  
  if (globalOri)
    cout<<"global ";
  std::cout<<"SURF matches: " << matches.size() << std::endl;
  
  return matchesVector;
}

#endif
