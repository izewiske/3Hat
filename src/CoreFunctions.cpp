#include "3Hat.h"
#include "3HatDefinitions.h"

float distance(int x1,int y1, int x2, int y2){
  return sqrt(pow(x2-x1,2) + pow(y2-y1,2));
}

cv::Mat locsToBool(vector<PixelLoc> contourPixels, cv::Mat img, int pixelBuffer){
  cv::Mat boolMat(img.rows,img.cols,CV_8UC3,cvScalar(0,0,0));


  for (unsigned int i=0; i<contourPixels.size(); i++){
    for(int k=-pixelBuffer; k<pixelBuffer; k++) {
      for(int h=-pixelBuffer; h<pixelBuffer;h++){
        for (int j=0; j<3; j++){
          if (contourPixels[i].y+h >= 0 && contourPixels[i].y+h < img.rows && contourPixels[i].x+k >= 0 && contourPixels[i].x+k < img.rows) {
            boolMat.at<cv::Vec3b>(contourPixels[i].y+h,contourPixels[i].x+k)[j]=1;
          }
        }
      }
    }
  }
  /*
   *  More preparation of the boolean contour array is possible but not implemented at this time.
   */ 

  /*
  // iterate over boolMat and find pits
  for(int i=x1-pixelBuffer; i<= x2+pixelBuffer; i++){
    for(int j = y1-pixelBuffer; j <= y2+pixelBuffer; j++){
      // if it's on the edge of matrix check for invalid memory
      if (boolMat.cols > i+1 && boolMat.rows > j+1 && i >= 1 && j >= 1) {
        // check if pixel is already part of contour
        if ( boolMat.at<cv::Vec3b>(i,j)[0]!=1 && boolMat.at<cv::Vec3b>(i,j)[1]!=1 && boolMat.at<cv::Vec3b>(i,j)[2]!=1) {
          // if it isn't then check that it has 2 or more neighbors
          if (pixelNeighbors(i,j,boolMat)) {
            // if it has 2 or more neighbors then it belongs in the contour
            boolMat.at<cv::Vec3b>(contourPixels[i].y,contourPixels[i].x)[j]=1;
          }
        }
      }
    }
  }
  */

  return boolMat;
}


/* mimg1 and mimg2 are slices of images within a bounding box
     around the contour; further processing optional
   globalOrientation indicates whether we will assign an orientation
     for the contour globally
  
   this will run basic OpenSURF on the two input images;
   orientation can be toggled from local to global;
   descriptors are not limited to the contours;
   there will be no partial features */
Plane matchStrengthsSimpleBounds(bool globalOri, cv::Mat mimg1, cv::Mat mimg2){
  // Copy the images to avoid screwing up other stuff
  cv::Mat mc1 = mimg1.clone();
  cv::Mat mc2 = mimg2.clone();
  
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
    OUT<<"global ";
  OUT<<"SURF matches: " << matches.size() << std::endl;
  
  return matchesVector;
}


/* mimg1 and mimg2 are slices of images within a bounding box
     around the contour; further processing optional
   bools1 and bools2 are boolean maps of the contour's location;
     we will only choose features within the contour
   globalOrientation indicates whether we will assign an orientation
     for the contour globally
  
   this will run basic SURF on the two input images, but
   descriptors will be confined to the contours;
   orientation can be toggled from local to global;
   there will be no partial features  */
Plane matchStrengthsSimpleBoundsInContour(bool globalOri, cv::Mat mimg1, cv::Mat mimg2, cv::Mat bools1, cv::Mat bools2){
  // Copy the images to avoid screwing up other stuff
  cv::Mat mc1 = mimg1.clone();
  cv::Mat mc2 = mimg2.clone();
  
  // Convert to IplImages for OpenSURF library use
  IplImage iimg1, iimg2;
  iimg1=mc1;
  iimg2=mc2;

  IplImage *img1, *img2;
  img1 = &iimg1;
  img2 = &iimg2;
  
  // Convert/initialize boolean contour maps
  // Ignore features outside the contour, but process the pixels 
  OUT<<"Contour-only ";
  cv::Mat bc1 = bools1.clone();
  cv::Mat bc2 = bools2.clone();

  IplImage bi1, bi2;
  bi1 = bc1;
  bi2 = bc2;

  IplImage *b1, *b2;
  b1 = &bi1;
  b2 = &bi2;

  // Vectors to fill with interest points from each image
  IpVec ipts1, ipts2;
 
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
    OUT<<"global ";
  OUT<<"SURF matches: " << matches.size() << std::endl;
  
  return matchesVector;
}


/* cimg1 and cimg2 are slices of images within a bounding box
     around the contour; pixels outside the contour are (0,0,0)
   bools1 and bools2 are boolean maps of the contour's location
     if assigned, we will only choose features within the contour
   globalOrientation indicates whether we will assign an orientation
     for the contour globally
  
   this will run contour-only SURF;
   pixels not in the contour will be ignored;
   orientation can be toggled from local to global;
   descriptors will be confined to the contours;
   there will be partial features.  */
Plane matchStrengthsContour(bool globalOri, cv::Mat cimg1, cv::Mat cimg2, cv::Mat bools1, cv::Mat bools2, bool partialFeatureAdjust=true){
  // Copy the images to avoid screwing up other stuff
  cv::Mat mc1 = cimg1.clone();
  cv::Mat mc2 = cimg2.clone();
  cv::Mat bc1 = bools1.clone();
  cv::Mat bc2 = bools2.clone();
  
  // Convert to IplImages for OpenSURF library use
  IplImage iimg1, iimg2, bi1, bi2;
  iimg1=mc1;
  iimg2=mc2;
  bi1 = bc1;
  bi2 = bc2;

  IplImage *img1, *img2, *b1, *b2;
  img1 = &iimg1;
  img2 = &iimg2;
  b1 = &bi1;
  b2 = &bi2;
  
  // Vectors to fill with interest points from each image
  IpVec ipts1, ipts2;
 
  surfDetDesContour(img1,ipts1,false,4,4,2,0.0001f,globalOri,b1);
  surfDetDesContour(img2,ipts2,false,4,4,2,0.0001f,globalOri,b2);

  MatchVec matches;
  getMatchesSymmetric(ipts1,ipts2,matches,partialFeatureAdjust);

  Plane matchesVector;
  for (unsigned int i = 0; i < matches.size(); ++i){
    matchesVector.leftImage.push_back(PixelLoc(matches[i].first.first.x,matches[i].first.first.y));
    matchesVector.rightImage.push_back(PixelLoc(matches[i].first.second.x,matches[i].first.second.y));
  }
  
  OUT<<"Contour-only ";

  if (globalOri)
    OUT<<"global ";
  OUT<<"SURF matches: " << matches.size() << std::endl;
  
  return matchesVector;
}


// function finds the best possible plane we could have computed with a series of points 
Plane bestPossibleComputedPlane(Plane computedPoints, Plane actualPoints){ 
  // determines the best fitting plane computed - best possible scenario
  Plane bestComputedPlane;
  for (int j=0; j< actualPoints.leftImage.size(); j++){
    // for each computed point search for a closest match to actuality
    int closestMatchIndex = 0;
    float closestDist = 10000.0;
    for(int i=0;i< computedPoints.leftImage.size(); i++){
      float distL = distance(computedPoints.leftImage[i].x,computedPoints.leftImage[i].y,actualPoints.leftImage[j].x,actualPoints.leftImage[j].y);
      float distR = distance(computedPoints.rightImage[i].x,computedPoints.rightImage[i].y,actualPoints.rightImage[j].x,actualPoints.rightImage[j].y);
      // if the left image point is darn close and the right image point is darn close
      //if ( distL <= closestDistL + APPROX_THRESHOLD && distR <= closestDistR + APPROX_THRESHOLD ) {
      // replaced threshold with a notion of average distance from desired
      if ( ((distL + distR) / 2) < closestDist){
        closestDist = ((distL + distR) / 2);
        closestMatchIndex = i;
      }
    }
    bestComputedPlane.leftImage.push_back(PixelLoc(computedPoints.leftImage[closestMatchIndex].x,computedPoints.leftImage[closestMatchIndex].y));
    bestComputedPlane.rightImage.push_back(PixelLoc(computedPoints.rightImage[closestMatchIndex].x,computedPoints.rightImage[closestMatchIndex].y));
  }
  return bestComputedPlane;
}


// function returns a summative overall quality of all matches 
double compareFeaturePoints(Plane computedPoints, Plane actualPoints){
  // This determines how well the features we have computed match those chosen by humans
  float sumOfDistances = 0;
  // TODO: consider using RANSAC here
  for(int i=0;i< computedPoints.leftImage.size(); i++){
    // for each computed point search for a closest match to actuality
    int closestMatchIndex = 0;
    float closestDist = 10000.0;
    for(int j=0; j< actualPoints.leftImage.size(); j++){
      float distL = distance(computedPoints.leftImage[i].x,computedPoints.leftImage[i].y,actualPoints.leftImage[j].x,actualPoints.leftImage[j].y);
      float distR = distance(computedPoints.rightImage[i].x,computedPoints.rightImage[i].y,actualPoints.rightImage[j].x,actualPoints.rightImage[j].y);
      // if the left image point is darn close and the right image point is darn close
      //if ( distL <= closestDistL + APPROX_THRESHOLD && distR <= closestDistR + APPROX_THRESHOLD ) {
      // replaced threshold with a notion of average distance from desired
      if ( ((distL + distR) / 2) < closestDist){
        closestDist = ((distL + distR) / 2);
        closestMatchIndex = i;
      }
    }
    sumOfDistances = sumOfDistances + closestDist;
  }

  // for now quality is determined by the sum of differences.
  // TODO: consider counting number of matches with distance of under threshold
  double quality = sumOfDistances;
  return quality;
}

Plane getUserDefinedPlane(std::string tileID,std::string imageID){
  // TODO: I don't know if these are ordered in a fashion such that coresponding points have same index
  std::string imageIDL = imageID + "L";
  std::string imageIDR = imageID + "R";
  std::vector<Coord> leftStandard = getFeaturePoints(tileID,imageIDL);
  std::vector<Coord> rightStandard = getFeaturePoints(tileID,imageIDR);
  if (leftStandard.size() != rightStandard.size()){
    ERR << "That's odd, the number of points in the right and left images of the gold standard differs. You broke science.\n";
    exit;
  }
  Plane actual;
  for (int i = 0; i < leftStandard.size() ; ++i) {
    actual.leftImage.push_back(PixelLoc(leftStandard[i].x,leftStandard[i].y));
    actual.rightImage.push_back(PixelLoc(rightStandard[i].x,rightStandard[i].y));
  }
  return actual;
}
