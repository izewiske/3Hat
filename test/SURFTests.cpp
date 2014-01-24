#include "surflib.h"
#include "kmeans.h"
#include <ctime>
#include <iostream>

//-------------------------------------------------------
// Define PROCEDURE as:
//  - 0 and supply image path to run on static image
//  - 1 to show matches between static images
//  - 2 to show matches between static images and print match strengths
#define PROCEDURE 2 

int mainImage(void)
{
  // Declare Ipoints and other stuff
  IpVec ipts;
  // Make image as a Mat; convert to IplImage for OpenSURF library actions
  cv::Mat mimg=cv::imread("OpenSURF/imgs/sf.jpg", CV_LOAD_IMAGE_COLOR);
  IplImage iimg=mimg;
  IplImage* img=&iimg;

  // Detect and describe interest points in the image
  clock_t start = clock();
  surfDetDes(img, ipts, false, 5, 4, 2, 0.0004f); 
  clock_t end = clock();

  std::cout<< "OpenSURF found: " << ipts.size() << " interest points" << std::endl;
  std::cout<< "OpenSURF took: " << float(end - start) / CLOCKS_PER_SEC  << " seconds" << std::endl;

  // Draw the detected points
  drawIpoints(img, ipts);
  
  // Display the result
  showImage(img);

  return 0;
}

//-------------------------------------------------------

int mainStaticMatch()
{
  // Make images as Mats; convert to IplImage for OpenSURF library actions
  cv::Mat mimg1, mimg2;
  mimg1=cv::imread("OpenSURF/imgs/img1.jpg", CV_LOAD_IMAGE_COLOR);
  mimg2=cv::imread("OpenSURF/imgs/img2.jpg", CV_LOAD_IMAGE_COLOR);

  IplImage iimg1, iimg2;
  iimg1=mimg1;
  iimg2=mimg2;

  IplImage *img1, *img2;
  img1 = &iimg1;
  img2 = &iimg2;

  IpVec ipts1, ipts2;
  surfDetDes(img1,ipts1,false,4,4,2,0.0001f);
  surfDetDes(img2,ipts2,false,4,4,2,0.0001f);

  IpPairVec matches;
  getMatches(ipts1,ipts2,matches);

  for (unsigned int i = 0; i < matches.size(); ++i)
  {
    drawPoint(img1,matches[i].first);
    drawPoint(img2,matches[i].second);
  
    const int & w = img1->width;
    cvLine(img1,cvPoint(matches[i].first.x,matches[i].first.y),cvPoint(matches[i].second.x+w,matches[i].second.y), cvScalar(255,255,255),1);
    cvLine(img2,cvPoint(matches[i].first.x-w,matches[i].first.y),cvPoint(matches[i].second.x,matches[i].second.y), cvScalar(255,255,255),1);
  }

  std::cout<< "Matches: " << matches.size();

  cvNamedWindow("1", CV_WINDOW_AUTOSIZE );
  cvNamedWindow("2", CV_WINDOW_AUTOSIZE );
  cvShowImage("1", img1);
  cvShowImage("2",img2);
  cvWaitKey(0);

  return 0;
}

//-------------------------------------------------------

int mainStaticMatchStrengths()
{
  bool matchGlobalOrientations = true;

  // Make images as Mats; convert to IplImage for OpenSURF library actions
  cv::Mat mimg1, mimg2;
  mimg1=cv::imread("OpenSURF/imgs/img1.jpg", CV_LOAD_IMAGE_COLOR);
  mimg2=cv::imread("OpenSURF/imgs/img2.jpg", CV_LOAD_IMAGE_COLOR);

  IplImage iimg1, iimg2;
  iimg1=mimg1;
  iimg2=mimg2;

  IplImage *img1, *img2;
  img1 = &iimg1;
  img2 = &iimg2;

  IpVec ipts1, ipts2;
  surfDetDes(img1,ipts1,false,4,4,2,0.0001f,matchGlobalOrientations);
  surfDetDes(img2,ipts2,false,4,4,2,0.0001f,matchGlobalOrientations);

  MatchVec matches;
  getMatchesSymmetric(ipts1,ipts2,matches);

  IpVec mpts1, mpts2;

  const int & w = img1->width;

  for (unsigned int i = 0; i < matches.size(); ++i)
  {
    float strengthOverThreshold = 1 - matches[i].second; // /MATCH_THRESHOLD;
    strengthOverThreshold*=255;
    CvScalar clr = cvScalar(strengthOverThreshold,strengthOverThreshold,strengthOverThreshold);
    clr = cvScalar(255,255,255);
    
    //drawPoint(img1,matches[i].first.first,clr);
    //drawPoint(img2,matches[i].first.second,clr),
    mpts1.push_back(matches[i].first.first);
    mpts2.push_back(matches[i].first.second);
  
    cvLine(img1,cvPoint(matches[i].first.first.x,matches[i].first.first.y),cvPoint(matches[i].first.second.x+w,matches[i].first.second.y), clr,1);
    cvLine(img2,cvPoint(matches[i].first.first.x-w,matches[i].first.first.y),cvPoint(matches[i].first.second.x,matches[i].first.second.y), clr,1);
  }

  drawIpoints(img1,mpts1);
  drawIpoints(img2,mpts2);

  std::cout<< "Matches: " << matches.size() << std::endl;

  cvNamedWindow("1", CV_WINDOW_AUTOSIZE );
  cvNamedWindow("2", CV_WINDOW_AUTOSIZE );
  cvShowImage("1", img1);
  cvShowImage("2",img2);
  cvWaitKey(0);

  // NOW DO IT AGAIN!
  cv::Mat mimg3, mimg4;
  mimg3=cv::imread("OpenSURF/imgs/img1.jpg", CV_LOAD_IMAGE_COLOR);
  mimg4=cv::imread("OpenSURF/imgs/img2.jpg", CV_LOAD_IMAGE_COLOR);

  IplImage iimg3, iimg4;
  iimg3=mimg3;
  iimg4=mimg4;

  IplImage *img3, *img4;
  img3 = &iimg3;
  img4 = &iimg4;

  IpVec ipts3, ipts4;
  surfDetDes(img3,ipts3,false,4,4,2,0.0001f,!matchGlobalOrientations);
  surfDetDes(img4,ipts4,false,4,4,2,0.0001f,!matchGlobalOrientations);

  matches.clear();
  getMatchesSymmetric(ipts3,ipts4,matches);

  IpVec mpts3, mpts4;

  for (unsigned int i = 0; i < matches.size(); ++i)
  {
    float strengthOverThreshold = 1 - matches[i].second; // /MATCH_THRESHOLD;
    strengthOverThreshold*=255;
    CvScalar clr = cvScalar(strengthOverThreshold,strengthOverThreshold,strengthOverThreshold);
    clr = cvScalar(255,255,255);
    
    //drawPoint(img1,matches[i].first.first,clr);
    //drawPoint(img2,matches[i].first.second,clr),
    mpts3.push_back(matches[i].first.first);
    mpts4.push_back(matches[i].first.second);
  
    cvLine(img3,cvPoint(matches[i].first.first.x,matches[i].first.first.y),cvPoint(matches[i].first.second.x+w,matches[i].first.second.y), clr,1);
    cvLine(img4,cvPoint(matches[i].first.first.x-w,matches[i].first.first.y),cvPoint(matches[i].first.second.x,matches[i].first.second.y), clr,1);
  }

  drawIpoints(img3,mpts3);
  drawIpoints(img4,mpts4);

  std::cout<< "Matches: " << matches.size() << std::endl;

  cvNamedWindow("3", CV_WINDOW_AUTOSIZE );
  cvNamedWindow("4", CV_WINDOW_AUTOSIZE );
  cvShowImage("3", img3);
  cvShowImage("4",img4);
  cvWaitKey(0);


  return 0;
}

//-------------------------------------------------------


int main(void) 
{
  if (PROCEDURE == 0) return mainImage();
  if (PROCEDURE == 1) return mainStaticMatch();
  if (PROCEDURE == 2) return mainStaticMatchStrengths();
}
