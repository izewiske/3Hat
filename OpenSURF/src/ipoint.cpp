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

#include <cv.h>
#include <vector>

#include "ipoint.h"

//! Populate IpPairVec with matched ipts 
void getMatches(IpVec &ipts1, IpVec &ipts2, IpPairVec &matches)
{
  float dist, d1, d2;
  Ipoint *match;

  matches.clear();

  for(unsigned int i = 0; i < ipts1.size(); i++) 
  {
    d1 = d2 = FLT_MAX;

    for(unsigned int j = 0; j < ipts2.size(); j++) 
    {
      dist = ipts1[i] - ipts2[j];  

      if(dist<d1) // if this feature matches better than current best
      {
        d2 = d1;
        d1 = dist;
        match = &ipts2[j];
      }
      else if(dist<d2) // this feature matches better than second best
      {
        d2 = dist;
      }
    }

    // If match has a d1:d2 ratio < MATCH_THRESHOLD ipoints are a match
    if(d1/d2 < MATCH_THRESHOLD) 
    { 
      // Store the change in position
      ipts1[i].dx = match->x - ipts1[i].x; 
      ipts1[i].dy = match->y - ipts1[i].y;
      matches.push_back(std::make_pair(ipts1[i], *match));
    }
  }
}

//-------------------------------------------------------

//! Populate IpPairVec with matched ipts AND strength of match (ratio of first over second match val)
void getMatches(IpVec &ipts1, IpVec &ipts2, MatchVec &matches){
  float dist, d1, d2;
  Ipoint *match;

  matches.clear();

  for(unsigned int i = 0; i < ipts1.size(); i++) 
  {
    d1 = d2 = FLT_MAX;

    for(unsigned int j = 0; j < ipts2.size(); j++) 
    {
      dist = ipts1[i] - ipts2[j];  

      if(dist<d1) // if this feature matches better than current best
      {
        d2 = d1;
        d1 = dist;
        match = &ipts2[j];
      }
      else if(dist<d2) // this feature matches better than second best
      {
        d2 = dist;
      }
    }

    // If match has a d1:d2 ratio < MATCH_THRESHOLD ipoints are a match
    if(d1/d2 < MATCH_THRESHOLD) 
    { 
      // Store the change in position, the match, and the match strength
      ipts1[i].dx = match->x - ipts1[i].x; 
      ipts1[i].dy = match->y - ipts1[i].y;
      matches.push_back(std::make_pair(std::make_pair(ipts1[i], *match), d1/d2));
    }
  }
}

//-------------------------------------------------------

//! Populate IpPairVec with matched ipts AND strength of match (ratio of first over second match val)
//! Invariant to ordering of &ipts1 vs &ipts2
void getMatchesSymmetric(IpVec &ipts1, IpVec &ipts2, MatchVec &matches, bool partial){
  float d1, d2;
  float** dists;
  Ipoint *match;

  dists = new float*[ipts1.size()];
  for (unsigned int i=0; i<ipts1.size(); i++)
    dists[i]= new float[ipts2.size()];

  // matches detected from each direction (comparing ipts1 to ipts2 vs. ipts2 to ipts1)
  MatchVec m1;
  MatchVec m2;

  matches.clear();

  for(unsigned int i = 0; i < ipts1.size(); i++) 
  {
    d1 = d2 = FLT_MAX;

    for(unsigned int j = 0; j < ipts2.size(); j++) 
    {
      if (partial)
        dists[i][j] = ipts1[i].partialDistance(ipts2[j]);
      else
        dists[i][j] = ipts1[i] - ipts2[j];  

      if(dists[i][j]<d1) // if this feature matches better than current best
      {
        d2 = d1;
        d1 = dists[i][j];
        match = &ipts2[j];
      }
      else if(dists[i][j]<d2) // this feature matches better than second best
      {
        d2 = dists[i][j];
      }
    }

    // If match has a d1:d2 ratio < MATCH_THRESHOLD ipoints are a match
    if(d1/d2 < MATCH_THRESHOLD) 
    { 
      // Store the change in position, the match, and the match strength
      ipts1[i].dx = match->x - ipts1[i].x; 
      ipts1[i].dy = match->y - ipts1[i].y;
      m1.push_back(std::make_pair(std::make_pair(ipts1[i], *match), d1/d2));
    }
  }

  for(unsigned int i = 0; i < ipts2.size(); i++) 
  {
    d1 = d2 = FLT_MAX;

    for(unsigned int j = 0; j < ipts1.size(); j++) 
    {
      if(dists[j][i]<d1) // if this feature matches better than current best
      {
        d2 = d1;
        d1 = dists[j][i];
        match = &ipts1[j];
      }
      else if(dists[j][i]<d2) // this feature matches better than second best
      {
        d2 = dists[j][i];
      }
    }

    // If match has a d1:d2 ratio < MATCH_THRESHOLD ipoints are a match
    if(d1/d2 < MATCH_THRESHOLD) 
    { 
      // Store the change in position, the match, and the match strength
      ipts2[i].dx = match->x - ipts2[i].x; 
      ipts2[i].dy = match->y - ipts2[i].y;
      m2.push_back(std::make_pair(std::make_pair(*match, ipts2[i]), d1/d2));
    }
  }
  
  for (unsigned int i = 0; i < m1.size(); i++)
    for (unsigned int j = 0; j< m2.size(); j++)
      if (m1[i].first==m2[j].first)
        matches.push_back(m1[i]);

}

//-------------------------------------------------------

//
// This function uses homography with CV_RANSAC (OpenCV 1.1)
// Won't compile on most linux distributions
//

//-------------------------------------------------------

//! Find homography between matched points and translate src_corners to dst_corners
int translateCorners(IpPairVec &matches, const CvPoint src_corners[4], CvPoint dst_corners[4])
{
#ifndef LINUX
  double h[9];
  CvMat _h = cvMat(3, 3, CV_64F, h);
  std::vector<CvPoint2D32f> pt1, pt2;
  CvMat _pt1, _pt2;
  
  int n = (int)matches.size();
  if( n < 4 ) return 0;

  // Set vectors to correct size
  pt1.resize(n);
  pt2.resize(n);

  // Copy Ipoints from match vector into cvPoint vectors
  for(int i = 0; i < n; i++ )
  {
    pt1[i] = cvPoint2D32f(matches[i].second.x, matches[i].second.y);
    pt2[i] = cvPoint2D32f(matches[i].first.x, matches[i].first.y);
  }
  _pt1 = cvMat(1, n, CV_32FC2, &pt1[0] );
  _pt2 = cvMat(1, n, CV_32FC2, &pt2[0] );

  // Find the homography (transformation) between the two sets of points
  if(!cvFindHomography(&_pt1, &_pt2, &_h, CV_RANSAC, 5))  // this line requires opencv 1.1
    return 0;

  // Translate src_corners to dst_corners using homography
  for(int i = 0; i < 4; i++ )
  {
    double x = src_corners[i].x, y = src_corners[i].y;
    double Z = 1./(h[6]*x + h[7]*y + h[8]);
    double X = (h[0]*x + h[1]*y + h[2])*Z;
    double Y = (h[3]*x + h[4]*y + h[5])*Z;
    dst_corners[i] = cvPoint(cvRound(X), cvRound(Y));
  }
#endif
  return 1;
}
