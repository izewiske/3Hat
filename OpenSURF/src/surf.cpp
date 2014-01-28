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

#include "utils.h"
#include "surf.h"

#include <iostream>
#include <cmath>

//-------------------------------------------------------
//! SURF priors (these need not be done at runtime)
const float pi = 3.14159f;

//! lookup table for 2d gaussian (sigma = 2.5) where (0,0) is top left and (6,6) is bottom right
const double gauss25 [7][7] = {
  0.02546481,	0.02350698,	0.01849125,	0.01239505,	0.00708017,	0.00344629,	0.00142946,
  0.02350698,	0.02169968,	0.01706957,	0.01144208,	0.00653582,	0.00318132,	0.00131956,
  0.01849125,	0.01706957,	0.01342740,	0.00900066,	0.00514126,	0.00250252,	0.00103800,
  0.01239505,	0.01144208,	0.00900066,	0.00603332,	0.00344629,	0.00167749,	0.00069579,
  0.00708017,	0.00653582,	0.00514126,	0.00344629,	0.00196855,	0.00095820,	0.00039744,
  0.00344629,	0.00318132,	0.00250252,	0.00167749,	0.00095820,	0.00046640,	0.00019346,
  0.00142946,	0.00131956,	0.00103800,	0.00069579,	0.00039744,	0.00019346,	0.00008024
};

//-------------------------------------------------------

//! Constructor
Surf::Surf(IplImage *img, IpVec &ipts)
: ipts(ipts)
{
  this->img = img;
}

//-------------------------------------------------------

//! Describe all features in the supplied vector
void Surf::getDescriptors(bool upright, IplImage* int_con)
{
  // Check there are Ipoints to be described
  if (!ipts.size()) return;

  // Get the size of the vector for fixed loop bounds
  int ipts_size = (int)ipts.size();

  if (upright)
  {
    // U-SURF loop just gets descriptors
    for (int i = 0; i < ipts_size; ++i)
    {
      // Set the Ipoint to be described
      index = i;

      // Extract upright (i.e. not rotation invariant) descriptors
      getDescriptor(true, int_con);
    }
  }
  else
  {
    // Main SURF-64 loop assigns orientations and gets descriptors
    for (int i = 0; i < ipts_size; ++i)
    {
      // Set the Ipoint to be described
      index = i;

      // Assign Orientations and extract rotation invariant descriptors
      getOrientation(int_con);
      getDescriptor(false, int_con);
    }
  }
}

//-------------------------------------------------------

//! Assign the supplied Ipoint an orientation
void Surf::getOrientation(IplImage* int_con)
{
  // the interest point to which we'd like to assign an orientation
  Ipoint *ipt = &ipts[index];
  float gauss = 0.f, scale = ipt->scale;
  const int s = fRound(scale), r = fRound(ipt->y), c = fRound(ipt->x);
  std::vector<float> resX(109), resY(109), Ang(109);

  // speed optimization over using absolute value evaluation
  const int id[] = {6,5,4,3,2,1,0,1,2,3,4,5,6};

  int idx = 0;

  // choose whether to use haar or haarContour
  if (int_con==NULL){
    // calculate haar responses for points within radius of 6*scale
    for(int i = -6; i <= 6; ++i) 
    {
      for(int j = -6; j <= 6; ++j) 
      {
        if(i*i + j*j < 36) // up to 6 away in any direction
        {
          // look up the value of the 2d gaussian at the given offset from center
          gauss = static_cast<float>(gauss25[id[i+6]][id[j+6]]);  // could use abs() but id lookup is faster
          // multiply distance away in i & j by scale for calculating wavelet responses
          resX[idx] = gauss * haarX(r+j*s, c+i*s, 4*s);
          resY[idx] = gauss * haarY(r+j*s, c+i*s, 4*s);
	  // angle of the gradient (up from x-axis)
          Ang[idx] = getAngle(resX[idx], resY[idx]);
          ++idx;
        }
      }
    }
  }
  else{
    // calculate haar responses for points within radius of 6*scale
    for(int i = -6; i <= 6; ++i) 
    {
      for(int j = -6; j <= 6; ++j) 
      {
        if(i*i + j*j < 36) // up to 6 away in any direction
        {
          // look up the value of the 2d gaussian at the given offset from center
          gauss = static_cast<float>(gauss25[id[i+6]][id[j+6]]);  // could use abs() but id lookup is faster
          // multiply distance away in i & j by scale for calculating wavelet responses
          resX[idx] = gauss * haarXContour(r+j*s, c+i*s, 4*s, int_con);
          resY[idx] = gauss * haarYContour(r+j*s, c+i*s, 4*s, int_con);
	  
	  // angle of the gradient (up from x-axis)
	  if (std::isfinite(resX[idx]) && std::isfinite(resY[idx]))
            Ang[idx] = getAngle(resX[idx], resY[idx]);
	  else{
	    Ang[idx] = FLT_MAX;
          }
          ++idx;
        }
      }
    }
  }

  // calculate the dominant direction 
  float sumX=0.f, sumY=0.f;
  float max=0.f, orientation = 0.f;
  float ang1=0.f, ang2=0.f;

  // loop slides pi/3 window from ang1 to ang2 around feature point
  // increments of .15 radians
  for(ang1 = 0; ang1 < 2*pi;  ang1+=0.15f) {
    ang2 = ( ang1+pi/3.0f > 2*pi ? ang1-5.0f*pi/3.0f : ang1+pi/3.0f);

    sumX = sumY = 0.f; 
    for(unsigned int k = 0; k < Ang.size(); ++k) 
    {
      // get angle from the x-axis of the sample point
      const float & ang = Ang[k];

      // determine whether the point's gradient is within the window
      // if so, add its x and y componenets to our sums 
      if (ang1 < ang2 && ang1 < ang && ang < ang2) 
      {
        sumX+=resX[k];  
        sumY+=resY[k];
      } 
      else if (ang2 < ang1 && 
        ((ang > 0 && ang < ang2) || (ang > ang1 && ang < 2*pi) )) 
      {
        sumX+=resX[k];  
        sumY+=resY[k];
      }
    }

    // if the vector produced from this window is longer than all 
    // previous vectors then this forms the new dominant direction
    if (sumX*sumX + sumY*sumY > max) 
    {
      // store largest orientation
      max = sumX*sumX + sumY*sumY;
      orientation = getAngle(sumX, sumY);
    }
  }

  // assign orientation of the dominant response vector
  ipt->orientation = orientation;
}

//-------------------------------------------------------

//! Get the modified descriptor. See Agrawal ECCV 08
//! Modified descriptor contributed by Pablo Fernandez
void Surf::getDescriptor(bool bUpright, IplImage* int_con)
{
  int y, x, sample_x, sample_y, count=0;
  int i = 0, ix = 0, j = 0, jx = 0, xs = 0, ys = 0;
  float scale, *desc, dx, dy, mdx, mdy, co, si;
  float gauss_s1 = 0.f, gauss_s2 = 0.f;
  float rx = 0.f, ry = 0.f, rrx = 0.f, rry = 0.f, len = 0.f;
  float cx = -0.5f, cy = 0.f; //Subregion centers for the 4x4 gaussian weighting

  Ipoint *ipt = &ipts[index];
  scale = ipt->scale;
  x = fRound(ipt->x);
  y = fRound(ipt->y);  
  desc = ipt->descriptor;

  if (bUpright)
  {
    co = 1;
    si = 0;
  }
  else
  {
    co = cos(ipt->orientation);
    si = sin(ipt->orientation);
  }

  i = -8;

  //choose whether to use haar or haarContour
  if (int_con==NULL){
    //Calculate descriptor for this interest point
    while(i < 12)
    {
      j = -8;
      i = i-4;
  
      cx += 1.f;
      cy = -0.5f;

      while(j < 12) 
      {
        dx=dy=mdx=mdy=0.f;
        cy += 1.f;

        j = j - 4;
  
        ix = i + 5;
        jx = j + 5;

        xs = fRound(x + ( -jx*scale*si + ix*scale*co));
        ys = fRound(y + ( jx*scale*co + ix*scale*si));

        for (int k = i; k < i + 9; ++k) 
        {
          for (int l = j; l < j + 9; ++l) 
          {
            //Get coords of sample point on the rotated axis
            sample_x = fRound(x + (-l*scale*si + k*scale*co));
            sample_y = fRound(y + ( l*scale*co + k*scale*si));
  
            //Get the gaussian weighted x and y responses
            gauss_s1 = gaussian(xs-sample_x,ys-sample_y,2.5f*scale);
            rx = haarX(sample_y, sample_x, 2*fRound(scale));
            ry = haarY(sample_y, sample_x, 2*fRound(scale));

            //Get the gaussian weighted x and y responses on rotated axis
            rrx = gauss_s1*(-rx*si + ry*co);
            rry = gauss_s1*(rx*co + ry*si);
  
            dx += rrx;
            dy += rry;
            mdx += fabs(rrx);
            mdy += fabs(rry);

          }  
        }

        //Add the values to the descriptor vector
        gauss_s2 = gaussian(cx-2.0f,cy-2.0f,1.5f);

        desc[count++] = dx*gauss_s2;
        desc[count++] = dy*gauss_s2;
        desc[count++] = mdx*gauss_s2;
        desc[count++] = mdy*gauss_s2;

        len += (dx*dx + dy*dy + mdx*mdx + mdy*mdy) * gauss_s2*gauss_s2;

        j += 9;
      }
      i += 9;
    }
  }
  else {
    //Calculate descriptor for this interest point
    while(i < 12)
    {
      j = -8;
      i = i-4;
  
      cx += 1.f;
      cy = -0.5f;

      while(j < 12) 
      {
        dx=dy=mdx=mdy=0.f;
        cy += 1.f;

        j = j - 4;
  
        ix = i + 5;
        jx = j + 5;

        xs = fRound(x + ( -jx*scale*si + ix*scale*co));
        ys = fRound(y + ( jx*scale*co + ix*scale*si));

        for (int k = i; k < i + 9; ++k) 
        {
          for (int l = j; l < j + 9; ++l) 
          {
            //Get coords of sample point on the rotated axis
            sample_x = fRound(x + (-l*scale*si + k*scale*co));
            sample_y = fRound(y + ( l*scale*co + k*scale*si));
  
            //Get the gaussian weighted x and y responses
            gauss_s1 = gaussian(xs-sample_x,ys-sample_y,2.5f*scale);
            rx = haarXContour(sample_y, sample_x, 2*fRound(scale), int_con);
            ry = haarYContour(sample_y, sample_x, 2*fRound(scale), int_con);

            //Get the gaussian weighted x and y responses on rotated axis
            rrx = gauss_s1*(-rx*si + ry*co);
            rry = gauss_s1*(rx*co + ry*si);
  
            dx += rrx;
            dy += rry;
            mdx += fabs(rrx);
            mdy += fabs(rry);

          }  
        }

        //Add the values to the descriptor vector
        gauss_s2 = gaussian(cx-2.0f,cy-2.0f,1.5f);

        desc[count++] = dx*gauss_s2;
        desc[count++] = dy*gauss_s2;
        desc[count++] = mdx*gauss_s2;
        desc[count++] = mdy*gauss_s2;

        len += (dx*dx + dy*dy + mdx*mdx + mdy*mdy) * gauss_s2*gauss_s2;

        j += 9;
      }
      i += 9;
    }
  }

  //Convert to Unit Vector
  len = sqrt(len);
  for(int i = 0; i < 64; ++i)
    desc[i] /= len;

}

//-------------------------------------------------------

//! Describe all features in the supplied vector
void Surf::getDescriptorsGlobal(bool upright, IplImage* int_con, const int init_sample)
{
  // Check there are Ipoints to be described
  if (!ipts.size()) return;

  // Get the size of the vector for fixed loop bounds
  int ipts_size = (int)ipts.size();

  if (upright)
  {
    // U-SURF loop just gets descriptors
    for (int i = 0; i < ipts_size; ++i)
    {
      // Set the Ipoint to be described
      index = i;

      // Extract upright (i.e. not rotation invariant) descriptors
      getDescriptor(true, int_con);
    }
  }
  else
  {
    // Main SURF-64 loop assigns orientations and gets descriptors
    for (int i = 0; i < ipts_size; ++i)
    {
      // Set the Ipoint to be described
      index = i;

      // Assign Orientations and extract rotation invariant descriptors
      getOrientationGlobal(int_con, init_sample);
      getDescriptorGlobal(false, int_con);
    }
  }
}

//-------------------------------------------------------

//! Determine the global orientation of the image at the given scale
void Surf::getOrientationGlobal(IplImage* int_con, const int init_sample)
{
  // setup
  Ipoint *ipt = &ipts[index];
  //int scale = fRound(ipt->scale);
  //std::cout<<scale<<std::endl;
  
  if (oris.count(-1)){
    ipt->orientation = oris[-1];
    return;
  }

  float sumScaleOris = 0.f;
  const int NUMSCALES = 10; 

  for (int scale = 2; scale < 2+NUMSCALES; scale++){

  

  // have we already calculated the global orientation at this scale?
  //else {
    // round scale to integer
    const int s = fRound(scale);
    // determine number of points in the image at this scale
    int width = img->width;
    int height = img->height;
    // what is init_sample? do we need this at all?
    int w = width/init_sample/s;
    int h = height/init_sample/s;

    // TEMP
    int sfactor = 1;
    w = width/(sfactor*s);
    h = height/(sfactor*s);
  
    const int numpoints = w*h;
  
    std::vector<float> resX(numpoints), resY(numpoints), Ang(numpoints);
    float totX=0.f, totY=0.f;

    // decide whether to use haar or haarContour
    if (int_con==NULL){
      // calculate haar response for the entire image at this scale
      for (int x=0; x<w; x++){
        for (int y=0; y<h; y++){
          // calculate wavelet response at this point and scale
          resX[x*h+y] = haarX(x*s*sfactor, y*s*sfactor, sfactor*s);
          totX+=resX[x*h+y];
          resY[x*h+y] = haarY(x*s*sfactor, y*s*sfactor, sfactor*s);
          totY+=resY[x*h+y];
          // angle of the gradient (up from x-axis)
          Ang[x*h+y] = getAngle(resX[x*h+y], resY[x*h+y]);

          //std::cout<<"x: "<<x*s*2<<"\ty: "<<y*s*2<<std::endl;
        }
      }
    }
    else {
      // calculate haar response for the entire image at this scale
      for (int x=0; x<w; x++){
        for (int y=0; y<h; y++){
          // calculate wavelet response at this point and scale
          resX[x*h+y] = haarXContour(x*s*sfactor, y*s*sfactor, sfactor*s, int_con);
	  if (std::isfinite(resX[x*h+y]))
            totX+=resX[x*h+y];
          resY[x*h+y] = haarYContour(x*s*sfactor, y*s*sfactor, sfactor*s, int_con);
	  if (std::isfinite(resY[x*h+y]))
            totY+=resY[x*h+y];
          // angle of the gradient (up from x-axis)
	  if (std::isfinite(resX[x*h+y]) && std::isfinite(resY[x*h+y]))
            Ang[x*h+y] =  getAngle(resX[x*h+y], resY[x*h+y]);
          else
            Ang[x*h+y] = FLT_MAX;

          //std::cout<<"x: "<<x*s*2<<"\ty: "<<y*s*2<<std::endl;
        }
      }
    }
    //std::cout<<" Angle from totals: "<<getAngle(totX,totY)<<" at scale "<<scale<<std::endl;
    
    float orientation = 0.f;

    // calculate the dominant direction 
    float sumX=0.f, sumY=0.f;
    float max=0.f;
    float ang1=0.f, ang2=0.f;

    // loop slides pi/3 window from ang1 to ang2
    // increments of .15 radians
    for(ang1 = 0; ang1 < 2*pi;  ang1+=0.15f) {
      ang2 = ( ang1+pi/3.0f > 2*pi ? ang1-5.0f*pi/3.0f : ang1+pi/3.0f);

      sumX = sumY = 0.f; 
      for(unsigned int k = 0; k < Ang.size(); ++k) 
      {
        // get angle from the x-axis of the sample point
        const float & ang = Ang[k];
  
        // determine whether the point's gradient is within the window
        // if so, add its x and y componenets to our sums 
        if (ang1 < ang2 && ang1 < ang && ang < ang2) 
        {
          sumX+=resX[k];  
          sumY+=resY[k];
        } 
        else if (ang2 < ang1 && 
          ((ang > 0 && ang < ang2) || (ang > ang1 && ang < 2*pi) )) 
        {
          sumX+=resX[k];  
          sumY+=resY[k];
        }
      }
  
      // if the vector produced from this window is longer than all 
      // previous vectors then this forms the new dominant direction
      if (sumX*sumX + sumY*sumY > max) 
      {
        // store largest orientation
        max = sumX*sumX + sumY*sumY;
        orientation = getAngle(sumX, sumY);
  
        //std::cout<<" New maximum: "<<max<<" at orientation "<<orientation<<std::endl;
      }
    }
    
    //std::cout<<" Image orientation is: "<<orientation<<" at scale "<<scale<<"\n"<<std::endl;

 //   orientation = getAngle(totX,totY);

    oris[scale] = orientation;

    //ipt->orientation = orientation;
    
    sumScaleOris += orientation;
  } //calc new orientation
  oris[-1] = sumScaleOris/(2+NUMSCALES);
  ipt->orientation = oris[-1];
}

//-------------------------------------------------------

//! Get the modified descriptor. See Agrawal ECCV 08
//! Modified descriptor contributed by Pablo Fernandez
//! Weighed, globally-oriented descriptor
void Surf::getDescriptorGlobal(bool bUpright, IplImage* int_con)
{
  int y, x, sample_x, sample_y, count=0;
  int i = 0, ix = 0, j = 0, jx = 0, xs = 0, ys = 0;
  float scale, *desc, dx, dy, mdx, mdy, co, si;
  float gauss_s1 = 0.f, gauss_s2 = 0.f;
  float rx = 0.f, ry = 0.f, rrx = 0.f, rry = 0.f, len = 0.f;
  float cx = -0.5f, cy = 0.f; //Subregion centers for the 4x4 gaussian weighting

  Ipoint *ipt = &ipts[index];
  scale = ipt->scale;
  x = fRound(ipt->x);
  y = fRound(ipt->y);  
  desc = ipt->descriptor;

  if (bUpright)
  {
    co = 1;
    si = 0;
  }
  else
  {
    co = cos(ipt->orientation);
    si = sin(ipt->orientation);
  }

  i = -8;

  //decide whether to use haar or haarContour
  if (int_con==NULL){
    //Calculate descriptor for this interest point
    while(i < 12)
    {
      j = -8;
      i = i-4;

      cx += 1.f;
      cy = -0.5f;

      while(j < 12) 
      {
        dx=dy=mdx=mdy=0.f;
        cy += 1.f;

        j = j - 4;

        ix = i + 5;
        jx = j + 5;

        xs = fRound(x + ( -jx*scale*si + ix*scale*co));
        ys = fRound(y + ( jx*scale*co + ix*scale*si));
  
        for (int k = i; k < i + 9; ++k) 
        {
          for (int l = j; l < j + 9; ++l) 
          {
            //Get coords of sample point on the rotated axis
            sample_x = fRound(x + (-l*scale*si + k*scale*co));
            sample_y = fRound(y + ( l*scale*co + k*scale*si));
  
            //Get the gaussian weighted x and y responses
            gauss_s1 = gaussian(xs-sample_x, ys-sample_y, 2.5f*scale);
            rx = haarX(sample_y, sample_x, 2*fRound(scale));
            ry = haarY(sample_y, sample_x, 2*fRound(scale));

            //Get the gaussian weighted x and y responses on rotated axis
            rrx = gauss_s1*(-rx*si + ry*co);
            rry = gauss_s1*(rx*co + ry*si);
  
            dx += rrx;
            dy += rry;
            mdx += fabs(rrx);
            mdy += fabs(rry);
          }
        }

        //Add the values to the descriptor vector
        gauss_s2 = gaussian(cx-2.0f,cy-2.0f,1.5f);

        desc[count++] = dx*gauss_s2;
        desc[count++] = dy*gauss_s2;
        desc[count++] = mdx*gauss_s2;
        desc[count++] = mdy*gauss_s2;

        len += (dx*dx + dy*dy + mdx*mdx + mdy*mdy) * gauss_s2*gauss_s2;

        j += 9;
      }
      i += 9;
    }
  }
  else {
    //Calculate descriptor for this interest point
    while(i < 12)
    {
      j = -8;
      i = i-4;

      cx += 1.f;
      cy = -0.5f;

      while(j < 12) 
      {
        dx=dy=mdx=mdy=0.f;
        cy += 1.f;

        j = j - 4;

        ix = i + 5;
        jx = j + 5;

        xs = fRound(x + ( -jx*scale*si + ix*scale*co));
        ys = fRound(y + ( jx*scale*co + ix*scale*si));
  
        for (int k = i; k < i + 9; ++k) 
        {
          for (int l = j; l < j + 9; ++l) 
          {
            //Get coords of sample point on the rotated axis
            sample_x = fRound(x + (-l*scale*si + k*scale*co));
            sample_y = fRound(y + ( l*scale*co + k*scale*si));
  
            //Get the gaussian weighted x and y responses
            gauss_s1 = gaussian(xs-sample_x, ys-sample_y, 2.5f*scale);
            rx = haarXContour(sample_y, sample_x, 2*fRound(scale), int_con);
            ry = haarYContour(sample_y, sample_x, 2*fRound(scale), int_con);

            //Get the gaussian weighted x and y responses on rotated axis
            rrx = gauss_s1*(-rx*si + ry*co);
            rry = gauss_s1*(rx*co + ry*si);
  
            dx += rrx;
            dy += rry;
            mdx += fabs(rrx);
            mdy += fabs(rry);
          }
        }

        //Add the values to the descriptor vector
        gauss_s2 = gaussian(cx-2.0f,cy-2.0f,1.5f);

        desc[count++] = dx*gauss_s2;
        desc[count++] = dy*gauss_s2;
        desc[count++] = mdx*gauss_s2;
        desc[count++] = mdy*gauss_s2;

        /*
        if (!std::isfinite(rx) || !std::isfinite(ry))
            std::cout<<dx*gauss_s2<<" -- "<<dy*gauss_s2<<" -- "<<mdx*gauss_s2<<" -- "<<mdy*gauss_s2<<std::endl;
        */

        if (std::isfinite(dx))
          len+=(dx*dx+mdx*mdx)*gauss_s2*gauss_s2;
        if (std::isfinite(dy))
          len+=(dy*dy+mdy*mdy)*gauss_s2*gauss_s2;

        j += 9;
      }
      i += 9;
    }
  }

  //Convert to Unit Vector
  len = sqrt(len);
  for(int i = 0; i < 64; ++i)
    desc[i] /= len;
}

//-------------------------------------------------------

//! Calculate the weight mask for the given offsets and orientation
float Surf::weightMask(int x, int y, float ori){
  float offsetAngle = getAngle((float) x,(float) y);
  float radDiff = ori-offsetAngle;
  //if the angle is within pi/8 of the orientation, first weight mask bin
  radDiff+=(0.125*pi);
  while (radDiff > 2*pi)
    radDiff-=2*pi;
  while (radDiff < 0)
    radDiff+=2*pi;
  if (radDiff<0.25*pi)      //bin 1
    return 1.0;
  else if (radDiff<0.5*pi)  //bin 2
    return 2.0;
  else if (radDiff<0.75*pi) //bin 3
    return 4.0;
  else if (radDiff<pi)      //bin 4
    return 8.0;
  else if (radDiff<1.25*pi) //bin 5
    return 16.0;
  else if (radDiff<1.5*pi)  //bin 6
    return 32.0;
  else if (radDiff<1.75*pi) //bin 7
    return 64.0;
  else                      //bin 8
    return 128.0;
}

//-------------------------------------------------------

//! Calculate the value of the 2d gaussian at x,y
inline float Surf::gaussian(int x, int y, float sig)
{
  return (1.0f/(2.0f*pi*sig*sig)) * exp( -(x*x+y*y)/(2.0f*sig*sig));
}

//-------------------------------------------------------

//! Calculate the value of the 2d gaussian at x,y
inline float Surf::gaussian(float x, float y, float sig)
{
  return 1.0f/(2.0f*pi*sig*sig) * exp( -(x*x+y*y)/(2.0f*sig*sig));
}

//-------------------------------------------------------

//! Calculate Haar wavelet responses in x direction
inline float Surf::haarX(int row, int column, int s)
{
  //recall integral starts at row col specified by 2nd, 3rd parameters, then runs # of rows and
  //columns specified by 4th and 5th parameters

  return BoxIntegral(img, row-s/2, column, s, s/2) 
    -1 * BoxIntegral(img, row-s/2, column-s/2, s, s/2);
}

//-------------------------------------------------------

//! Calculate Haar wavelet responses in y direction
inline float Surf::haarY(int row, int column, int s)
{
  return BoxIntegral(img, row, column-s/2, s/2, s) 
    -1 * BoxIntegral(img, row-s/2, column-s/2, s/2, s);
}

//-------------------------------------------------------

//! Calculate Haar wavelet responses in x direction for a contour
inline float Surf::haarXContour(int row, int column, int s, IplImage* int_con)
{
  float inverse_plus = 1.f/BoxIntegral(int_con, row-s/2, column, s, s/2);
  float inverse_minus = 1.f/BoxIntegral(int_con, row-s/2, column-s/2, s, s/2);
  return BoxIntegral(img, row-s/2, column, s, s/2) * inverse_plus
         -1 * BoxIntegral(img, row-s/2, column-s/2, s, s/2) * inverse_minus;
}

//-------------------------------------------------------

//! Calculate Haar wavelet responses in y direction for contour
inline float Surf::haarYContour(int row, int column, int s, IplImage* int_con)
{
  float inverse_plus = 1.f/BoxIntegral(int_con, row, column-s/2, s/2, s);
  float inverse_minus = 1.f/BoxIntegral(int_con, row-s/2, column-s/2, s/2, s);
  return BoxIntegral(img, row, column-s/2, s/2, s) * inverse_plus
         -1 * BoxIntegral(img, row-s/2, column-s/2, s/2, s) * inverse_minus;
}


//-------------------------------------------------------

//! Get the angle from the +ve x-axis of the vector given by (X Y)
float Surf::getAngle(float X, float Y)
{
  if(X > 0 && Y >= 0)
    return atan(Y/X);

  if(X < 0 && Y >= 0)
    return pi - atan(-Y/X);

  if(X < 0 && Y < 0)
    return pi + atan(Y/X);

  if(X > 0 && Y < 0)
    return 2*pi - atan(-Y/X);

  return 0;
}
