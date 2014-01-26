#include "ContourMat.h"
#include <cv.h>

//(included in ContourMat.h)
//#include "responselayer.h"
//#include <vector>

bool ContourMat::inContour(int row, int col){
  //placeholder (obviously)
  return true;

  /*
  return contour[row*obb.cols+col];`
  */
}

bool ContourMat::inContour(int row, int col, ResponseLayer* src){
  //placeholder (obviously)
  return true;

  /*
  int scale = obb.cols / src->width;
  return contour[row*scale*obb.cols+col*scale];
  */
}
