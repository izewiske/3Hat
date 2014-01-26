#ifndef CONTOURMAT_H
#define CONTOURMAT_H
#include <cv.h>
#include <vector>
#include "responselayer.h"
//#include "Contour.h"

class ContourMat{
  //0 - not in contour
  //1 - in contour
  std::vector<bool> contour;

 public:
  //oriented bounding box Mat (holds image data)
  cv::Mat obb;

  bool inContour(int row, int col);
  bool inContour(int row, int col, ResponseLayer* src);
};

#endif
