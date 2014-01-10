//Image.h extracted from Eriol Header by Olaf Hall-Holt, 2007-2011

#include "PixelLoc.h"
#include "Coord.h"
#include "Color.h"
#include "BBox.h"
#include <iostream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
using std::vector;
using std::ostream;
using std::istream;
using std::string;
using std::stringstream;
using std::ws;
#include <cassert>
#include <math.h>

#ifndef __IMAGE__
#define __IMAGE__

class Image {  // matrix we can draw on the screen
public:
  enum ChannelType { NONE_CHAN, BOOL_CHAN, CHAR_CHAN, SHORT_CHAN, INT_CHAN, FLOAT_CHAN, DOUBLE_CHAN, MMAP_CHAN };
  string name;
  bool isSimpleFormat;  // homemade format for a tiny segmentation
  bool readOnly, writeOnly;  // for mmap
  int mmapOffset, numBytes, fd;  // for mmap
  float imageScale;  // for display
  Coord windowUL[4];    // for display
  float perPixelScale;  // for display

private:
  unsigned width, height, numChannels;
  ChannelType channelType;
  vector<bool> boolData;
  vector<unsigned char> charData;
  vector<unsigned short> shortData;
  vector<int> intData;
  vector<float> floatData;
  vector<double> doubleData;
  unsigned char *mmapRoot, *mmapData;

  void assertInImage(unsigned x, unsigned y) const
          { assert(x < width); assert(y < height); }
public:
  static unsigned char clamp255(double d)
          { return (d<0) ? 0 : (d>255) ? 255 : round(d); }
  static unsigned char clamp65535(double d)
          { return (d<0) ? 0 : (d>65535) ? 65535 : round(d); }
  void initData();
  void initZoomParameters();
  void init(istream &, bool useMmap=false);
  Image(): isSimpleFormat(false), width(0), height(0),
     numChannels(0), channelType(NONE_CHAN) {}
  Image(unsigned w, unsigned h, ChannelType t, unsigned d):
     isSimpleFormat(false), width(w), height(h),
     numChannels(d), channelType(t) { initData(); initZoomParameters(); }
  Image(istream &);
  Image(const char *fname, bool useMmap=false);
  Image(unsigned w, unsigned h, unsigned d, const char *fname); // mmap write
  ~Image();

  unsigned getWidth() const { return width; }
  unsigned getHeight() const { return height; }
  unsigned getNumChannels() const { return numChannels; }
  ChannelType getChannelType() const { return channelType; }
  string basicImageInfo() const {
     stringstream ss;
     ss << name << " " << width << "x" << height;
     return ss.str(); }
  bool hasType(int n, ChannelType t) {
     return (n == static_cast<int>(numChannels)) && (t == channelType); }
  void *getData() {  // for low-level access
    //  low level access doesn't work correctly with bool, due to packing?
    // if ( BOOL_CHAN == channelType ) return &(boolData[0]); else
    if ( CHAR_CHAN == channelType ) return &(charData[0]);
    else if ( SHORT_CHAN == channelType ) return &(shortData[0]);
    else if ( INT_CHAN == channelType ) return &(intData[0]);
    else if ( FLOAT_CHAN == channelType ) return &(floatData[0]);
    else if ( DOUBLE_CHAN == channelType ) return &(doubleData[0]);
    else if ( MMAP_CHAN == channelType ) return mmapData;
    else assert(false);
  }
  const void *getData() const {
    //  low level access doesn't work correctly with bool, due to packing?
    // if ( BOOL_CHAN == channelType ) return &(boolData[0]); else
    if ( CHAR_CHAN == channelType ) return &(charData[0]);
    else if ( SHORT_CHAN == channelType ) return &(shortData[0]);
    else if ( INT_CHAN == channelType ) return &(intData[0]);
    else if ( FLOAT_CHAN == channelType ) return &(floatData[0]);
    else if ( DOUBLE_CHAN == channelType ) return &(doubleData[0]);
    else if ( MMAP_CHAN == channelType ) return mmapData;
    else assert(false);
  }
  bool getBoolPixel(PixelLoc p, unsigned channel=0) const {
    assert( BOOL_CHAN == channelType );
    assert( channel < numChannels );
    assertInImage(p.x,p.y);
    return boolData[numChannels*(p.y*width+p.x)+channel];
  }
  unsigned char getCharPixel(PixelLoc p, unsigned channel=0) const {
    assert( channel < numChannels );
    assertInImage(p.x,p.y);
    if ( CHAR_CHAN == channelType )
      return charData[numChannels*(p.y*width+p.x)+channel];
    else {
      assert( MMAP_CHAN == channelType);
      assert(readOnly);
      return mmapData[numChannels*(p.y*width+p.x)+channel];
    }
  }
  Color getPixel(PixelLoc p) const {
    assert( 3 == numChannels );
    assertInImage(p.x,p.y);
    if ( CHAR_CHAN == channelType )
      return Color(charData[numChannels*(p.y*width+p.x)+0],
                   charData[numChannels*(p.y*width+p.x)+1],
                   charData[numChannels*(p.y*width+p.x)+2] );
    else if ( SHORT_CHAN == channelType )
      return Color(shortData[numChannels*(p.y*width+p.x)+0],
                   shortData[numChannels*(p.y*width+p.x)+1],
                   shortData[numChannels*(p.y*width+p.x)+2] );
    else {
      assert( MMAP_CHAN == channelType);
      assert(readOnly);
      return Color(mmapData[numChannels*(p.y*width+p.x)+0],
                   mmapData[numChannels*(p.y*width+p.x)+1],
                   mmapData[numChannels*(p.y*width+p.x)+2] );
    }
  }
  unsigned short getShortPixel(PixelLoc p, unsigned channel=0) const {
    assert( SHORT_CHAN == channelType );
    assert( channel < numChannels );
    assertInImage(p.x,p.y);
    return shortData[numChannels*(p.y*width+p.x)+channel];
  }
  int getIntPixel(PixelLoc p, unsigned channel=0) const {
    assert( INT_CHAN == channelType );
    assert( channel < numChannels );
    assertInImage(p.x,p.y);
    return intData[numChannels*(p.y*width+p.x)+channel];
  }
  float getFloatPixel(PixelLoc p, unsigned channel=0) const {
    assert( FLOAT_CHAN == channelType );
    assert( channel < numChannels );
    assertInImage(p.x,p.y);
    return floatData[numChannels*(p.y*width+p.x)+channel];
  }
  double getDoublePixel(PixelLoc p, unsigned channel=0) const {
    assert( DOUBLE_CHAN == channelType );
    assert( channel < numChannels );
    assertInImage(p.x,p.y);
    return doubleData[numChannels*(p.y*width+p.x)+channel];
  }
  void setPixel(PixelLoc p, bool col, unsigned chan=0 ) {
    assert( BOOL_CHAN == channelType );
    assert( chan < numChannels );
    assertInImage(p.x,p.y);
    boolData[numChannels*(p.y*width+p.x)+chan] = col; }
  void setPixel(PixelLoc p, unsigned char col, unsigned chan=0 ) {
    assert( chan < numChannels );
    assertInImage(p.x,p.y);
    if ( CHAR_CHAN == channelType )
      charData[numChannels*(p.y*width+p.x)+chan] = col;
    else {
      assert( MMAP_CHAN == channelType);
      assert(writeOnly);
      mmapData[numChannels*(p.y*width+p.x)+chan] = col; } }
  void setPixel(PixelLoc p, unsigned short col, unsigned chan=0 ) {
    assert( SHORT_CHAN == channelType );
    assert( chan < numChannels );
    assertInImage(p.x,p.y);
    intData[numChannels*(p.y*width+p.x)+chan] = col; }
  void setPixel(PixelLoc p, int col, unsigned chan=0 ) {
    assert( INT_CHAN == channelType );
    assert( chan < numChannels );
    assertInImage(p.x,p.y);
    intData[numChannels*(p.y*width+p.x)+chan] = col; }
  void setPixel(PixelLoc p, float col, unsigned chan=0 ) {
    assert( FLOAT_CHAN == channelType );
    assert( chan < numChannels );
    assertInImage(p.x,p.y);
    floatData[numChannels*(p.y*width+p.x)+chan] = col; }
  void setPixel(PixelLoc p, double col, unsigned chan=0 ) {
    assert( DOUBLE_CHAN == channelType );
    assert( chan < numChannels );
    assertInImage(p.x,p.y);
    doubleData[numChannels*(p.y*width+p.x)+chan] = col; }
  void setPixel(PixelLoc p, Color col) {
    assert( 3 == numChannels );
    assertInImage(p.x,p.y);
    if ( CHAR_CHAN == channelType ) {
      charData[numChannels*(p.y*width+p.x)+0] = clamp255(col.r);
      charData[numChannels*(p.y*width+p.x)+1] = clamp255(col.g);
      charData[numChannels*(p.y*width+p.x)+2] = clamp255(col.b);
    } else if ( SHORT_CHAN == channelType ) {
      shortData[numChannels*(p.y*width+p.x)+0] = clamp65535(col.r);
      shortData[numChannels*(p.y*width+p.x)+1] = clamp65535(col.g);
      shortData[numChannels*(p.y*width+p.x)+2] = clamp65535(col.b);
    } else {
      assert(MMAP_CHAN == channelType);
      assert(writeOnly);
      mmapData[numChannels*(p.y*width+p.x)+0] = clamp255(col.r);
      mmapData[numChannels*(p.y*width+p.x)+1] = clamp255(col.g);
      mmapData[numChannels*(p.y*width+p.x)+2] = clamp255(col.b);
    }
  }
  void flipX();
  void flipY();
  void flipXY();
  void print(ostream &os, bool justHeader=false, bool forceChar=false) const;
  void print(const char *) const;
  void draw(Coord offset=noCoord, int imageIndex=0) const;
  void readFromFrameBuffer();
  Color getWindowAvg(const BBox &);
  double getWindowAvgGray(const BBox &);
  void lowPassFilter();
  void edgeDetect(bool vert);
  void subSample(bool saveFile=true);
};


#endif //__IMAGE__