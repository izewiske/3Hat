// by Olaf Hall-Holt, 2007-2014
#ifndef _ERIOL_
#define _ERIOL_

#include<iostream>
#include<algorithm>
#include<vector>
#include<string>
#include<sstream>
using std::vector;
using std::ostream;
using std::istream;
using std::string;
using std::stringstream;
using std::ws;
#include<cassert>
#include<math.h>

template<typename T>
void makeUnique(vector<T> &v)
{
  sort(v.begin(), v.end());
  typename vector<T>::iterator iter = unique(v.begin(), v.end());
  v.erase(iter, v.end());
}
template<class T>
bool memberOf( T j, const vector<T> &possible )
{
  return find(possible.begin(), possible.end(), j) != possible.end();
}
template<class T>
ostream& operator<<(ostream &os, const vector<T> &all)
{
  for(unsigned i=0, i_end=all.size(); i<i_end; ++i)
    os << all[i] << ' ';
  return os;
}
template<class T, class T2>
ostream& operator<<(ostream &os, const pair<T, T2> &both)
{
  os << both.first << ' ' << both.second;
  return os;
}
template<class T>
void rotateVector(vector<T> &pt, unsigned n)
{
  assert( n < pt.size() );
  vector<T> tmp(pt);
  pt.clear();
  pt.insert(pt.end(), tmp.begin()+n, tmp.end());
  pt.insert(pt.end(), tmp.begin(), tmp.begin()+n);
}
template<class T>
bool removeIndicesFromVector(const vector<int> &toRemove, vector<T> &v)
{
  // assumes that toRemove is sorted
  if ( ! toRemove.empty() ) {
    int j = 0, j_end = v.size();
    vector<T> tmp;
    tmp.reserve(v.size() - toRemove.size());
    for ( unsigned i=0,i_end=toRemove.size(); i<i_end; ++i ) {
      if ( i+1<i_end ) assert( toRemove[i] < toRemove[i+1] );
      while ( j < toRemove[i] ) tmp.push_back(v[j++]);
      ++j;
    }
    while ( j < j_end ) tmp.push_back(v[j++]);
    v = tmp;
    return true;
  }
  return false;
}

class Color;
double hueDiff(const Color &a, const Color &b);
double hueDiff2(Color a, Color b);
double robustHueDiff(const Color &a, const Color &b);
struct Color {
  double r,g,b;
  Color() { r = 0.; g = 0.; b = 0.; }
  Color(double _r, double _g, double _b) : r(_r), g(_g), b(_b) {}
  Color(istream &is) { is >> r >> g >> b; }
  Color& operator=(const int& i) { 
    r = g = b = i;
    return *this; }
  Color operator+(const Color & p) const {
    return Color(r + p.r, g + p.g, b + p.b); }
  Color operator-(const Color & p) const {
    return Color(r - p.r, g - p.g, b - p.b); }
  Color operator/(double denom) const {
    assert(denom != 0);
    return Color(r/denom,g/denom,b/denom); }
  Color & operator+=(const Color &p) {
    r += p.r; g += p.g; b += p.b;
    return *this; }
  Color & operator/=(double denom) {
    assert(denom != 0);
    r /= denom; g /= denom; b /= denom;
    return *this; }
  double operator*(const Color &p) const {
    return r*p.r + g*p.g + b*p.b; }
  double mag() const {
    return sqrt(r*r + g*g + b*b); }
  double mag2() const {
    return r*r + g*g + b*b; }
  void normalize() {
    double m = mag();
    if ( m != 0 ) { r /= m; g /= m; b /= m; } }
  double& operator[](int channel) {
    if (channel == 0) return r;
    else if (channel == 1) return g;
    else return b; }
  double operator[](int channel) const {
    if (channel == 0) return r;
    else if (channel == 1) return g;
    else return b; }
  friend bool isZero(const Color &c) {
    return 0==c.r && 0==c.g && 0==c.b; }
  friend Color operator*(double s, const Color &c) {
    return Color(s*c.r, s*c.g, s*c.b); }
  friend double dist(const Color &a, const Color &b) {
    // return (a-b).mag();
    return hueDiff(a,b);
  }
  friend double dist2(const Color &a, const Color &b) {
    // return (a-b).mag2();
    return hueDiff2(a,b);
  }
  friend double absDiff(const Color &a, const Color &b) {
    return fabs(a.r-b.r) + fabs(a.g-b.g) + fabs(a.b-b.b); }
  friend Color colorDiff(Color a, Color b) {
    return Color(128,128,128) + (b-a); }
  friend bool isNoColor(Color c) {
    return (c.r==-1e10) && (c.g==-1e10) && (c.b==-1e10); } // must equal noColor
  friend bool approxEqual(const Color &a, const Color &b, double thresh=.01) {
    return fabs(a.r - b.r) < thresh &&
           fabs(a.g - b.g) < thresh &&
           fabs(a.b - b.b) < thresh; }
  friend ostream& operator<<(ostream &os, const Color &p) {
    os << p.r << " " << p.g << " " << p.b;
    return os; }
  friend istream & operator>>(istream& istr, Color &c) {
    c = Color(istr);  return istr; }
};
void sendColor(double, double, double);
void sendColor(Color);
static const Color WHITE(255,255,255);
static const Color BLACK(0,0,0);
static const Color RED(255,0,0);
static const Color ORANGE(255,125,0);
static const Color YELLOW(255,255,0);
static const Color GREEN(0,255,0);
static const Color CYAN(0,255,255);
static const Color BLUE(0,0,255);
static const Color PURPLE(255,0,255);
static const Color DARKRED(125,0,0);
static const Color DARKGREEN(0,125,0);
static const Color DARKBLUE(0,0,125);
static const Color DARKYELLOW(125,125,0);
static const Color DARKPURPLE(125,0,125);
static const Color GREY(125,125,125);
static const Color noColor(-1e10,-1e10,-1e10);

double variance(const vector<Color> &col);

struct PixelLoc {
  int x,y;
  PixelLoc() : x(-1), y(-1) {}
  PixelLoc(int _x, int _y) : x(_x), y(_y) {}
  PixelLoc(istream &is) { is >> x; assert(',' == is.get()); is >> y; }
  int &operator[](int i) {
    return (0==i) ? x : y; }
  int operator[](int i) const {
    return (0==i) ? x : y; }
  PixelLoc &operator+=(const PixelLoc &b) {
    x += b.x; y += b.y; return *this; }
  PixelLoc &operator-=(const PixelLoc &b) {
    x -= b.x; y -= b.y; return *this; }
  PixelLoc &operator*=(int s) {
    x *= s; y *= s; return *this; }
  bool operator<(const PixelLoc &b) const {
    if (x!=b.x) return x<b.x; return y<b.y; }
  bool operator==(const PixelLoc &b) const {
    return x==b.x && y==b.y; }
  bool operator!=(const PixelLoc &b) const {
    return !(*this==b); }
  friend PixelLoc operator+(const PixelLoc &a, const PixelLoc &b) {
    return PixelLoc(a.x + b.x, a.y + b.y); }
  friend PixelLoc operator-(const PixelLoc &a, const PixelLoc &b) {
    return PixelLoc(a.x - b.x, a.y - b.y); }
  friend PixelLoc operator*(int s, const PixelLoc &a) {
    return PixelLoc(s*a.x, s*a.y); }
  friend double dist2(PixelLoc a, PixelLoc b) {
    double dx = a.x - b.x, dy = a.y - b.y;
    return dx*dx + dy*dy; }
  friend double dist(PixelLoc a, PixelLoc b) {
    return sqrt(dist2(a,b)); }
  friend bool isNilLoc(PixelLoc a) {
    return (a.x==-1000000000) && (a.y==-1000000000); }  // must equal nilLoc
  friend ostream& operator<<(ostream &os, const PixelLoc &p) {
    os << p.x << "," << p.y;
    return os; }
  friend istream & operator>>(istream& istr, PixelLoc &p);
};
static const PixelLoc UP(0, -1);
static const PixelLoc LEFT(-1, 0);
static const PixelLoc DOWN(0, 1);
static const PixelLoc RIGHT(1, 0);
static const PixelLoc nilLoc(-1000000000,-1000000000);

vector<PixelLoc> getNbr8(PixelLoc pos);

struct Coord {
  double x,y;
  Coord(): x(0), y(0) {}
  Coord(double _x, double _y): x(_x), y(_y) {}
  Coord(istream &is);

  double dot(Coord c) const {
    return x*c.x+y*c.y; }
  double operator*(Coord c) const {
    return dot(c); }
  Coord& operator+=(const Coord& c) {
    x += c.x;
    y += c.y;
    return *this; }
  Coord& operator-=(const Coord& c) {
    x -= c.x;
    y -= c.y;
    return *this; }
  Coord& operator*=(double s) {
    x *= s;
    y *= s;
    return *this; }
  double mag() const {
    return sqrt(x*x + y*y); }
  void normalize() {
    double len = mag();
    if (len) { x /= len; y /= len; } }
  void moveToward(Coord intersection, double howMuch);
  friend bool closeTo(double a, double b, double thresh=0.00001) {
    return fabs(a-b) < thresh; }
  bool closeTo(Coord b, double thresh=0.00001) const {
    return dist2(*this, b) < thresh; }
  void draw() const;

  friend Coord operator+(Coord a, Coord b) {
    return Coord(a.x + b.x, a.y + b.y); }
  friend Coord operator-(Coord a, Coord b) {
    return Coord(a.x - b.x, a.y - b.y); }
  friend Coord operator*(double s, Coord c) {
    return Coord(s*c.x, s*c.y); }
  friend Coord operator/(Coord c,double s) {
    assert(s != 0);
    return Coord(c.x/s, c.y/s); }
  friend double dist2(Coord a, Coord b) {
    double dx = a.x - b.x, dy = a.y - b.y;
    return dx*dx + dy*dy; }
  friend double dist(Coord a, Coord b) {
    return sqrt(dist2(a,b)); }
  friend Coord quarterTurn(Coord c, bool inv=false) {
    return inv ? Coord(c.y, -c.x) : Coord(-c.y, c.x); }
  friend PixelLoc asPixelLoc(Coord c) {
    return PixelLoc( floor(c.x),  floor(c.y) ); }
  friend Coord asCoord(PixelLoc a) {
    return Coord(a.x, a.y); }
  friend Coord asShiftedCoord(const PixelLoc &a) {
    return Coord(a.x + .5, a.y + .5); }
  friend bool isNoCoord(Coord c) {
    return (c.x==1e10) && (c.y==1e10); }  // must equal noCoord
// friend Coord warpPoint(Coord p, double *mat);
  friend ostream & operator<<(ostream& ostr, const Coord& c) {
    ostr << c.x << "," << c.y;
    return ostr; }
  friend istream & operator>>(istream& is, Coord &c) { is >> ws;
    // if ( EOF == is.peek() ) is.setstate(ios::badbit); else
    c = Coord(is);
    return is; }
};
static const Coord noCoord(1e+10, 1e+10); // see isNoCoord above

static const double EPSILON = 1e-10;

double angleOf(Coord dir);
double angleOf(Coord a, Coord b, Coord c);
double shorterAngleDiff(double, double);
double fullAngleDiff(double, double);
double jitter(double width=0.00001);
Coord asCoord(PixelLoc a);
Coord asShiftedCoord(const PixelLoc &a);
bool closeTo(double a, double b, double thresh);
bool containsPoint(const vector<Coord> &bdy, Coord query);

struct BBox {
  PixelLoc ul;
  int w,h;
  BBox(): ul(0,0), w(-1), h(-1) {}
  BBox(PixelLoc pos): ul(pos), w(1), h(1) {}
  BBox(PixelLoc pos, int _w, int _h): ul(pos), w(_w), h(_h) {}
  BBox(PixelLoc pos, PixelLoc pos2): ul(pos), w(pos2.x-pos.x+1), h(pos2.y-pos.y+1) {}
  BBox(const vector<Coord> &pt): ul(0,0), w(-1), h(-1) {
    if ( !pt.empty() ) {
      ul = asPixelLoc(pt[0]); w = h = 1;
      for ( unsigned i=1,i_end=pt.size(); i<i_end; ++i ) expand(pt[i]);
    } }
  BBox(const BBox &a, const BBox &b): ul(a.ul), w(a.w), h(a.h) {
    expand(b.ul);
    expand(b.ul + PixelLoc(b.w-1,b.h-1)); }
  int minX() const { return ul.x; }
  int minY() const { return ul.y; }
  int maxX() const { return ul.x+w-1; }
  int maxY() const { return ul.y+h-1; }
  PixelLoc minXY() const { return ul; }
  PixelLoc maxXY() const { return PixelLoc( maxX(), maxY() ); }
  int xLim() const { return ul.x + w; }
  int yLim() const { return ul.y + h; }
  bool contains( PixelLoc p ) const {
    return p.x >= ul.x && p.x < xLim() &&
           p.y >= ul.y && p.y < yLim(); }
  bool contains( Coord p ) const {
    return p.x >= ul.x && p.x < xLim() &&
           p.y >= ul.y && p.y < yLim(); }
  void incrX() { ++w; }
  void incrY() { ++h; }
  void decrX() { --ul.x; ++w; }
  void decrY() { --ul.y; ++h; }
  void expand(Coord c) {
    if ( c.x < ul.x ) {
      w += (ul.x - floor(c.x));
      ul.x = floor(c.x);
    }
    if ( c.y < ul.y ) {
      h += (ul.y - floor(c.y));
      ul.y = floor(c.y);
    }
    if ( c.x >= ul.x+w ) w = ceil(c.x) - ul.x;
    if ( c.y >= ul.y+h ) h = ceil(c.y) - ul.y;
  }
  void expand(PixelLoc p) { expand( asShiftedCoord(p) ); }
  void expand(BBox b) { expand(b.ul);
                        expand(b.ul + PixelLoc(b.w-1,b.h-1)); }
  void expandByOne() { --ul.x; --ul.y; w+=2; h+=2; }
  bool overlaps(int amin, int amax, int bmin, int bmax) const
    { return (amin <= bmax) && (bmin <= amax); }
  bool overlaps(const BBox &bb) const {
    return overlaps(bb.minX(), bb.maxX(), minX(), maxX()) &&
           overlaps(bb.minY(), bb.maxY(), minY(), maxY()); }
  bool contains(const BBox &bb) const {
    return contains(bb.minXY()) && contains(bb.maxXY()); }
  vector<PixelLoc> pixelsInside() const { vector<PixelLoc> res;
    PixelLoc pos;
    for ( pos.y = minY(); pos.y <= maxY(); ++pos.y )
    for ( pos.x = minX(); pos.x <= maxX(); ++pos.x )
      res.push_back(pos);
    return res; }
  void draw() const;
  friend ostream& operator<<(ostream &os, const BBox &b) {
    os << b.w << "x" << b.h << "+" << b.ul.x << "+" << b.ul.y;
    return os; }
};

class Segment;
class Arc;
class Curvelet;
BBox computeBBox(const Segment &, double offset=0);
BBox computeBBox(const Arc &, double offset=0);
BBox computeBBox(const Curvelet &);
bool overlaps(const BBox &, const Segment &);

class Image {  // matrix we can draw on the screen
public:
  enum ChannelType { NONE_CHAN, BOOL_CHAN, CHAR_CHAN, SHORT_CHAN, INT_CHAN, FLOAT_CHAN, DOUBLE_CHAN, MMAP_CHAN };
  string name;
  bool isSimpleFormat;  // homemade format for a tiny segmentation
  bool active;  // if not active, we have the name, but we still need to load the data
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
  bool init(istream &, bool useMmap=false, bool startActive=true);
  Image(): isSimpleFormat(false), active(false), width(0), height(0),
     numChannels(0), channelType(NONE_CHAN) {}
  Image(unsigned w, unsigned h, ChannelType t, unsigned d):
     isSimpleFormat(false), active(true), width(w), height(h),
     numChannels(d), channelType(t) { initData(); initZoomParameters(); }
  Image(istream &);
  Image(const char *fname, bool useMmap=false, bool startActive=true);
  Image(unsigned w, unsigned h, unsigned d, const char *fname); // mmap write
  ~Image();
  void attemptToLoad(bool useMmap, bool startActive);
  void makeActive();

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
    assert(active);
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
    assert(active);
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
    assert(active);
    assert( BOOL_CHAN == channelType );
    assert( channel < numChannels );
    assertInImage(p.x,p.y);
    return boolData[numChannels*(p.y*width+p.x)+channel];
  }
  unsigned char getCharPixel(PixelLoc p, unsigned channel=0) const {
    assert(active);
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
    assert(active);
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
    assert(active);
    assert( SHORT_CHAN == channelType );
    assert( channel < numChannels );
    assertInImage(p.x,p.y);
    return shortData[numChannels*(p.y*width+p.x)+channel];
  }
  int getIntPixel(PixelLoc p, unsigned channel=0) const {
    assert(active);
    assert( INT_CHAN == channelType );
    assert( channel < numChannels );
    assertInImage(p.x,p.y);
    return intData[numChannels*(p.y*width+p.x)+channel];
  }
  float getFloatPixel(PixelLoc p, unsigned channel=0) const {
    assert(active);
    assert( FLOAT_CHAN == channelType );
    assert( channel < numChannels );
    assertInImage(p.x,p.y);
    return floatData[numChannels*(p.y*width+p.x)+channel];
  }
  double getDoublePixel(PixelLoc p, unsigned channel=0) const {
    assert(active);
    assert( DOUBLE_CHAN == channelType );
    assert( channel < numChannels );
    assertInImage(p.x,p.y);
    return doubleData[numChannels*(p.y*width+p.x)+channel];
  }
  void setPixel(PixelLoc p, bool col, unsigned chan=0 ) {
    assert(active);
    assert( BOOL_CHAN == channelType );
    assert( chan < numChannels );
    assertInImage(p.x,p.y);
    boolData[numChannels*(p.y*width+p.x)+chan] = col; }
  void setPixel(PixelLoc p, unsigned char col, unsigned chan=0 ) {
    assert(active);
    assert( chan < numChannels );
    assertInImage(p.x,p.y);
    if ( CHAR_CHAN == channelType )
      charData[numChannels*(p.y*width+p.x)+chan] = col;
    else {
      assert( MMAP_CHAN == channelType);
      assert(writeOnly);
      mmapData[numChannels*(p.y*width+p.x)+chan] = col; } }
  void setPixel(PixelLoc p, unsigned short col, unsigned chan=0 ) {
    assert(active);
    assert( SHORT_CHAN == channelType );
    assert( chan < numChannels );
    assertInImage(p.x,p.y);
    intData[numChannels*(p.y*width+p.x)+chan] = col; }
  void setPixel(PixelLoc p, int col, unsigned chan=0 ) {
    assert(active);
    assert( INT_CHAN == channelType );
    assert( chan < numChannels );
    assertInImage(p.x,p.y);
    intData[numChannels*(p.y*width+p.x)+chan] = col; }
  void setPixel(PixelLoc p, float col, unsigned chan=0 ) {
    assert(active);
    assert( FLOAT_CHAN == channelType );
    assert( chan < numChannels );
    assertInImage(p.x,p.y);
    floatData[numChannels*(p.y*width+p.x)+chan] = col; }
  void setPixel(PixelLoc p, double col, unsigned chan=0 ) {
    assert(active);
    assert( DOUBLE_CHAN == channelType );
    assert( chan < numChannels );
    assertInImage(p.x,p.y);
    doubleData[numChannels*(p.y*width+p.x)+chan] = col; }
  void setPixel(PixelLoc p, Color col) {
    assert(active);
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
  void zeroOtherChannels(unsigned chan);
  void addChromaticAberration();
  void flipX();
  void flipY();
  void flipXY();
  void print(ostream &os, bool justHeader=false, bool forceChar=false) const;
  void print(const char *) const;
  void draw(Coord offset=noCoord, int windowIndex=0, float scale=-1) const;
  void readFromFrameBuffer();
  Color getWindowAvg(const BBox &);
  double getWindowAvgGray(const BBox &);
  void lowPassFilter();
  void edgeDetect(bool vert);
  void subSample(bool saveFile=true);
  void squeeze(int squeezeFactor, bool vert, bool saveFile=true);
  Image *createCrop(const BBox &bbox);
};

namespace ImageUI {

  extern vector<Image *> allImages;
  extern int currImageIndex[4];
  extern Image *faceMapImage, *faceMapImageR, *sourceColorImage;
  // extern float imageScale0;  // there may also exist others in Lattice
  // extern Coord windowUL0;    //   to control other images separately
  extern double windowPixelWidth[4], windowPixelHeight[4];
  extern double origWindowPixelWidth[4], origWindowPixelHeight[4];
  extern double snapRange;
  extern bool tetherMotion, dualView, quadView, isFullScreen[4];
  extern BBox bbox;

  enum UIMode {NORMAL, ADJUST, ADD_VERTEX, SELECT, CORRESP, LATTICE, PIXELSEG};

  float &imageScale(int windowIndex);
  Coord &windowUL(int windowIndex);
  double getWindowPixelWidth(int windowIndex);
  double getWindowPixelHeight(int windowIndex);
  void addImage(const char *f, bool startActive=true);
  Image*& currImage(int windowIndex=0);
  int currImageWidth(int windowIndex=0);
  int currImageHeight(int windowIndex=0);
  Image *getFaceMap();
  void drawImage(int windowIndex);
  void drawImgBBox();
  bool isValidIndex(int);
  void setImageIndexChar(char, int windowIndex);
  void setImageIndex(int, bool toTitle=false, int windowIndex=0);
  void putImageNameInWindowTitle(int windowIndex);
  void setImage(Image *img, int windowIndex=0);
  void setImageByPrefix(const string &fname, int windowIndex=0);
  void incrImageIndex(bool fwd=true, int windowIndex=0);
  Coord asImageCoord(int x, int y, int windowIndex=0);
  void cropImage();
  void frameGrab(const char *fname=0);
}

class PixelEdge;
bool inImage(PixelLoc);
bool inImage(const Image *, PixelLoc);
bool inImage(Coord);
bool inImage(const Image *, Coord);
bool inImage(const PixelEdge &e);
bool onImageBorder(const PixelEdge &e);

extern int myWindowID[];
extern bool NO_DISPLAY, SILENT;
extern bool leftMouseButtonIsDown;
extern bool middleMouseButtonIsDown;
extern bool rightMouseButtonIsDown;
extern bool shiftPressed;
extern bool ctrlPressed;
extern PixelLoc currXY;
extern Coord startImagePt;

class Color;
Color asColor(PixelLoc pos, Image *img=0);
vector<Color> asColor(const vector<PixelLoc> &loc);
Color asColorOrGray(PixelLoc pos, Image *img=0);
Color asBoundedColor(PixelLoc pos, Image *img=0);
Color asInterpolatedColor(Coord pt, Image *img=0);
double asGray(PixelLoc pos, Image *img=0);
double asInterpolatedGray(Coord pt, Image *img=0);

class Segment;
class Arc;
class Curvelet;
class Curve;

struct CPoly {
  vector<Coord> bdy;
  string text;
  Color bdyColor, fillColor;
  bool widen, thick;
  CPoly() : widen(false), thick(false) {}
  void draw() const;

  // globals
  static vector<CPoly> allCPoly[4], mouseInfo;
  static bool showCPoly;
  static void drawAll(int windowIndex=0);
};
void addCPoly(Coord, Color col=RED, int windowIndex=0);
void addCPoly(Coord, Coord, Color col=WHITE, bool w=false, int windowIndex=0);
void addCPoly(const vector<Coord> &pt, Color col=WHITE, bool w=false, int windowIndex=0);
void addCPoly(Coord, double, Color col=RED, int windowIndex=0);
void addCPoly(Coord, const vector<int> &, Color col=RED, int windowIndex=0);
void addCPoly(Coord, const vector<double> &, Color col=RED, int windowIndex=0);
void addCPoly(Coord, const char *, Color col=RED, int windowIndex=0);

void addCPoly(PixelEdge, Color col=WHITE);
void addCPoly(const Segment &, Color col=WHITE);
void addCPoly(const Arc &, Color col=WHITE);
void addCPoly(const Curvelet &, Color col=WHITE, bool left=true);
void addCPoly(const Curve &, Color col=WHITE);

void moveImgBBox(PixelLoc delta);
int windowFocus();
void setWindowFocus(int);
void requestWindowRedraw(int);
void zoomInBy(double s, int x, int y, int windowIndex);
void zoomInBy(double s, Coord fixed, int windowIndex);
void resizeWindow(double s, int windowIndex); // in gl.cpp
void fullScreenWindow(int windowIndex); // in gl.cpp
void unFullScreenWindow(int windowIndex); // in gl.cpp
void translateWindowBy(int dx, int dy); // in gl.cpp
void doubleWindow();         // in gl.cpp
void quadrupleWindow();         // in gl.cpp
void resetZoom(int windowIndex);
void reshape(int windowIndex);  // in ui.cpp
void reshape(int w, int h, int windowIndex);  // in ui.cpp
void estimateImageScale(int w, int h, float &scale);
void fitWindowToImage(bool updateScale, int windowIndex);
void drawMouseInfo(int x, int y, bool shiftPressed, bool ctrlPressed, int windowIndex);
void undrawMouseInfo();
void drawWindow(int windowIndex=0);
void drawNumber(Coord pt, double d);
void sendWidth(float w);
void sendPointSize(float s);
void printTitle(const char *s);
void printTitle(bool t);
void printTitle(const char *s, double d);
void printTitle(const char *s, int i);
void printTitle(const char *s, bool t);
void printTitle(double t);
void printTitle(Coord c);
void printTitle(const Curvelet &c);
void printTitle(const Curve &c);
void displayAllConnections(bool sameBdy, int i);
void displayConnectionInfo(bool sameBdy, int i=0, int j=0, bool init=false);
void mysrand(long val);
string getFilenameExtension(const string &filename);
void removeFilenameExtension(string &filename);
string getFinalDirectory(const string &query);
string getFilename(const string &query);
string getFilenameBase(const string &query);
bool fileExists(const string &query);
bool removeFile(const string &query);
void prependDotToFilename(string &s);
void removeSpacesAndTabs(istream &is);
#include<vector>
using std::vector;
void drawNumbers(Coord pt, const vector<int> &d);
void drawNumbers(Coord pt, const vector<double> &d);

struct Matrix3x3
{
  double m[9];
  double &operator()(int i, int j) { return m[3*j+i]; }
  const double &operator()(int i, int j) const { return m[3*j+i]; }
  Matrix3x3() { for(int i=0;i<9;++i) m[i] = 0; }
  Matrix3x3(bool t) { for(int i=0;i<9;++i) m[i] = 0; m[0]=m[4]=m[8]=t; }
  Matrix3x3(double angle) { for(int i=0;i<9;++i) m[i] = 0;
      m[0] = m[4] = cos(angle);  m[1] = m[3] = sin(angle); m[1]*=-1; m[8] = 1; }
  Matrix3x3(Coord trans) { for(int i=0;i<9;++i) m[i] = 0;
      m[2] = trans.x; m[5] = trans.y; m[8] = 1; }
  Matrix3x3(istream &is) { for(int i=0;i<9;++i) is >> m[i]; }
  Matrix3x3 &operator*=(double scale)
    { for(int i=0;i<9;++i) m[i] *= scale; return *this; }
};
ostream & operator<<(ostream & ostr,const Matrix3x3 &);
static const Matrix3x3 zero3x3, identity3x3(true);
Matrix3x3 computeHomography(const vector<Coord> &p, const vector<Coord> &q);

static const double DEFAULT_LINE_WIDTH = 1.5;

namespace Wright {
  extern Image *leftImg, *rightImg;
  extern bool showOverlay, oneTileMode;
  extern string selectedTileID;
};

class Contour {
  string camName;
  vector<int> superPixelID;
public:
  vector<PixelLoc> bdy; // cached set of locations of the contour boundary
  vector< vector<Coord> > approxBdy; // cached simplified version of the boundary
  Coord centroid;
  Color avgCol, highlightCol;

  Contour(const string &cam) : camName(cam), centroid(noCoord) {}

  string getCamName() const { return camName; }
  vector<int> getSuperPixels() const { return superPixelID; }
  void addSuperPixel(const string &tID, int i);
  void removeSuperPixel(const string &tID, int i);
  
  void drawOnOverlayImg(bool onRight);
  void computeAvgCol(bool onRight);
  void computeHighlightCol(bool onRight);
  void computeSimplifiedBdy(bool onRight, bool justRaw=false);
};
extern Contour noContour;
ostream &operator<<(ostream &os, const Contour &m);

string nextTileID();
class Tile {
  string id;  // four digits followed by username
  vector<Contour> contour;  // actual contours for this tile
  vector<string> cornerID;  // ids of corners
  string name;

public:
  int poly3Index;  // index of the associated poly3

  Tile() : id(nextTileID()), poly3Index(-1) {}
  Tile(const string &_id) : id(_id), poly3Index(-1) {}
  Tile(double,double,double) : id("noTile"), poly3Index(-1) {}

  string getID() const { return id; }
  static bool isValidID(const string &s);
  bool hasCamera(const string &camName) const;
  vector<string> getCameras() const;
  Contour &getContour(bool onRight);  // for ui
  Contour &getContour(const string &camName);
  const Contour &getContour(const string &camName) const;
  vector<int> getSuperPixels(bool onRight) const;
  vector<int> getSuperPixels(const string &camName) const;
  void addFiber(const string &camName, int spID);
  void addFibers(const string &camName, const vector<int> &spID);
  void removeFiber(const string &camName, int spID);
  void removeFibers(const string &camName, const vector<int> &spID);
  void removeAllFiber(const string &camName);
  bool hasCorner(const string &cornerID) const;
  vector<string> getCorners() const { return cornerID; }
  void addCorner(const string &cornerID);
  void removeCorner(const string &cornerID);
  string getName() const { return name; }
  void setName(const string &n) { name = n; }
  friend vector<string> outputAllTiles();
  Color getHighlightColor() const;

  static int getIntegerPart(const string &s);
};
extern vector< Tile > allTiles;
extern Tile noTile;
ostream &operator<<(ostream &os, const Tile &t);
bool noTileSelected();
bool noCornerSelected();
Tile &getTile(const string &tileID);
Tile inputTile(const string &, bool standalone=false, bool newPoly3=true);
string outputTile(const Tile &);

void initAllSides();
void init_gl_window();
bool shiftIsPressed();
bool ctrlIsPressed();
bool isDigit(char c);
vector<PixelLoc> getContour(const string &tileID, const string &imgName);
Image getContourColors(const string &tileID, const string &imgName);
#endif // _ERIOL_
