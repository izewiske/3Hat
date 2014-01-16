// BBox.h by Olaf Hall-Holt, 2007-2011

#include "PixelLoc.h"
#include "Coord.h"
#include <vector>
#include <iostream>


#ifndef __BBOX__
#define __BBOX__

struct BBox {
  PixelLoc ul;
  int w,h;
  BBox(): ul(0,0), w(-1), h(-1) {}
  BBox(PixelLoc pos): ul(pos), w(1), h(1) {}
  BBox(PixelLoc pos, int _w, int _h): ul(pos), w(_w), h(_h) {}
  BBox(PixelLoc pos, PixelLoc pos2): ul(pos), w(pos2.x-pos.x+1), h(pos2.y-pos.y+1) {}
  BBox(const std::vector<Coord> &pt): ul(0,0), w(-1), h(-1) {
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
  std::vector<PixelLoc> pixelsInside() const { std::vector<PixelLoc> res;
    PixelLoc pos;
    for ( pos.y = minY(); pos.y <= maxY(); ++pos.y )
    for ( pos.x = minX(); pos.x <= maxX(); ++pos.x )
      res.push_back(pos);
    return res; }
  void draw() const;
  friend std::ostream& operator<<(std::ostream &os, const BBox &b) {
    os << b.w << "x" << b.h << "+" << b.ul.x << "+" << b.ul.y;
    return os; }
};


#endif //__BBOX__