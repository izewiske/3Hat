//PixelLoc.h by Olaf Hall-Holt, 2007-2011

#include <iostream>
#include <cassert>
#include <math.h>

#ifndef __PIXEL_LOC__
#define __PIXEL_LOC__


struct PixelLoc {
  int x,y;
  PixelLoc() : x(-1), y(-1) {}
  PixelLoc(int _x, int _y) : x(_x), y(_y) {}
  PixelLoc(std::istream &is) { is >> x; assert(',' == is.get()); is >> y; }
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
  friend std::ostream& operator<<(std::ostream &os, const PixelLoc &p) {
    os << p.x << "," << p.y;
    return os; }
  friend std::istream & operator>>(std::istream& istr, PixelLoc &p) {
    p = PixelLoc(istr);  return istr; }
};

#endif //__PIXEL_LOC__