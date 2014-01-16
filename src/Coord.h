// Coord.h by Olaf Hall-Holt, 2007-2011

#include "PixelLoc.h"
#include <iostream>

#ifndef __COORD__
#define __COORD__

struct Coord {
  double x,y;
  Coord(): x(0), y(0) {}
  Coord(double _x, double _y): x(_x), y(_y) {}
  Coord(std::istream &is) { is >> x; is >> std::ws; assert(is.get() == ','); is >> y; }

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
  friend std::ostream & operator<<(std::ostream& ostr, const Coord& c) {
    ostr << c.x << "," << c.y;
    return ostr; }
  friend std::istream & operator>>(std::istream& is, Coord &c) { is >> std::ws;
    // if ( EOF == is.peek() ) is.setstate(ios::badbit); else
    c = Coord(is);
    return is; }
};
static const Coord noCoord(1e+10, 1e+10);
Coord asShiftedCoord(const PixelLoc &a);


#endif //__COORD__