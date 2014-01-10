// Color.h by Olaf Hall-Holt, 2007-2011
#include <iostream>
#include <string>
#include <cassert>
#include <vector>
#include <math.h>


#ifndef __COLOR__
#define __COLOR__

class Color;
double hueDiff(const Color &a, const Color &b);
double hueDiff2(Color a, Color b);
double robustHueDiff(const Color &a, const Color &b);
struct Color {
  double r,g,b;
  Color() { r = 0.; g = 0.; b = 0.; }
  Color(double _r, double _g, double _b) : r(_r), g(_g), b(_b) {}
  Color(std::istream &is) { is >> r >> g >> b; }
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
  friend std::ostream& operator<<(std::ostream &os, const Color &p) {
    os << p.r << " " << p.g << " " << p.b;
    return os; }
  friend std::istream & operator>>(std::istream& istr, Color &c) {
    c = Color(istr);  return istr; }
};
void sendColor(double, double, double);
void sendColor(Color);


#endif //__COLOR__