#ifndef Vector3_h
#define Vector3_h

#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include "Log.h"


struct Vector3 {
  struct { double x,y,z; };

  Vector3() { }
  Vector3(double x, double y, double z) : x(x),y(y),z(z) { }

//  float x(){ return x; }
//  float y(){ return y; }
//  float z(){ return z; }

  Vector3 operator+(const Vector3& b) const { return Vector3(x+b.x, y+b.y, z+b.z); }
  Vector3 operator-(const Vector3& b) const { return Vector3(x-b.x, y-b.y, z-b.z); }
  Vector3 operator*(double b) const { return Vector3(x*b, y*b, z*b); }
  Vector3 operator/(double b) const { return Vector3(x/b, y/b, z/b); }

  // Component-wise multiply and divide
  Vector3 cmul(const Vector3& b) const { return Vector3(x*b.x, y*b.y, z*b.z); }
  Vector3 cdiv(const Vector3& b) const { return Vector3(x/b.x, y/b.y, z/b.z); }

  // dot (inner) product
  double operator*(const Vector3& b) const {
    return x*b.x + y*b.y + z*b.z;
  }

  // Cross Product
  Vector3 operator^(const Vector3& b) const {
    return Vector3(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x);
  }

  // sample-wise div
  Vector3 operator/(const Vector3& b) const { return Vector3(x/b.x, y/b.y, z/b.z); }

  // Handy component indexing
  double& operator[](const unsigned int i) { return (&x)[i]; }
  const double& operator[](const unsigned int i) const { return (&x)[i]; }
};

inline Vector3 operator*(double a, const Vector3&b)  { return Vector3(a*b.x,a*b.y,a*b.z); }

// Component-wise min
inline Vector3 min(const Vector3& a, const Vector3& b) {
  return Vector3(std::min(a.x,b.x), std::min(a.y,b.y), std::min(a.z,b.z));
}

// Component-wise max
inline Vector3 max(const Vector3& a, const Vector3& b) {
  return Vector3(std::max(a.x,b.x), std::max(a.y,b.y), std::max(a.z,b.z));
}

// Length of a vector
inline double length(const Vector3& a) {
  return sqrtf(a*a);
}

// Make a vector unit length
inline Vector3 normalize(const Vector3& in) {
  double norm = sqrt(in.x*in.x + in.y*in.y + in.z*in.z);
  return Vector3(in.x/norm, in.y/norm, in.z/norm);
}

#endif
