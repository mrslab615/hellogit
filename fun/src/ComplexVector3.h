//
//  ComplexVector3.h
//  test_vector3complex
//
//  Created by Steve Chiang on 6/2/17.
//  Copyright (c) 2017 Steve Chiang. All rights reserved.
//

#ifndef ComplexVector3_h
#define ComplexVector3_h

#include <iostream>
#include <iomanip>
#include <complex>
#include <string>
#include <cmath>
#include <ostream>
#include "Vector3_new.h"


using namespace std;

typedef complex<double> CXF;

class CXV3 {
public:
	CXF x,y,z;
public:
	CXV3() { }
	CXV3(CXF x, CXF y, CXF z) : x(x),y(y),z(z) { }
	CXV3 operator+(const CXV3& b) const { return CXV3(x+b.x, y+b.y, z+b.z); }
	CXV3 operator-(const CXV3& b) const { return CXV3(x-b.x, y-b.y, z-b.z); }
	CXV3 operator*(CXF b)         const { return CXV3(x*b, y*b, z*b); }
	CXV3 operator/(CXF b)         const { return CXV3(x/b, y/b, z/b); }
	CXF  operator*(const CXV3& b) const { return x*b.x + y*b.y + z*b.z; }	// dot (inner) product
	CXV3 operator^(const CXV3& b) const { return CXV3(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.z); } // Cross Product
	CXV3 operator/(const CXV3& b) const { return CXV3(x/b.x, y/b.y, z/b.z); }	// sample-wise div
	//new function
	CXV3 operator+(const Vector3 &b) const { return CXV3(x+b.x, y+b.y, z+b.z); }
	CXF  operator*(const Vector3 &b) const { return (x*b.x + y*b.y+z*b.z);}	// dot (inner) product
	CXV3 operator*(const double &b) const { return CXV3(x*b, y*b, z*b); }




	friend ostream& operator<<(ostream& os, const CXV3& a){
		os<<"[ "<<a.x<<", "<<a.y<<", "<<a.z<<" ]";
		return os;
	}
};

inline CXV3 operator*(CXF a, const CXV3&b)  { return CXV3(a*b.x,a*b.y,a*b.z); }
//new
//inline CXF operator*(CXF a, const float&b)  { return CXF(a*b); }







#endif
