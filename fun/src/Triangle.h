#ifndef Triangle_h_
#define Triangle_h_

#include <cmath>
#include "Object.h"
#include <iostream>
using namespace std;


// 取三浮點數最大值
float Max3( float a1,  float a2,  float a3){
	float out = a1;
	if(a2 > out){ out = a2; }
	if(a3 > out){ out = a3; }
	return out;
}

// 取三浮點數最小值
float Min3(float a1, float a2, float a3){
	float out = a1;
	if(a2 < out){ out = a2; }
	if(a3 < out){ out = a3; }
	return out;
}

//! For the purposes of demonstrating the BVH, a simple triangle
struct Triangle : public Object {

// 物件細節
public:
	Vector3 v0,v1,v2;	// vertex of triangle
	Vector3 center;		// center of triangle
	Vector3 N;			// Normal vector

// 物件成員方程式(實作方法)
public:
	// 無定義初始化物件
	Triangle(){}
	
	// 有定義初始化物件
	Triangle(const Vector3& V0, const Vector3& V1, const Vector3& V2){
		v0 = V0;
		v1 = V1;
		v2 = V2;
		center = (v0+v1+v2)/3.0;

		Vector3 v0v1 = v1 - v0;
		Vector3 v0v2 = v2 - v0;

		N = normalize(v0v1 ^ v0v2);
	}


	// 取Bounding box定義
	BBox getBBox() const {

		float x_min = Min3(v0.x, v1.x, v2.x);
		float y_min = Min3(v0.y, v1.y, v2.y);
		float z_min = Min3(v0.z, v1.z, v2.z);

		float x_max = Max3(v0.x, v1.x, v2.x);
		float y_max = Max3(v0.y, v1.y, v2.y);
		float z_max = Max3(v0.z, v1.z, v2.z);
		
		return BBox(Vector3(x_min,y_min,z_min), Vector3(x_max,y_max,z_max));
	}

	// 波束與三角形交點
	bool getIntersection(const Ray& ray, IntersectionInfo* I) const {
		// using MOLLER TRUMBORE algorithm
		// Ref: https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection
		Vector3 v0v1 = v1 - v0;
		Vector3 v0v2 = v2 - v0;
		Vector3 pvec = ray.d ^ v0v2;	// cross product
		float det = v0v1 * pvec;	// inner product

		// if the determinant is negative the triangle is backfacing
		// if the determinant is close to 0, the ray misses the triangle
		if (det < 1E-10) return false;





		float invDet = 1 / det;

		Vector3 tvec = ray.o - v0;
		float u = (tvec * pvec) * invDet;
		if (u < 0 || u > 1) return false;

		Vector3 qvec = tvec ^ v0v1;
		float v = (ray.d * qvec) * invDet;
		if (v < 0 || u + v > 1) return false;

		I->object = this;

		I->t = (v0v2 * qvec) * invDet;

		if(I->t<=0){
			return false;
		}

		I->hit = ray.o + I->t*ray.d;


		return true; 
	}
    
	// 取三角形法線向量
    Vector3 getNormal(IntersectionInfo& I) const {
    	return N;
    }

    // 取三角形中心點座標
	Vector3 getCentroid() const {
		return center;
	}



};

#endif


