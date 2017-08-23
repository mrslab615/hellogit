#ifndef IntersectionInfo_h_
#define IntersectionInfo_h_

class Object;

struct IntersectionInfo {
  double t; // Intersection distance along the ray
  const Object* object; // Object that was hit
  Vector3 hit; // Location of the intersection
};

#endif
