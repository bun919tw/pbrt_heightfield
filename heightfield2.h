#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_SHAPES_HEIGHTFIELD2_H
#define PBRT_SHAPES_HEIGHTFIELD2_H

// shapes/heightfield2.h*
#include "shape.h"
#include "shapes/trianglemesh.h"

// Heightfield2 Declarations
class Heightfield2 : public Shape {
public:
    // Heightfield2 Public Methods
    Heightfield2(const Transform *o2, const Transform *w2o, bool ro, int nu, int nv, const float *zs);
    ~Heightfield2();
    bool CanIntersect() const;
    bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                   DifferentialGeometry *dg) const;
    bool IntersectP(const Ray &ray) const;
    bool TriangleInterscetion(const Ray &objRay, float *tHit, DifferentialGeometry *dg, Point p1a, Point p2a, Point p3a, const float *uvs, int index) const;
    bool TriangleInterscetionP(const Ray &objRay, Point p1a, Point p2a, Point p3a) const;
    void Refine(vector<Reference<Shape> > &refined) const;
    BBox ObjectBound() const;
private:
    // Heightfield2 Private Data
    float *z;
    float **grid;
    int *verts;
    int nverts;
    Point *P;
    float *uvs;
    int nx, ny;
};


Heightfield2 *CreateHeightfield2Shape(const Transform *o2w, const Transform *w2o,
                                      bool reverseOrientation, const ParamSet &params);

#endif // PBRT_SHAPES_HEIGHTFIELD2_H
