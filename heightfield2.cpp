// shapes/heightfield2.cpp*
#include "stdafx.h"
#include "shapes/heightfield2.h"
#include "shapes/trianglemesh.h"
#include "paramset.h"

// Heightfield2 Method Definitions
Heightfield2::Heightfield2(const Transform *o2w, const Transform *w2o,
                           bool ro, int x, int y, const float *zs)
    : Shape(o2w, w2o, ro) {
    nx = x;
    ny = y;
    z = new float[nx*ny];
    memcpy(z, zs, nx*ny*sizeof(float));

    grid = new float*[nx];
    for(int i = 0; i < nx; i++)
    {
        grid[i] = new float[ny];
        for(int j = 0; j < ny; j++)
        {
            grid[i][j] = z[i + j*ny];
        }
    }

    int ntris = 2*(nx-1)*(ny-1);
    verts = new int[3*ntris];
    P = new Point[nx*ny];
    uvs = new float[2*nx*ny];
    nverts = nx*ny;
    int i, j;
    // Compute heightfield2 vertex positions
    int pos = 0;
    for (j = 0; j < ny; ++j) {
        for (i = 0; i < nx; ++i) {
            P[pos].x = uvs[2*pos]   = (float)x / (float)(nx-1);
            P[pos].y = uvs[2*pos+1] = (float)y / (float)(ny-1);
            P[pos].z = z[pos];
            ++pos;
        }
    }

    // Fill in heightfield2 vertex offset array
    int *vp = verts;
    for (j = 0; j < ny-1; ++j) {
        for (i = 0; i < nx-1; ++i) {
#define VERT(x,y) ((x)+(y)*nx)
            *vp++ = VERT(i, j);
            *vp++ = VERT(i+1, j);
            *vp++ = VERT(i+1, j+1);

            *vp++ = VERT(i, j);
            *vp++ = VERT(i+1, j+1);
            *vp++ = VERT(i, j+1);
        }
#undef VERT
    }
//    ParamSet paramSet;
//    paramSet.AddInt("indices", verts, 3*ntris);
//    paramSet.AddFloat("uv", uvs, 2 * nverts);
//    paramSet.AddPoint("P", P, nverts);
//    shape.push_back(CreateTriangleMeshShape(ObjectToWorld, WorldToObject, ReverseOrientation, paramSet));
//    delete[] P;
//    delete[] uvs;
//    delete[] verts;
}


Heightfield2::~Heightfield2() {
    delete[] z;
}


BBox Heightfield2::ObjectBound() const {
    float minz = z[0], maxz = z[0];
    for (int i = 1; i < nx*ny; ++i) {
        if (z[i] < minz) minz = z[i];
        if (z[i] > maxz) maxz = z[i];
    }
    return BBox(Point(0,0,minz), Point(1,1,maxz));
}


bool Heightfield2::CanIntersect() const {
    return true;
}


bool Heightfield2::TriangleInterscetion(const Ray &objRay, float *tHit, DifferentialGeometry *dg, Point p1a, Point p2a, Point p3a, const float *uvs, int index) const
{
    // Get triangle vertices in _p1_, _p2_, and _p3_
    const Point &p1 = p1a;
    const Point &p2 = p2a;
    const Point &p3 = p3a;

    Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(objRay.d, e2);
    float divisor = Dot(s1, e1);
    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;
    // Compute first barycentric coordinate
    Vector d = objRay.o - p1;
    float b1 = Dot(d, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;
    // Compute second barycentric coordinate
    Vector s2 = Cross(d, e1);
    float b2 = Dot(objRay.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;
    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < objRay.mint || t > objRay.maxt)
        return false;

    // Compute triangle partial derivatives
    Vector dpdu, dpdv;
    float _uvs[3][2];
    _uvs[0][0] = uvs[2 * index];
    _uvs[0][1] = uvs[2 * index + 1];
    _uvs[1][0] = uvs[2 * index + 2];
    _uvs[1][1] = uvs[2 * index + 3];
    _uvs[2][0] = uvs[2 * index + 4];
    _uvs[2][1] = uvs[2 * index + 5];

    // Compute deltas for triangle partial derivatives
    float du1 = _uvs[0][0] - _uvs[2][0];
    float du2 = _uvs[1][0] - _uvs[2][0];
    float dv1 = _uvs[0][1] - _uvs[2][1];
    float dv2 = _uvs[1][1] - _uvs[2][1];
    Vector dp1 = p1 - p3, dp2 = p2 - p3;
    float determinant = du1 * dv2 - dv1 * du2;
    if (determinant == 0.f) {
        // Handle zero determinant for triangle partial derivative matrix
        CoordinateSystem(Normalize(Cross(e2, e1)), &dpdu, &dpdv);
    }
    else {
        float invdet = 1.f / determinant;
        dpdu = ( dv2 * dp1 - dv1 * dp2) * invdet;
        dpdv = (-du2 * dp1 + du1 * dp2) * invdet;
    }

    // Interpolate $(u,v)$ triangle parametric coordinates
    float b0 = 1 - b1 - b2;
    float tu = b0 * _uvs[0][0] + b1 * _uvs[1][0] + b2 * _uvs[2][0];
    float tv = b0 * _uvs[0][1] + b1 * _uvs[1][1] + b2 * _uvs[2][1];

    const Transform &o2w = *ObjectToWorld;
    *dg = DifferentialGeometry(o2w(objRay(t)), o2w(dpdu), o2w(dpdv),
                               Normal(0,0,0), Normal(0,0,0),
                               tu, tv, this);
    *tHit = t;
    return true;
}


bool Heightfield2::TriangleInterscetionP(const Ray &objRay, Point p1a, Point p2a, Point p3a) const
{
    const Point &p1 = p1a;
    const Point &p2 = p2a;
    const Point &p3 = p3a;

    Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(objRay.d, e2);
    float divisor = Dot(s1, e1);
    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;
    // Compute first barycentric coordinate
    Vector d = objRay.o - p1;
    float b1 = Dot(d, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;
    // Compute second barycentric coordinate
    Vector s2 = Cross(d, e1);
    float b2 = Dot(objRay.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;
    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < objRay.mint || t > objRay.maxt)
        return false;
    
    return true;
}


bool Heightfield2::Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                             DifferentialGeometry *dg) const {
    Ray r;
    (*WorldToObject)(ray, &r);

    // Check ray against overall grid bounds
    float rayT;
    BBox bounds = this->ObjectBound();

    if (bounds.Inside(r(r.mint)))
        rayT = r.mint;
    else if (!bounds.IntersectP(r, &rayT))
        return false;

    Point gridIntersect = r(rayT);

    // Set up 2D DDA for ray
    float NextCrossingT[2], DeltaT[2];
    int Step[2], Out[2], Pos[2];
    float Width[2] = {1., 1.};
    int NVoxels[2] = {nx-1, ny-1};
    for (int axis = 0; axis < 2; ++axis) {
        // Compute current voxel for axis
        gridIntersect[axis] = gridIntersect[axis] * NVoxels[axis];
        Pos[axis] = Clamp(Float2Int(gridIntersect[axis]), 0, NVoxels[axis]-1);
        if (r.d[axis] >= 0) {
            // Handle ray with positive direction for voxel stepping
            NextCrossingT[axis] = rayT + (Pos[axis]+1 - gridIntersect[axis]) / (r.d[axis] * NVoxels[axis]);
            DeltaT[axis] = Width[axis] / (r.d[axis] * NVoxels[axis]);
            Step[axis] = 1;
            Out[axis] = NVoxels[axis];
        }
        else {
            // Handle ray with negative direction for voxel stepping
            NextCrossingT[axis] = rayT + (Pos[axis] - gridIntersect[axis]) / (r.d[axis] * NVoxels[axis]);
            DeltaT[axis] = -Width[axis] / (r.d[axis] * NVoxels[axis]);
            Step[axis] = -1;
            Out[axis] = -1;
        }
    }

    // Walk ray through voxel grid
    for(;;)
    {
        const Point &p1 = Point((float) Pos[0] / NVoxels[0], (float) Pos[1] / NVoxels[1], grid[Pos[0]][Pos[1]]);
        const Point &p2 = Point((float) (Pos[0] + 1) / NVoxels[0], (float) Pos[1] / NVoxels[1],	grid[Pos[0] + 1][Pos[1]]);
        const Point &p3 = Point((float) (Pos[0] + 1) / NVoxels[0], (float) (Pos[1] + 1) / NVoxels[1],	grid[Pos[0] + 1][Pos[1] + 1]);
        const Point &p4 = Point((float) Pos[0] / NVoxels[0], (float) (Pos[1] + 1) / NVoxels[1],	grid[Pos[0]][Pos[1] + 1]);

        // Get the index!!!
        int index = verts[(Pos[1] * nx + Pos[0]) * 3];
        // printf("%f %d\n", P[index].x, Pos[0]);
        // printf("%d %d\n", Pos[0], Pos[1]);
        // printf("index = %d\n", index);

        float thit = 0.0;

        if (TriangleInterscetion(r, &thit, dg, p1, p2, p3, uvs, index)) {
            *tHit = thit;
            *rayEpsilon = 1e-3F * thit;
            return true;
        }

        index = verts[(Pos[1] * nx + Pos[0]) * 2 + 1];

        if (TriangleInterscetion(r, &thit, dg, p1, p4, p3, uvs, index)) {
            *tHit = thit;
            *rayEpsilon = 1e-3F * thit;
            return true;
        }

//        if (TriangleInterscetion(r, &thit1, dg, p1, p2, p3, uvs, index) || TriangleInterscetion(r, &thit2, dg, p1, p4, p3, uvs, index))
//        {
//            *tHit = min(thit1, thit2);
//            return true;
//        }
        // Advance to next voxel
        // Find _stepAxis_ for stepping to next voxel
        int stepAxis = NextCrossingT[0] < NextCrossingT[1] ? 0 : 1;
        if (r.maxt < NextCrossingT[stepAxis])
            break;
        Pos[stepAxis] += Step[stepAxis];
        if (Pos[stepAxis] == Out[stepAxis])
            break;
        NextCrossingT[stepAxis] += DeltaT[stepAxis];
    }

    return false;
}


bool Heightfield2::IntersectP(const Ray &ray) const {
    Ray r;
    (*WorldToObject)(ray, &r);

    // Check ray against overall grid bounds
    float rayT;
    BBox bounds = this->ObjectBound();

    if (bounds.Inside(r(r.mint)))
        rayT = r.mint;
    else if (!bounds.IntersectP(r, &rayT))
        return false;

    Point gridIntersect = r(rayT);

    // Set up 2D DDA for ray
    float NextCrossingT[2], DeltaT[2];
    int Step[2], Out[2], Pos[2];
    float Width[2] = {1., 1.};
    int NVoxels[2] = {nx-1, ny-1};
    for (int axis = 0; axis < 2; ++axis) {
        // Compute current voxel for axis
        gridIntersect[axis] = gridIntersect[axis] * NVoxels[axis];
        Pos[axis] = Clamp(Float2Int(gridIntersect[axis]), 0, NVoxels[axis]-1);
        if (r.d[axis] >= 0) {
            // Handle ray with positive direction for voxel stepping
            NextCrossingT[axis] = rayT + (Pos[axis]+1 - gridIntersect[axis]) / (r.d[axis] * NVoxels[axis]);
            DeltaT[axis] = Width[axis] / (r.d[axis] * NVoxels[axis]);
            Step[axis] = 1;
            Out[axis] = NVoxels[axis];
        }
        else {
            // Handle ray with negative direction for voxel stepping
            NextCrossingT[axis] = rayT + (Pos[axis] - gridIntersect[axis]) / (r.d[axis] * NVoxels[axis]);
            DeltaT[axis] = -Width[axis] / (r.d[axis] * NVoxels[axis]);
            Step[axis] = -1;
            Out[axis] = -1;
        }
    }

    // Walk ray through voxel grid
    for(;;)
    {
        const Point &p1 = Point((float) Pos[0] / NVoxels[0], (float) Pos[1] / NVoxels[1], grid[Pos[0]][Pos[1]]);
        const Point &p2 = Point((float) (Pos[0] + 1) / NVoxels[0], (float) Pos[1] / NVoxels[1],	grid[Pos[0] + 1][Pos[1]]);
        const Point &p3 = Point((float) (Pos[0] + 1) / NVoxels[0], (float) (Pos[1] + 1) / NVoxels[1],	grid[Pos[0] + 1][Pos[1] + 1]);
        const Point &p4 = Point((float) Pos[0] / NVoxels[0], (float) (Pos[1] + 1) / NVoxels[1],	grid[Pos[0]][Pos[1] + 1]);

        if (TriangleInterscetionP(r, p1, p2, p3) || TriangleInterscetionP(r, p1, p4, p3))
        {
            return true;
        }
        // Advance to next voxel
        // Find _stepAxis_ for stepping to next voxel
        int stepAxis = NextCrossingT[0] < NextCrossingT[1] ? 0 : 1;
        if (r.maxt < NextCrossingT[stepAxis])
            break;
        Pos[stepAxis] += Step[stepAxis];
        if (Pos[stepAxis] == Out[stepAxis])
            break;
        NextCrossingT[stepAxis] += DeltaT[stepAxis];
    }

    return false;
}


void Heightfield2::Refine(vector<Reference<Shape> > &refined) const {
    int ntris = 2*(nx-1)*(ny-1);
    refined.reserve(ntris);
    int *verts = new int[3*ntris];
    Point *P = new Point[nx*ny];
    float *uvs = new float[2*nx*ny];
    int nverts = nx*ny;
    int x, y;
    // Compute heightfield2 vertex positions
    int pos = 0;
    for (y = 0; y < ny; ++y) {
        for (x = 0; x < nx; ++x) {
            P[pos].x = uvs[2*pos]   = (float)x / (float)(nx-1);
            P[pos].y = uvs[2*pos+1] = (float)y / (float)(ny-1);
            P[pos].z = z[pos];
            ++pos;
        }
    }

    // Fill in heightfield2 vertex offset array
    int *vp = verts;
    for (y = 0; y < ny-1; ++y) {
        for (x = 0; x < nx-1; ++x) {
#define VERT(x,y) ((x)+(y)*nx)
            *vp++ = VERT(x, y);
            *vp++ = VERT(x+1, y);
            *vp++ = VERT(x+1, y+1);

            *vp++ = VERT(x, y);
            *vp++ = VERT(x+1, y+1);
            *vp++ = VERT(x, y+1);
        }
#undef VERT
    }
    ParamSet paramSet;
    paramSet.AddInt("indices", verts, 3*ntris);
    paramSet.AddFloat("uv", uvs, 2 * nverts);
    paramSet.AddPoint("P", P, nverts);
    refined.push_back(CreateTriangleMeshShape(ObjectToWorld, WorldToObject, ReverseOrientation, paramSet));
    delete[] P;
    delete[] uvs;
    delete[] verts;
}


Heightfield2 *CreateHeightfield2Shape(const Transform *o2w, const Transform *w2o,
                                      bool reverseOrientation, const ParamSet &params) {
    int nu = params.FindOneInt("nu", -1);
    int nv = params.FindOneInt("nv", -1);
    int nitems;
    const float *Pz = params.FindFloat("Pz", &nitems);
    Assert(nitems == nu*nv);
    Assert(nu != -1 && nv != -1 && Pz != NULL);
    return new Heightfield2(o2w, w2o, reverseOrientation, nu, nv, Pz);
}