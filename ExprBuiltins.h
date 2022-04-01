/*
* Copyright Disney Enterprises, Inc.  All rights reserved.
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License
* and the following modification to it: Section 6 Trademarks.
* deleted and replaced with:
*
* 6. Trademarks. This License does not grant permission to use the
* trade names, trademarks, service marks, or product names of the
* Licensor and its affiliates, except as required for reproducing
* the content of the NOTICE file.
*
* You may obtain a copy of the License at
* http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef ExprBuiltins_h
#define ExprBuiltins_h

#include "Vec.h"

namespace SeExpr2 {

// trig
inline double deg(double angle) { return angle * (180 / M_PI); }
inline double rad(double angle) { return angle * (M_PI / 180); }
inline double cosd(double x) { return cos(rad(x)); }
inline double sind(double x) { return sin(rad(x)); }
inline double tand(double x) { return tan(rad(x)); }
inline double acosd(double x) { return deg(acos(x)); }
inline double asind(double x) { return deg(asin(x)); }
inline double atand(double x) { return deg(atan(x)); }
inline double atan2d(double y, double x) { return deg(atan2(y, x)); }

// clamping
inline double clamp(double x, double lo, double hi) { return x < lo ? lo : x > hi ? hi : x; }
inline double round(double x) { return x < 0 ? ceil(x - 0.5) : floor(x + 0.5); }
inline double max(double x, double y) { return x > y ? x : y; }
inline double min(double x, double y) { return x < y ? x : y; }

// blending / remapping
inline double invert(double x) { return 1 - x; }
double compress(double x, double lo, double hi);
double expand(double x, double lo, double hi);
double fit(double x, double a1, double b1, double a2, double b2);
double gamma(double x, double g);
double bias(double x, double b);
double contrast(double x, double c);
double boxstep(double x, double a);
double linearstep(double x, double a, double b);
double smoothstep(double x, double a, double b);
double gaussstep(double x, double a, double b);
double remap(double x, double s, double r, double f, double interp);
double mix(double x, double y, double alpha);
Vec3d hsi(Vec3d const& rgb, double h, double s, double i, double m=1);
Vec3d midhsi(Vec3d const& rgb, double h, double s, double i, double m, double falloff=1, int interp=0);
Vec3d rgbtohsl(const Vec3d& rgb);
Vec3d hsltorgb(const Vec3d& hsl);
Vec3d saturate(const Vec3d& rgb, double amt);

// noise
double hash(int n, double const* args);
double noise(Vec3d const& p);
double noise(double x);
double noise(double x, double y);
double noise(double x, double y, double z);
double noise(double x, double y, double z, double w);
double snoise(const Vec3d& p);
double snoise(const Vec3d& p, double w);
Vec3d cnoise(const Vec3d& p);
Vec3d cnoise(const Vec3d& p, double w);
Vec3d vnoise(const Vec3d& p);
Vec3d vnoise(const Vec3d& p, double w);
double turbulence(const Vec3d& p, int octaves = 6, double lacunarity = 2, double gain = 0.5);
Vec3d vturbulence(const Vec3d& p, int octaves = 6, double lacunarity = 2, double gain = 0.5);
Vec3d cturbulence(const Vec3d& p, int octaves = 6, double lacunarity = 2, double gain = 0.5);
double fbm(const Vec3d& p, int octaves = 6, double lacunarity = 2, double gain = 0.5);
Vec3d vfbm(const Vec3d& p, int octaves = 6, double lacunarity = 2, double gain = 0.5);
Vec3d cfbm(const Vec3d& p, int octaves = 6, double lacunarity = 2, double gain = 0.5);
double fbm4(const Vec3d& p, double t, int octaves = 6, double lacunarity = 2, double gain = 0.5);
Vec3d vfbm4(const Vec3d& p, double t, int octaves = 6, double lacunarity = 2, double gain = 0.5);
Vec3d cfbm4(const Vec3d& p, double t, int octaves = 6, double lacunarity = 2, double gain = 0.5);
double cellnoise(const Vec3d& p);
Vec3d ccellnoise(const Vec3d& p);
double pnoise(const Vec3d& p, const Vec3d& period);

// vectors
double dist(double ax, double ay, double az, double bx, double by, double bz);
double length(const Vec3d& v);
double hypot(double x, double y);
double dot(const Vec3d& a, const Vec3d& b);
Vec3d norm(const Vec3d& a);
Vec3d cross(const Vec3d& a, const Vec3d& b);
double angle(const Vec3d& a, const Vec3d& b);
Vec3d ortho(const Vec3d& a, const Vec3d& b);
Vec3d up(const Vec3d& vec, const Vec3d& upvec);

// variations
double cycle(double index, double loRange, double hiRange);
double pick(int n, double* params);
double choose(int n, double* params);
double wchoose(int n, double* params);
double spline(int n, double* params);

// voronoi
struct VoronoiPointData {
    Vec3d points[27];
    Vec3d cell;
    double jitter;
    VoronoiPointData() : jitter(-1) {}
};
double voronoi(VoronoiPointData& cache, Vec3d p, int type=1,float jitter=0.5, float fbmScale=0, int fbmOctaves=4,float fbmLacunarity=2, float fbmGain=.5);
Vec3d cvoronoi(VoronoiPointData& cache, Vec3d p, int type=1,float jitter=0.5, float fbmScale=0, int fbmOctaves=4,float fbmLacunarity=2, float fbmGain=.5);
Vec3d pvoronoi(VoronoiPointData& cache, Vec3d p, int type=1,float jitter=0.5, float fbmScale=0, int fbmOctaves=4,float fbmLacunarity=2, float fbmGain=.5);

}

#endif
