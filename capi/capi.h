#pragma once

#ifdef _MSC_VER
#if defined(EXPORT_DLL) || defined(IMPORT_DLL)
#ifdef EXPORT_DLL
#define DLL __declspec(dllexport)
#else
#define DLL __declspec(dllimport)
#endif
#else  // EXPORT_DLL or IMPORT_DLL
#define DLL
#endif // EXPORT_DLL or IMPORT_DLL
#else  // _MSC_VER
#define DLL
#endif // _MSC_VER

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Vec3d { double x,y,z; } Vec3d;

// clamping
DLL double clamp(double x, double lo, double hi);
//DLL double round(double x);

// blending / remapping
DLL double invert(double x);
DLL double compress(double x, double lo, double hi);
DLL double expand(double x, double lo, double hi);
DLL double fit(double x, double a1, double b1, double a2, double b2);
DLL double gamma(double x, double g);
DLL double bias(double x, double b);
DLL double contrast(double x, double c);
DLL double boxstep(double x, double a);
DLL double linearstep(double x, double a, double b);
DLL double smoothstep(double x, double a, double b);
DLL double gaussstep(double x, double a, double b);
DLL double remap(double x, double s, double r, double f, double interp);
DLL double mix(double x, double y, double alpha);
DLL Vec3d  hsi(Vec3d rgb, double h, double s, double i, double m);
DLL Vec3d  midhsi(Vec3d rgb, double h, double s, double i, double m, double falloff, int interp);
DLL Vec3d  rgbtohsl(Vec3d rgb);
DLL Vec3d  hsltorgb(Vec3d hsl);
DLL Vec3d  saturate(Vec3d rgb, double amt);

// noise
DLL double hash(int n, double const* args);
DLL double noisep(Vec3d p);
DLL double noise1(double x);
DLL double noise2(double x, double y);
DLL double noise3(double x, double y, double z);
DLL double noise4(double x, double y, double z, double w);
DLL double snoise(Vec3d p);
DLL double snoise4(Vec3d p, double w);
DLL Vec3d  cnoise(Vec3d p);
DLL Vec3d  vnoise(Vec3d p);
DLL Vec3d  cnoise4(Vec3d p, double w);
DLL Vec3d  vnoise4(Vec3d p, double w);
DLL double turbulence(Vec3d p, int octaves, double lacunarity, double gain);
DLL Vec3d  vturbulence(Vec3d p, int octaves, double lacunarity, double gain);
DLL Vec3d  cturbulence(Vec3d p, int octaves, double lacunarity, double gain);
DLL double fbm(Vec3d p, int octaves, double lacunarity, double gain);
DLL Vec3d  vfbm(Vec3d p, int octaves, double lacunarity, double gain);
DLL Vec3d  cfbm(Vec3d p, int octaves, double lacunarity, double gain);
DLL double fbm4(Vec3d p, double t, int octaves, double lacunarity, double gain);
DLL Vec3d  vfbm4(Vec3d p, double t, int octaves, double lacunarity, double gain);
DLL Vec3d  cfbm4(Vec3d p, double t, int octaves, double lacunarity, double gain);
DLL double cellnoise(Vec3d p);
DLL Vec3d  ccellnoise(Vec3d p);
DLL double pnoise(Vec3d p, Vec3d period);

// vectors
DLL double dist(double ax, double ay, double az, double bx, double by, double bz);
DLL double length(Vec3d v);

//DLL double hypot(double x, double y);
DLL double dot(Vec3d a, Vec3d b);
DLL Vec3d  norm(Vec3d a);
DLL Vec3d  cross(Vec3d a, Vec3d b);
DLL double angle(Vec3d a, Vec3d b);
DLL Vec3d  ortho(Vec3d a, Vec3d b);
DLL Vec3d  up(Vec3d vec, Vec3d upvec);

// variations
DLL double cycle(double index, double loRange, double hiRange);
DLL double pick(int n, double* params);
DLL double choose(int n, double* params);
DLL double wchoose(int n, double* params);
DLL double spline(int n, double* params);

typedef struct VoronoiPointData VoronoiPointData;
DLL void   initVoronoi(VoronoiPointData** cache);
DLL void   freeVoronoi(VoronoiPointData*  cache);
DLL double voronoi (VoronoiPointData* cache, Vec3d p, int type,float jitter, float fbmScale, int fbmOctaves,float fbmLacunarity, float fbmGain);
DLL Vec3d  cvoronoi(VoronoiPointData* cache, Vec3d p, int type,float jitter, float fbmScale, int fbmOctaves,float fbmLacunarity, float fbmGain);
DLL Vec3d  pvoronoi(VoronoiPointData* cache, Vec3d p, int type,float jitter, float fbmScale, int fbmOctaves,float fbmLacunarity, float fbmGain);

typedef struct ColorCurve ColorCurve;
typedef struct NumCurve   NumCurve;
DLL void initColorCurve(ColorCurve** outCurve, int cnt, double const* points, Vec3d const* colors, int const* interps);
DLL void initNumCurve(NumCurve** outCurve, int cnt, double const* points, double const* values, int const* interps);
DLL Vec3d  sampleColorCurve(ColorCurve const* curve, double p);
DLL double sampleNumCurve(NumCurve const* curve, double p);
DLL void freeColorCurve(ColorCurve*  curve);
DLL void freeNumCurve(NumCurve*  curve);

#ifdef __cplusplus
} // extern "C"
#endif

