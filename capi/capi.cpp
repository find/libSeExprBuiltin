#include "../ExprBuiltins.h"
#include "../Curve.h"
#include "capi.h"

namespace S = SeExpr2;

static inline Vec3d V(S::Vec3d const& v) { return *reinterpret_cast<Vec3d const*>(&v); }
static inline S::Vec3d V(Vec3d const& v) { return *reinterpret_cast<S::Vec3d const*>(&v); }

// clamping
extern "C" double clamp(double x, double lo, double hi) { return S::clamp(x,lo,hi); }

//extern "C" double round(double x) { return S::round(x); } // CRT has this

// blending / remapping
extern "C" double invert(double x) { return S::invert(x); }
extern "C" double compress(double x, double lo, double hi) { return S::compress(x,lo,hi); }
extern "C" double expand(double x, double lo, double hi) { return S::expand(x,lo,hi); }
extern "C" double fit(double x, double a1, double b1, double a2, double b2) { return S::fit(x,a1,b1,a2,b2); }
extern "C" double gamma(double x, double g) { return S::gamma(x, g); }
extern "C" double bias(double x, double b) { return S::bias(x, b); }
extern "C" double contrast(double x, double c) { return S::contrast(x, c); }
extern "C" double boxstep(double x, double a) { return S::boxstep(x, a); }
extern "C" double linearstep(double x, double a, double b) { return S::linearstep(x, a, b); }
extern "C" double smoothstep(double x, double a, double b) { return S::smoothstep(x, a, b); }
extern "C" double gaussstep(double x, double a, double b) { return S::gaussstep(x, a, b); }
extern "C" double remap(double x, double s, double r, double f, double interp) { return S::remap(x, s, r, f, interp); }
extern "C" double mix(double x, double y, double alpha) { return S::mix(x,y,alpha); }
extern "C" Vec3d  hsi(Vec3d rgb, double h, double s, double i, double m) { return V(S::hsi(V(rgb), h, s, i, m)); }
extern "C" Vec3d  midhsi(Vec3d rgb, double h, double s, double i, double m, double falloff, int interp) { return V(S::midhsi(V(rgb), h, s, i, m, falloff, interp)); }
extern "C" Vec3d  rgbtohsl(Vec3d rgb) { return V(S::rgbtohsl(V(rgb))); }
extern "C" Vec3d  hsltorgb(Vec3d hsl) { return V(S::hsltorgb(V(hsl))); }
extern "C" Vec3d  saturate(Vec3d rgb, double amt) { return V(S::saturate(V(rgb), amt)); }

// noise
extern "C" double hash(int n, double const* args) { return S::hash(n, args); }
extern "C" double noisep(Vec3d p) { return S::noise(V(p)); }
extern "C" double noise1(double x) { return S::noise(x); }
extern "C" double noise2(double x, double y) { return S::noise(x,y); }
extern "C" double noise3(double x, double y, double z) { return S::noise(x,y,z); }
extern "C" double noise4(double x, double y, double z, double w) { return S::noise(x,y,z,w); }
extern "C" double snoise(Vec3d p) { return S::snoise(V(p)); }
extern "C" double snoise4(Vec3d p, double w) { return S::snoise(V(p), w); }
extern "C" Vec3d  cnoise(Vec3d p) { return V(S::cnoise(V(p))); }
extern "C" Vec3d  vnoise(Vec3d p) { return V(S::vnoise(V(p))); }
extern "C" Vec3d  cnoise4(Vec3d p, double w) { return V(S::cnoise(V(p), w)); }
extern "C" Vec3d  vnoise4(Vec3d p, double w) { return V(S::vnoise(V(p), w)); }
extern "C" double turbulence(Vec3d p, int octaves, double lacunarity, double gain) { return S::turbulence(V(p), octaves, lacunarity, gain); }
extern "C" Vec3d  vturbulence(Vec3d p, int octaves, double lacunarity, double gain) { return V(S::vturbulence(V(p), octaves, lacunarity, gain)); }
extern "C" Vec3d  cturbulence(Vec3d p, int octaves, double lacunarity, double gain) { return V(S::cturbulence(V(p), octaves, lacunarity, gain)); }
extern "C" double fbm(Vec3d p, int octaves, double lacunarity, double gain) { return S::fbm(V(p), octaves, lacunarity, gain); }
extern "C" Vec3d  vfbm(Vec3d p, int octaves, double lacunarity, double gain) { return V(S::vfbm(V(p), octaves, lacunarity, gain)); }
extern "C" Vec3d  cfbm(Vec3d p, int octaves, double lacunarity, double gain) { return V(S::cfbm(V(p), octaves, lacunarity, gain)); }
extern "C" double fbm4(Vec3d p, double t, int octaves, double lacunarity, double gain) { return S::fbm4(V(p), t, octaves, lacunarity, gain); }
extern "C" Vec3d  vfbm4(Vec3d p, double t, int octaves, double lacunarity, double gain) { return V(S::vfbm4(V(p), t, octaves, lacunarity, gain)); }
extern "C" Vec3d  cfbm4(Vec3d p, double t, int octaves, double lacunarity, double gain) { return V(S::cfbm4(V(p), t, octaves, lacunarity, gain)); }
extern "C" double cellnoise(Vec3d p) { return S::cellnoise(V(p)); }
extern "C" Vec3d  ccellnoise(Vec3d p) { return V(S::ccellnoise(V(p))); }
extern "C" double pnoise(Vec3d p, Vec3d period) { return S::pnoise(V(p), V(period)); }

// vectors
extern "C" double dist(double ax, double ay, double az, double bx, double by, double bz) { return S::dist(ax,ay,az,bx,by,bz); }
extern "C" double length(Vec3d v) { return S::length(V(v)); }
//extern "C" double hypot(double x, double y) { return S::hypot(x,y); } // crt has this
extern "C" double dot(Vec3d a, Vec3d b) { return S::dot(V(a), V(b)); }
extern "C" Vec3d  norm(Vec3d a) { return V(S::norm(V(a))); }
extern "C" Vec3d  cross(Vec3d a, Vec3d b) { return V(S::cross(V(a), V(b))); }
extern "C" double angle(Vec3d a, Vec3d b) { return S::angle(V(a), V(b)); }
extern "C" Vec3d  ortho(Vec3d a, Vec3d b) { return V(S::ortho(V(a), V(b))); }
extern "C" Vec3d  up(Vec3d vec, Vec3d upvec) { return V(S::up(V(vec), V(upvec))); }

// variations
extern "C" double cycle(double index, double loRange, double hiRange) { return S::cycle(index, loRange, hiRange); }
extern "C" double pick(int n, double* params) { return S::pick(n, params); }
extern "C" double choose(int n, double* params) { return S::choose(n, params); }
extern "C" double wchoose(int n, double* params) { return S::wchoose(n, params); }
extern "C" double spline(int n, double* params) { return S::spline(n, params); }

struct VoronoiPointData {
  S::VoronoiPointData data;
};

extern "C" void   initVoronoi(VoronoiPointData** cache) { *cache = new VoronoiPointData; }
extern "C" void   freeVoronoi(VoronoiPointData*  cache) { delete cache; }

extern "C" double voronoi (VoronoiPointData* cache, Vec3d p, int type,float jitter, float fbmScale, int fbmOctaves,float fbmLacunarity, float fbmGain) {
  return S::voronoi(cache->data, V(p), type, jitter, fbmScale, fbmOctaves, fbmLacunarity, fbmGain);
}
extern "C" Vec3d  cvoronoi(VoronoiPointData* cache, Vec3d p, int type,float jitter, float fbmScale, int fbmOctaves,float fbmLacunarity, float fbmGain) {
  return V(S::cvoronoi(cache->data, V(p), type, jitter, fbmScale, fbmOctaves, fbmLacunarity, fbmGain));
}
extern "C" Vec3d  pvoronoi(VoronoiPointData* cache, Vec3d p, int type,float jitter, float fbmScale, int fbmOctaves,float fbmLacunarity, float fbmGain) {
  return V(S::pvoronoi(cache->data, V(p), type, jitter, fbmScale, fbmOctaves, fbmLacunarity, fbmGain));
}

struct ColorCurve {
  S::Curve<S::Vec3d> curve;
};
struct NumCurve {
  S::Curve<double> curve;
};
extern "C" void initColorCurve(ColorCurve** outCurve, int cnt, double const* points, Vec3d const* colors, int const* interps) {
  *outCurve = new ColorCurve;
  for (int i=0; i<cnt; ++i)
    S::makecurve_helper(outCurve[0]->curve, points[i], V(colors[i]), interps[i]);
  outCurve[0]->curve.preparePoints();
}
extern "C" void initNumCurve(NumCurve** outCurve, int cnt, double const* points, double const* values, int const* interps) {
  *outCurve = new NumCurve;
  for (int i=0; i<cnt; ++i)
    S::makecurve_helper(outCurve[0]->curve, points[i], values[i], interps[i]);
  outCurve[0]->curve.preparePoints();
}
extern "C" Vec3d  sampleColorCurve(ColorCurve const* curve, double p) {
  return V(curve->curve.getValue(p));
}
extern "C" double sampleNumCurve(NumCurve const* curve, double p) {
  return curve->curve.getValue(p);
}
extern "C" void freeColorCurve(ColorCurve* curve) {
  delete curve;
}
extern "C" void freeNumCurve(NumCurve* curve) {
  delete curve;
}


