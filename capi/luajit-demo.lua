local ffi = require('ffi')

ffi.cdef([[
typedef struct Vec3d { double x,y,z; } Vec3d;

// clamping
double clamp(double x, double lo, double hi);
//double round(double x);

// blending / remapping
double invert(double x);
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
Vec3d  hsi(Vec3d rgb, double h, double s, double i, double m);
Vec3d  midhsi(Vec3d rgb, double h, double s, double i, double m, double falloff, int interp);
Vec3d  rgbtohsl(Vec3d rgb);
Vec3d  hsltorgb(Vec3d hsl);
Vec3d  saturate(Vec3d rgb, double amt);

// noise
double hash(int n, double const* args);
double noisep(Vec3d p);
double noise1(double x);
double noise2(double x, double y);
double noise3(double x, double y, double z);
double noise4(double x, double y, double z, double w);
double snoise(Vec3d p);
double snoise4(Vec3d p, double w);
Vec3d  cnoise(Vec3d p);
Vec3d  vnoise(Vec3d p);
Vec3d  cnoise4(Vec3d p, double w);
Vec3d  vnoise4(Vec3d p, double w);
double turbulence(Vec3d p, int octaves, double lacunarity, double gain);
Vec3d  vturbulence(Vec3d p, int octaves, double lacunarity, double gain);
Vec3d  cturbulence(Vec3d p, int octaves, double lacunarity, double gain);
double fbm(Vec3d p, int octaves, double lacunarity, double gain);
Vec3d  vfbm(Vec3d p, int octaves, double lacunarity, double gain);
Vec3d  cfbm(Vec3d p, int octaves, double lacunarity, double gain);
double fbm4(Vec3d p, double t, int octaves, double lacunarity, double gain);
Vec3d  vfbm4(Vec3d p, double t, int octaves, double lacunarity, double gain);
Vec3d  cfbm4(Vec3d p, double t, int octaves, double lacunarity, double gain);
double cellnoise(Vec3d p);
Vec3d  ccellnoise(Vec3d p);
double pnoise(Vec3d p, Vec3d period);

// vectors
double dist(double ax, double ay, double az, double bx, double by, double bz);
double length(Vec3d v);

//double hypot(double x, double y);
double dot(Vec3d a, Vec3d b);
Vec3d  norm(Vec3d a);
Vec3d  cross(Vec3d a, Vec3d b);
double angle(Vec3d a, Vec3d b);
Vec3d  ortho(Vec3d a, Vec3d b);
Vec3d  up(Vec3d vec, Vec3d upvec);

// variations
double cycle(double index, double loRange, double hiRange);
double pick(int n, double* params);
double choose(int n, double* params);
double wchoose(int n, double* params);
double spline(int n, double* params);

typedef struct VoronoiPointData VoronoiPointData;
void   initVoronoi(VoronoiPointData** cache);
void   freeVoronoi(VoronoiPointData*  cache);
double voronoi (VoronoiPointData* cache, Vec3d p, int type,float jitter, float fbmScale, int fbmOctaves,float fbmLacunarity, float fbmGain);
Vec3d  cvoronoi(VoronoiPointData* cache, Vec3d p, int type,float jitter, float fbmScale, int fbmOctaves,float fbmLacunarity, float fbmGain);
Vec3d  pvoronoi(VoronoiPointData* cache, Vec3d p, int type,float jitter, float fbmScale, int fbmOctaves,float fbmLacunarity, float fbmGain);

typedef struct ColorCurve ColorCurve;
typedef struct NumCurve   NumCurve;
void initColorCurve(ColorCurve** outCurve, int cnt, double const* points, Vec3d const* colors, int const* interps);
void initNumCurve(NumCurve** outCurve, int cnt, double const* points, double const* values, int const* interps);
Vec3d  sampleColorCurve(ColorCurve const* curve, double p);
double sampleNumCurve(NumCurve const* curve, double p);
void freeColorCurve(ColorCurve*  curve);
void freeNumCurve(NumCurve*  curve);

int savepng(char const* filename, uint8_t const* image, int w, int h, int nchannels);
]])

local vecmt = {
  __tostring = function(self) return string.format('[%f, %f, %f]', self.x, self.y, self.z) end
}

local se = ffi.load('libse')
local vec = ffi.metatype('Vec3d', vecmt)

--print(se.vturbulence(vec(0.2,0.1,0), 6, 2, 0.5))
local function color(u, v, curve)
  local _uv = vec(u, v, 1)
  local _offset = vec(0.552941,0.709804,0.262745)
  local _oct = 7
  local _gain = 0.6652
  local _lac = 2.7928
  local _size = 1.486

  _uv = vec(_uv.x*_size+_offset.x, _uv.y*_size+_offset.y, _uv.z*_size+_offset.z)
  _uv = se.vturbulence(_uv, _oct, _lac, _gain)
  local n = se.fbm(_uv, _oct, _lac, _gain)
  _uv = se.sampleColorCurve(curve, n)
  return _uv
end

local makecurve = function(...)
  local cnt = select('#',...)/3
  local pts = ffi.new('double[?]', cnt)
  local vls = ffi.new('Vec3d[?]', cnt)
  local itp = ffi.new('int[?]', cnt)
  local args = {...}
  for i=0,cnt-1 do
    pts[i] = args[i*3+1]
    vls[i] = args[i*3+2]
    itp[i] = args[i*3+3]
  end

  local pcurve = ffi.new('ColorCurve*[1]')
  se.initColorCurve(pcurve, cnt, pts, vls, itp)
  return pcurve[0]
end

local curve = makecurve(0.280255,vec(0.341176,0.0901961,0.372549),4,0.609428,vec(0.823529,0.670588,0.541176),4,0,vec(0,0,0),4,1,vec(1,1,1),4);
--print(color(0.5, 0.5, curve))
local w,h = 1024,1024
local buf = ffi.new('uint8_t[?]', w*h*3)
local now = os.clock()
for i=0,h-1 do
  for j=0,w-1 do
    local pix = color(j/w, i/h, curve)
    local off = (i*w+j)*3
    buf[off]   = pix.x*255
    buf[off+1] = pix.y*255
    buf[off+2] = pix.z*255
  end
end
print(string.format('calculation takes %f seconds', os.clock()-now))
se.savepng("output.png", buf, w, h, 3)
print('done.')

se.freeColorCurve(curve)


