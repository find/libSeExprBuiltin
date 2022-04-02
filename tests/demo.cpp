#include "../ExprBuiltins.h"
#include "../Curve.h"
#include "fpng.h"
#include <chrono>
#include <iostream>

using namespace SeExpr2;


Vec3d color(double u, double v, Curve<Vec3d> const& curve) {
  auto _uv = Vec3d(u, v, 1);
  auto _offset = Vec3d(0.552941,0.709804,0.262745);
  auto _oct = 7;
  auto _gain = 0.6652;
  auto _lac = 2.7928;
  auto _size = 1.486;

  _uv *= _size;
  _uv += _offset;
  _uv = vturbulence(_uv, _oct, _lac, _gain);
  auto n = fbm(_uv, _oct, _lac, _gain);
  _uv = curve.getValue(n);
  return _uv;
}

int main(int c, char** v) {
  int w = 4096, h = 4096;
  if (c==2) {
    if (w = atoi(v[1]); w==0) {
      std::cerr<<"usage: "<<v[0]<<" [resolution]"<<std::endl;
      return 1;
    } else {
      h = w;
    }
  }

  auto curve = makecurve(0.280255,Vec3d(0.341176,0.0901961,0.372549),4,0.609428,Vec3d(0.823529,0.670588,0.541176),4,0,Vec3d(0,0,0),4,1,Vec3d(1,1,1),4);
  auto now = std::chrono::high_resolution_clock::now();
  uint8_t* img = new uint8_t[w*h*3];
#pragma omp parallel for
  for (int i=0; i<h; ++i) {
    for (int j=0; j<w; ++j) {
      auto pix = color(double(j)/w, double(i)/h, curve);
      size_t offset = i*w*3+j*3;
      img[offset  ] = uint8_t(255*pix[0]);
      img[offset+1] = uint8_t(255*pix[1]);
      img[offset+2] = uint8_t(255*pix[2]);
    }
  }
  auto dur = std::chrono::high_resolution_clock::now()-now;
  std::cout<<"generation : "<<std::chrono::nanoseconds(dur).count()*1e-9<<"s\n";
  now = std::chrono::high_resolution_clock::now();
  fpng::fpng_encode_image_to_file("output.png", img, w,h,3);
  dur = std::chrono::high_resolution_clock::now()-now;
  std::cout<<"writing png: "<<std::chrono::nanoseconds(dur).count()*1e-9<<"s\n";
  delete[] img;
  return 0;
}
