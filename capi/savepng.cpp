#include "../tests/fpng.h"
#include "../tests/fpng.cpp"

extern "C"
#ifdef _MSC_VER
__declspec(dllexport)
#endif
int savepng(char const* filename, uint8_t const* image, int w, int h, int nchannels)
{
  return fpng::fpng_encode_image_to_file(filename, image, w, h, nchannels);
}
