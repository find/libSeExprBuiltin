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

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <cfloat>

#include "Vec.h"
#include "Curve.h"
#include "ExprBuiltins.h"
#include "Noise.h"

namespace SeExpr2 {

static const char* fabs_docstring = "float abs(float x)\nabsolute value of x";

// angle conversion functions
static const char* deg_docstring = "float deg(float angle)\nradians to degrees";
static const char* rad_docstring = "float deg(float angle)\ndegrees to radians";
// trig in degrees
static const char* cosd_docstring = "float cosd(float angle)\ncosine in degrees";
static const char* sind_docstring = "float sind(float angle)\nsine in degrees";
static const char* tand_docstring = "float tand(float angle)\ntangent in degrees";
static const char* acosd_docstring = "float acosd(float angle)\narc cosine in degrees";
static const char* asind_docstring = "float asind(float angle)\narc sine in degrees";
static const char* atand_docstring = "float atand(float angle)\narc tangent in degrees";
static const char* atan2d_docstring =
    "float atan2d(float y,float x)\narc tangent in degrees of y/x between -180 and 180";
// trig in radians
static const char* cos_docstring = "float cos(float angle)\ncosine in radians";
static const char* sin_docstring = "float sin(float angle)\nsine in radians";
static const char* tan_docstring = "float tan(float angle)\ntangent in radians";
static const char* acos_docstring = "float acos(float angle)\narc cosine in radians";
static const char* asin_docstring = "float asin(float angle)\narc sine in radians";
static const char* atan_docstring = "float atan(float angle)\narc tangent in radians";
static const char* atan2_docstring = "float atan2(float y,float x)\narc tangent in radians of y/x between -PI and PI";
// hyperbolic trig
static const char* cosh_docstring = "float cosh(float angle)\nhyperbolic cosine in radians";
static const char* sinh_docstring = "float sinh(float angle)\nhyperbolic sine in radians";
static const char* tanh_docstring = "float tanh(float angle)\nhyperbolic tangent in radians";
static const char* acosh_docstring = "float acosh(float angle)\nhyperbolic arc cosine in radians";
static const char* asinh_docstring = "float asinh(float angle)\nhyperbolic arc sine in radians";
static const char* atanh_docstring = "float atanh(float angle)\nhyperbolic arc tangent in radians";
// clamping/rounding
static const char* clamp_docstring = "float clamp(float x,float lo,float hi)\nconstrain x to range [lo,hi]";
static const char* round_docstring = "float round(float x)\nconstrain x to range [lo,hi]";
static const char* max_docstring = "float max(float a,float b)\ngreater of a and b";
static const char* min_docstring = "float min(float a,float b)\nlesser of a and b";
static const char* trunc_docstring = "float trunc(float a)\nnearest integer towards zero";
static const char* floor_docstring = "float floor(float a)\nnext lower integer";
static const char* ceil_docstring = "float ceil(float a)\nnext higher integer";
// misc math
static const char* invert_docstring = "float invert(float a)\nDefined as 1-x";
static const char* cbrt_docstring = "float cbrt(float x)\ncube root";
static const char* sqrt_docstring = "float sqrt(float x)\nsquare root";
static const char* exp_docstring = "float exp(float x)\nE raised to the x power";
static const char* pow_docstring = "float pow(float x)\nx to the y power, also available as ^";
static const char* log_docstring = "float log(float x)\nNatural logarithm";
static const char* log10_docstring = "float log10(float x)\nBase 10 logarithm";
static const char* fmod_docstring = "float fmod(float x,float y)\nremainder of x/y (also available as % operator)";
static const char* turbulence_docstring =
    "float turbulence(vector v,int octaves=6,float lacunarity=2,float gain=.5)\nAbsolute value of each noise term is "
    "taken. This gives billowy appearance";
static const char* cturbulence_docstring =
    "color cturbulence(vector v,int octaves=6,float lacunarity=2,float gain=.5)\nAbsolute value of each noise term is "
    "taken. This gives billowy appearance";
static const char* vturbulence_docstring =
    "vector vturbulence(vector v,int octaves=6,float lacunarity=2,float gain=.5)\nAbsolute value of each noise term is "
    "taken. This gives billowy appearance";

double compress(double x, double lo, double hi) { return (hi - lo) * x + lo; }
static const char* compress_docstring = "float compress(float x,float lo,float hi)\nRemaps x in [0,1] to [lo,hi]";

double expand(double x, double lo, double hi) {
    if (lo == hi) return x < lo ? 0 : 1;
    return (x - lo) / (hi - lo);
}
static const char* expand_docstring = "float expand(float x,float lo,float hi)\nRemaps x in [lo,hi] to [0,1]";

double fit(double x, double a1, double b1, double a2, double b2) {
    return (x * (b2 - a2) - a1 * b2 + b1 * a2) / (b1 - a1);
}
static const char* fit_docstring =
    "float fit(float x,float a1,float b1,float a2,float b2)\nLinearly remaps x in [a1,b1] to [a2,b2]";

double gamma(double x, double g) { return pow(x, 1 / g); }
static const char* gamma_docstring = "float gamma(float x, float g)\nGamma correction of x with gamma factor g";

double bias(double x, double b) {
    static double C = 1 / log(0.5);
    return pow(x, log(b) * C);
}
static const char* bias_docstring =
    "float bias(float x, float g)\nVariation of gamma where values less than 0.5 pull the curve down\nand values "
    "greater than 0.5 pull the curve up\npow(x,log(b)/log(0.5))";

double contrast(double x, double c) {
    if (x < 0.5)
        return 0.5 * bias(1 - c, 2 * x);
    else
        return 1 - 0.5 * bias(1 - c, 2 - 2 * x);
}
static const char* contrast_docstring =
    "float contrast(float x,float x)\nAdjust the contrast.&nbsp; For c from 0 to 0.5, the contrast is decreased.&nbsp; "
    "For c &gt; 0.5, the contrast is increased.";

double boxstep(double x, double a) { return x < a ? 0.0 : 1.0; }
static const char* boxstep_docstring = "float boxstep(float x,float a)\n if x < a then 0 otherwise 1";

double linearstep(double x, double a, double b) {
    if (a < b) {
        return x < a ? 0 : (x > b ? 1 : (x - a) / (b - a));
    } else if (a > b) {
        return 1 - (x < b ? 0 : (x > a ? 1 : (x - b) / (a - b)));
    }
    return boxstep(x, a);
}
static const char* linearstep_docstring =
    "float linearstep(float x,float a,float b)\n if x &lt; a then 0, if x &gt; b then 1, and\nx transitions linearly "
    "when &lt; x &lt; b ";

double smoothstep(double x, double a, double b) {
    if (a < b) {
        if (x < a) return 0;
        if (x >= b) return 1;
        x = (x - a) / (b - a);
    } else if (a > b) {
        if (x <= b) return 1;
        if (x > a) return 0;
        x = 1 - (x - b) / (a - b);
    } else
        return boxstep(x, a);
    return x * x * (3 - 2 * x);
}
static const char* smoothstep_docstring =
    "float smoothstep(float x,float a,float b)\n if x &lt; a then 0, if x &gt; b then 1, and\nx transitions smoothly "
    "(cubic) when &lt; x &lt; b";

double gaussstep(double x, double a, double b) {
    if (a < b) {
        if (x < a) return 0;
        if (x >= b) return 1;
        x = 1 - (x - a) / (b - a);
    } else if (a > b) {
        if (x <= b) return 1;
        if (x > a) return 0;
        x = (x - b) / (a - b);
    } else
        return boxstep(x, a);
    return pow(2, -8 * x * x);
}
static const char* gaussstep_docstring =
    "float gasussstep(float x,float a,float b)\n if x &lt; a then 0, if x &gt; b then 1, and\nx transitions smoothly "
    "(exponentially) when &lt; x &lt; b";

double remap(double x, double source, double range, double falloff, double interp) {
    range = fabs(range);
    falloff = fabs(falloff);

    if (falloff == 0) return fabs(x - source) < range;

    double a, b;
    if (x > source) {
        a = source + range;
        b = a + falloff;
    } else {
        a = source - range;
        b = a - falloff;
    }

    switch (int(interp)) {
        case 0:
            return linearstep(x, b, a);
        case 1:
            return smoothstep(x, b, a);
        default:
            return gaussstep(x, b, a);
    }
}
static const char* remap_docstring =
    "remap(float x, float\n"
    "source, float range, float falloff, float interp)\nGeneral remapping function.\n"
    "When x is within +/- <i>range</i> of source, the result is one.\n"
    "The result falls to zero beyond that range over <i>falloff</i> distance.\n"
    "The falloff shape is controlled by <i>interp</i>. Numeric values\n"
    "or named constants may be used:\n"
    "&nbsp;&nbsp;&nbsp;&nbsp;int <b>linear</b>\n"
    "= 0\n"
    "&nbsp;&nbsp;&nbsp;&nbsp;int <b>smooth</b> = 1\n"
    "&nbsp;&nbsp;&nbsp;&nbsp;int <b>gaussian</b> = 2\n";

double mix(double x, double y, double alpha) { return x * (1 - alpha) + y * alpha; }
static const char* mix_docstring = "mix(float a,float b,float alpha)\nBlend of a and b according to alpha.";

Vec3d hsiAdjust(const Vec3d& rgb, double h, double s, double i) {
    Vec3d hsl = rgbtohsl(rgb);
    hsl[0] += h * (1.0 / 360);
    hsl[1] *= s;
    return hsltorgb(hsl) * i;
}

Vec3d hsi(const Vec3d& rgb, double h, double s, double i, double m) {
    h *= m;
    s = (s - 1) * m + 1;
    i = (i - 1) * m + 1;
    return hsiAdjust(rgb, h, s, i);
}
static const char* hsi_docstring =
    "color  hsi(color x, float h, float s, float i, float map=1)\n"
    "The hsi function shifts the hue by h\n"
    "(in degrees) and scales the saturation and intensity by s and i\n"
    "respectively.&nbsp; An map may be supplied which will control the shift\n"
    "- the full shift will happen when the map is one and no shift will\n"
    "happen when the map is zero.&nbsp; The shift will be scaled back for\n"
    "values between zero and one.";

Vec3d midhsi(Vec3d const& rgb, double h, double s, double i, double m, double falloff, int interp) {
    // remap from [0..1] to [-1..1]
    m = m * 2 - 1;
    // add falloff (if specified)
    if (m < 0)
        m = -remap(-m, 1, 0, falloff, interp);
    else
        m = remap(m, 1, 0, falloff, interp);

    // scale hsi values according to mask (both directions)
    h *= m;
    float absm = fabs(static_cast<float>(m));
    s = s * absm + 1 - absm;
    i = i * absm + 1 - absm;
    if (m < 0) {
        s = 1 / s;
        i = 1 / i;
    }
    return hsiAdjust(rgb, h, s, i);
}
static const char* midhsi_docstring =
    "color midhsi(color x, float h, float s, float i, float map, float falloff=1, int interp=0)\n"
    "The midhsi function is just like the hsi function except that\n"
    "the control map is centered around the mid point (value of 0.5)\n"
    "and can scale the shift in both directions.";

Vec3d rgbtohsl(const Vec3d& rgb) {
    // RGB to HSL color space conversion
    // This is based on Foley, Van Dam (2nd ed; p. 595)
    // but extended to allow rgb values outside of 0..1
    double R, G, B, H, S, L, x, y, sum, diff;
    R = rgb[0];
    G = rgb[1];
    B = rgb[2];
    x = R < G ? (R < B ? R : B) : (G < B ? G : B);  // min(R,G,B)
    y = R > G ? (R > B ? R : B) : (G > B ? G : B);  // max(R,G,B)

    // compute lightness = avg of min and max rgb vals
    sum = x + y;
    diff = y - x;
    L = sum / 2;
    if (diff < 1e-6)  // achromatic
        return Vec3d(0, 0, L);

    // compute saturation
    if (L <= .5) {
        if (x < 0)
            S = 1 - x;
        else
            S = diff / sum;
    } else {
        if (y > 1)
            S = y;
        else
            S = diff / (2 - sum);
    }

    // compute hue
    if (R == y)
        H = (G - B) / diff;
    else if (G == y)
        H = (B - R) / diff + 2;
    else
        H = (R - G) / diff + 4;
    H *= 1 / 6.;
    H -= floor(H);  // make sure hue is in range 0..1

    return Vec3d(H, S, L);
}
static const char* rgbtohsl_docstring =
    "color rgbtohsl(color rgb)\n"
    "RGB to HSL color space conversion.\n"
    "HSL is Hue, Saturation, Lightness (all in range [0..1] )\n"
    "These functions have also been extended to support rgb and hsl values\n"
    "outside of the range [0..1] in a reasonable way.&nbsp; For any rgb or\n"
    "hsl value (except for negative s values), the conversion is\n"
    "well-defined and reversible.";

static double hslvalue(double x, double y, double H) {
    H -= floor(H);  // make sure hue is in range 0..1

    if (H < 1 / 6.)
        return x + (y - x) * H * 6;
    else if (H < 3 / 6.)
        return y;
    else if (H < 4 / 6.)
        return x + (y - x) * (4 / 6. - H) * 6;
    else
        return x;
}

Vec3d hsltorgb(const Vec3d& hsl) {
    // HSL to RGB color space conversion
    // This is based on Foley, Van Dam (2nd ed; p. 596)
    // but extended to allow rgb values outside of 0..1
    double H, S, L, R, G, B, x, y;
    H = hsl[0];
    S = hsl[1];
    L = hsl[2];
    if (S <= 0)  // achromatic
        return Vec3d(L, L, L);

    // find min/max rgb values
    if (L < 0.5) {
        if (S > 1)
            y = 2 * L + S - 1;
        else
            y = L + L * S;
    } else {
        if (S > 1)
            y = S;
        else
            y = L + S - L * S;
    }
    x = 2 * L - y;

    // reconstruct rgb from min,max,hue
    R = hslvalue(x, y, H + (1 / 3.));
    G = hslvalue(x, y, H);
    B = hslvalue(x, y, H - (1 / 3.));
    return Vec3d(R, G, B);
}
static const char* hsltorgb_docstring =
    "color hsltorgb(color hsl)\n"
    "RGB to HSL color space conversion.\n"
    "HSL is Hue, Saturation, Lightness (all in range [0..1] )\n"
    "These functions have also been extended to support rgb and hsl values\n"
    "outside of the range [0..1] in a reasonable way.&nbsp; For any rgb or\n"
    "hsl value (except for negative s values), the conversion is\n"
    "well-defined and reversible.";

Vec3d saturate(const Vec3d& Cin, double amt) {
    const Vec3d lum(.2126, .7152, .0722);  // rec709 luminance
    Vec3d result = Vec3d(Cin.dot(lum) * (1 - amt)) + Cin * amt;
    if (result[0] < 0) result[0] = 0;
    if (result[1] < 0) result[1] = 0;
    if (result[2] < 0) result[2] = 0;
    return result;
}
static const char* saturate_docstring =
    "color saturate(color val, float amt)\n"
    "Scale saturation of color by amt.\n"
    "The color is scaled around the rec709 luminance value,\n"
    "and negative results are clamped at zero.\n";

double hash(int n, double const* args) {
    // combine args into a single seed
    uint32_t seed = 0;
    for (int i = 0; i < n; i++) {
        // make irrational to generate fraction and combine xor into 32 bits
        int exp = 0;
        double frac = frexp(args[i] * double(M_E * M_PI), &exp);
        uint32_t s = (uint32_t)(frac * UINT32_MAX) ^ (uint32_t)exp;

        // blend with seed (constants from Numerical Recipes, attrib. from Knuth)
        static const uint32_t M = 1664525, C = 1013904223;
        seed = seed * M + s + C;
    }

    // tempering (from Matsumoto)
    seed ^= (seed >> 11);
    seed ^= (seed << 7) & 0x9d2c5680UL;
    seed ^= (seed << 15) & 0xefc60000UL;
    seed ^= (seed >> 18);

    // permute
    static unsigned char p[256] = {
        148, 201, 203, 34,  85,  225, 163, 200, 174, 137, 51,  24,  19,  252, 107, 173, 110, 251, 149, 69,  180, 152,
        141, 132, 22,  20,  147, 219, 37,  46,  154, 114, 59,  49,  155, 161, 239, 77,  47,  10,  70,  227, 53,  235,
        30,  188, 143, 73,  88,  193, 214, 194, 18,  120, 176, 36,  212, 84,  211, 142, 167, 57,  153, 71,  159, 151,
        126, 115, 229, 124, 172, 101, 79,  183, 32,  38,  68,  11,  67,  109, 221, 3,   4,   61,  122, 94,  72,  117,
        12,  240, 199, 76,  118, 5,   48,  197, 128, 62,  119, 89,  14,  45,  226, 195, 80,  50,  40,  192, 60,  65,
        166, 106, 90,  215, 213, 232, 250, 207, 104, 52,  182, 29,  157, 103, 242, 97,  111, 17,  8,   175, 254, 108,
        208, 224, 191, 112, 105, 187, 43,  56,  185, 243, 196, 156, 246, 249, 184, 7,   135, 6,   158, 82,  130, 234,
        206, 255, 160, 236, 171, 230, 42,  98,  54,  74,  209, 205, 33,  177, 15,  138, 178, 44,  116, 96,  140, 253,
        233, 125, 21,  133, 136, 86,  245, 58,  23,  1,   75,  165, 92,  217, 39,  0,   218, 91,  179, 55,  238, 170,
        134, 83,  25,  189, 216, 100, 129, 150, 241, 210, 123, 99,  2,   164, 16,  220, 121, 139, 168, 64,  190, 9,
        31,  228, 95,  247, 244, 81,  102, 145, 204, 146, 26,  87,  113, 198, 181, 127, 237, 169, 28,  93,  27,  41,
        231, 248, 78,  162, 13,  186, 63,  66,  131, 202, 35,  144, 222, 223};
    union {
        uint32_t i;
        unsigned char c[4];
    } u1, u2;
    u1.i = seed;
    u2.c[3] = p[u1.c[0]];
    u2.c[2] = p[(u1.c[1] + u2.c[3]) & 0xff];
    u2.c[1] = p[(u1.c[2] + u2.c[2]) & 0xff];
    u2.c[0] = p[(u1.c[3] + u2.c[1]) & 0xff];

    // scale to [0.0 .. 1.0]
    return u2.i * (1.0 / UINT32_MAX);
}
static const char* hash_docstring =
    "float hash(float seed1,[float seed2, ...])\n"
    "Like rand, but with no internal seeds. Any number of seeds may be given\n"
    "and the result will be a random function based on all the seeds.";

double noise(Vec3d const& p) {
  double result;
  Noise<3,1>(&p[0], &result);
  return result*0.5+0.5;
}
double noise(double x) {
  double result;
  Noise<1,1>(&x, &result);
  return result*0.5+0.5;
}
double noise(double x, double y) {
  double result;
  double in[] = {x,y};
  Noise<2,1>(in, &result);
  return result*0.5+0.5;
}
double noise(double x, double y, double z) {
  double result;
  double in[] = {x,y,z};
  Noise<3,1>(in, &result);
  return result*0.5+0.5;
}
double noise(double x, double y, double z, double w) {
  double result;
  double in[] = {x,y,z,w};
  Noise<4,1>(in, &result);
  return result*0.5+0.5;
}
static const char* noise_docstring =
    "float noise ( vector v ) <br>\n"
    "float noise ( float x, float y )\n"
    "float noise ( float x, float y, float z )\n"
    "float noise ( float x, float y, float z, float w )\n"
    "Original perlin noise at location (C2 interpolant)";

double snoise(const Vec3d& p) {
    double result;
    double args[3] = {p[0], p[1], p[2]};
    Noise<3, 1>(args, &result);
    return result;
}
static const char* snoise_docstring =
    "float snoise ( vector v)\n"
    "signed noise w/ range -1 to 1 formed with original perlin noise at location (C2 interpolant)";

Vec3d vnoise(const Vec3d& p) {
    Vec3d result;
    double args[3] = {p[0], p[1], p[2]};
    Noise<3, 3>(args, &result[0]);
    return result;
}
static const char* vnoise_docstring =
    "vector vnoise ( vector v)\n"
    "vector noise formed with original perlin noise at location (C2 interpolant)";

Vec3d cnoise(const Vec3d& p) { return .5 * vnoise(p) + Vec3d(.5); }
static const char* cnoise_docstring =
    "color cnoise ( vector v)\n"
    "color noise formed with original perlin noise at location (C2 interpolant)";

double snoise(const Vec3d& p, double w) {
    double result;
    double args[4] = {p[0], p[1], p[2], w};
    Noise<4, 1>(args, &result);
    return result;
}
static const char* snoise4_docstring =
    "float snoise4 ( vector v,float t)\n"
    "4D signed noise w/ range -1 to 1 formed with original perlin noise at location (C2 interpolant)";

Vec3d vnoise(const Vec3d& p, double w) {
    Vec3d result;
    double args[4] = {p[0], p[1], p[2], w};
    Noise<4, 3>(args, &result[0]);
    return result;
}
static const char* vnoise4_docstring =
    "vector vnoise4 ( vector v,float t)\n"
    "4D vector noise formed with original perlin noise at location (C2 interpolant)";

Vec3d cnoise(const Vec3d& p, double w) { return .5 * vnoise(p, w) + Vec3d(.5); }
static const char* cnoise4_docstring =
    "color cnoise4 ( vector v,float t)\n"
    "4D color noise formed with original perlin noise at location (C2 interpolant)";

double turbulence(const Vec3d& p, int octaves, double lacunarity, double gain) {
    double result = 0;
    double P[3] = {p[0], p[1], p[2]};
    FBM<3, 1, true>(P, &result, octaves, lacunarity, gain);
    return .5 * result + .5;
}

Vec3d vturbulence(const Vec3d& p, int octaves, double lacunarity, double gain) {
    Vec3d result;
    double P[3] = {p[0], p[1], p[2]};
    FBM<3, 3, true>(P, &result[0], octaves, lacunarity, gain);
    return result;
}

Vec3d cturbulence(const Vec3d& p, int octaves, double lacunarity, double gain) {
  return vturbulence(p, octaves, lacunarity, gain) * .5 + Vec3d(.5);
}

double fbm(const Vec3d& p, int octaves, double lacunarity, double gain) {
    double result = 0.0;
    double P[3] = {p[0], p[1], p[2]};
    FBM<3, 1, false>(P, &result, octaves, lacunarity, gain);
    return .5 * result + .5;
}
static const char* fbm_docstring =
    "float fbm(vector v,int octaves=6,float lacunarity=2,float gain=.5)\n"
    "fbm (Fractal Brownian Motion) is a multi-frequency noise function. \n"
    "The base frequency is the same as the \"noise\" function. The total \n"
    "number of frequencies is controlled by octaves. The lacunarity is the \n"
    "spacing between the frequencies - a value of 2 means each octave is \n"
    "twice the previous frequency. The gain< controls how much each \n"
    "frequency is scaled relative to the previous frequency.";

Vec3d vfbm(const Vec3d& p, int octaves, double lacunarity, double gain) {
    Vec3d result = 0.0;
    double P[3] = {p[0], p[1], p[2]};
    FBM<3, 3, false>(P, &result[0], octaves, lacunarity, gain);
    return result;
}
static const char* vfbm_docstring = "vector vfbm(vector vint octaves=6,float lacunarity=2,float gain=.5)";

double fbm4(const Vec3d& p, double time, int octaves, double lacunarity, double gain) {
    double result = 0.0;
    double P[4] = {p[0], p[1], p[2], time};
    FBM<4, 1, false>(P, &result, octaves, lacunarity, gain);
    return .5 * result + .5;
}
static const char* fbm4_docstring =
    "float fbm4(vector v,float time,int octaves=6,float lacunarity=2,float gain=.5)\n"
    "fbm (Fractal Brownian Motion) is a multi-frequency noise function. \n"
    "The base frequency is the same as the \"noise\" function. The total \n"
    "number of frequencies is controlled by octaves. The lacunarity is the \n"
    "spacing between the frequencies - a value of 2 means each octave is \n"
    "twice the previous frequency. The gain< controls how much each \n"
    "frequency is scaled relative to the previous frequency.";

Vec3d vfbm4(const Vec3d& p, double time, int octaves, double lacunarity, double gain) {
    Vec3d result = 0.0;
    double P[4] = {p[0], p[1], p[2], time};
    FBM<4, 3, false>(P, &result[0], octaves, lacunarity, gain);
    return result;
}
static const char* vfbm4_docstring = "vector vfbm4(vector v,float time,int octaves=6,float lacunarity=2,float gain=.5)";

Vec3d cfbm(const Vec3d& p, int octaves, double lacunarity, double gain) {
  return vfbm(p, octaves, lacunarity, gain) * .5 + Vec3d(.5);
}
static const char* cfbm_docstring = "color cfbm(vector vint octaves=6,float lacunarity=2,float gain=.5)";

Vec3d cfbm4(const Vec3d& p, double time, int octaves, double lacunarity, double gain) {
  return vfbm4(p, time, octaves, lacunarity, gain) * .5 + Vec3d(.5);
}
static const char* cfbm4_docstring = "color cfbm4(vector v,float time,int octaves=6,float lacunarity=2,float gain=.5)";

double cellnoise(const Vec3d& p) {
    double result;
    double args[3] = {p[0], p[1], p[2]};
    CellNoise<3, 1>(args, &result);
    return result;
}
static const char* cellnoise_docstring =
    "float cellnoise(vector v)\n"
    "cellnoise generates a field of constant colored cubes based on the integer location.\n"
    "This is the same as the prman cellnoise function.";

Vec3d ccellnoise(const Vec3d& p) {
    Vec3d result;
    double args[3] = {p[0], p[1], p[2]};
    CellNoise<3, 3>(args, &result[0]);
    return result;
}
static const char* ccellnoise_docstring =
    "color cellnoise(vector v)\n"
    "cellnoise generates a field of constant colored cubes based on the integer location.\n"
    "This is the same as the prman cellnoise function.";

double pnoise(const Vec3d& p, const Vec3d& period) {
    double result;
    double args[3] = {p[0], p[1], p[2]};
    int pargs[3] = {std::max((int)1, (int)period[0]),
                    std::max((int)1, (int)period[1]),
                    std::max((int)1, (int)period[2])};
    PNoise<3, 1>(args, pargs, &result);
    return result;
}
static const char* pnoise_docstring =
    "float pnoise ( vector v, vector period )\n"
    "periodic noise";

static Vec3d* voronoi_points(VoronoiPointData& data, const Vec3d& cell, double jitter) {
    if (cell == data.cell && jitter == data.jitter) return data.points;
    data.cell = cell;
    data.jitter = jitter;

    int n = 0;
    for (double i = -1; i <= 1; i++) {
        for (double j = -1; j <= 1; j++) {
            for (double k = -1; k <= 1; k++, n++) {
                Vec3d testcell = cell + Vec3d(i, j, k);
                data.points[n] = testcell + jitter * (ccellnoise(testcell) - Vec3d(.5));
            }
        }
    }
    return data.points;
}

static void voronoi_f1_3d(VoronoiPointData& data, const Vec3d& p, double jitter, double& f1, Vec3d& pos1) {
    // from Advanced Renderman, page 257
    Vec3d thiscell(floor(p[0]) + 0.5, floor(p[1]) + 0.5, floor(p[2]) + 0.5);

    f1 = 1000;
    Vec3d* pos = voronoi_points(data, thiscell, jitter);
    Vec3d* end = pos + 27;

    for (; pos != end; pos++) {
        Vec3d offset = *pos - p;
        double dist = offset.dot(offset);
        if (dist < f1) {
            f1 = dist;
            pos1 = *pos;
        }
    }
    f1 = sqrt(f1);
}

static void voronoi_f1f2_3d(VoronoiPointData& data,
                            const Vec3d& p,
                            double jitter,
                            double& f1,
                            Vec3d& pos1,
                            double& f2,
                            Vec3d& pos2) {
    // from Advanced Renderman, page 258
    Vec3d thiscell(floor(p[0]) + 0.5, floor(p[1]) + 0.5, floor(p[2]) + 0.5);
    f1 = f2 = 1000;
    Vec3d* pos = voronoi_points(data, thiscell, jitter);
    Vec3d* end = pos + 27;

    for (; pos != end; pos++) {
        Vec3d offset = *pos - p;
        double dist = offset.dot(offset);
        if (dist < f1) {
            f2 = f1;
            pos2 = pos1;
            f1 = dist;
            pos1 = *pos;
        } else if (dist < f2) {
            f2 = dist;
            pos2 = *pos;
        }
    }
    f1 = sqrt(f1);
    f2 = sqrt(f2);
}

double voronoi(VoronoiPointData& data, Vec3d p, int type, float jitter, float fbmScale, int fbmOctaves ,float fbmLacunarity, float fbmGain) {
    if (fbmScale > 0) {
        p += fbmScale * vfbm(2*p, fbmOctaves, fbmLacunarity, fbmGain);
    }

    double f1, f2;
    Vec3d pos1, pos2;
    if (type >= 3)
        voronoi_f1f2_3d(data, p, jitter, f1, pos1, f2, pos2);
    else
        voronoi_f1_3d(data, p, jitter, f1, pos1);

    switch (type) {
        case 1:
            pos1[0] += 10;
            return cellnoise(pos1);
        case 2:
            return f1;
        case 3:
            return f2;
        case 4:
            return f2 - f1;
        case 5: {
            float scalefactor = static_cast<float>((pos2 - pos1).length() / ((pos1 - p).length() + (pos2 - p).length()));
            return smoothstep(f2 - f1, 0, 0.1 * scalefactor);
        }
    }

    return 0.0;
}
const static char* voronoi_docstring =
    "float voronoi(vector v, int type=1,float jitter=0.5, float fbmScale=0, int fbmOctaves=4,float fbmLacunarity=2, "
    "float fbmGain=.5)\n"
    "voronoi is a cellular noise pattern. It is a jittered variant of cellnoise.";

Vec3d cvoronoi(VoronoiPointData& data, Vec3d p, int type, float jitter, float fbmScale, int fbmOctaves ,float fbmLacunarity, float fbmGain) {
    if (fbmScale > 0) {
        p += fbmScale * vfbm(2*p, fbmOctaves, fbmLacunarity, fbmGain);
    }

    double f1, f2;
    Vec3d pos1, pos2;
    if (type >= 3)
        voronoi_f1f2_3d(data, p, jitter, f1, pos1, f2, pos2);
    else
        voronoi_f1_3d(data, p, jitter, f1, pos1);

    Vec3d color = ccellnoise(pos1);
    switch (type) {
        case 1:
            pos1[0] += 10;
            return color;
        case 2:
            return f1 * color;
        case 3:
            return f2 * color;
        case 4:
            return (f2 - f1) * color;
        case 5: {
            float scalefactor = static_cast<float>((pos2 - pos1).length() / ((pos1 - p).length() + (pos2 - p).length()));
            return smoothstep(f2 - f1, 0, 0.1 * scalefactor) * color;
        }
    }

    return 0.0;
}
const static char* cvoronoi_docstring =
    "color cvoronoi(vector v, int type=1,float jitter=0.5, float fbmScale=0, int fbmOctaves=4,float fbmLacunarity=2, "
    "float fbmGain=.5)\n"
    "returns color in cellular pattern. It is a jittered variant of cellnoise.";

Vec3d pvoronoi(VoronoiPointData& data, Vec3d p, int type, float jitter, float fbmScale, int fbmOctaves ,float fbmLacunarity, float fbmGain) {
    if (fbmScale > 0) {
        p += fbmScale * vfbm(2*p, fbmOctaves, fbmLacunarity, fbmGain);
    }

    double f1;
    Vec3d pos1;
    voronoi_f1_3d(data, p, jitter, f1, pos1);
    return pos1;
}
const static char* pvoronoi_docstring =
    "color pvoronoi(vector v, int type=1,float jitter=0.5, float fbmScale=0, int fbmOctaves=4,float fbmLacunarity=2, "
    "float fbmGain=.5)\n"
    "returns center of voronoi cell.";

double dist(double ax, double ay, double az, double bx, double by, double bz) {
    double x = ax - bx;
    double y = ay - by;
    double z = az - bz;
    return sqrt(x * x + y * y + z * z);
}
static const char* dist_docstring =
    "float dist(vector a, vector b)\n"
    "distance between two points";

double length(const Vec3d& v) { return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); }
static const char* length_docstring =
    "float length(vector v)\n"
    "length of vector";

double hypot(double x, double y) { return sqrt(x * x + y * y); }
static const char* hypot_docstring =
    "float hypot(vector v)\n"
    "length of 2d vector [x,y]";

double dot(const Vec3d& a, const Vec3d& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
static const char* dot_docstring =
    "float dot(vector a,vector b)\n"
    "vector dot product";

Vec3d norm(const Vec3d& a) {
    double len = length(a);
    if (len == 0)
        return 0.0;
    else
        return a / len;
}
static const char* norm_docstring =
    "vector norm(vector v)\n"
    "vector scaled to unit length";

Vec3d cross(const Vec3d& a, const Vec3d& b) {
    return Vec3d(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
static const char* cross_docstring =
    "vector cross(vector a,vector b)\n"
    "vector cross product";

double angle(const Vec3d& a, const Vec3d& b) {
    double len = length(a) * length(b);
    if (len == 0) return 0;
    return acos(dot(a, b) / len);
}
static const char* angle_docstring =
    "float angle(vector a,vector b)\n"
    "angle between two vectors (in radians)";

Vec3d ortho(const Vec3d& a, const Vec3d& b) { return norm(cross(a, b)); }
static const char* ortho_docstring =
    "vector angle(vector a,vector b)\n"
    "normalized vector orthogonal to a and b scaled to unit length";

Vec3d rotate(int n, const Vec3d* args) {
    if (n != 3) return 0.0;
    const Vec3d& P = args[0];
    const Vec3d& axis = args[1];
    float angle = static_cast<float>(args[2][0]);
    double len = axis.length();
    if (!len) return P;
    return P.rotateBy(axis / len, angle);
}
static const char* rotate_docstring =
    "vector rotate(vector v,vector axis,float angle)\n"
    "rotates v around axis by given angle (in radians)";

Vec3d up(const Vec3d& P, const Vec3d& upvec) {
    // rotate vec so y-axis points to upvec
    Vec3d yAxis(0, 1, 0);
    return P.rotateBy(ortho(upvec, yAxis), angle(upvec, yAxis));
}
static const char* up_docstring =
    "vector up(vector P,vector upvec)\n"
    "rotates v such that the Y axis points in the given up direction";

double cycle(double index, double loRange, double hiRange) {
    int lo = int(loRange);
    int hi = int(hiRange);
    int range = hi - lo + 1;
    if (range <= 0) return lo;
    int result = int(index) % range;
    if (result < 0) result += range;
    return lo + result;
}
static const char* cycle_docstring =
    "int cycle(int index, int loRange, int hiRange )\n"
    "Cycles through values between loRange and hiRange based on supplied index.\n"
    "This is an offset \"mod\" function. The result is rotates v such that the\n"
    "Y axis points in the given up direction";

double pick(int n, double* params) {
    if (n < 3) return 0;
    double index = hash(1, &params[0]);
    int loRange = int(params[1]);
    int hiRange = int(params[2]);
    int range = hiRange - loRange + 1;
    if (range <= 0) return loRange;
    int numWeights = n - 3;
    if (numWeights > range) numWeights = range;

    // build cutoff points based on weights
    double* cutoffs = (double*)alloca(sizeof(double) * range);
    double* weights = (double*)alloca(sizeof(double) * range);
    double total = 0;
    for (int i = 0; i < range; i++) {
        double weight = i < numWeights ? params[i + 3] : 1;
        total += weight;
        cutoffs[i] = total;
        weights[i] = weight;
    }

    if (total == 0) return loRange;

    // scale value from [0..1] to [0..total] range
    index *= total;

    // bsearch cutoff table to find index that spans value
    int lo = 0, hi = range - 1;
    while (lo < hi) {
        int m = (lo + hi) / 2;
        if (index <= cutoffs[m])
            hi = m;
        else
            lo = m + 1;
    }

    // skip zero-length intervals
    if (weights[lo] == 0) {
        if (lo > 0 && cutoffs[lo] > 0)  // scan backward if possible
            while (--lo > 0 && weights[lo] == 0)
                ;
        else if (lo < range - 1)  // else scan forward if possible
            while (++lo < range - 1 && weights[lo] == 0)
                ;
    }

    // add offset and return result
    return loRange + lo;
}
static const char* pick_docstring =
    "int pick(float index, int loRange, int hiRange, [float weights, ...] )\n"
    "Picks values randomly between loRange and hiRange based on supplied index (which is\n"
    "automatically hashed).&nbsp; The values will be distributed according\n"
    "to the supplied weights.&nbsp; Any weights not supplied are assumed to\n"
    "be 1.0.";

double swatch(int n, double* params) { return choose(n, params); }
static const char* swatch_docstring =
    "color swatch(float index, color choice0, color choice1, color choice2, [...])\n"
    "Chooses one of the supplied color choices based on the index (assumed to be in range [0..1]).";

double choose(int n, double* params) {
    if (n < 3) return 0;
    double key = params[0];
    // NaN protection
    if (key != key) return 0;
    int nvals = n - 1;
    return params[1 + int(clamp(key * nvals, 0, nvals - 1))];
}
static const char* choose_docstring =
    "float choose(float index,float choice1, float choice2, [...])\n"
    "Chooses one of the supplied choices based on the index (assumed to be in range [0..1]).";

double wchoose(int n, double* params) {
    if (n < 5) return 0;
    double key = params[0];
    // NaN protection
    if (key != key) return 0;
    int nvals = (n - 1) / 2;  // nweights = nvals

    // build cutoff points based on weights
    double* cutoffs = (double*)alloca(sizeof(double) * nvals);
    double* weights = (double*)alloca(sizeof(double) * nvals);
    double total = 0;
    for (int i = 0; i < nvals; i++) {
        double weight = params[i * 2 + 2];
        total += weight;
        cutoffs[i] = total;
        weights[i] = weight;
    }

    if (total == 0) return params[1];

    // scale value from [0..1] to [0..total] range
    key *= total;

    // bsearch cutoff table to find index that spans value
    int lo = 0, hi = nvals - 1;
    while (lo < hi) {
        int m = (lo + hi) / 2;
        if (key <= cutoffs[m])
            hi = m;
        else
            lo = m + 1;
    }

    // skip zero-length intervals
    if (weights[lo] == 0) {
        if (lo > 0 && cutoffs[lo] > 0)  // scan backward if possible
            while (--lo > 0 && weights[lo] == 0)
                ;
        else if (lo < nvals - 1)  // else scan forward if possible
            while (++lo < nvals - 1 && weights[lo] == 0)
                ;
    }

    // return corresponding value
    return params[lo * 2 + 1];
}
static const char* wchoose_docstring =
    "float wchoose(float index,float choice1, float weight1, float choice2, float weight2, [...] )\n"
    "Chooses one of the supplied choices based on the index (assumed to be in range[0..1]).\n"
    "The values will be distributed according to the supplied weights.";

double spline(int n, double* params) {
    if (n < 5) return 0;
    double u = clamp(params[0], 0, 1);
    if (u == 0) return params[2];
    if (u == 1) return params[n - 2];
    int nsegs = n - 4;
    double seg;
    u = modf(u * nsegs, &seg);
    double* p = &params[int(seg) + 1];
    double u2 = u * u;
    double u3 = u2 * u;
    return 0.5 * (p[0] * (-u3 + 2 * u2 - u) + p[1] * (3 * u3 - 5 * u2 + 2) + p[2] * (-3 * u3 + 4 * u2 + u) +
                  p[3] * (u3 - u2));
}
static const char* spline_docstring =
    "float spline(float param,float y1,float y2,float y3,float y4,[...])\n\n"
    "Interpolates a set of values to the parameter specified where y1, ..., yn are\n"
    "distributed evenly from [0...1]";

}
