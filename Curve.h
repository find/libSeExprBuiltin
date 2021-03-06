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
#pragma once

#include "Vec.h"
#include <vector>

#include <cfloat>
#include <cassert>

namespace SeExpr2 {

//! Interpolation curve class for double->double and double->Vec3D
/**   Interpolation curve class for mapping from double -> double or double -> Vec3D
      Subject to some interpolation points.

      Each segment is interpolated according to the interpolation type specified on the
      left control point. Interpolation types supported are members of InterpType
      below.

      Valid instantiation types for this are double, or Vec3D
*/
template <class T>
class Curve {
    mutable int cacheCV;

  public:
    //! Supported interpolation types
    enum InterpType {
        kNone = 0,
        kLinear,
        kSmooth,
        kSpline,
        kMonotoneSpline
    };
    struct CV {
        CV(double pos, const T& val, InterpType type) : _pos(pos), _val(val), _interp(type) {}

        double _pos;
        T _val, _deriv;
        InterpType _interp;
    };

  private:
    std::vector<CV> _cvData;
    bool prepared;

  public:
    Curve();

    //! Adds a point to the curve
    void addPoint(double position, const T& val, InterpType type);

    //! Prepares points for evaluation (sorts and computes boundaries, clamps extrema)
    void preparePoints();

    //! Evaluates curve and returns full value
    T getValue(const double param) const;

    //! Evaluates curve for a sub-component of the interpolation values
    //! must call preparePoints() before this is ok to call
    double getChannelValue(const double param, int channel) const;

    //! Returns the control point that is less than the parameter, unless there is no
    //! point, in which case it returns the right point or nothing
    CV getLowerBoundCV(const double param) const;

    //! Returns whether the given interpolation type is supported
    static bool interpTypeValid(InterpType interp);

    //! CV Parameter ordering (cv1._pos < cv2._pos)
    static bool cvLessThan(const CV& cv1, const CV& cv2);

  private:
    //! Performs hermite derivative clamping in canonical space
    void clampCurveSegment(const T& delta, T& d1, T& d2);

    //! Returns a component of the given value
    static double comp(const T& val, const int i);
};

static inline void makecurve_helper(Curve<Vec3d>& curve, double p, Vec3d const& c, int interp)
{
  curve.addPoint(p, c, static_cast<Curve<Vec3d>::InterpType>(interp));
}
template <typename... T>
static inline void makecurve_helper(Curve<Vec3d>& curve, double p, Vec3d const& c, int interp, T... rest)
{
  makecurve_helper(curve, p, c, interp);
  makecurve_helper(curve, rest...);
}
template <typename... T>
static inline Curve<Vec3d> makecurve(double p, Vec3d const& c, int interp, T... rest)
{
  Curve<Vec3d> curve;
  makecurve_helper(curve, p, c, interp);
  makecurve_helper(curve, rest...);
  curve.preparePoints();
  return curve;
}

static inline void makecurve_helper(Curve<double>& curve, double p, double v, int interp)
{
  curve.addPoint(p, v, static_cast<Curve<double>::InterpType>(interp));
}
template <typename... T>
static inline void makecurve_helper(Curve<double>& curve, double p, double v, int interp, T... rest)
{
  makecurve_helper(curve, p, v, interp);
  makecurve_helper(curve, rest...);
}
template <typename... T>
static inline Curve<double> makecurve(double p, double v, int interp, T... rest)
{
  Curve<double> curve;
  makecurve_helper(curve, p, v, interp);
  makecurve_helper(curve, rest...);
  curve.preparePoints();
  return curve;
}

// FYI. here's a little VIM script to covert SeExpr vectors to Vec3d:
// %s/\v\[([^]]*)\]/Vec3d(\1)/g
// after this, ccurve(...) expressions can be translated to makecurve call easily

} // namespace SeExpr2

