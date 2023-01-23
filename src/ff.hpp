/*
 * Copyright (c) 2023 Jordi Pereira, Marcus Ritt
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#pragma once

#include "util.hpp"

#include <cmath>

////////////////////////
//// Float handling ////
////////////////////////

// Structure to model a float that compares with epsilons
// Example usage:
// ff(x) < ff(y), or ffuple(x, y, 5.0) < ffuple(10.0, 0.0, z)
struct ff {
  static constexpr double EPS = 1e-5;
  static bool fEq(double a, double b) { return std::abs(a - b) < EPS; }
  ff() = default;
  ff(double v) : v(v) {}
  bool operator==(const ff& o) const { return fEq(v, o.v); }
  bool operator!=(const ff& o) const { return not fEq(v, o.v); }
  bool operator<=(const ff& o) const { return v < o.v or fEq(v, o.v); }
  bool operator>=(const ff& o) const { return v > o.v or fEq(v, o.v); }
  bool operator<(const ff& o) const { return v + EPS < o.v; }
  bool operator>(const ff& o) const { return o < *this; }
  operator double() const { return v; }
  static void fix_zero(double& v) { if (ff(v) == ff(0.0)) v = 0.0; }
  static bool isnull(double v) { return ff(v) == ff(0.0); }
  double v;
};
