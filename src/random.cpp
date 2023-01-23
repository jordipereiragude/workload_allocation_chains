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
#include "random.hpp"

#include <fstream>
using namespace std;

namespace wap {
  mt19937 rng;

  unsigned setupRandom(unsigned seed) {
    if (seed==0) {
      seed = time(0);
      ifstream f("/dev/urandom");
      if (f.good()) {
	f.read((char *)(&seed), sizeof(unsigned int));
      }
    }
    rng.seed(seed);
    srand48(seed);
    return seed;
  }

  int getRandom(int min, int max) {
#if defined(JORDI)
    return lrand48()%(max-min+1)+min;
#else
    uniform_int_distribution<> U(min,max);
    return U(rng);
#endif
  }

  double getRandom() {
#if defined(JORDI)
    return drand48();
#else
    uniform_real_distribution<> U;
    return U(rng);
#endif
  }
  double getRandom(double min,double max) {
    uniform_real_distribution<> U(min,max);
    return U(rng);
  }
}
