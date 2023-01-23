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
#include <utility>

#include "util.hpp"
#include "hs.hpp"
#include "status.hpp"

namespace sadpF {
  std::pair<Solution,Solver_Status> sadpFeasibility(const Instance& I, Solution &S, bool=false); 
}

namespace sadpO {

  struct State {
    int tasks;
    int numImplicitWorkers;
    unsigned *explicitWorkers;
  };

  std::pair<Solution,bool> sadpOptimization(const Instance& I, Solution &S); 
  unsigned long checkLowerBoundSADP(const Instance& I, unsigned long current, int t);
  void freeMemoryCBFSO();
}

// starting after `lastTaskPerformed` return maximum task that worker `w` can do staying below `currentC`
// (helper, also used in beamSearch.cpp)
inline int maxTask(const Instance& I, unsigned w, int lastTaskPerformed, long int currentC) {
  long unsigned int maxValue=I.d[lastTaskPerformed][w]+currentC;
  if (I.d[I.n-1][w]<maxValue) {
    return I.n-1;
  } else {
    int li = lastTaskPerformed; // <  maxValue
    int ui = I.n-1;		// >= maxValue
    while (li+1<ui) {
      int mi = (li+ui)/2;
      if (I.d[mi][w]<maxValue)
	li = mi;
      else
	ui = mi;
    }
    return li;
  }
}
