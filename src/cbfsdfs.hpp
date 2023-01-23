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
#include "hs.hpp"

using namespace std;

namespace cbfsDfsF {

struct FState {
  unsigned long workers;
  int lastTask;
  vector<unsigned int> workersInUse; //which worker works at workstation 

  inline bool operator<(const FState& other) const {
    return lastTask<other.lastTask;
  }
};

pair<Solution,bool> cbfsDfsFeasibility(const Instance& I, Solution &S); 
}

namespace cbfsDfsO {

struct OState{
  unsigned long workers;
  int lastTask;
  unsigned long c;
  vector<unsigned int> workersInUse;

  inline bool operator<(const OState& other) const {
    if(lastTask<other.lastTask) return true;
    if((lastTask==other.lastTask)&&(c<other.c)) return true;
    return false;
  }
};


pair<Solution,bool> cbfsDfsOptimization(const Instance& I, Solution &S); 
}
