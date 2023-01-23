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

#include <vector>
#include <ostream>
#include <unordered_map>
#include <utility>

#include "instance.hpp"
#include "solution.hpp"

// lb1: minimum task times over number of workers
double minSum(const Instance&);

// lb2: optimal cost with worker repetition by DP
double optRep(const Instance&);

namespace DP {
  struct State {
    State() : i(0), j(0), S(0) {}
    State(uint32_t i, uint16_t j, uint16_t S) : i(i), j(j), S(S) {}
    uint32_t i; // next task
    uint16_t j; // # of free workers
    uint16_t S; // set of available controlled workers
    bool operator==(State other) const {
      return i==other.i && j==other.j && S==other.S;
    }
    friend std::ostream& operator<<(std::ostream& o, const State& S) {
      return o << "[" << S.i << ":" << S.j << " " << S.S << "]";
    }
  };
}

namespace std {
  template<>
  struct hash<DP::State> {
    inline std::size_t operator()(const DP::State& S) const {
      return (size_t(S.i) << 32) | (size_t(S.j) << 16) | size_t(S.S);
    }
  };
}

namespace DP {
  struct optnRepDP {
    const Instance& I;
    std::unordered_map<State,double> c; // c[i,j,S] is the optimal cost doing i,i+1,...,n with j workers from [w]\W, and workers S
    std::vector<unsigned> W;
    std::vector<int> ctr; // ctr[w]<0: not controlled, otherwise: bit position

    optnRepDP(const Instance& I);

    unsigned getW(const State& S) const;

    // cost of segment, worker used, next state
    std::tuple<double,unsigned,State> bestSegment(const State& S);
    // best value without repeating workers in W
    double optnRep(State S);
    // set controlled workers
    void setW(const std::vector<unsigned>& W_);
    State getInit() const;
    double opt();
    // extract solution, and an array with workers frequencies
    std::pair<Solution,std::vector<unsigned>> rep();
    // produce a lower bound for a partial worker sequence starting with `pi` at task `i`
    Time getLB(unsigned i, std::vector<unsigned>& pi);
    Time getLB(unsigned i, unsigned long B);
  };
}

// lb3: optimal cost avoiding to repeat up to w workers, determined dynamically
double optnRep(const Instance& I, unsigned wmax, bool=false);
// ditto, with external storage
double optnRep(const Instance& I, DP::optnRepDP& O, unsigned wmax, bool output);
