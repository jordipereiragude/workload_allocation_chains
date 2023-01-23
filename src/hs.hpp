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

#include <boost/multi_array.hpp>

#include "util.hpp"
#include "instance.hpp"
#include "solution.hpp"

Solution optW(const Instance& I, const std::vector<unsigned>& w); // optimal for fixed permutation
std::vector<unsigned> optT(const Instance& I, const std::vector<Time>& t); // optimal for constant times
Solution sample(const Instance& I, int n); // sample n permutations

template<typename Construct, typename Improve>
Solution multi_start(const Instance& I, int n, Construct construct, Improve improve) {
  Solution S(I);

  std::vector<unsigned> w(I.m);
  iota(w.begin(),w.end(),0);
  while (n-->0) {
    construct(I,w);
    Solution c = optW(I,w);
    improve(I,c);
    vprint(4,"Sample {}\n",c.to_string());
    if (c.value<S.value)
      S = c;
  }
  return S;
}

// optimal solution for given segments, by matching them to workers
// 1..m: segments, m+1..m+m: workers
struct optS {
  const Instance& I;
  boost::multi_array<Time,2> cost; // cost for doing segment i by worker j

  std::vector<int> peer,dist;	// matching
  bool free(int v) { return peer[v]==-1; }

  optS(const Instance& I) : I(I), cost(boost::extents[I.m][I.m]), peer(2*I.m), dist(2*I.m) {};

  bool bfs(Time c);
  bool dfs(unsigned v, Time c);
  unsigned hk(Time c);
  Solution match(const std::vector<unsigned>& l);
};

//// standard neighborhoods

// adjacent swap, swap, shift neighborhoods, starting from `w0` (for round-robin); return next worker to continue search
int nb_adj_swap(const Instance&,Solution&,int w0);
int nb_swap(const Instance&,Solution&,int w0);
int nb_shift(const Instance&,Solution&,int w0);

// local search
template<typename N>
void ls(const Instance& I, Solution& S, const N& n) {
  int w = 0;
  do {
    w = n(I,S,w);
  } while (w>-1);
}

// solution polishing: optimal segments (optW) plus local searches
void polish(const Instance&, Solution &);
// only local searches
void all_nb(const Instance&, Solution &);
// iterated local search
struct ILS {
  const Instance& I;
  unsigned iter = 0, bdist = 0, bcount = 0, btotal = 0;
  Solution B;
  Time tbest;
  Time lb = 0;
  
  ILS(const Instance& I) : I(I), B(I) {}

  void setLB(Time);
  void run();
};
