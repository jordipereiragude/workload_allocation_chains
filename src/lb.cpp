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
#include "lb.hpp"

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <limits>
using namespace std;

#include <boost/multi_array.hpp>
#include <boost/container_hash/hash.hpp>
using namespace boost;

#include "util.hpp"

// lb1: minimum task times over number of workers
double minSum(const Instance& I) {
  double tsum = 0.0;
  for(unsigned i=0; i!=I.n; i++) {
    double t = I.t[i][0];
    for(unsigned w=1; w!=I.m; w++)
      t = min(t,I.t[i][w]);
    tsum += t;
  }
  return tsum/I.m;
}

// lb2: optimal cost with worker repetition by DP
double optRep(const Instance& I) {

  multi_array<double,2> c(extents[I.n+1][I.m+1]); // c[i,j] is the optimal cost for doing i,...,n-1, i=0,...,n-1 with j workers

  // single worker: cheapest total segment cost over all workers
  vector<double> seg(I.m,0);
  for(int i=I.n-1; i>=0; i--) {
    seg[0] += I.t[i][0];
    c[i][1]=seg[0];
    for(unsigned w=1; w!=I.m; w++) {
      seg[w]+=I.t[i][w];
      c[i][1]=min(c[i][1],seg[w]);
    }
  }
  // no tasks: no cost
  for(unsigned j=1; j!=I.m+1; j++)
    c[I.n][j]=0;
  c[I.n][0]=numeric_limits<double>::max();

  // all other segments and workers: cheapest cost over all workers
  for(int i=I.n-1; i>=0; i--) {
    fill(seg.begin(),seg.end(),0);
    for(unsigned j=2; j!=I.m+1; j++)
      c[i][j]=c[i][j-1];
    for(unsigned k=i; k!=I.n; k++) {
      // segment [i,k]
      for(unsigned w=0; w!=I.m; w++)
	seg[w] += I.t[k][w];
      double best = *min_element(seg.begin(),seg.end());
      for(unsigned j=2; j!=I.m+1; j++)
	c[i][j]=min(c[i][j],max(best,c[k+1][j-1]));
    }
  }

  return c[0][I.m];
}

namespace DP {

  optnRepDP::optnRepDP(const Instance& I) : I(I), ctr(I.m,-1) {}

  unsigned optnRepDP::getW(const State& S) const {
    unsigned nw = S.j;
    for(unsigned w=0; w!=I.m; w++)
      if (ctr[w]>=0 && checkbit(S.S,ctr[w]))
	nw++;
    return nw;
  }

  // cost of segment, worker used, next state
  tuple<double,unsigned,State> optnRepDP::bestSegment(const State& S) {
    const unsigned nw = getW(S);

    double bc = numeric_limits<double>::max();
    unsigned bw = 0;
    State bS;
    for(unsigned w=0; w!=I.m; w++) {
      // (1) select w, if possible
      State T(S);
      if (ctr[w]<0 && S.j>0)
	T.j--;
      else if (ctr[w]>=0 && checkbit(S.S,ctr[w])) {
	unsetbit(T.S,ctr[w]);
      } else continue;

      // (2) sole worker: assign all
      if (nw==1) {
	T.i=I.n; double c1 = I.segCost(S.i,I.n-1,w); //max(I.segCost(S.i,I.n-1,w),optnRep(T));
	if (c1<bc) {
	  bS=T; bw=w; bc=c1;
	}
	continue;
      }

      // (3) binary search for best segment split
      unsigned l = S.i, u = I.n;
      // inv: for l: seg<res; for u: seg>=res
      unsigned k = S.i+(I.n-S.i)/nw;
      while (l+1<u) {
	double seg = I.segCost(S.i,k-1,w);
	T.i=k;
	double res = optnRep(T);
	if (seg<res)
	  l = k;
	else
	  u = k;
	k = (l+u)/2;
      }
      T.i=l; double c1 = max(I.segCost(S.i,l-1,w),optnRep(T));
      if (c1<bc) {
	bS=T; bw=w; bc=c1;
      }
      T.i=u; double c2 = max(I.segCost(S.i,u-1,w),optnRep(T));
      if (c2<bc) {
	bS=T; bw=w; bc=c2;
      }
    }
    return {bc,bw,bS};
  }

  // best value without repeating workers in W
  inline double optnRepDP::optnRep(State S) {
    assert(W.size()<=16);
    double r;
    auto p = c.find(S);
    if (p == c.end()) {
      // base: done
      if (S.i>=I.n)
	r=0;
      else if (S.j==0 && S.S==0)
	r=numeric_limits<double>::max();
      else {
	auto s = bestSegment(S);
	assert(get<0>(s)<=c[get<2>(s)]);
	r=max(get<0>(s),c[get<2>(s)]);
      }
      c[S]=r;
      return r;
    }
    return p->second;
  }

  // set controlled workers
  void optnRepDP::setW(const vector<unsigned>& W_) {
    assert(W_.size()<16); // restriction in State to uint16
    W=W_;
    fill(ctr.begin(),ctr.end(),-1);
    unsigned b = 0;
    for(auto w : W)
      ctr[w]=b++;
    c.clear();
  }

  State optnRepDP::getInit() const {
    State S(0,I.m-W.size(),0);
    for(auto w : W)
      setbit(S.S,ctr[w]);
    return S;
  }

  double optnRepDP::opt() {
    return optnRep(getInit());
  }

  // extract solution, and an array with workers frequencies
  pair<Solution,vector<unsigned>> optnRepDP::rep() {
    Solution sol(I);
    sol.value=0;
    State S(getInit());
    vector<unsigned> N(I.m,0);
    while (S.i<I.n) {
      auto s = bestSegment(S);
      N[get<1>(s)]++;
      sol.w.push_back(get<1>(s));
      sol.l.push_back(get<2>(s).i-S.i);
      sol.value=max(sol.value,I.segCost(S.i,get<2>(s).i-1,get<1>(s)));
      S = get<2>(s);
    }
    return {sol,N};
  }

  // produce a lower bound for a partial worker sequence starting with `pi` at task `i`
  Time optnRepDP::getLB(unsigned i, vector<unsigned>& pi) {
    State S(i,I.m-W.size(),0);
    for(auto w : W)
      setbit(S.S,ctr[w]);
    for(auto w : pi)
      if (ctr[w]>=0)
	unsetbit(S.S,ctr[w]); // controlled worker: mark used
      else
	S.j--; // one free worker less
    return optnRep(S);
  }
  // ditto with a bit state
  Time optnRepDP::getLB(unsigned i, unsigned long B) {
    State S(i,I.m-W.size(),0);
    for(auto w : W)
      setbit(S.S,ctr[w]);
    for(auto w=0u; w!=I.m; w++) {
      if (checkbit(B,w)==0)
	continue;
      if (ctr[w]>=0)
	unsetbit(S.S,ctr[w]); // controlled worker: mark used
      else
	S.j--; // one free worker less
    }
    return optnRep(S);
  }
}

// lb3: optimal cost avoiding to repeat up to w workers, determined dynamically
double optnRep(const Instance& I, unsigned wmax, bool output) {
  DP::optnRepDP O(I);
  return optnRep(I,O,wmax,output);
}

// ditto, with external storage
double optnRep(const Instance& I, DP::optnRepDP& O, unsigned wmax, bool output) {
  vector<unsigned> W;

  double lb;
  vector<unsigned> N;
  Solution s(I);

  do {
    timer lb3;
    O.setW(W);
    lb = O.opt();
    tie(s,N) = O.rep();
    vprint(1,"optnRep: lower bound {:7.4f} in time {:.2f} with workers {{{}}} in {} states.\n{}\n", lb, lb3.elapsed_secs(), fmt::join(W,","),O.c.size(),s.to_string()); // view::transform(W,[](auto w) { return w+1; }) with C++20

    if (output) {
      fmt::print("{} {:.2f} ",lb,lb3.elapsed_secs());
      fflush(stdout);
    }

    // find index of max. freq. worker `mw` with highest number of tasks
    unsigned nw = *max_element(N.begin(),N.end());
    unsigned mw = max_element(N.begin(),N.end())-N.begin();
    vector<unsigned> tt(I.m,0);
    for(unsigned w=0; w!=I.m; w++)
      tt[s.w[w]]+=s.l[w];
    unsigned mt = tt[mw];
    for(unsigned w=0; w!=I.m; w++)
      if (N[w]==nw && tt[w]>mt) {
	mt=tt[w];
	mw=w;
      }

    assert(N[mw]>0);
    if (N[mw]==1) {
      if (output) {
	while (W.size()<wmax) {
	  W.push_back(0);
	  fmt::print("{} NA ",lb);
	}
      }
      break;
    }
    W.push_back(mw);
  } while (W.size()<=wmax);
  return lb;
}
