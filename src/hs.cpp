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
#include "hs.hpp"

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <algorithm>
using namespace std;

#include <boost/multi_array.hpp>
using namespace boost;

#include "ff.hpp"
#include "random.hpp"
#include "options.hpp"
using namespace wap;

inline void setMinMax(const Instance& I, int i, int j, unsigned w, boost::multi_array<Time,2>& C, boost::multi_array<int,2>& l) {
  if (ff(I.segCost(i,I.n-1,w))<=ff(C[I.n][j+1])) {
    C[i][j]=0;
    l[i][j]=I.n-i;
  } else {
    int li = i;     // prefix[i..li-1]<=suffix[li,n]
    int ui = I.n;   // prefix[i..ui-1]> suffix[ui,n]
    while (li+1<ui) {
      int mi = (li+ui)/2;
      if (ff(I.segCost(i,mi-1,w))<=ff(C[mi][j+1]))
	li = mi;
      else
	ui = mi;
    }
    if (ff(I.segCost(i,ui-1,w))<=ff(C[li][j+1])) {
      C[i][j]=I.segCost(i,ui-1,w);
      l[i][j]=ui-i;
    } else {
      C[i][j]=C[li][j+1];
      l[i][j]=li-i;
    }
  }
}

inline void setMinMax_old(const Instance& I, int i, int j, unsigned w, boost::multi_array<Time,2>& C, boost::multi_array<int,2>& l) {
  Time c = 0.0;
  for(unsigned k=i; k!=I.n+1; k++) {
    Time c_ = max(c,C[k][j+1]);
    if (c_<C[i][j]) {
      C[i][j]=c_;
      l[i][j]=k-i;
    }
    if (k<I.n)
      c += I.t[k][w]; // I.d[k+1][w[j]]-I.d[i-1][w[j]]
  }
}

// optimal solution for fixed worker permutation `w`
Solution optW(const Instance& I, const vector<unsigned>& w) {
  Solution result(I,I.m);
  result.w = w;

  boost::multi_array<Time,2> C(extents[I.n+1][I.m+1]); // C[i..n-1][j..m-1]: best subproblem cycle time
  boost::multi_array<int,2> l(extents[I.n+1][I.m+1]);  // l[i..n-1][j..m-1]: corr. interval length

  for(unsigned i=0; i!=I.n; i++)
    C[i][I.m]=inf_time;
  for(unsigned j=0; j!=I.m+1; j++)
    C[I.n][j]=0.0;

  for(int j=I.m-1; j>=0; j--) {
    for(int i=I.n-1; i>=0; i--) {
      C[i][j]=inf_time;
      setMinMax(I,i,j,w[j],C,l);
    }
  }
  unsigned i = 0, j = 0;
  while (i<I.n) {
    result.l[j] = l[i][j];
    if (ff(I.segCost(i,i+l[i][j]-1,w[j]))==ff(C[0][0]))
      result.critical[j]=true;
    i += l[i][j];
    j += 1;
  }
  result.value = C[0][0];
  return result;
}

// optimal solution when all workers are the same with task times `t`
// returns: segment lengths
vector<unsigned> optT(const Instance& I, const vector<Time>& t) {
  assert(t.size()==I.n);

  boost::multi_array<Time,2> C(extents[I.n+1][I.m+1]); // C[i][j]: cycle time for tasks i..n-1, with j workers
  boost::multi_array<int,2> l(extents[I.n+1][I.m+1]);  // corresponding length

  for(unsigned i=0; i!=I.n; i++)
    C[i][I.m]=inf_time;
  for(unsigned j=0; j!=I.m+1; j++)
    C[I.n][j]=0;

  for(int j=I.m-1; j>=0; j--) {
    for(int i=I.n-1; i>=0; i--) {
      C[i][j]=inf_time;
      Time c = 0.0;
      for(unsigned k=i; k!=I.n+1; k++) {
	Time c_ = max(c,C[k][j+1]);
	if (c_<C[i][j]) {
	  C[i][j]=c_;
	  l[i][j]=k-i;
	}
	if (k<I.n)
	  c += t[k];
      }
    }
  }

  vector<unsigned> result(I.m);
  unsigned i = 0, j = 0;
  while (i<I.n) {
    result[j] = l[i][j];
    i += l[i][j];
    j += 1;
  }
  return result;
}

bool optS::bfs(Time c) {
  // put free vertices in vertex queue
  fill(dist.begin(),dist.end(),-1);
  vector<unsigned> Q;
  for(unsigned i=0; i!=I.m; i++)
    if (free(i)) {
      dist[i]=0;
      Q.push_back(i);
    }

  // process vertex queue
  bool found = false;
  while (Q.size()>0 && !found) {
    vector<unsigned> R;
    for(unsigned u : Q) {
      assert(u<I.m);
      for(unsigned v=I.m, ve=2*I.m; v!=ve; v++)  {
	if (cost[u][v-I.m]>c)
	  continue;
	if (dist[v]==-1) {
	  dist[v]=dist[u]+1;
	  if (free(v))
	    found = true;
	  else {
	    dist[peer[v]]=dist[v]+1;
	    R.push_back(peer[v]);
	  }
	}
      }
    }
    Q.swap(R);
  }
  return found;
}

// dfs with bottleneck c
bool optS::dfs(unsigned v, Time c) {
  assert(v>=I.m);

  for(unsigned u=0; u!=I.m; u++)  {
    if (cost[u][v-I.m]>c)
      continue;
    if (dist[u]+1 == dist[v]) {
      if (!free(u)) {
	dist[u]=-1;
	if (dfs(peer[u],c)) {
	  peer[v]=u;
	  peer[u]=v;
	  return true;
	}
	return false;
      } else {
	peer[v]=u;
	peer[u]=v;
	dist[u]=-1;
	return true;
      }
    }
  }
  return false;
}

// hopcroft-karp algorithm with bottleneck c excluding worker w
//   sets globals: peer,dist
unsigned optS::hk(Time c) {
  unsigned result = 0;
  fill(peer.begin(),peer.end(),-1);
  while (bfs(c)) {
    for(unsigned s=I.m; s!=2*I.m; s++)
      if (free(s) && dfs(s,c))
	result++;
  }
  return result;
}

// optimal solution for given segments: bottleneck matching
Solution optS::match(const std::vector<unsigned>& l) {
  assert(l.size()==I.m);
  // compute costs
  for(unsigned w=0; w!=I.m; w++) {
    int j = 0;
    for(unsigned i=0; i!=I.m; i++) {
      cost[i][w]=0.0;
      for(int je=j+l[i]; j!=je; j++)
	cost[i][w]+=I.t[j][w];
    }
  }

  Time t_min = inf_time, t_sum = 0.0;
  for(unsigned w=0; w!=I.m; w++)
    for(unsigned i=0; i!=I.n; i++)
      t_min = min(t_min,I.t[i][w]);
  for(unsigned i=0; i!=I.n; i++)
    t_sum += I.t[i][0];

  // binary search
  Time lt = t_min, ut = t_sum;
  while (lt+t_min<ut) {
    Time t = (lt+ut)/2;
    unsigned m = hk(t);
    if (m==I.m)
      ut = t;
    else
      lt = t+t_min;
  }

  // solution
  hk(ut);
  Solution S(I,I.m);
  for(unsigned i=0; i!=I.m; i++) {
    assert(peer[i]!=-1);
    S.w[i]=peer[i]-I.m;
  }
  S.l=l;
  S.value=ut;
  return S;
}

Solution sample(const Instance& I, int n) {
  Solution S(I);

  vector<unsigned> w(I.m);
  iota(w.begin(),w.end(),0);
  while (n-->0) {
    shuffle(w.begin(),w.end(),wap::rng); // JPG. random order for workers
    Solution c = optW(I,w);
    vprint(4,"Sample {}\n",c.to_string());
    if (c.value<S.value)
      S = c;
  }
  return S;
}

// does worker permutation `wp` permit a feasible solution of time less than `ub`?
// time O(m log n)
  bool feasible(const Instance& I, vector<unsigned>& wp, Time ub) {
  unsigned long C = ub*precision;
  int t = -1; // last done
  for(unsigned i=0u, ie=wp.size(); i!=ie; ++i) { // go over workers
    auto w = wp[i];
    // find largest u in [t,n-1] s.t. I.d[u][w]-I.d[t][w]<C, i.e. I.d[u][w]<C+I.d[t][w]
    unsigned long C_ = C + I.d[t][w];
    if (I.d[I.n-1][w]<C_)
      return true;
    int li = t;     // good
    int ui = I.n-1; // bad
    while (li+1<ui) {
      int mi = (li+ui)/2;
      if (I.d[mi][w]<C_)
	li = mi;
      else
	ui = mi;
    }
    t = li;
  }
  return false;
}

// check if worker permutation `wp` improves `S`, update
bool improved(const Instance& I, Solution& S, vector<unsigned>& wp) {
  if (feasible(I,wp,S.value)) {
    Solution S_ = optW(I,wp);
    assert((ff(S_.value)<ff(S.value)));
    S=S_;
    return true;
  }
  return false;
}

// check if worker permutation `wp` improves `S`, update
bool slow_improved(const Instance& I, Solution& S, vector<unsigned>& wp) {
  Solution S_ = optW(I,wp);
  if (ff(S_.value)<ff(S.value)) {
    S=S_;
    return true;
  } else
    return false;
}

template <typename Fwd>
void rotate_left(Fwd begin, Fwd end) {
  rotate(begin,std::next(begin),end);
}

template <typename Fwd>
void rotate_right(Fwd begin, Fwd end) {
  rotate(begin,prev(end),end);
}

int nb_swap(const Instance& I, Solution& S, int w0) {
  vector<unsigned> wp(S.w);
  for(unsigned w_=w0, w=w0, we=w0+I.m; w_!=we; w_++, w++) {
    if (w==I.m)
      w=0;
    if (!S.critical[w])
      continue;
    for(unsigned x_=w+1, x=w+1, xe=w+I.m; x_!=xe; x_++, x++) { // swap w and x
      if (x==I.m)
	x=0;
      swap(wp[w],wp[x]);
      if (improved(I,S,wp)) return w+1;
      swap(wp[w],wp[x]);
    }
  }
  return -1;
}

int nb_shift(const Instance& I, Solution& S, int w0) {
  vector<unsigned> wp(S.w);
  for(unsigned w_=w0, w=w0, we=w0+I.m; w_!=we; w_++, w++) {
    if (w==I.m)
      w=0;
    if (!S.critical[w])
      continue;
    for(unsigned d=0; d!=w; ++d) { // shift left to d
      rotate_right(wp.begin()+d,wp.begin()+w+1);
      if (improved(I,S,wp)) return w+1;
      rotate_left(wp.begin()+d,wp.begin()+w+1);
    }
    for(unsigned d=w+1; d!=I.m; d++) { // shift right to d
      rotate_left(wp.begin()+w,wp.begin()+d+1);
      if (improved(I,S,wp)) return w+1;
      rotate_right(wp.begin()+w,wp.begin()+d+1);
    }
  }
  return -1;
}

int nb_adj_swap(const Instance& I, Solution& S, int w0) {
  vector<unsigned> wp(S.w);
  for(unsigned w_=w0, w=w0, we=w0+I.m; w_!=we; w_++, w++) {
    if (w==I.m)
      w=0;
    if (!S.critical[w])
      continue;

    unsigned pred = w==0?I.m-1:w-1;
    if (!S.critical[pred]) {
      swap(wp[w],wp[pred]);
      if (improved(I,S,wp)) return w+1;
      swap(wp[w],wp[pred]);
    }

    unsigned succ = w+1==I.m?0:w+1;
    if (!S.critical[succ]) {
      swap(wp[w],wp[succ]);
      if (improved(I,S,wp)) return w+1;
      swap(wp[w],wp[succ]);
    }

  }
  return -1;
}

void pt_swap(const Instance& I, Solution& S) {
  unsigned w = getRandom(0,I.m-2);
  unsigned x = getRandom(int(w+1),I.m-1);
  swap(S.w[w],S.w[x]);
}

void pt_shift(const Instance& I, Solution& S) {
  unsigned w = getRandom(0,I.m-1);
  unsigned d = getRandom(0,I.m-2);
  if (d<w)
      rotate_right(S.w.begin()+d,S.w.begin()+w+1);
  else {
    d++;
    rotate_left(S.w.begin()+w,S.w.begin()+d+1);
  }
}

void pt_adj_swap(const Instance& I, Solution& S) {
  unsigned w = getRandom(0,I.m-1);
  if (w<I.m-1)
    swap(S.w[w],S.w[w+1]);
  else
    swap(S.w[w],S.w[0]);
}

void polish(const Instance& I, Solution &S) {
  vprint(2,"polish: init {:7.6} ",S.value);
  S = optW(I,S.w);
  vprint(2,"opt {:7.6} ",S.value);
  all_nb(I,S);
}

void all_nb(const Instance& I, Solution &S) {
  vprint(2,"all_nb: init {:7.6} ",S.value);
  ls(I,S,nb_adj_swap);
  vprint(2,"adj {:7.6} ",S.value);
  ls(I,S,nb_shift);
  vprint(2,"shift {:7.6} ",S.value);
  ls(I,S,nb_swap);
  vprint(2,"swap {:7.6} ",S.value);
}

void ILS::run() {
  Solution S = sample(I,1);
  B = S;
  tbest = ::S.elapsed();
  all_nb(I,S);
  if (ff(S.value)<ff(B.value)) B=S;
  vprint(1," {:7.6}\n",B.value);

  while (::S.elapsed()<opt.tlim && ff(B.value)>ff(lb)) {
    Solution P = S;
    for(auto i=0u; i!=opt.pstrength; ++i)
      switch (getRandom(0,2)) {
        case 0: pt_swap(I,P); break;
        case 1: pt_shift(I,P); break;
        case 2: pt_adj_swap(I,P); break;
      }
    polish(I,P);
    iter++; bdist++;
    if (ff(P.value)<ff(B.value)) {
      S=B=P;
      btotal+=bdist;
      bdist=0;
      bcount++;
      tbest = ::S.elapsed();
      vprint(1,"* {:4} {:6.2f} {:6.2} {:6.2f} {:7.6f} \n",iter,::S.elapsed(),double(btotal)/bcount,tbest,B.value);
    } else if (ff(P.value)<ff(B.value*1.05)) {
      S=P;
      vprint(2,"✔️ {:4} {:6.2f} {:6.2} {:6.2f} {:7.6f} \n",iter,::S.elapsed(),double(btotal)/bcount,tbest,B.value);
    } else if (getRandom()<0.015) {
      S=B;
      bdist=0;
      vprint(2,"< {:4} {:6.2f} {:6.2} {:6.2f} {:7.6f} \n",iter,::S.elapsed(),double(btotal)/bcount,tbest,B.value);
    } else
      vprint(2,"x {:4} {:6.2f} {:6.2} {:6.2f} {:7.6f} \n",iter,::S.elapsed(),double(btotal)/bcount,tbest,B.value);
  }
}

void ILS::setLB(Time nlb) {
  lb = nlb;
}
