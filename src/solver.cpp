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
#include "solver.hpp"

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <limits>
using namespace std;

#include <boost/heap/pairing_heap.hpp>
using namespace boost;

#include "hs.hpp"
#include "options.hpp"
#include "beamSearch.hpp"

const unsigned int inf_ubound = numeric_limits<unsigned int>::max();
const size_t no_pred = numeric_limits<size_t>::max();

// state data structures
typedef unsigned long WState;	// used workers (bitfield)
struct CState {
  CState() : w(0), t(0), l(0), u(inf_ubound), pred(no_pred) {}
  CState(unsigned int lb, unsigned int ub) : w(0), t(0), l(lb), u(ub), pred(no_pred) {}

  WState w;			// used workers
  unsigned int t;               // number of processed tasks
  unsigned int l, u;		// lower and upper bound
  size_t pred;			// index of predecessor state (for extracting solutions)
};

vector<CState> state; // global state table: there can be multiple states for the same workers, with different lower bounds
typedef pair<WState,unsigned int> BState; // (worker state x lower bound)
template <> struct fmt::formatter<BState> {
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) { return ctx.end(); }
  template <typename FormatContext>
  auto format(const BState& s, FormatContext& ctx) -> decltype(ctx.out()) {
    return format_to(ctx.out(), "[{}:{}]", s.first, s.second);
  }
};
map<BState,unsigned long> sindex; // maps a state to its index

struct BState_heap_order {
  inline bool operator()(const BState& a, const BState& b) const {
    return a.second > b.second;
  }
};
typedef heap::pairing_heap<BState,boost::heap::compare<BState_heap_order>> PQ;

const unsigned no_solution = numeric_limits<unsigned>::max();
unsigned bsi = no_solution; // best solution index

vector<PQ> q;
map<BState,PQ::handle_type> handle;

string state_string(const Instance& I, const BState& current) {
  string d;
  for(unsigned w=0; w!=I.m; ++w)
    d+=checkbit(current.first,w)==0?"0":"1";
  d+=" ["+to_string(time_double(current.second));
  if (sindex.find(current)!=sindex.end()) {
    d+="/"+to_string(time_double(state[sindex[current]].l));
    d+=","+to_string(time_double(state[sindex[current]].u));
  }
  d+="]";
  return d;
}

void process_load(const Instance& I, unsigned ci, BState current, unsigned w, unsigned t, unsigned int l, unsigned depth) {
  assert(sindex.find(current)!=sindex.end() && sindex[current]==ci);

  // (1) create and search the new state
  assert(!checkbit(current.first,w));
  setbit(current.first,w);
  l=max(l,current.second); // final loads may be less
  current.second=l;
  auto it=handle.find(current);
  auto si=sindex.find(current);

  // (2) full state
  bool store = true;
  if (t==I.n) {
    unsigned csi = si==sindex.end()?state.size():si->second; // first case relies on step (4) to create this state
    if (bsi==no_solution || state[bsi].l>l) {
      bsi = csi;
      vprint(1,"* {}\n", time_double(l));
    } else {
      store = false;
    }
  // (3) partial state
  } else if (it==handle.end()) { // not active
    if (si == sindex.end() || state[si->second].t<t) { // // new or improving: activate
      handle[current]=q[depth+1].push(current);
    } else { // dominated: skip
    }
  } else { // is active
    assert(si!=sindex.end());
  }

  // (3) update global state table
  if (!store)
    return;

  unsigned delta = t<I.n?I.t[t][w]*precision:0;
  if (si==sindex.end()) {
    // new global state
    CState cs = state[ci];
    cs.w = current.first;
    cs.t = t;
    cs.l = l;
    cs.u = min(cs.u,l+delta);
    cs.pred = ci;
    sindex[current]=state.size();
    state.push_back(cs);
  } else if (state[si->second].t<t) {
    state[si->second].t = t;
    state[si->second].u = min(state[ci].u,state[si->second].l+delta); // coming via another path, u may be inconsistent
    state[si->second].pred = ci;
  }
}

void process_state(const Instance& I, const BState& current, unsigned depth) {
  // (1) find the full state
  assert(sindex.find(current)!=sindex.end());
  auto s = sindex[current];

  if (bsi!=no_solution && state[s].l>state[bsi].l)
    return;

  // (2) create all loads in [l,u)
  for(unsigned w=0; w!=I.m; w++) { // check all workers
    unsigned int t = state[s].t; // t is the next task state to be processed
    unsigned int l = 0;
    if (checkbit(current.first,w))
      continue;

    // (2.1) find the minimal feasible load
    while (t<I.n && l+I.t[t][w]*precision<state[s].l) {
      l += I.t[t][w]*precision;
      t++;
    }
    if (t!=I.n && l+I.t[t][w]*precision<state[s].l) {
      continue;
    }
    if (l>=state[s].u || t==state[s].t) {
      continue;
    }
    process_load(I, s, current, w, t, l, depth);

    // (2.) find all remaining feasible loads
    while (t<I.n) {
      l += I.t[t][w]*precision;
      t++;
      if (l>=state[s].u)
	break;
      process_load(I, s, current, w, t, l, depth);
    }
  }
}

/*
  Combined solver for best cycle time. We keep partial solutions and explore them in a CBFS order.

  Each load has a current state (workers used W, last done task t) plus a range [l,u) that limits the possible cycle times. When branching on such a state, for worker w:
  - We branch on all numbers of tasks n, such that l<=L<u where L=load(t+1,t+n,w).
  - We update l=L, u=min{u,load(t+1,t+n+1,w)}.

  Further:
  - We make probes (i.e. greedy forward allocation within bounds) to detect trivial branches

  Open:
  - We need a better, possible compatible lower bound for choosing candidates.
*/
pair<Time,Solution> csolve(const Instance& I, unsigned int lb, const Solution& S) {
  // (1) create the initial state
  state.clear();
  CState s0(lb,S.value*precision);
  state.push_back(s0);
  sindex[{0,lb}]=0;
  // store the states of the initial solution
  unsigned t = 0;
  for(auto i=0u; i!=I.m; i++) {
    setbit(s0.w,S.w[i]);
    t+=S.l[i];
    s0.l=max(s0.l,unsigned(I.d[t][S.w[i]]-I.d[s0.t][S.w[i]]));
    s0.t=t;
    s0.pred=state.size()-1;
    sindex[{s0.w,s0.l}]=state.size();
    state.push_back(s0);
  }

  // (2) setup priority queues, handles, and push index of initial state
  q.clear();
  q.resize(I.m+1);
  handle.clear();
  handle[{0,lb}]=q[0].push({0,lb});

  // (3) setup expansion
  unsigned mindepth = 0;
  unsigned passes = 0, expansions = 0;
  bool expanded = true;

  // (4) expand
  while (expanded && passes<opt.maxpasses) {
    expanded = false;
    for(unsigned depth=mindepth; depth!=I.m+1; depth++) {
      // take best from current depth and expand
      if (q[depth].size()==0) {
        if (depth==mindepth)
          mindepth++;
        continue;
      }

      expanded = true;
      expansions++;
      auto current = q[depth].top();
      q[depth].pop();
      assert(handle.find(current)!=handle.end());
      handle.erase(current);


      process_state(I,current,depth);
    } // depth loop
    passes++;
    if (passes % 10==0) {
      vprint(2,"{:6}; {:5} passes: [{}] ", state.size(), passes, mindepth);
      for(unsigned depth=mindepth; depth!=I.m+1; ++depth)
	vprint(2,"{} ", q[depth].size());
      vprint(2,"\n");
    }
  } // passes

  // (5) extract the best solution and lower bound
  Solution S_ = S;
  if (bsi!=no_solution) {
    vector<unsigned> wperm;
    // (5.1) get workers from best state
    size_t csi = bsi;
    while (csi != no_pred && state[csi].pred != no_pred) {
      assert(state[csi].pred != csi);
      for(unsigned w=0; w!=I.m; w++)
	if (checkbit(state[csi].w,w) && !checkbit(state[state[csi].pred].w,w)) {
	  wperm.push_back(w);
	  break;
	}
      csi = state[csi].pred;
    }
    reverse(wperm.begin(),wperm.end());
    // (5.2) complete with free workers
    for(unsigned w=0; w!=I.m; w++)
      if (checkbit(state[bsi].w,w)==0) // free
	wperm.push_back(w);
    S_ = optW(I,wperm);
  }
  // set lower bound to ub, or least value of any open node
  Time lb_ = S_.value;
  for(auto q : handle)
    lb_ = min(lb_,time_double(q.first.second));
  
  return {lb_,S_};
}

Solution combined_solver(const Instance& I, Time lb, Solution& S, DP::optnRepDP& O) {
  int verbose = 0;

  // (1) print initial solution
  vprint(1, "Solver starts with bounds [{:7.4},{:7.4}].\n", lb, S.value);

  // (2) run CBFS
  swap(verbose,verbosec.count);
  Solver_Status status;
  tie(S,status)=CBFS(I,S,O);
  swap(verbose,verbosec.count);
  vprint(1, "After CBFS we have bounds [{:7.4},{:7.4}].\n", lb, S.value);

  // the above need: some time management, and an alternating search between good upper and lower bounds, until we conclude that progress is unlikely; then we switch to the method below

  // (3) here comes the rest of the complete solving approach
  Time lb_;
  tie(lb_,S)=csolve(I,lb*precision,S);
  vprint(1, "After csolve we have bounds [{:7.4},{:7.4}].\n", lb_, S.value);

  return S;
}
