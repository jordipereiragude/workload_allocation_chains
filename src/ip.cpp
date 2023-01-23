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
#include "ip.hpp"

#include <string>
#include <queue>
#include <unordered_set>
using namespace std;

#include <boost/functional/hash.hpp>
using namespace boost;

#include "util.hpp"
#include "ff.hpp"
#include "options.hpp"

#include <tuple>

namespace std {
  template<typename... T>
  struct hash<tuple<T...>> {
    size_t operator()(tuple<T...> const& arg) const noexcept {
        return boost::hash_value(arg);
    }
  };
}

map<IloAlgorithm::Status,string>
status_string = {
		 { IloAlgorithm::Unknown, "Unknown" },
		 { IloAlgorithm::Feasible, "Feasible" },
		 { IloAlgorithm::Optimal, "Optimal" },
		 { IloAlgorithm::Infeasible, "Infeasible" },
		 { IloAlgorithm::Unbounded, "Unbounded" },
		 { IloAlgorithm::InfeasibleOrUnbounded, "InfeasibleOrUnbounded" },
		 { IloAlgorithm::Error, "Error" }
};

string to_string(IloAlgorithm::Status status) {
  return status_string[status];
}

// build a list of candidate intervals
// check lower and upper bounds, and Battarra's rule (Prop 4.1)
vector<Interval> get_candidate_intervals(const Instance& I, double ub, double lb=0.0) {
  vector<Interval> ci;
  for(unsigned w=0; w!=I.m; w++) {
    double tc = I.segCost(0,I.n-1); // minimal total cost
    for(unsigned i=0; i!=I.n; i++) {
      double sc = 0.0;
      for(unsigned j=i; j!=I.n; j++) {
	sc += I.t[j][w];
	if (ff(lb)<=ff(sc) && ff(sc)<ff(ub) && ff((tc-sc)/(I.m-1))<ff(ub))
	  ci.push_back({i,j,w});
      }
    }
  }
  return ci;
}

void M2I::addVars(bool ups) {
  for(unsigned w=0; w!=I.m; w++)
    for(unsigned i=0; i!=I.n; i++) {
	x.add(IloNumVar(env, 0.0, 1.0, ILOINT));
      x[x.getSize()-1].setName(("x["+to_string(i)+","+to_string(w)+"]").c_str());
    }

  x.add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT));
  obj.setLinearCoef(x[x.getSize()-1], 1.0);
  if (ups)
    for(unsigned w=0; w!=I.m; w++)
      for(unsigned i=0; i!=I.n; i++) {
	u.add(IloNumVar(env, 0.0, 1.0, ILOFLOAT));
	u[u.getSize()-1].setName(("u["+to_string(i)+","+to_string(w)+"]").c_str());
      }
}

void M2I::addCycle() {
  IloRangeArray cycle(env);
  for(unsigned w=0; w!=I.m; w++) {
    cycle.add(IloRange(env, -IloInfinity, 0.0, ("cycle"+to_string(w)).c_str()));
    cycle[w].setLinearCoef(x[I.m*I.n],-1.0);
    for(unsigned i=0; i!=I.n; i++)
      cycle[w].setLinearCoef(x[w*I.n+i],I.t[i][w]);
  }
  model.add(cycle);
}

void M2I::addAssignment() {
  IloRangeArray assign(env);
  for(unsigned i=0; i!=I.n; i++) {
    assign.add(IloRange(env, 1.0, 1.0, ("assign"+to_string(i)).c_str()));
    for(unsigned w=0; w!=I.m; w++)
      assign[i].setLinearCoef(x[w*I.n+i],1.0);
  }
  model.add(assign);
}

void M2I::addContinuity() {
  IloRangeArray continuity(env);
  int back = -1;

  for(unsigned i=0; i!=I.n; i++)
    for(unsigned h=i+1; h!=I.n; h++)
      for(unsigned j=h+1; j!=I.n; j++)
	for(unsigned w=0; w!=I.m; w++) {
	  continuity.add(IloRange(env, -IloInfinity, 1.0, ("cont"+to_string(i)+","+to_string(h)+","+to_string(j)+"-"+to_string(w)).c_str()));
	  back++;
	  continuity[back].setLinearCoef(x[w*I.n+i],1.0);
	  continuity[back].setLinearCoef(x[w*I.n+j],1.0);
	  continuity[back].setLinearCoef(x[w*I.n+h],-1.0);
	}
  model.add(continuity);
}

void M2I::addContinuity2() {
  IloRangeArray continuity(env);
  int back = -1;
  for(unsigned w=0; w!=I.m; w++) {
    for(unsigned i=0; i!=I.n; i++) {
      continuity.add(IloRange(env, 0.0, IloInfinity, ("cont"+to_string(i)+","+to_string(w)).c_str()));
      back++;
      continuity[back].setLinearCoef(x[w*I.n+i  ],-1.0);
      if (i>0)
	continuity[back].setLinearCoef(x[w*I.n+i-1], 1.0);
      continuity[back].setLinearCoef(u[w*I.n+i],   1.0);
    }
    continuity.add(IloRange(env, -IloInfinity, 1.0, ("cont"+to_string(w)).c_str()));
    back++;
    for(unsigned i=0; i!=I.n; i++)
      continuity[back].setLinearCoef(u[w*I.n+i], 1.0);
  }
  model.add(continuity);
}

void setupSolver(IloCplex solver) {
  IloEnv env = solver.getEnv();
  if (verbosec.count<3) {
    solver.setOut(env.getNullStream());
    solver.setWarning(env.getNullStream());
  }
  solver.setParam(IloCplex::Param::Threads, opt.nthreads);
  solver.setParam(IloCplex::Param::TimeLimit, opt.tlim);
}

IloCplex solveRelaxed(IloModel model, IloNumVarArray x) {
  IloEnv env = model.getEnv();
  IloModel relax(env);
  relax.add(model);
  relax.add(IloConversion(env, x, ILOFLOAT));

  IloCplex solver(relax);
  setupSolver(solver);
  if (false)
    solver.exportModel("relax.lp");
  solver.solve();

  return solver;
}

double M2I::solveRelaxed() {
  IloCplex solver = ::solveRelaxed(model,x);

  if (verbosec.count>1)
    showSolution(solver);

  return solver.getObjValue();
}

void M2I::showSolution(IloCplex solver) {
  IloNumArray value(env);
  solver.getValues(x, value);
  for(unsigned w=0; w!=I.m; w++) {
    vprint(1,"W{:2} ",w);
    for(unsigned i=0; i!=I.n; i++) {
      vprint(1," ");
      if (ff::isnull(value[w*I.n+i]))
	vprint(1,"     .");
      else
	vprint(1,"{:6.2}",value[w*I.n+i]);
    }
    vprint(1,"\n");
  }
}

Solution M2I::getSolution(IloCplex solver) {
  IloNumArray value(env);
  solver.getValues(x, value);
  Solution S(I,I.m);
  vector<int> pi(I.m); iota(pi.begin(),pi.end(),0); // position in worker permutation; handle empty workers by permuting them to the end
  unsigned wi = 0;
  for(unsigned i=0; i!=I.n; i++) {
    for(unsigned w=0; w!=I.m; w++) {
      if (!ff::isnull(value[w*I.n+i])) {
	if (wi == 0 || S.w[wi-1]!=w) {
	  swap(S.w[wi],S.w[pi[w]]);
	  pi[S.w[pi[w]]]=pi[w]; pi[w]=wi;
	  S.l[wi]=1;
	  wi++;
	} else
	  S.l[wi-1]++;
	break;
      }
    }
  }
  S.value = solver.getObjValue();
  return S;
}

void M2I::buildFull(bool ups) {
  buildBase(ups);
  if (!ups)
    addContinuity();
}

void M2I::buildBase(bool ups) {
  addVars(ups);
  addCycle();
  addAssignment();
  if (ups)
    addContinuity2();
}

void M2I::setSolution(IloCplex solver, const Solution& S) {
  IloNumArray v(env,x.getSize());
  for(unsigned i=0, ie=x.getSize(); i!=ie; i++)
    v[i]=0.0;
  int ct = 0;
  for(unsigned i=0; i!=I.m; i++) {
    for(unsigned j=0; j!=S.l[i]; j++)
      v[S.w[i]*I.n+ct+j]=1.0;
    ct += S.l[i];
  }
  solver.addMIPStart(x,v);
  v.end();
}

pair<Solution,IloAlgorithm::Status> M2I::solve(const Solution& S) {
  IloCplex solver(model);
  setupSolver(solver);
  if (false)
    solver.exportModel("m2i.lp");

  if (S.isValid())
    setSolution(solver,S);

  solver.solve();

  Solution T(I);
  IloAlgorithm::Status status = solver.getStatus();
  if (status==IloAlgorithm::Feasible || status==IloAlgorithm::Optimal)
    T = getSolution(solver);
  else
    vprint(1,"M2I: No feasible solution.\n");
  solver.end();
  return {T,status};
}

ILOLAZYCONSTRAINTCALLBACK3(lazyContiguity, IloNumVarArray, x, const Instance&, I, int&, ncuts) {
  IloEnv env = getEnv();
  IloNumArray v(env,x.getSize());
  getValues(v, x);
  for(unsigned w=0; w!=I.m; w++) {
    unsigned i=0, k=I.n-1;
    while (i<k && ff::isnull(v[w*I.n+i])) i++;
    while (k>i && ff::isnull(v[w*I.n+k])) k--;
    for(unsigned j=i+1; j<k; j++) {
      add(x[w*I.n+i]+x[w*I.n+k]-1 <= x[w*I.n+j]).end();
      ncuts++;
    }
  }
  v.end();
  return;
}

// do the `mcuts` deepest cuts
ILOUSERCUTCALLBACK3(cutContiguity, IloNumVarArray, x, const Instance&, I, int&, ncuts) {
  if (!isAfterCutLoop())
    return;
  int lcuts = 0;
  typedef tuple<double,int,int,int,int> Cut;
  priority_queue<Cut,vector<Cut>,greater<Cut>> Q;
  IloEnv env = getEnv();
  IloNumArray v(env,x.getSize());
  getValues(v, x);
  for(unsigned w=0; w!=I.m; w++)
    for(unsigned i=0; i!=I.n; i++)
      for(unsigned h=i+1; h!=I.n; h++)
	for(unsigned j=h+1; j!=I.n; j++) {
	  double e = v[w*I.n+i]+v[w*I.n+j]-1-v[w*I.n+h];
	  if (ff(e)>ff(0.0)) {
	    Q.push({e,w,i,h,j});
	  }
	  if (Q.size()>opt.mcuts)
	    Q.pop();
	}
  while (Q.size()>0) {
    auto c = Q.top(); Q.pop();
    auto [d,w,i,h,j] = c;
    add(x[w*I.n+i]+x[w*I.n+j]-1 <= x[w*I.n+h]).end();
    ncuts++; lcuts++;
  }
  v.end();
}

// do the deepest cut for each interval & worker
ILOUSERCUTCALLBACK3(cutContiguity2, IloNumVarArray, x, const Instance&, I, int&, ncuts) {
  if (!isAfterCutLoop())
    return;

  IloEnv env = getEnv();
  IloNumArray v(env,x.getSize());
  getValues(v, x);
  for(unsigned w=0; w!=I.m; w++)
    for(unsigned i=0; i!=I.n; i++)
      for(unsigned j=i+2; j!=I.n; j++) {
	// deepest separator
	double de = 0.0; int dh = -1;
	for(unsigned h=i+1; h<j; h++) {
	  double e = v[w*I.n+i]+v[w*I.n+j]-1-v[w*I.n+h];
	  if (ff(e)>ff(de)) {
	    de = e; dh = h;
	  }
	}
	if (dh>-1) {
	  ncuts++;
	  add(x[w*I.n+i]+x[w*I.n+j]-1 <= x[w*I.n+dh]).end();
	}
      }
  v.end();
}

// do the central cut for each (pseudo-)convex region
// idea here: for a fixed worker, segment the values into convex regions
// for each such region, the deepest cut is between the borders and the minimum in the interior
ILOUSERCUTCALLBACK3(cutContiguity3, IloNumVarArray, x, const Instance&, I, int&, ncuts) {
  if (!isAfterCutLoop())
    return;

  IloEnv env = getEnv();
  IloNumArray v(env,x.getSize());
  getValues(v, x);
  for(unsigned w=0; w!=I.m; w++) {

    unsigned i=0;
    while (ff(v[w*I.n+0])==ff(v[w*I.n+1]) && i!=I.n)
      i++;
    if (i==I.n)
      continue;
    int up=-1, dn=-1;
    bool onup;
    if (ff(v[w*I.n+i])<ff(v[w*I.n+i+1])) {
      up = i; onup = true;
    } else {
      dn = i; onup = false;
    }
    i++;
    while (i+1!=I.n) {
      if (ff(v[w*I.n+i])<ff(v[w*I.n+i+1]) && !onup) {
	up = i; onup = true;
      }
      if (ff(v[w*I.n+i])>ff(v[w*I.n+i+1]) && onup) {
	if (dn>-1) {
	  // here (dn,up,i) is deepest
	  ncuts++;
	  double e = v[w*I.n+dn]+v[w*I.n+i]-1-v[w*I.n+up];
	  if (ff(e)>ff(0.0)) {
	    add(x[w*I.n+dn]+x[w*I.n+i]-1 <= x[w*I.n+up]).end();
	  }
	}
	dn = i; onup = false;
      }
      i++;
    }
    if (onup && dn>-1) {
      // here (dn,up,i) is deepest
      ncuts++;
      double e = v[w*I.n+dn]+v[w*I.n+i]-1-v[w*I.n+up];
      if (ff(e)>ff(0.0)) {
	add(x[w*I.n+dn]+x[w*I.n+i]-1 <= x[w*I.n+up]).end();
      }
    }
  }
  v.end();
}

// other strategies
// per worker: only the widest, narrowest, deepest cut?

pair<Solution,IloAlgorithm::Status> M2I::solveBC(const Solution& S) {
  IloCplex solver(model);
  setupSolver(solver);
  if (false)
    solver.exportModel("m2i-bc.lp");

  int ncuts = 0;
  solver.use(lazyContiguity(env,x,I,ncuts));
  solver.use(cutContiguity3(env,x,I,ncuts));

  if (S.isValid())
    setSolution(solver,S);

  solver.solve();
  vprint(1,"Added {} cuts.\n",ncuts);

  Solution T(I);
  IloAlgorithm::Status status = solver.getStatus();
  if (status==IloAlgorithm::Feasible || status==IloAlgorithm::Optimal)
    T = getSolution(solver);
  else
    vprint(1,"M2I: No feasible solution.\n");
  solver.end();
  return {T,status};
}

double M3I::addVars(double ub, bool addC) {
  const bool only_maximal = false;
  if (only_maximal)
    vprint(2,"M3I: considers only maximal segments with costs at most {}.\n",ub);
  
  int nvars = 0;
  double maxcost = 0.0;
  for(unsigned w=0; w!=I.m; w++)
    for(unsigned i=0; i!=I.n; i++)
      for(unsigned j=i; j!=I.n; j++) {
	double seg = I.segCost(i,j,w);
	bool canGrowLeft = i>0 && ff(I.segCost(i-1,j,w))<=ff(ub);
	bool canGrowRight = j+1<I.n && ff(I.segCost(i,j+1,w))<=ff(ub);
	bool infeas = ff(seg)>ff(ub) || (only_maximal && (canGrowLeft || canGrowRight));
	if (!infeas) {
	  maxcost=max(maxcost,seg);
	  nvars++;
	  idx[{i,j,w}]=y.getSize();
	  y.add(IloNumVar(env, 0.0, 1.0, ILOINT));
	  y[y.getSize()-1].setName(("y["+to_string(i)+","+to_string(j)+","+to_string(w)+"]").c_str());
	  if (!addC)
	    obj.setLinearCoef(y[y.getSize()-1], seg);
	}
      }
  // cycle time
  if (addC) {
    idx_c = y.getSize();
    y.add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT));
    y[idx_c].setName("c");
    obj.setLinearCoef(y[idx_c], 1.0);
  }

  vprint(2,"M3I: added {} out of {} possible variables.\n",nvars,I.m*I.n*(I.n+1)/2);
  return maxcost;
}

void M3I::addCycle() {
  IloExprArray cexpr(env,I.m);
  for(unsigned w=0; w!=I.m; w++) {
    cexpr[w]=IloExpr(env);
    cexpr[w] -= y[idx_c];
  }

  for(auto ijw : idx) {
    auto [i,j,w] = ijw.first;
    cexpr[w] += I.segCost(i,j,w)*y[ijw.second];
  }

  for(unsigned w=0; w!=I.m; w++) {
    cycle.add(IloRange(env, -IloInfinity, cexpr[w], 0.0, ("cycle"+to_string(w)).c_str()));
    cexpr[w].end();
  }

  model.add(cycle);
}

void M3I::addAssignment() {
  IloExprArray assexpr(env,I.n);
  for(unsigned k=0; k!=I.n; k++)
    assexpr[k]=IloExpr(env);
  for(auto ijw : idx)
    for(unsigned k=get<0>(ijw.first), ke=get<1>(ijw.first)+1; k!=ke; k++)
      assexpr[k] += y[ijw.second];
  for(unsigned k=0; k!=I.n; k++) {
    assign.add(IloRange(env, 1.0, assexpr[k], IloInfinity, ("assign"+to_string(k)).c_str()));
    assexpr[k].end();
  }
  model.add(assign);
}

void M3I::addWorker() {
  IloExprArray wexpr(env,I.m);
  for(unsigned w=0; w!=I.m; w++)
    wexpr[w]=IloExpr(env);

  for(auto ijw : idx) {
    int w = get<2>(ijw.first);
    wexpr[w] += y[ijw.second];
  }

  for(unsigned w=0; w!=I.m; w++) {
    worker.add(IloRange(env, -IloInfinity, wexpr[w], 1.0, ("worker"+to_string(w)).c_str())); // we relax here
    wexpr[w].end();
  }
  model.add(worker);
}

void M3I::buildFull(const Solution& S, bool feasibility) {
  double ub = S.isValid()? S.value : inf_time;
  addVars(ub,!feasibility);
  if (!feasibility)
    addCycle();
  addAssignment();
  addWorker();
}

double M3I::buildFull(double ub, bool feasibility) {
  double maxcost = addVars(ub,!feasibility);
  if (!feasibility)
    addCycle();
  addAssignment();
  addWorker();
  return maxcost;
}

void M3I::buildCore(const Solution& S) {
  // take vars from solution
  int ct = 0;
  for(unsigned i=0; i!=I.m; i++)
    if (S.l[i]>0) {
      idx[{ct,ct+S.l[i]-1,S.w[i]}]=y.getSize();
      y.add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT));
      y[y.getSize()-1].setName(fmt::format("y[{},{},{}]",ct,ct+S.l[i]-1,S.w[i]).c_str());
      ct += S.l[i];
    }
  idx_c = y.getSize();
  y.add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT)); // cycle time
  y[idx_c].setName("c");
  obj.setLinearCoef(y[idx_c], 1.0);

  // build initial model with these
  addCycle();
  addAssignment();
  addWorker();
}

void M3I::addColumn(unsigned i, unsigned j, unsigned w) {
  // add var
  int idx_new = y.getSize();
  idx[{i,j,w}]=idx_new;

  IloNumColumn col = obj(0.0);
  col += cycle[w](I.segCost(i,j,w));
  for(unsigned k=i; k!=j+1; k++)
    col += assign[k](1.0);
  col += worker[w](1.0);
  
  y.add(IloNumVar(col, 0.0, IloInfinity, ILOFLOAT));
  y[y.getSize()-1].setName(fmt::format("y[{},{},{}]",i,j,w).c_str());
}

void M3I::addColumn(Interval iv) {
  addColumn(get<0>(iv),get<1>(iv),get<2>(iv));
}

double M3I::solveRelaxedCG(double ub, double lb) {
  IloCplex solver(model);
  setupSolver(solver);

  // for checking multiple columns
  unordered_set<tuple<int,int,int>> colpool;

  // build a list of candidate intervals
  vector<Interval> ci = get_candidate_intervals(I,ub,lb);
  vprint(1,"There are {} candidate intervals (out of {}) for bounds [{},{}].\n",ci.size(),I.m*I.n*(I.n+1)/2,lb,ub);

  sort(ci.begin(),ci.end(),[](Interval a, Interval b) { return (get<1>(a)-get<0>(a))>(get<1>(b)-get<0>(b)); });

  int iter = 0;
  double ret = 0.0;
  auto next_ci = ci.begin();

  do {
    vprint(1,"{:5.1f} {:4} {:6.3f} {:5}",::S.elapsed(),iter,ret,y.getSize()); fflush(stdout);
    solver.solve();
    ret = solver.getObjValue();
    vprint(1,".\n"); fflush(stdout);

    // extract duals
    IloNumArray dc(env,cycle.getSize()), da(env,assign.getSize()), dw(env,worker.getSize());
    solver.getDuals(dc,cycle);
    solver.getDuals(da,assign);
    solver.getDuals(dw,worker);

    vector<double> dap(I.n+1,0.0); // partial sum: dap[i]=da[0]+...+da[i-1]
    for(unsigned k=0; k!=I.n; k++)
      dap[k+1]=dap[k]+da[k];

    // search for columns of small reduced cost
    vector<pair<Interval,double>> cols;
    unsigned checked = 0;
    while (cols.size()<opt.ncols && checked<ci.size()) {
      auto const& [i,j,w] = *next_ci++;
      if (next_ci==ci.end())
	next_ci=ci.begin();
      checked++;

      double rc = dap[j+1]-dap[i]+dw[w]+I.segCost(i,j,w)*dc[w];
      if (ff(rc)>ff(0.0))
	cols.push_back({{i,j,w},rc});
    }
    if (cols.size()==0)
      break;
    vprint(1,"({}/{}) ",checked,ci.size());
    for(auto const& [iv,rc] : cols) {
      if (colpool.find(iv)!=colpool.end()) {
	vprint(1,fg(fmt::color::red),"Column {} seen twice, reduced cost {}.\n",iv,-rc);
	throw;
      } else
	colpool.insert(iv);
      addColumn(iv);
    }
    iter++;
  } while (true || iter<50);  /* new columns added */
  solver.end();
  return ret;
}

pair<Solution,IloAlgorithm::Status> M3I::solvedTenseCG() {
  model.add(IloConversion(env, y, ILOINT));
  IloCplex solver(model);
  setupSolver(solver);
  solver.solve();

  Solution T(I);
  IloAlgorithm::Status status = solver.getStatus();
  if (status==IloAlgorithm::Feasible || status==IloAlgorithm::Optimal)
    T = getSolution(solver);
  else
    vprint(1,"M3I/tenseCG: No feasible solution.\n");
  solver.end();
  return {T,status};
}

void M3I::setSolution(IloCplex solver, const Solution& S) {
  IloNumArray v(env,y.getSize());
  for(unsigned i=0, ie=y.getSize(); i!=ie; i++)
    v[i]=0.0;
  unsigned ct = 0;
  for(unsigned i=0; i!=I.m; i++) {
    if (S.l[i]>0) {
      auto e = idx.find({ct,ct+S.l[i]-1,S.w[i]});
      if (e != idx.end()) {
	v[e->second]=1.0;
      } else {
      // find an interval that contains it
	bool found = false;
	for(auto ijw : idx) {
	  auto [i_,j,w] = ijw.first;
	  if (i_<=ct && ct+S.l[i]-1<=j && w==S.w[i]) {
	    v[ijw.second]=1.0;
	    found = true;
	    break;
	  }
	}
	if (!found) {
	  fmt::print(fg(fmt::color::yellow),"M3I::setSolution: Warning: no interval contains [{},{},{}].\n",ct,ct+S.l[i]-1,S.w[i]);
	}
      }
      ct += S.l[i];
    }
  }
  solver.addMIPStart(y,v);
  v.end();
}

Solution M3I::getSolution(IloCplex solver) {
  IloNumArray value(env);
  solver.getValues(y, value);
  Solution S(I,I.m);
  vector<int> pi(I.m); iota(pi.begin(),pi.end(),0); // position in worker permutation; handle empty workers by permuting them to the end

  vector<tuple<int,int,int>> interval;

  for (auto const& [iv, val] : idx)
    if (!ff::isnull(value[val]))
      interval.push_back(iv);

  sort(interval.begin(),interval.end());

  unsigned ibegin = 0;
  Time C = 0.0;
  for(unsigned i=0, ie=interval.size(); i!=ie; i++) {
    int w = get<2>(interval[i]);
    swap(S.w[i],S.w[pi[w]]);
    pi[S.w[pi[w]]]=pi[w]; pi[w]=i;
    unsigned iend = get<1>(interval[i]);
    S.l[i]=(ibegin<=iend)?(iend-ibegin+1):0;
    C=max(C,I.segCost(ibegin,ibegin+S.l[i]-1,w));
    ibegin=get<1>(interval[i])+1;
  }
  S.value = C;
  return S;
}

double M3I::solveRelaxed() {
  IloCplex solver = ::solveRelaxed(model,y);

  if (verbosec.count>3)
    showSolution(solver);

  double value = solver.getObjValue();
  solver.end();

  return value;
}

bool M3I::solveRelaxedFeasibility() {
  IloCplex solver = ::solveRelaxed(model,y);

  if (verbosec.count>3)
    showSolution(solver);

  bool feasible = solver.getStatus()==IloAlgorithm::Feasible || solver.getStatus()==IloAlgorithm::Optimal;

  // check for integrality
  if (feasible) {
    IloNumArray value(env);
    solver.getValues(y, value);
    bool integral = true;
    for(auto ijw : idx) {
      if (ff(value[ijw.second])!=ff(0.0) && ff(value[ijw.second])!=ff(1.0)) {
	integral = false;
	fmt::print("Witness {}.\n",value[ijw.second]);
	showSolution(solver);
	solver.exportModel("notintegral.lp");
	break;
      }
    }
    fmt::print("Solution is {} integral.\n",(integral?"":"not"));
  }
  solver.end();
  return feasible;
}

pair<Solution,IloAlgorithm::Status> M3I::solve(const Solution& S) {
  IloCplex solver(model);
  setupSolver(solver);
  if (false)
    solver.exportModel("m3i.lp");

  if (S.isValid())
    setSolution(solver,S);

  solver.solve();
  Solution T(I);
  IloAlgorithm::Status status = solver.getStatus();
  if (status==IloAlgorithm::Feasible || status==IloAlgorithm::Optimal)
    T = getSolution(solver);
  else
    vprint(1,"M3I: No feasible solution.\n");
  solver.end();
  return {T,status};
}

void M3I::showSolution(IloCplex solver) {
  IloNumArray value(env);
  solver.getValues(y, value);
  for(unsigned w=0; w!=I.m; w++) {
    bool first = true;
    for(unsigned i=0; i!=I.n; i++)
      for(unsigned j=i; j!=I.n; j++) {
	if (idx.find({i,j,w}) == idx.end())
	  continue;
	if (!ff::isnull(value[idx[{i,j,w}]])) {
	  if (first)
	    vprint(1,"W{:2} ",w);
	  else
	    vprint(1,"    ");

	  first=false;
	  if (I.n<50) {
	    for(unsigned k=0; k!=I.n; k++) {
	      if (k<i || k>j)
		vprint(1,"      .");
	      else
		vprint(1," {:6.2}",value[idx[{i,j,w}]]);
	    }
	  } else
	    vprint(1,"[Skipping {} values.]\n",j-i+1);
	  vprint(1,"\n");
	}
      }
  }
}

void Corf::addVars() {
  const unsigned o = I.s0.size();
  for(unsigned w=0; w!=I.m; w++)
    for(unsigned i=0; i!=o; i++) {
      p.add(IloNumVar(env, 0.0, 1.0, ILOINT));

      q.add(IloNumVar(env, 0.0, 1.0, ILOINT));

      r.add(IloNumVar(env, 0.0, 1.0-ff::EPS, ILOFLOAT));
    }
  // cycle time
  c = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
  obj.setLinearCoef(c, 1.0);
}

void Corf::addCycle() {
  const unsigned o = I.s0.size();
  for(unsigned w=0; w!=I.m; w++) {
    IloExpr expr(env);
    for(unsigned i=0; i!=o; i++)
      expr += I.B*I.t[i*I.B][w]*(p[i+w*o]+r[i+w*o]);
    expr -= c;
    cycle.add(IloRange(env, -IloInfinity, expr, 0.0, fmt::format("cycle{}",w).c_str()));
    expr.end();
  }
  model.add(cycle);
}

void Corf::addAssignment() {
  const unsigned o = I.s0.size();
  IloExprArray oexpr(env,o);
  for(unsigned i=0; i!=o; i++) {
    oexpr[i]=IloExpr(env);
    for(unsigned w=0; w!=I.m; w++)
      oexpr[i] += p[i+w*o]+r[i+w*o];
  }
  for(unsigned i=0; i!=o; i++) {
    assign.add(IloRange(env, 1.0, oexpr[i], 1.0, fmt::format("assign{}",i).c_str()));
    oexpr[i].end();
  }
  model.add(assign);
}

void Corf::addContinuity() {
  const unsigned o = I.s0.size();
  for(unsigned i=0; i!=o; i++)
    for(unsigned k=i+2; k<o; k++)
      for(unsigned j=i+1; j<k; j++)
	for(unsigned w=0; w!=I.m; w++) {
	  IloExpr expr(env);
	  expr += q[i+w*o]+q[k+w*o]+p[i+w*o]+p[k+w*o]-p[j+w*o];
	  continuity.add(IloRange(env, -IloInfinity, expr, 1.0, fmt::format("cont{},{},{},{}",i,j,k,w).c_str()));
	  expr.end();
	}
  model.add(continuity);
}

void Corf::addPacking() {
  const unsigned o = I.s0.size();
  for(unsigned i=0; i!=o; i++)
    for(unsigned w=0; w!=I.m; w++) {
      IloExpr expr(env);
      expr += p[i+w*o]+q[i+w*o];
      packing.add(IloRange(env, -IloInfinity, expr, 1.0, fmt::format("pack{},{}",i,w).c_str()));
      expr.end();
    }
  model.add(packing);
}

void Corf::addLinking() {
  const unsigned o = I.s0.size();
  for(unsigned i=0; i!=o; i++)
    for(unsigned w=0; w!=I.m; w++) {
      IloExpr expr(env);
      expr += r[i+w*o]-q[i+w*o];
      linking.add(IloRange(env, -IloInfinity, expr, 0.0, fmt::format("pack{},{}",i,w).c_str()));
      expr.end();
    }
  model.add(linking);
}

void Corf::addFringe() {
  const unsigned o = I.s0.size();
  for(unsigned w=0; w!=I.m; w++) {
      IloExpr expr(env);
      for(unsigned i=0; i!=o; i++)
	expr += q[i+w*o];
      fringe.add(IloRange(env, -IloInfinity, expr, 2.0, fmt::format("fringe{}",w).c_str()));
      expr.end();
  }
  model.add(fringe);
}

void Corf::build() {
  addVars();
  addCycle();
  addAssignment();
  addContinuity();
  addPacking();
  addLinking();
  addFringe();
}

pair<double,IloAlgorithm::Status> Corf::solve() {
  IloCplex solver(model);
  setupSolver(solver);
  solver.solve();

  double C = 0.0;
  IloAlgorithm::Status status = solver.getStatus();
  if (status==IloAlgorithm::Feasible || status==IloAlgorithm::Optimal)
    C = solver.getObjValue();
  else
    vprint(1,"Corf: No feasible solution.\n");
  solver.end();
  return {C,status};
}
