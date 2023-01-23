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
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
using namespace std;

#include "util.hpp"
#include "options.hpp"
#include "statistics.hpp"
#include "random.hpp"
#include "instance.hpp"
#include "ip.hpp"
#include "hs.hpp"
#include "beamSearch.hpp"
#include "reducedIP.hpp"
#include "lb.hpp"
#include "dp.hpp"
#include "sadp.hpp"
#include "ff.hpp"
#include "solver.hpp"
#include "cbfs.hpp"
#include "cbfsdfs.hpp"
#include "memory.hpp"
#include "solution.hpp"
using namespace wap;

using namespace boost;

//// globals
string iname;

void test_heuristics(const Instance& I, unsigned seed) {
  vprint(1,fg(fmt::color::green),"Test heuristics\n");
  const int nsamples = 100;
  Solution S = sample(I,nsamples);
  vprint(1,"Best of {} samples ", nsamples);
  polish(I,S);
  vprint(1,"\n{}\n",S.to_string());

  Solution Sms = multi_start(I,nsamples,[](const Instance& I,vector<unsigned>& w) { shuffle(w.begin(),w.end(),wap::rng);},all_nb);
  vprint(1,"Best of {} multi-start samples ", nsamples);
  vprint(1,"{}\n",Sms.to_string());

  // uniform segments
  optS seg(I);
  vector<unsigned> l(I.m);
  for(auto i=0u; i!=I.m; i++)
    l[i]=(i+1)*I.n/I.m - i*I.n/I.m;
  Solution Su = seg.match(l);
  vprint(1,"Uniform time ");
  polish(I,Su);
  vprint(1,"\n{}\n",Su.to_string());

  // segments for mean time
  vector<Time> m(I.n,0.0);
  for(auto i=0u; i!=I.n; i++) {
    for(auto w=0u; w!=I.m; w++)
      m[i]+=I.t[i][w];
    m[i]/=I.m;
  }
  Solution Sm = seg.match(optT(I,m));
  vprint(1,"Mean time ");
  polish(I,Sm);
  vprint(1,"\n{}\n",Sm.to_string());

  // segments from standard time
  Solution Ss = seg.match(optT(I,I.s));
  vprint(1,"Standard time ");
  polish(I,Ss);
  vprint(1,"\n{}\n",Ss.to_string());
  
  // segments from minimum time
  fill(m.begin(),m.end(),inf_time);
  for(auto i=0u; i!=I.n; i++)
    for(auto w=0u; w!=I.m; w++)
      m[i]=min(m[i],I.t[i][w]);
  Solution Sm_ = seg.match(optT(I,m));
  vprint(1,"Minimum times ");
  polish(I,Sm_);
  vprint(1,"\n{}\n",Sm_.to_string());

  fmt::print("BASEHEUR {} {} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f}\n", iname, seed, S.value, Sms.value, Su.value, Sm.value, Ss.value, Sm_.value);
}

Solution initial_segments(const Instance& I) {
  vtprint(1,fg(fmt::color::green),"Run initial heuristics.\n");
  Solution S(I);

  // uniform segments
  optS seg(I);
  vector<unsigned> l(I.m);
  for(auto i=0u; i!=I.m; i++)
    l[i]=(i+1)*I.n/I.m - i*I.n/I.m;
  Solution Su = seg.match(l);
  vprint(1,"Uniform time ");
  polish(I,Su);
  vprint(1,"\n{}\n",Su.to_string());
  if (Su.value<S.value)
    S = Su;

  // segments for mean time
  vector<Time> m(I.n,0.0);
  for(auto i=0u; i!=I.n; i++) {
    for(auto w=0u; w!=I.m; w++)
      m[i]+=I.t[i][w];
    m[i]/=I.m;
  }
  Solution Sm = seg.match(optT(I,m));
  vprint(1,"Mean time ");
  polish(I,Sm);
  vprint(1,"\n{}\n",Sm.to_string());
  if (Sm.value<S.value)
    S = Sm;

  // segments from standard time
  Solution Ss = seg.match(optT(I,I.s));
  vprint(1,"Standard time ");
  polish(I,Ss);
  vprint(1,"\n{}\n",Ss.to_string());
  if (Ss.value<S.value)
    S = Ss;

  // segments from minimum time
  fill(m.begin(),m.end(),inf_time);
  for(auto i=0u; i!=I.n; i++)
    for(auto w=0u; w!=I.m; w++)
      m[i]=min(m[i],I.t[i][w]);
  Solution Sm_ = seg.match(optT(I,m));
  vprint(1,"Minimum times ");
  polish(I,Sm_);
  vprint(1,"\n{}\n",Sm_.to_string());
  if (Sm_.value<S.value)
    S = Sm_;

  return S;
}

void read_instance(Instance& I) {
  vtprint(1,"Read input {}.\n",iname);
  ifstream in(opt.instance);
  if (!in.is_open()) {
    cerr << "Cannot open " << iname << "." << endl;
    exit(1);
  }
  if (opt.format=="guess") {
    if (fs::path(opt.instance).extension()==".txt")
      opt.format="std";
    else
      opt.format="alwabp";
  }
  vprint(2,"Input format \"{}\".\n",opt.format);
  try {
    if (opt.format=="std")
      I.readSTD(in);
    else if (opt.format=="alwabp")
      I.readALWABP(in);
    else
      throw "Unknown format.";
  } catch (const char *s) {
    cerr << "Failed to read instance: " << s << endl;
    exit(1);
  } catch (const string s) {
    cerr << "Failed to read instance: " << s << endl;
    exit(1);
  } catch (...) {
    cerr << "Failed to read instance." << endl;
    exit(1);
  }
  in.close();

  vtprint(1,fg(fmt::color::green),"Instance read.\n");
  vprint(1,"Instance has {} tasks, and {} workers. Batch size is {}.\n",I.n,I.m,I.B);
}

int main(int argc, char *argv[]) {
  // (0) process commandline
  S.start = Clock::now();
  S.iter = 0;

  po::variables_map vm;
  if (!process_options(argc,argv,opt,vm))
    return 1;

  // (1) random seed and basics
  unsigned seed = setupRandom(vm["seed"].as<unsigned>());

  iname = string(fs::path(opt.instance).parent_path().filename());
  if (iname!="")
    iname = iname+"/";
  iname = iname+string(fs::path(opt.instance).stem());

  // (2) read the instance
  Instance I;
  read_instance(I);

  Time lb1 = minSum(I);
  Time lb2 = optRep(I);
  fmt::print("LB {} {} {:.6} {:.6} ", iname, seed, lb1, lb2);
  DP::optnRepDP O(I);
  Time lb3 = optnRep(I,O,opt.lb3w,opt.showlb3);
  if (!opt.showlb3)
    fmt::print("{:.6}",lb3);
  fmt::print("\n");

  if (opt.onlylb)
    return 0;

  Time lb = max(lb1,max(lb2,lb3));

  // (3) compute a seed solution
  Solution S(I);
  if (opt.initial=="sample")
    S = sample(I,opt.nsamples);
  else if (opt.initial=="rnd")
    S = sample(I,1);
  else if (opt.initial=="seg")
    S = initial_segments(I);
  else if (opt.initial=="test") {
    test_heuristics(I,seed); return 0;
  } else if (opt.initial=="ils") {
    ILS ils(I);
    ils.setLB(lb);
    ils.run();
    fmt::print("ILS {} {} {:5.2f} {} {:5.2f} {:.5f} {}\n", iname, seed, ::S.elapsed(), ils.iter, ils.tbest, ils.B.value, ff(ils.B.value)==ff(lb)?"opt":"feas");
    return 0;
  } else {
    vprint(1,fg(fmt::color::green),"Invalid initial heuristic \"{}\".\n",opt.initial);
    return 1;
  }

  Time ini = S.value;

  timer theur;
  Time heur;

  // (4) run main heuristic
  bool optimal=ff(S.value)==ff(lb);

  map<string,string> hname { { "bs", "Beam search" }, { "cbfs", "CBFS" } , { "cmb", "Combined heuristic" } };
  const bool hknown = hname.find(opt.heuristic)!=hname.end();
  Solver_Status status = FEASIBLE;
  if (!hknown || optimal)
    vtprint(1,fg(fmt::color::green),"Skipped heuristic solution ({}).\n", optimal?"optimal":"unknown");
  else {
    vtprint(1,fg(fmt::color::green),"{}\n",hname[opt.heuristic]);
    if (opt.heuristic=="bs")
      S=beamSearch(I,S,opt.bsmin,opt.bsmax);
    else if (opt.heuristic=="cbfs")
      tie(S,status)=CBFS(I,S,O);
    else if (opt.heuristic=="cmb") {
      S=combined_solver(I,lb,S,O);
      return 0;
    }
    polish(I,S);
    vprint(1,"{}\n",S.to_string());
  }

  heur = S.value;
  double htime = theur.elapsed_secs();

  // (5) run exact solver
  optimal=ff(S.value)==ff(lb) || status==OPTIMAL;
  timer texact;

  map<string,string> ename { { "ip", "Reduced IP" }, { "dp", "Dynamic Programming" } , { "sadpF", "Successive DP approximations" }, { "sadpFS", "Successive tail DP approximations" }, { "sadpO", "Successive DP approximations, optimization" }, {"cbfsF", "cbfs, feasibility"}, {"cbfsO", "cbfs, optimization"}, {"cbfsDfsF", "cbfs/dfs hybrid, feasibility"}, {"cbfsDfsO", "cbfs/dfs hybrid, optimization"} };
  const bool eknown = ename.find(opt.exact)!=ename.end();
  if (!eknown || optimal)
    vtprint(1,fg(fmt::color::green),"Skipped exact solution ({}).\n", optimal?"optimal":"unknown");
  else {
    optimal=true; // by default methods solve to optimality, only IP solver respects time limit
    vtprint(1,fg(fmt::color::green),"{}.\n",ename[opt.exact]);
    if (opt.exact=="ip")
      optimal = reducedIP(I,S) ;
    else if (opt.exact=="dp")
      S=dynamicProgramming(I,S);
    else if (opt.exact=="sadpF") {
      optimal = false;
      tie(S,status)=sadpF::sadpFeasibility(I,S);
    } else if (opt.exact=="sadpFS") {
      optimal = false;
      tie(S,status)=sadpF::sadpFeasibility(I,S,true);
    } else if (opt.exact=="sadpO")
      tie(S,optimal)=sadpO::sadpOptimization(I,S);
    else if (opt.exact=="cbfsF")
      tie(S,optimal)=cbfsF::cbfsFeasibility(I,S);
    else if (opt.exact=="cbfsO")
      tie(S,optimal)=cbfsO::cbfsOptimization(I,S);
    else if (opt.exact=="cbfsDfsF")
      tie(S,optimal)=cbfsDfsF::cbfsDfsFeasibility(I,S);
    else if (opt.exact=="cbfsDfsO")
      tie(S,optimal)=cbfsDfsO::cbfsDfsOptimization(I,S);
    vprint(1,"{}\n",S.to_string());
  }

  Time exact = S.value;
  double etime = texact.elapsed_secs();
  if (optimal)
    status = OPTIMAL;

  fmt::print("SUMM {} {} {:5.2f} {:.5f} {:.5f} {:.2f} {:.5f} {} {:.2f}\n", iname, seed, ::S.elapsed(), ini, heur, htime, exact, status_string.at(status), etime);
  fmt::print("STAT {} {} {} {} {}\n", iname, seed, usedVMproc(), stat.expansions, opt.sadpmemory);

  // (6) run models
  Solution e2a(I);
  if (opt.m2i) {
    if (I.n*I.n*I.n*I.m<100000 || opt.m2iup) {
      timer tm2i;
      M2I m2i(I);
      m2i.buildFull(opt.m2iup);
      double lb4 = 0;//m2i.solveRelaxed();
      IloAlgorithm::Status status = IloAlgorithm::Unknown;
      tie(e2a,status) = m2i.solve(S);
      vprint(1,"Exact solution M2I {}\n",e2a.to_string());
      fmt::print("MODEL m2i{} {} {} {} {} {}\n",opt.m2iup?"'":"", iname, lb4, e2a.value, to_string(status), tm2i.elapsed_secs());
    } else
      vprint(1,fg(fmt::color::red),"Skipping M2I model, we have about {} continuity constraints.\n",I.n*I.n*I.n*I.m);
  }

  Solution e2b(I);
  if (opt.m2ibc) {
    timer tm2ibc;
    M2I m2ibc(I);
    m2ibc.buildBase(opt.m2iup);
    IloAlgorithm::Status status;
    tie(e2b,status) = m2ibc.solveBC(S);
    vprint(1,"Exact solution via B&C {}\n",e2b.to_string());
    fmt::print("MODEL m2i{} {} {} {} {}\n",opt.m2iup?"'":"", iname, e2b.value, to_string(status), tm2ibc.elapsed_secs());
  } else
    vprint(1,fg(fmt::color::red),"Skipping M2I B&C model.\n");

  Solution e3a(I);
  if (opt.m3i) {
    timer tm3i;
    M3I m3i(I);
    m3i.buildFull(S,false);
    double lb5 = m3i.solveRelaxed();
    vprint(1,"Lower bound M3I {}\n",lb5);
    IloAlgorithm::Status status;
    tie(e3a,status)  = m3i.solve(S);
    fmt::print("MODEL m3i{} {} {} {} {} {}\n",opt.m2iup?"'":"", iname, lb5, e3a.value, to_string(status), tm3i.elapsed_secs());
    vprint(1,"Exact solution M3I {}\n",e3a.to_string());
  } else
    vprint(1,fg(fmt::color::red),"Skipping full M3I model.\n");

  if (opt.m3if) {
    timer tm3if;
    double mt = I.getMinTime();
    double lb = mt, ub = S.value;
    vprint(1,"Initial interval [{:.6f},{:.6f}], minimum time {:.5f}\n",lb,ub,mt);
    while (ff(lb)<ff(ub)) {
      double m = (lb+ub)/2;
      M3I m3i(I);
      double maxcost = m3i.buildFull(m,true);
      bool feasible = m3i.solveRelaxedFeasibility();
      vprint(1,"Trial {:.6f} is {}feasible, max cost {:.6f}",m,(feasible?"  ":"in"),maxcost);
      if (feasible) {
	if (ff(ub)<=ff(maxcost)) {
	  vprint(1,", is optimal within tolerance.\n");
	  break;
	}
	ub = min(maxcost,m);
      } else
	lb = max(maxcost,m);
      vprint(1," next interval [{:.6f},{:.6f}].\n",lb,ub);
    }
    M3I m3i(I);
    double maxcost = m3i.buildFull(ub,true);
    fmt::print("Post-search integer solution for ub {:.6f} maxcost {:.6f}.\n",ub,maxcost);
    IloAlgorithm::Status status;
    Solution is(I);    
    tie(is,status) = m3i.solve(Solution(I));
    if (is.isValid())
      vprint(1,"Integer solution {}.\n", is.value);
    else
      vprint(1,"No valid integer solution.\n");
    fmt::print("MODEL m3if {} {} {} {}\n", iname, is.value, to_string(status), tm3if.elapsed_secs());
  } else
    vprint(1,fg(fmt::color::red),"Skipping M3I repeated feasibility model.\n");
  
  if (opt.m3icg) {
    Solution e3b(I);
    timer tm3icg;
    double lb6 = 0.0;
    M3I m3i(I);
    m3i.buildCore(S);
    lb6 = m3i.solveRelaxedCG(S.value,0.0);
    IloAlgorithm::Status status = IloAlgorithm::Unknown;
    //tie(e3b,status) = m3i.solvedTenseCG();
    fmt::print("MODEL m3icg/relaxed {} {} {} {} {}\n", iname, lb6, e3b.value, to_string(status), tm3icg.elapsed_secs());
    vprint(1,"Heuristic solution M3I {}\n",e3b.to_string());
  } else
    vprint(1,fg(fmt::color::red),"Skipping M3I CG model.\n");

  if (opt.corf) {
    timer tcorf;
    double lb7 = 0.0;
    Corf corf(I);
    corf.build();
    IloAlgorithm::Status status = IloAlgorithm::Unknown;
    tie(lb7,status) = corf.solve();
    fmt::print("MODEL CORF {} {} {} {}\n", iname, lb7, to_string(status), tcorf.elapsed_secs());
    vprint(1,"Lower bound via CORF {}\n",lb7);
  } else
    vprint(1,fg(fmt::color::red),"Skipping CORF model.\n");
}
