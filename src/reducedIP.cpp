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
#include <boost/multi_array.hpp>

#include <ilcplex/ilocplex.h>

#include "reducedIP.hpp"
#include "util.hpp"
#include "hs.hpp"
#include "options.hpp"

using namespace std;
using namespace boost;

struct ipVariable {
  int worker;
  int st; //starting task. It is included, or equivalent to [ in the range
  int et; //ending task. It is not included, or equivalent to  ) in the range
  int t;  //sum of task times [st,et)
};

bool tryIP(const Instance &I, Solution &S, unsigned long currentC) {
  int tt;
  long unsigned maxValue;

  //(1) generation of variables
  vector<int> usefulTasks(I.n+1,false);
  vector<bool> stillActive(I.m+1,true);
  vector<ipVariable> ipVariables;
  usefulTasks[0]=true;
  for(int t=0; t!=int(I.n); t++) {
    if (usefulTasks[t]) {
      for(unsigned w=0;w<I.m;w++) {
	if (!stillActive[w])
	  continue;
	ipVariable v;
	v.worker=w;
	v.st=t;
	//find the first task that does not fit with currentC
	maxValue=I.d[t-1][w]+currentC;
	for(tt=t; tt!=int(I.n) && I.d[tt][w]<maxValue; tt++) {}
	v.et=tt;
	v.t=I.d[tt-1][w]-I.d[t-1][w];
	assert(v.t<currentC);
	ipVariables.push_back(v);
	if (tt<int(I.n)) usefulTasks[tt]=true;
	else stillActive[w]=false;
      }
    }
  }
  vprint(2,"{} useful tasks, {} variables, I.n {}, I.m {}.\n",count(usefulTasks.begin(),usefulTasks.end(),true),ipVariables.size(),I.n,I.m);
  if (verbosec.count>3) {
    for(unsigned t=0;t<ipVariables.size();t++)
      printf("%d %d %d %d\n",ipVariables[t].worker,ipVariables[t].st,ipVariables[t].et,ipVariables[t].t);
  }

  //(2) cplex initialization
  int remain_time = opt.tlim-Duration(Clock::now()-::S.start).count();
  if (remain_time<5)
    return false;
  vprint(2,"Solving with time limit {}s.\n",remain_time);
  
  IloEnv env;
  IloModel model(env);
  IloCplex solver(model);
  solver.setParam(IloCplex::Param::Threads, 1);
  solver.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 1E-9);
  solver.setParam(IloCplex::Param::TimeLimit, remain_time);
  solver.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1);	//remove heuristic from b&b
  solver.setParam(IloCplex::Param::MIP::Limits::Solutions,1); //any feasible solution does the trick for us
  if (verbosec.count<3) {
    solver.setOut(env.getNullStream());
    solver.setWarning(env.getNullStream());
  }

  //(3) cplex model
  IloRangeArray constraintTasks(env,I.n,1.0,IloInfinity);
  IloRangeArray constraintWorkers(env,I.m,1.0,1.0);
  model.add(constraintTasks);
  model.add(constraintWorkers);
  IloNumVarArray variables(env);
  IloObjective obj=IloAdd(model,IloMinimize(env)); //any objective function serves, may be worth testing alternatives
  for(unsigned i=0;i<ipVariables.size();i++) {
    IloNumColumn col = obj(ipVariables[i].t);
    for(int t=ipVariables[i].st;t<ipVariables[i].et;t++) col += constraintTasks[t](1.0);
    col += constraintWorkers[ipVariables[i].worker](1.0);
    variables.add(IloNumVar(col, 0.0, 1.0, ILOINT));
    col.end();
  }
  solver.solve();
  vprint(1,"cplex reports: {}\n",solver.getStatus());

  if (solver.getStatus()==IloAlgorithm::Feasible || solver.getStatus()==IloAlgorithm::Optimal) {
    //(4) cplex solution
    vector<unsigned> vWorkers(I.m);
    unsigned numWorkers=0;
    for(unsigned t=0;t<ipVariables.size();t++) {
      if (solver.getValue(variables[t])>0.9) {
        vWorkers[numWorkers]=ipVariables[t].worker;
        numWorkers++;
        assert(numWorkers<=I.m);
      }
    }
    assert(numWorkers==I.m);
    S=optW(I,vWorkers);
    env.end();
    return false;
  }
  bool status = solver.getStatus()==IloAlgorithm::Infeasible;
  env.end();
  return status;
}

// reducedIP: repeatedly solve IP on maximum segments
bool reducedIP(const Instance& I, Solution &S) {
  unsigned int currentC=S.value*precision; //current cycle time (bad practice, sorry. I'll make it more general if we move these data to the instance definition)

  vprint(1,"Starting reduced IP currentC {} numTasks: {} numWorkers: {}\n",time_double(currentC),I.n,I.m);

  while(Duration(Clock::now()-::S.start).count()<opt.tlim) { //keep trying while solution improves
    vprint(1,"Try to improve {}.\n",time_double(currentC));
    bool opt = tryIP(I,S,currentC);
    if (S.value*precision<currentC) currentC=S.value*precision;
    else return opt;
  }
  return false;
}
