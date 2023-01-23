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
#include <assert.h>
#include <queue>
#include <vector>
#include <map>

#include "util.hpp"
#include "hs.hpp"
#include "sadp.hpp"
#include "ap.hpp"
#include "options.hpp"
#include "memory.hpp"

// helper for controlling resources
struct Control {
  unsigned memMB, memLimit;
  double timeUsed;

  Control(unsigned memLimit) : memLimit(memLimit) {
    update();
  }

  void update() {
    memMB = usedVMproc()/(1l<<20);
    timeUsed = ::S.elapsed();
  }
  
  bool terminate() {
    update();
    return oom() || oot();
  }

  bool oom() {
    return memMB>memLimit;
  }

  bool oot() {
    return timeUsed>opt.tlim;
  }
};

using namespace std;

// There are two alternative successive approximation approaches
// A feasibility based approximation which solves feasibility versions of the problem for decreasing values of c
// An optimization based approach which solves an optimization version of the problem
// They differ on the state space representation and approach to solve the recurrence.
// Also, the optimization based approach may
// be adapted to provide lower bounds which can then be used in other methods
// Only the feasibility version has been implemented

namespace sadpF { //feasibility version

  /******************************************************
   * Overall scheme:
   *
   * The method considers multiple feasibility problems (much like the dp approach) until no improving solution is possible.
   * Instead of solving a complete dp formulation, it solves multiple simplified formulations.
   *
   * Each simplified formulation corresponds to a DP model where states online a subset of workers explicitly and aggregates the remaining workers.
   * After solving each simplified formulation, one of the following cases apply:
   *
   * 1 No feasible solution exists. Then, there is no solution with the desired cycle time and we can stop.
   * 2 A feasible solutions exists. The solutions may be feasible for the complete problem (after dissagregation each worker is used once).
   * 2.1 If the solution uses each worker once, we found a feasible solution for the problem and we can return it.
   * 2.2 Otherwise, one or more workers have been used two or more times. Select one of these workers and add it to the relaxation until 2.1 applies.
   *
   *******************************************************/

#define SOL_INFEASIBLE 0	 //No solution is possible with the cycle time
#define SOL_FEASIBLE 1	 //The relaxation is feasible

map<vector<int>,int> stateSpace;   // for each state (a vector) maximum number of tasks
vector<int> inUse;    //inUse[i] -> which one is the i-th explicit worker
vector<bool> explicitWorkers; // explicitWorkers[i] -> is worker i an explicit worker?

vector<int> topState; // for tail relaxation

inline int old_maxTask(const Instance& I, unsigned w, int lastTaskPerformed, long int currentC) {
  long unsigned int maxValue=I.d[lastTaskPerformed][w]+currentC;
  int t;
  for(t=lastTaskPerformed+1;t<int(I.n)&&(I.d[t][w]<maxValue);t++) {}
  return t-1;
}

// solves state space relaxation with workers tracked according to `inUse` and `explicitWorkers`
int internalFeasibilityRelaxation(const Instance& I, long int currentC) {
  int tp,t;
  unsigned numExplicitWorkers=inUse.size();
  priority_queue<vector<int>> oldStates;
  vector<vector<int>> newStates;
  stateSpace.clear();
  assert(stateSpace.empty());

  // store initial state
  vector<int> state(numExplicitWorkers+1,0); // {0,1} for explicit workers, one extra value for "used remaining workers"
  oldStates.push(state);
  stateSpace[state]=-1;

  for(unsigned depth=0; depth!=I.m; ++depth) {
    assert(newStates.empty());
    while (!oldStates.empty()) {
      state=oldStates.top();
      int lastTaskPerformed=stateSpace[state];
      oldStates.pop();
      for(unsigned w=0; w!=numExplicitWorkers; ++w) { // check all workers considered within the relaxation
	if (state[w]==1) continue;
	state[w]=1;
	t = maxTask(I,inUse[w],lastTaskPerformed,currentC);
	
	auto it=stateSpace.find(state);
	if (it==stateSpace.end()) { // state it does not exist, save
	  stateSpace[state]=t;
	  newStates.push_back(state);
	} else if (it->second<t) { // update
	  it->second=t;
	}
	state[w]=0;
      }
      
      if (state[numExplicitWorkers]+numExplicitWorkers<I.m) { // an implicit worker is possible?
        state[numExplicitWorkers]++;
        t=-1;
        for(auto w=0u; w!=I.m; ++w) {
	  if (explicitWorkers[w]) continue; // skip explicit workers
	  tp = maxTask(I,w,lastTaskPerformed,currentC);
	  if (t<tp) t=tp;
        }
        assert(t>=0);
        auto it=stateSpace.find(state);
        if (it==stateSpace.end()) { // state does not exist, save it
          stateSpace[state]=t;
          newStates.push_back(state);
        } else if (it->second<t) { // state improves, update it
	  it->second=t;
        }
        //state[numExplicitWorkers]--; no need to restore since we loop & pop
      }
    } // !oldStates.empty()

    // copy newStates to oldStates
    for(auto& state : newStates) oldStates.push(state);
    newStates.clear();
  } // max depth = num workers

  while (!oldStates.empty()) {
    state=oldStates.top();
    if (stateSpace[state]<int(I.n-1))
      return SOL_INFEASIBLE;
    return SOL_FEASIBLE;
  }
  assert(false);
  return SOL_INFEASIBLE;
}

// solves state space relaxation with tracking the tail of the last q workers
int tailinternalFeasibilityRelaxation(const Instance& I, long int currentC, unsigned q) {
  priority_queue<vector<int>> oldStates;
  vector<vector<int>> newStates;
  stateSpace.clear();

  //store initial state
  vector<int> state(q+1,-1); // last q workers, plus number of workers
  state[q]=0;
  oldStates.push(state);
  stateSpace[state]=-1;

  for(auto depth=0u; depth!=I.m; ++depth) {
    assert(newStates.empty());
    while (!oldStates.empty()) {
      state=oldStates.top();
      int lastTaskPerformed=stateSpace[state];
      oldStates.pop();
       
      // check all workers
      for(auto w=0u; w!=I.m; ++w) {
	if (find(state.begin(),state.begin()+q,int(w))!=state.begin()+q)
	  continue;

	auto t = maxTask(I,w,lastTaskPerformed,currentC);
	
	vector<int> sstate(state);
	rotate(sstate.begin(), sstate.begin()+1, sstate.begin()+q);
	sstate[q-1]=w;
	sstate[q]=state[q]+1;

	auto it=stateSpace.find(sstate);
	if (it==stateSpace.end()) { //state does not exist, save it
	  stateSpace[sstate]=t;
	  newStates.push_back(sstate);
	} else {
	  if (it->second<t) {
	    it->second=t;
	  } else {
	  }
	}
      }
    } // !oldStates.empty()

    // copy newStates to oldStates
    for(auto& state : newStates) oldStates.push(state);
    newStates.clear();
  } // max depth = num workers
  
  assert(!oldStates.empty());
  while (!oldStates.empty()) {
    topState=oldStates.top();
    if (!(stateSpace[topState]<int(I.n-1))) return SOL_FEASIBLE;
    oldStates.pop();
  }
  return SOL_INFEASIBLE;
}

// returns currentC, if it does not improve, it returns original currentC
// Due to the use of a Pascal-translated function (vectors and arrays go from 1 to n) there are a few +1 and -1 in the code
int bottleneckHeuristic(const Instance& I, Solution& S, const vector<int>& cumTasks, long int currentC) {
  //1. Populate array with costs (I.d[numTasks][worker]
  //first row does not need substraction
  vprint(4,"checking bottleneckHeuristic ");
  for(auto i=0u;i<I.m;i++) { //rows, worker count
    ap::c[i+1][1]=I.d[cumTasks[0]][i];
    for(auto j=1u;j<I.m;j++) { //columns, station count
      ap::c[i+1][j+1]=I.d[cumTasks[j]][i]-I.d[cumTasks[j-1]][i];
    }
  }
  //2. Solve bottleneck AP
  auto newC=ap::bottleneckAP(I.m);
  vector<unsigned int> vWorkers(I.m,I.m+1);  //I.m+1 means unassigned (new solution assignment)
  for(auto i=0u;i<I.m;i++) {
    vWorkers[i]=ap::y[i+1]-1;
  }
  Solution S2=optW(I,vWorkers);
  //3. If improves, change solution
  if (S2.value<S.value) {
    //update solution
    //S=optW(I,vWorkers);
    S=S2;
    vprint(4,"current: {} new: {}\n",time_double(currentC),S.value);
    currentC=S.value*precision;
  } else {
    vprint(4,"does not improve: fixed value {} for current value {}.\n",time_double(newC),time_double(currentC));
  }
  return currentC;
}

//updates state space relaxation
//  We have to ``backtrack'' to obtain the solution. We have aggregated workers and solutions where multiple workers may do the trick
//  specially when there are dummy workstations. This procedure follows the ``decode'' routine found in beamSearch.cpp
int checkInternalFeasibilitySolution(const Instance& I, Solution& S, long int* currentC) {
  const unsigned numExplicitWorkers=inUse.size();

  //fmt::print("checkInt: state space {}, {} explicit workers.\n", stateSpace.size(), numExplicitWorkers);

  // representation of extracted solution: vWorkers[i], i=0,...,m-1 is the ith worker; it does tasks cumTasks[i-1]+1,...,cumTasks[i], where cumTasks[-1]=0 by convention, cumTasks[m-1]=n-1 by definition.
  vector<int> cumTasks(I.m,0);
  const unsigned int FREE = I.m+1; //I.m+1 means unassigned
  vector<unsigned int> vWorkers(I.m,FREE);
  
  vector<int> state(numExplicitWorkers+1,1);
  state[numExplicitWorkers]=I.m-numExplicitWorkers;
  int t=I.n-1;
  for(int depth=I.m-1;depth>0;depth--) {
    //there are two types of workers, explicit and aggregated (implicit)
    for(unsigned w=0;w<numExplicitWorkers;w++) {
      if (state[w]==1) { //explicit worker
        state[w]=0;
        auto it=stateSpace.find(state);
        if (it!=stateSpace.end()) { //if it exists, it may be the worker in depth position (otherwise, it is not the case)
          if ((I.d[t][inUse[w]]-I.d[it->second][inUse[w]])<(*currentC)) { //requires less than current c, possible to construct a path
            vWorkers[depth]=inUse[w];
            cumTasks[depth]=t;
            t=it->second;
            break;
          }
        }
        state[w]=1; //it is not this explicit worker
      }
    }
    
    if (vWorkers[depth]==FREE) {
      assert(state[numExplicitWorkers]>0);
      //try aggregated (implicit) workers
      state[numExplicitWorkers]--;
      auto it=stateSpace.find(state);
      assert(it!=stateSpace.end());
      for(auto w=0u; w<I.m; w++) {
	if (explicitWorkers[w]) continue; // skip explicit workers
	if ((I.d[t][w]-I.d[it->second][w])<(*currentC)) {
	  vWorkers[depth]=w;
	  cumTasks[depth]=t;
	  t=it->second;
	  break;
	}
      }
      //if(vWorkers[depth]>I.m) state[numExplicitWorkers]++; unneeded
    }
    assert(vWorkers[depth]<I.m);
  }
  
  //last (actually first) worker
  for(unsigned w=0u;w<numExplicitWorkers;w++) {
    if (state[w]==1) { //explicit worker
      state[w]=0;
      vWorkers[0]=inUse[w];
      cumTasks[0]=t;
      break;
    }
  }
  
  if (vWorkers[0]>I.m) { //aggregated worker
    assert(state[numExplicitWorkers]==1);
    for(unsigned w=0;w<I.m;w++) {
      if(!explicitWorkers[w]) {
        if(I.d[t][w]<(*currentC)) {
          vWorkers[0]=w;
          cumTasks[0]=t;
          break;
        }
      }
    }
  }
  assert(vWorkers[0]<I.m);

  //assignment of tasks to workers available
  //build solution

  vector<unsigned> whist(I.m,0);
  vector<double> ctime(I.m,0.0);
  unsigned whmax = 0;
  ctime[0]=time_double(I.d[cumTasks[0]][vWorkers[0]]);
  double crit = ctime[0];
  for(auto i=0u; i!=I.m; ++i) {
    whist[vWorkers[i]]++;
    if (i>0) {
      ctime[i]=time_double(I.d[cumTasks[i]][vWorkers[i]]-I.d[cumTasks[i-1]][vWorkers[i]]);
      crit=max(crit,ctime[i]);
    }
  }
  vprint(3,"[");
  for(auto i=0u; i!=I.m; ++i) {
    if (ctime[i]==crit)
      vprint(3,fg(fmt::color::red),"{:2}{}",vWorkers[i],i+1==I.m?"":" ");
    else if (whist[vWorkers[i]]>1)
      vprint(3,fg(fmt::color::yellow),"{:2}{}",vWorkers[i],i+1==I.m?"":" ");
    else
      vprint(3,"{:2}{}",vWorkers[i],i+1==I.m?"":" ");
  }
  vprint(3,"]: ");
  for(auto i=0u; i!=I.m; ++i)
    if (whist[i]>1) {
      vprint(3,"{}x{} ",i,whist[i]);
      if (whist[whmax]<whist[i])
	whmax=i;
    }
  vprint(3,"; {}\n",crit);
  for(auto i=0u; i!=I.m; ++i) {
    if (ctime[i]==crit)
      vprint(3,fg(fmt::color::red),"{:.6f} ",ctime[i]);
    else if (whist[vWorkers[i]]>1)
      vprint(3,fg(fmt::color::yellow),"{:.6f} ",ctime[i]);
    else
      vprint(3,"{:.6f} ",ctime[i]);
  }
  vprint(3,"\n");

  // test: highest frequency
  if (false && whist[whmax]>1) {
    explicitWorkers[whmax]=true;
    inUse.push_back(whmax);
    *currentC=bottleneckHeuristic(I,S,cumTasks,*currentC);
    vprint(3,"Tracking worker {} now, c={:.6f}\n", whmax, time_double(*currentC));
    return SOL_INFEASIBLE;
  }

  // test: cheapest of highest frequency
  if (whist[whmax]>1) {
    int cmin = -1;
    for(auto i=0u; i!=I.m; ++i)
      if (whist[vWorkers[i]]==whist[whmax])
	if (cmin==-1 || ctime[i]>ctime[cmin])
	  cmin=i;
  
    explicitWorkers[vWorkers[cmin]]=true;
    inUse.push_back(vWorkers[cmin]);
    long int currentC_old = *currentC;
    *currentC=bottleneckHeuristic(I,S,cumTasks,*currentC);
    vprint(3,"Tracking worker {} now, c={:.6f}\n", vWorkers[cmin], time_double(*currentC));
    if (*currentC<currentC_old)
      return SOL_FEASIBLE;
    else
      return SOL_INFEASIBLE;
  }
  
  vector<bool> inSolution(I.m,false);
  for(unsigned i=0;i<I.m;i++) {
    //add to inSolution
    if (!inSolution[vWorkers[i]]) inSolution[vWorkers[i]]=true;
    else { //vWorkers[i] should be added to the formulation
      assert(!explicitWorkers[vWorkers[i]]);
      explicitWorkers[vWorkers[i]]=true;
      inUse.push_back(vWorkers[i]);
      vprint(3,"Tracking worker {} now.\n", vWorkers[i]);
      *currentC=bottleneckHeuristic(I,S,cumTasks,*currentC);
      return SOL_INFEASIBLE;
    }
    if (cumTasks[i]==int(I.n-1)) { //solution found
      for(unsigned ii=i+1;ii<I.m;ii++) { //if i is smaller than I.m, random robots suffice to obtain a solution. We are forced to change for valid robots
        vWorkers[ii]=(-1);
        for(unsigned a=0u;(a<I.m)&&(vWorkers[ii]==(-1));a++) {
          if (!inSolution[a]) {
            vWorkers[ii]=a;
            inSolution[a]=true;
          }
        }
        assert(vWorkers[ii]>=0);
      }
      //valid assignment of workers to stations
      S=optW(I,vWorkers);
      return SOL_FEASIBLE;
    }
  }
  //never here
  assert(false);
  exit(0);
}

//updates state space relaxation
//  We have to ``backtrack'' to obtain the solution. We have aggregated workers and solutions where multiple workers may do the trick
//  specially when there are dummy workstations. This procedure follows the ``decode'' routine found in beamSearch.cpp
int tailcheckInternalFeasibilitySolution(const Instance& I, Solution& S, long int* currentC, unsigned q) {
  vector<int> cumTasks(I.m,0);
  const unsigned int FREE = I.m+1; //I.m+1 means unassigned
  vector<unsigned int> vWorkers(I.m,FREE);
  
  vector<int> state(topState);
  int t=I.n-1;
  
  for(int depth=I.m-1; depth>=0; depth--) {
    unsigned lw = state[q-1];
    assert(depth+1==state[q]);
    vWorkers[depth]=lw;
    cumTasks[depth]=t;
    
    // now we need to find the predecessor state
    for(auto w=-1; w!=int(I.m); ++w) {
      vector<int> pstate(state);
      rotate(pstate.begin(), pstate.begin()+q-1, pstate.begin()+q);
      pstate[0]=w;
      pstate[q]=state[q]-1;

      auto it=stateSpace.find(pstate);
      if (it!=stateSpace.end()) { //if it exists, it may be the worker in depth position (otherwise, it is not the case)
	if (w==-1 || (I.d[t][lw]-I.d[it->second][lw])<*currentC) { //requires less than current c, possible to construct a path
	  state = pstate;
	  t=it->second;
	  break;
	}
      }
    }
  }
  
  //assignment of tasks to workers available
  //build solution

  vector<unsigned> whist(I.m,0);
  vprint(3,"[");
  for(auto i=0u; i!=I.m; ++i) {
    whist[vWorkers[i]]++;
    vprint(3,"{:2}{}",vWorkers[i],i+1==I.m?"":" ");
  }
  vprint(3,"]: ");
  for(auto i=0u; i!=I.m; ++i)
    if (whist[i]>1)
      vprint(3,"{}x{} ",i,whist[i]);
  vprint(3,"\n");
  
  vector<bool> inSolution(I.m,false);
  for(auto i=0u; i!=I.m; ++i) {
    //add to inSolution
    if (inSolution[vWorkers[i]]) {
      *currentC=bottleneckHeuristic(I,S,cumTasks,*currentC);
      return SOL_INFEASIBLE;
    }
    inSolution[vWorkers[i]]=true;
    
    if (cumTasks[i]==int(I.n-1)) { //solution found
      for(unsigned ii=i+1;ii<I.m;ii++) { //if i is smaller than I.m, random robots suffice to obtain a solution. We are forced to change for valid robots
        vWorkers[ii]=FREE;
        for(unsigned a=0u; a<I.m && vWorkers[ii]==FREE; a++) {
          if (!inSolution[a]) {
            vWorkers[ii]=a;
            inSolution[a]=true;
          }
        }
        assert(vWorkers[ii]>=0);
      }
      //valid assignment of workers to stations
      S=optW(I,vWorkers);
      return SOL_FEASIBLE;
    }
  }
  assert(false);
  return SOL_INFEASIBLE;
}

// solve feasibility problem via successive relaxation
  Solution externalFeasibilityRelaxation(const Instance& I, Solution& S, long int currentC, Control& con) {
  inUse.clear();
  for(auto i=0u;i<I.m;i++) explicitWorkers[i]=false;
  
  while (!con.terminate()) {
    vprint(3,"Current relaxation has {} (of {}) uniquely identified workers (memory {} MB).\n", inUse.size(), I.m, usedVMproc()/(1l<<20));
    int retV=internalFeasibilityRelaxation(I,currentC); //solve relaxed dp model
    if (retV==SOL_INFEASIBLE) break;
    
    retV=checkInternalFeasibilitySolution(I,S,&currentC); //add worker or update solution
    if (retV==SOL_FEASIBLE) break;
  }
  vprint(2,fg(fmt::color::blue),"Last relaxation has {} (of {}) uniquely identified workers.\n", inUse.size(), I.m);
  return S;
}

Solution tailexternalFeasibilityRelaxation(const Instance& I, Solution& S, long int currentC) {
  unsigned q = 1;
  int retV;
  
  while (true) {
    vprint(3,"Current relaxation has a tail of {} (of {}) workers.\n", q, I.m);
    retV=tailinternalFeasibilityRelaxation(I,currentC,q); //solve relaxed dp model
    if (retV==SOL_INFEASIBLE) break;
    
    retV=tailcheckInternalFeasibilitySolution(I,S,&currentC,q); //add worker or update solution
    if (retV==SOL_FEASIBLE) break;
    
    q++;
  }
  vprint(2,"Last relaxation has a tail of {} (of {}) workers.\n", q, I.m);
  vprint(2,"It is {}feasible.\n",retV==SOL_FEASIBLE?"":"in");
  return S;
}

// successive approximative DP, feasibility version
//   call successive approximation with successively smaller cycle time
std::pair<Solution,Solver_Status> sadpFeasibility(const Instance& I, Solution &S, bool tail) {
  long unsigned int currentC=S.value*precision; //current cycle time (bad practice).
  Solver_Status status = FEASIBLE;
  
  ap::initAP(I.m);
  assert(explicitWorkers.size()==0);
  for(auto i=0u;i<I.m;i++) explicitWorkers.push_back(false);
  vprint(2,"Starting successive approximations initialC: {} numTasks: {} numWorkers: {}\n", time_double(currentC), I.n, I.m);
  Control con(opt.sadpFmemory);
  while (!con.terminate()) {
    vprint(2,"currentBest: {} with {} explicit workers, memory {} MB.\n", time_double(currentC), inUse.size(), con.memMB);
    if (tail)
      S=tailexternalFeasibilityRelaxation(I,S,currentC);
    else
      S=externalFeasibilityRelaxation(I,S,currentC,con);
    if (S.value*precision<currentC)
      currentC=S.value*precision;
    else break;
  }
  if (con.oom())
    status = MEMLIMIT;
  else if (con.oot())
    status = TIMELIMIT;
  else
    status = OPTIMAL;
  vprint(2, "sadp found optimal objective {} (memory {} MB).\n", time_double(currentC), usedVMproc()/(1l<<20));
  ap::freeAP(I.m);
  return {S,status};
}

} // namespace sadpF

namespace sadpO { //optimization version

int debugInf=0;
char letra;


unsigned numExplicitWorkers;
int numImplicitWorkers;
vector<bool> explicitWorkers; // explicitWorkers[i] -> is worker i an explicit worker?
vector<unsigned int> mapExplicitWorkers;
unsigned int* states = 0;
unsigned int* trace = 0;
int id;
int taskStep; //value that needs to be substracted from a state to reach a state with one less task
long int currentC;

int constantTime;
int constantImplicitWorkers;
State st,stOld,stNew; //contains an explicit representation of the state
State stTest;

void decodeState(const Instance& I,State& st,int id) {
  st.tasks=id/constantTime;
  id=id-st.tasks*constantTime;
  st.numImplicitWorkers=id/constantImplicitWorkers;
  id=id-st.numImplicitWorkers*constantImplicitWorkers;

  for(auto i=int(numExplicitWorkers)-1;i>=0;i--) {
    if(id%2==1) {
      st.explicitWorkers[i]=1;
    } else {
      st.explicitWorkers[i]=0;
    }
    id=id/2 ;
  }
  return;
}

int encodeState(const Instance& I,State& st) {
  vprint(4,"Encoding state: {} {} .",st.tasks, st.numImplicitWorkers);
  for(auto i=0u; i!=numExplicitWorkers; ++i) vprint(4,"{} ",st.explicitWorkers[i]);
  vprint(4,"\nconstantes: {} {}\n", constantTime, constantImplicitWorkers);
  int id=st.tasks*constantTime;
  id += st.numImplicitWorkers*constantImplicitWorkers;
  int idWorkers=0;
  for(auto i=0u;i<numExplicitWorkers;i++) {
    idWorkers = idWorkers*2;
    if (st.explicitWorkers[i]==1)
      idWorkers++;
  }
  id += idWorkers;
  vprint(4," equates to : {} tasks, numImplicitWorkers: {}", st.tasks, st.numImplicitWorkers);
  if (numExplicitWorkers>0) {
    vprint(4," and explicit workers ");
    for(auto i=0u; i!=numExplicitWorkers; ++i) vprint(4,"{} ", st.explicitWorkers[i]);
  }
  vprint(4, " . id: {}\n", id);
  return id;
}

//this a bad practice. It repeats the same code than sadpF with a change in cumTasks (here goes from 1 to n)
int bottleneckHeuristic(const Instance& I, Solution& S, const vector<int>& cumTasks, long int currentC) {
  //1. Populate array with costs (I.d[numTasks][worker]
  //first row does not need substraction
  vprint(4,"\nchecking bottleneckHeuristic\n");
  for(auto i=0u;i<I.m;i++) { //rows, worker count
    ap::c[i+1][1]=I.d[cumTasks[1]-1][i];
    for(auto j=1u;j<I.m;j++) { //columns, station count
      ap::c[i+1][j+1]=I.d[cumTasks[j+1]-1][i]-I.d[cumTasks[j]-1][i];
    }
  }
  //2. Solve bottleneck AP
  auto newC=ap::bottleneckAP(I.m);
  vector<unsigned int> vWorkers(I.m,I.m+1);  //I.m+1 means unassigned (new solution assignment)
  //check if solution is located in x or y
  for(auto i=0u;i<I.m;i++) {
    vWorkers[i]=ap::y[i+1]-1;
  }
  Solution S2=optW(I,vWorkers);
  //3. If improves, change solution
  if(S2.value<S.value) {
    //update solution
    S=optW(I,vWorkers);
    vprint(4,"current: {} new: {}\n",time_double(currentC),S.value);
    currentC=S.value*precision;
  }
  else {
    vprint(4,"does not improve: fixed value {} for current value {}.\n",time_double(newC),time_double(currentC));
  }
  return currentC;
}

unsigned int optimalityRelaxation(const Instance& I,Solution& S,int stateSize) {
  timer optRelax;
  
  unsigned int returnValue=I.m;
  //1. recurrence
  const unsigned int inf=RAND_MAX/10; //still thinking about it.
  constantImplicitWorkers=pow(2,numExplicitWorkers);
  constantTime=constantImplicitWorkers*(numImplicitWorkers+1);
  vprint(1,"Optimality relaxation with {} states and {} explicit workers (", stateSize, numExplicitWorkers);
  for(auto i=0u; i!=numExplicitWorkers; ++i)
    vprint(1,"{} ", mapExplicitWorkers[i]+1);
  vprint(1,") time {}.\n",::S.elapsed());
  for(auto i=0u;i<numExplicitWorkers;i++) st.explicitWorkers[i]=1;
  st.tasks=I.n;
  st.numImplicitWorkers=numImplicitWorkers;

  for(auto i=0;i<stateSize-1;i++) states[i]=inf; //initialize cycle times
  states[stateSize-1]=0;
  //dp recurrence
  for(auto id=stateSize-1;id>0;id--) { //loop through all states
    decodeState(I,stTest,id);
    if(states[id]==inf) {
      continue; //avoids some calculations for unreachable states
    }
    decodeState(I,st,id); //the state contains an integer that codes the complete information about tasks performed, workers available, etc.
    if(st.numImplicitWorkers>0) { //loop through one less implicit worker
      st.numImplicitWorkers--;
      auto i=encodeState(I,st);
      for(auto t=st.tasks-1;t>=0;t--) {
        i=i-constantTime;
        assert(i>=0);
        auto best=inf;
        for(auto w=0u;w<I.m;w++) { 
          if(explicitWorkers[w]==false) {
            auto accTime=I.d[st.tasks-1][w];
            if(t>0) accTime -= I.d[t-1][w];
            if(best>accTime) best=accTime;
          }
        }
        if(states[i]>max(best,states[id])) {
          states[i]=max(best,states[id]);
          trace[i]=id;
        }
        if(best>currentC) break; //avoid unnecessary work (TO DO. The code allows for solution identical to optimum, maybe change easy to detect by setting a condition to test states[0]=inf
      }
      st.numImplicitWorkers++;
    }
    for(auto w=0u; w<numExplicitWorkers; w++) { //loop through explicit workers and check their use
      if(st.explicitWorkers[w]==1) {
        st.explicitWorkers[w]=0;
        auto i=encodeState(I,st);
        for(auto t=st.tasks-1;t>=0;t--) {
          i=i-constantTime;
          assert(i>=0);
          unsigned int accTime=I.d[st.tasks-1][mapExplicitWorkers[w]];
          if(t>0) accTime -= I.d[t-1][mapExplicitWorkers[w]];
          if(states[i]>max(accTime,states[id])) {
            states[i]=max(accTime,states[id]);
            trace[i]=id;
          }
          if(accTime>currentC) break; //avoid unnecessary work (TO DO. The code allows for solution identical to optimum, maybe change easy to detect by setting a condition to test states[0]=inf
        }
        st.explicitWorkers[w]=1;
      }
    }
  }
  //state  0 contains optimal cycle time
  assert((double)(states[0])/(double)(precision)-1.0/(double)(precision)<S.value); //a lower bound cannot be better than the incumbent
  //3. rebuild solution
  vector<int> cumTasks(I.m+1,0);
  vector<unsigned int> workersInUse(I.m+1,I.m);
  vector<int> inSolution(I.m,0);
  int depth=0;
  auto oldState=0;
  auto currentState=0;
  while(currentState!=(stateSize-1))   {
    decodeState(I,stOld,currentState);
    oldState=currentState;
    currentState=trace[oldState];
    decodeState(I,stNew,currentState);
    depth++;
    cumTasks[depth]=currentState/constantTime;
    //check if explicit worker
    bool isExplicitWorker=false;
    for(auto i=0u;i<numExplicitWorkers;i++) {
      if(stOld.explicitWorkers[i]!=stNew.explicitWorkers[i]) {
        //assert((stOld.explicitWorkers[i]==0)&&(stNew.explicitWorkers[i]==1));
        workersInUse[depth]=mapExplicitWorkers[i]; //TO DO
        isExplicitWorker=true;
        break;
      }
    }
    //else, first worker that satisfies cycle time constraint (but annotate all workers that may perform the task)
    if(isExplicitWorker==false) {
      for(auto i=0u;i<I.m;i++) {
        if((explicitWorkers[i]==false)&&((I.d[stNew.tasks-1][i]-I.d[stOld.tasks-1][i])<=states[oldState])) {
          if(workersInUse[depth]==I.m) workersInUse[depth]=i;
          inSolution[i]++;
          //break;
        }
      }
    }
    assert(workersInUse[depth]!=I.m);
  }
  vprint(4," State: {} ... \t{}\n",currentState,cumTasks[depth],states[currentState]);
  // Select the one with minimum sum of task times among those that repeat the most
  returnValue=I.m;
  for(auto i=0u;i<I.m;i++) {
    if(inSolution[i]>1) {
      if(returnValue==I.m) {
        returnValue=i;
      } else {
        if(inSolution[i]>inSolution[returnValue]) {
          returnValue=i;
        } else {
          if((inSolution[i]==inSolution[returnValue])&&(I.T[I.n][i]<I.T[I.n][returnValue])) {
            returnValue=i;
          }
        }
      }
    }
  }
      
  if (returnValue!=I.m) { //Apply heuristic 
    //we need to call bottleneckHeuristic using cumTasks (goes from 0 to m-1)
    bottleneckHeuristic(I,S,cumTasks,currentC);
    currentC=S.value*precision;
    if(states[0]==S.value*precision) returnValue=I.m; //also means optimality
    vprint(3,"  Bound {:.4f}, solution {:.4f}, time {} returnValue {}\n", time_double(states[0]), S.value, optRelax.elapsed_secs(),returnValue);
  } else { //evaluate solution
    vector<unsigned int> vWorkers(I.m,I.m);
    for(auto i=0u;i<I.m;i++) vWorkers[i]=workersInUse[i+1];
    S=optW(I,vWorkers);
    currentC=S.value*precision;
    vprint(3,"New feasible solution with cost: {}\n",currentC);
  }
  return returnValue;
}

Solution initOptimalityRelaxation(const Instance& I, Solution& S) {
  assert(explicitWorkers.size()==0);
  for(auto i=0u;i<I.m;i++) explicitWorkers.push_back(false);
  assert(mapExplicitWorkers.size()==0);
  numExplicitWorkers=0;
  numImplicitWorkers=I.m;
  //resize dp tables
  int stateSize=(I.n+1)*(numImplicitWorkers+1);
  states=(unsigned int *) malloc(stateSize*sizeof(unsigned int));
  trace=(unsigned int *) malloc(stateSize*sizeof(unsigned int));
  debugInf=0;
  while(1) {
    //solve relaxation and return worker to add to relaxation
    auto whichWorker=optimalityRelaxation(I,S,stateSize);
    if(whichWorker>=I.m) break; //means feasible
    else {
      if(::S.elapsed()>opt.tlim) break;
      //update relaxation and resize dp tables
      assert(explicitWorkers[whichWorker]==false);
      if(((I.n+1)*(numImplicitWorkers)*pow(2,numExplicitWorkers+1))>opt.sadpmemory) break; //no more memory
      explicitWorkers[whichWorker]=true; //add explicit worker
      mapExplicitWorkers.push_back(whichWorker);
      numExplicitWorkers++;
      numImplicitWorkers--;
      stateSize=(I.n+1)*(numImplicitWorkers+1)*pow(2,numExplicitWorkers);
      states=(unsigned int *) realloc(states,stateSize*sizeof(unsigned int)); //realloc states
      trace=(unsigned int *) realloc(trace,stateSize*sizeof(unsigned int)); //realloc trace
    }
  }
  return(S);
}

pair<Solution,bool> sadpOptimization(const Instance& I, Solution& S) {
  //memory allocation
  bool optimal=true;
  ap::initAP(I.m);
  st.explicitWorkers=(unsigned *)malloc(I.m*sizeof(unsigned));
  stOld.explicitWorkers=(unsigned *)malloc(I.m*sizeof(unsigned));
  stNew.explicitWorkers=(unsigned *)malloc(I.m*sizeof(unsigned));
  stTest.explicitWorkers=(unsigned *)malloc(I.m*sizeof(unsigned));
  //initialize
  currentC=S.value*precision;
  vprint(2,"Starting successive approximations initialC: {} numTasks: {} numWorkers: {}\n",S.value,I.n,I.m);
  //call relaxation
  S=initOptimalityRelaxation(I,S);
  vprint(2,"sadp found optimal objective {}\n",S.value);
  //free memory
  if(abs(S.value-time_double(states[0]))<0.0001) {
    vprint(2,"sadp found optimal objective {}\n",S.value);
  } else {
    vprint(2,"sadp found objective {} lower bound is: {}\n",S.value,time_double(states[0]));
    optimal=false;
  }
  //free memory
  ap::freeAP(I.m);
  if(opt.sadpmemory==UINT_MAX) freeMemoryCBFSO();
  return {S,optimal};
}

void freeMemoryCBFSO() {
  free(st.explicitWorkers);
  free(stOld.explicitWorkers);
  free(stNew.explicitWorkers);
  free(stTest.explicitWorkers);
  if (states) free(states);
  if (trace) free(trace);
}

unsigned long checkLowerBoundSADP(const Instance& I,unsigned long current,int t) {
  //1. encode state
  stTest.tasks=t+1;
  stTest.numImplicitWorkers=0;
  for(auto i=0u;i<I.m;i++) { //implicit workers
    if((explicitWorkers[i]==false)&&(checkbit(current,i)>0))
      stTest.numImplicitWorkers++;
  }
  for(auto i=0u;i!=numExplicitWorkers;i++) { //explicitWorkers
    if(checkbit(current,mapExplicitWorkers[i])>0) {
      stTest.explicitWorkers[i]=true;
    } else {
      stTest.explicitWorkers[i]=false;
    }
  }
  vprint(3,"({},{}) ",stTest.tasks,stTest.numImplicitWorkers);
  //2. find the state id
  auto id=encodeState(I,stTest);
  //3. return value in states
  return(states[id]);
}
}
