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
#include "beamSearch.hpp"

#include <queue>
#include <map>
#include <cassert>
#include <unordered_map>
using namespace std;

#include <boost/multi_array.hpp>
#include <boost/heap/pairing_heap.hpp>
using namespace boost;

#include "util.hpp"
#include "statistics.hpp"
#include "options.hpp"
#include "hs.hpp"
#include "sadp.hpp"

struct TPW {

  // worker `j` starting from task `i` can do tpw[i][j], for fixed candidate cycle time.
  int **tpw;

  TPW() : tpw(nullptr) {}

  TPW(const Instance& I) {
    allocate(I);
  }

  ~TPW() {
    free();
  }

  void allocate(const Instance& I) {
    tpw=(int**)malloc(sizeof(int) * (I.n+1) * (I.m+1) + sizeof(int *) * (I.n+1) );
    for(unsigned i=0; i<=I.n; i++)  tpw[i]=(int *)(tpw+I.n+1)+i*(I.m+1);
  }

  void free() {
    if (tpw)
      ::free(tpw);
  }

  void compute(const Instance& I, unsigned cc) {
    //first compute number of tasks per worker per operation-id
    for(unsigned j=0; j<I.m; j++) {
      unsigned last=0;
      for(unsigned i=0; i<I.n; i++) {
	while (last<I.n && (I.d[last][j]-I.d[i][j])<cc)
	  last++;
	tpw[i][j]=last-i;
      }
    }
    for(unsigned j=0; j<I.m; j++)
      tpw[I.n][j]=0;


    //second select largest form i to n (reverse operation is better and avoids keeping second array.
    for(unsigned j=0; j<I.m; j++)
      for(int i=I.n-1; i>=0; i--)
	if(tpw[i][j]<tpw[i+1][j])
	  tpw[i][j]=tpw[i+1][j];
  }

  void computeBackward(const Instance& I, unsigned cc) {
    //first compute number of tasks per worker per operation-id
    for(unsigned j=0; j<I.m; j++) {
      unsigned last=0;
      for(unsigned i=0; i<I.n; i++) {
	while (last<I.n && (I.db[last][j]-I.db[i][j])<cc)
	  last++;
	tpw[i][j]=last-i;
      }
    }
    for(unsigned j=0; j<=I.m; j++)
      tpw[I.n][j]=0;

    //second select largest form i to n (reverse operation is better and avoids keeping second array.
    for(unsigned j=0; j<I.m; j++)
      for(int i=I.n-1; i>=0; i--)
	if(tpw[i][j]<tpw[i+1][j])
	  tpw[i][j]=tpw[i+1][j];
  }

  // access
  int operator()(int t, int w) const {
    return tpw[t][w];
  }
};

// return the sum of the tasks all unassigned workers can do, all starting from task `t`
//   tpw[t][w] contains the maximum number of consecutive tasks that worker w can perform starting from task t (t included)
int taskBasedBound(unsigned numWorkers, const TPW& tpw, unsigned long partialSolution, int t) {
  int ret=0;
  for(unsigned w=0; w<numWorkers; w++)
    if (checkbit(partialSolution,w)==0)
      ret += tpw(t,w);
  return ret;
}

// The key is an integer bitset representing a DP state in the exploration, value is the number of performed tasks in this state.
//   Shared map to use in pqOrder (should pass as a parameter to the compare function (define as a struct, etcetera).
map<unsigned long,int> avoidRepetition;
// state predecessors, used in CBFS
unordered_map<unsigned long,int> pred;

bool pqOrder(unsigned long a, unsigned long b) {
  return avoidRepetition[a]>avoidRepetition[b];
}

// Decode solution `current` by following the chain of stored predecessors `pred`.
Solution decode_pred(const Instance& I, unsigned long current, unsigned long currentC) {
  vprint(3, "Decode solution for {}\n", currentC);
  vector<unsigned> vWorkers;
  int t = I.n-1;
  while (pred[current]>=0) {
    unsigned w = pred[current];
    vWorkers.push_back(w);
    unsetbit(current,w);
    int tn = avoidRepetition[current];
    t = tn;
  }
  assert(current==0);
  reverse(vWorkers.begin(),vWorkers.end());
  // return the optimal solution for this worker sequence (assert: opw(I,vWorkers).value<=currentC)
  return optW(I,vWorkers);
}

// decode solution `current` by checking possible predecessors
Solution decode(const Instance& I, unsigned long current, unsigned long currentC)
{
  vprint(3, "Decode solution for {}\n", currentC);

  vector<unsigned> vWorkers(I.m); // the extracted worker sequence

  // (1) go backwards over tasks, extract workers
  int t = I.n-1;
  for(int depth=I.m-1; depth>0; depth--)
  {
    for(unsigned w=0; w<I.m; w++)
    {
      if(checkbit(current,w)>0)   //if assigned, it may be the worker assigned to the said position
      {
        unsetbit(current,w);
        //look if it exists in avoidRepetition
        auto it=avoidRepetition.find(current);
        if (it!=avoidRepetition.end())   //if it exists, it may be the worker in depth position (otherwise, we are sure that it is not the case
        {
          if ((I.d[t][w]-I.d[it->second][w])<currentC)   // requires less than current c, means it is possible to construct the path
          {
            vWorkers[depth]=w;
            t=it->second;
            //l is missing
            break; //stop loop
          }
        }
        setbit(current,w); //if it arrives here, it means it is not this worker
      }
    }
  }

  // (2) extract the last worker
  for(unsigned i=0; i<I.m; i++)
  {
    if (checkbit(current,i)>0)
    {
      vWorkers[0]=i;
      unsetbit(current,i);
      break;
    }
  }
  assert(current==0);

  // sequence is ready. Now I want to double-check cycle time
  return optW(I,vWorkers);
}

//not the best way to program backward decoded. Kept separately from decode not to mess with DP and CBFS code. Also it makes additional step to "swap" sequence of workers
Solution decodeBackward(const Instance& I, unsigned long current, unsigned long currentC)
{
  vprint(3,"Decode solution for {}\n",currentC);
  vector<unsigned> vWorkers(I.m);
  vector<unsigned> vWorkersBackward(I.m);
  int t = I.n-1;
  for(int depth=I.m-1; depth>0; depth--)
  {
    for(unsigned w=0; w<I.m; w++)
    {
      if(checkbit(current,w)>0)   //if assigned, it may be the worker assigned to the said position
      {
        unsetbit(current,w);
        //look if it exists in avoidRepetition
        auto it=avoidRepetition.find(current);
        if (it!=avoidRepetition.end())   //if it exists, it may be the worker in depth position (otherwise, we are sure that it is not the case
        {
          if ((I.db[t][w]-I.db[it->second][w])<currentC)   // requires less than current c, means it is possible to construct the path
          {
            vWorkers[depth]=w;
            t=it->second;
            //l is missing
            break; //stop loop
          }
        }
        setbit(current,w); //if it arrives here, it means it is not this worker
      }
    }
  }
  for(unsigned i=0; i<I.m; i++)
  {
    if (checkbit(current,i)>0)
    {
      vWorkers[0]=i;
      unsetbit(current,i);
      break;
    }
  }
  assert(current==0);
  for(unsigned i=0; i<I.m; i++)
    vWorkersBackward[i]=vWorkers[I.m-1-i];
  //sequence is ready. Now I want to double-check cycle time
  return optW(I,vWorkersBackward);
}


Solution tryBS(const Instance &I, Solution &S, unsigned long currentC, unsigned beamSize, TPW& tpw, int isForward)
{
  //memory structures for beam search
  map<unsigned long,int>::iterator it;  //key is an integer bitset to consider only one copy of a DP state in the exploration.
  priority_queue<unsigned long,std::vector<unsigned long>, std::function<bool(unsigned long,unsigned long)>> oldStates(pqOrder); //must order from worst to best element according to value recorded in map
  vector<unsigned long> newStates;
  unsigned long current=0;	// current state (allocated workers)
  unsigned int maxValue, t;
  int lastTaskPerformed;		// number of last performed task

  //empty for repeated use
  avoidRepetition.clear();

  if (isForward)
    tpw.compute(I,currentC);
  else
    tpw.computeBackward(I,currentC);

  vprint(2,"BS testing currentC: {} beamSize: {}\n",time_double(currentC),beamSize);

  oldStates.push(current);
  avoidRepetition[current]=(-1);
  for(unsigned depth=0; depth<I.m; depth++)
  {
    assert(newStates.empty());
    while(!(oldStates.empty()))
    {
      current=oldStates.top();
      lastTaskPerformed=avoidRepetition[current];
      oldStates.pop();
      for(unsigned w=0; w<I.m; w++) //check all workers
      {
        if(checkbit(current,w)==0) //unassigned
        {
          setbit(current,w);
          //we must look at the number of tasks that the worker can perform for a given c
          if(isForward==true)
          {
            maxValue=I.d[lastTaskPerformed][w]+currentC;
            for(t=lastTaskPerformed+1; (t<I.n)&&(I.d[t][w]<maxValue); t++) {  }
          }
          else
          {
            maxValue=I.db[lastTaskPerformed][w]+currentC;
            for(t=lastTaskPerformed+1; (t<I.n)&&(I.db[t][w]<maxValue); t++) {  }
          }
          //check if it may improve (only implemented for forward, hence we always pass the test if we work in backward direction
          if((t+taskBasedBound(I.m,tpw,current,t))>=I.n)
          {
            t=t-1;
	    stat.expansions++;

            //check if already exists
            it=avoidRepetition.find(current);
            if(it==avoidRepetition.end()) //it it does not exist, save
            {
              avoidRepetition[current]=t;
              newStates.push_back(current);
            }
            else //if it exists, save changes to maxTasks only
            {
              if(it->second<int(t))
              {
                it->second=t;
              }
            }
          }
          unsetbit(current,w);
        }
      }
    }
    vprint(3,"finished with size: {}\n",newStates.size());
    //now we have to move newStates into oldStates
    for (unsigned i=0; i<newStates.size(); i++) //for each generated stage
    {
      if(oldStates.size()<beamSize)  //if not filled, save
      {
        oldStates.push(newStates[i]);
      }
      else //if filled only save if it improves
      {
        if(avoidRepetition[oldStates.top()]<avoidRepetition[newStates[i]])
        {
          oldStates.pop();
          oldStates.push(newStates[i]);
        }
      }
    }
    newStates.clear();
  }
  if(oldStates.size()<1)
  {
    return S;
  }
  assert(oldStates.size()==1); //check all workers have been assigned
  current=oldStates.top();
  t=avoidRepetition.find(current)->second; //check all tasks have been assigned (otherwise there is no improvement)
  assert((t+1)<=I.n);
  if((t+1)!=I.n) return S; //otherwise no improvement, hence return
#ifndef NDEBUG
  current=0;
  for(unsigned i=0; i<I.m; i++) setbit(current,i);
  assert(oldStates.top()==current);
#endif

  current=oldStates.top();
  if(isForward==true) S=decode(I,current,currentC);
  else S=decodeBackward(I,current,currentC);
  //sequence is ready. Now I want to double-check cycle time
  assert(currentC>(S.value*precision));
  return S;
}

typedef tuple<unsigned long,int,double> BState; // (workers,last done task)
namespace std {
  inline bool operator<(const BState& a, const BState& b) {
    return get<2>(a)>get<2>(b);
  }
}

typedef heap::pairing_heap<BState> PQ;

bool tasks_consistent(int lastTaskPerformed, map<unsigned long,int>& avoidRepetition, unsigned long current) {
  if (lastTaskPerformed!=avoidRepetition[current])
    fmt::print("For {} : {} versus {}\n", current, lastTaskPerformed, avoidRepetition[current]);
  return lastTaskPerformed<=avoidRepetition[current];
}

// cyclic best-first search for fixed cycle time `currentC`
//   (maybe to be revised: we maintain global state (avoidRepetition etc.) as well as local state (priority queues); the first for extraction, the second for processing.)
pair<Solution,Solver_Status> tryCBFS(const Instance& I, Solution& S, unsigned long currentC, TPW& tpw, DP::optnRepDP& O) {
  // setup priority queues and state
  vector<PQ> q(I.m);
  unordered_map<unsigned long,PQ::handle_type> handle;
  handle[0]=q[0].push({0,-1,O.getLB(0,0)*precision}); // initial state
  avoidRepetition.clear();
  avoidRepetition[0]=-1;
  pred.clear();
  pred[0]=-1;
  double lastReport = -1;

  // setup expansion
  unsigned mindepth = 0;
  unsigned passes = 0, expansions = 0;
  bool expanded = true;

  //tpw.compute(I,currentC);
  // repeated expansion passes
  vprint(3,"cycle: {}\n", currentC);

  while (expanded && passes<opt.maxpasses) {
    expanded = false;
    for(unsigned depth=mindepth; depth!=I.m; depth++) {
      // take best from current depth and expand
      if (q[depth].size()==0) {
        if (depth==mindepth)
          mindepth++;
        continue;
      }

      expanded = true;
      expansions++;
      auto [current, lastTaskPerformed, bound] = q[depth].top();
      q[depth].pop();
      handle.erase(current);

      assert(tasks_consistent(lastTaskPerformed,avoidRepetition,current));

      if (lastTaskPerformed==avoidRepetition[current]) { // otherwise there is a better method to reach the state, hence no need to explore
        for(unsigned w=0; w!=I.m; w++) { // check all workers
          if (checkbit(current,w)!=0) // assigned
            continue;
          setbit(current,w);

	  unsigned t = maxTask(I,w,lastTaskPerformed,currentC);
	  t++; // make t: # of done tasks = first undone task

	  // check LB3
	  unsigned long lb3 = O.getLB(t,current)*precision;
	  
	  if (lb3>=currentC) {
	    unsetbit(current,w);
	    continue;
	  }

          // process
	  t=t-1; // make t: last task done
	  stat.expansions++;

	  // full solution
	  if (depth+1==I.m) {
	    if (t+1==I.n) {
	      // full solution found; decode and return
	      avoidRepetition[current]=t;
	      pred[current]=w;
	      S=decode_pred(I,current,currentC);
	      vprint(3,"{:6}; {:5} passes: [{}] ", avoidRepetition.size(), passes, mindepth);
	      for(unsigned depth=mindepth; depth!=I.m; depth++)
		vprint(3,"{} ",q[depth].size());
	      vprint(3,"\n");
	      return {S,IMPROVED};
	    }
	    continue;
	  }

	  // update handles

	  //check if already exists
	  auto it=handle.find(current);
	  if (it==handle.end() && (avoidRepetition.find(current)==avoidRepetition.end() || avoidRepetition[current]<int(t))) { // not active and never seen or better than seen: make active
	    handle[current]=q[depth+1].push({current,t,O.getLB(t+1,current)});
	  } else if (it!=handle.end()) {  // active
	    if (get<1>(*handle[current])<int(t)) { // and better: update
	      q[depth+1].update(handle[current],{current,t,O.getLB(t+1,current)});
	    } // and worse: do nothing
	  } // not active, but worse that seen: do nothing

	    // update best found
	  if (avoidRepetition.find(current)==avoidRepetition.end() || avoidRepetition[current]<int(t)) {
	    avoidRepetition[current]=t;
	    pred[current]=w;
	  }
	  unsetbit(current,w);
	}
      } // lastTaskPerformend == avoidRepetition
    } // depth loop
    passes++;
    if (verbosec.count>=2 && ::S.elapsed()-lastReport>0.5) {
      vprint(2,"{:6}; {:5} passes: [{}] ", avoidRepetition.size(), passes, mindepth);
      for(unsigned depth=mindepth; depth!=I.m; ++depth)
	vprint(2,"{} ", q[depth].size());
      vprint(2,"\n");
      lastReport = ::S.elapsed();
    }
  } // passes


  return {S,passes<opt.maxpasses?OPTIMAL:LIMIT};
}

//simple iterative beam search method (needs initial solution)
Solution beamSearch(const Instance& I, Solution &S,int minBeamSize,int maxBeamSize) {
  unsigned int currentC=S.value*precision; //current cycle time (bad practice, sorry. I'll make it more general if we move these data to the instance definition)

  TPW tpw(I);

  vprint(2,"Starting beam search currentC {} numTasks: {} numWorkers: {}\n",time_double(currentC),I.n,I.m);

  int beamSize=minBeamSize;
  while(::S.elapsed()<opt.tlim) { //keep trying while solution improves
    vprint(1,"{} ",time_double(currentC)); fflush(stdout);
    S=tryBS(I,S,currentC,beamSize,tpw,true);
    if (S.value*precision<currentC)
      currentC=S.value*precision;
    else {
      S=tryBS(I,S,currentC,beamSize,tpw,false);
      if (S.value*precision<currentC)
	currentC=S.value*precision; //tries to improve reverse order
      else {
        beamSize=beamSize*2;
        if (beamSize>maxBeamSize) break;
      }
    }
  }
  vprint(1,"\n");

  return S;
}

// cyclic best-first search starting from initial solution `S`
pair<Solution,Solver_Status> CBFS(const Instance& I, Solution &S, DP::optnRepDP& O) {
  unsigned int currentC=S.value*precision; //current cycle time (bad practice, sorry. I'll make it more general if we move these data to the instance definition)

  TPW tpw;

  while (::S.elapsed()<opt.tlim) { // keep trying while solution improves
    vprint(1,"{} ",time_double(currentC)); fflush(stdout);
    Solver_Status status;
    tie(S,status)=tryCBFS(I,S,currentC,tpw,O);
    if (S.value*precision<currentC) {
      assert(status==IMPROVED);
      currentC=S.value*precision;
    } else {
      vprint(1,"\n");
      return {S,status};
    }
  }
  assert(false);
  return {S,TIMELIMIT};
}
