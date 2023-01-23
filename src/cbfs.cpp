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
#include <utility>
#include "util.hpp"
#include "hs.hpp"
#include "sadp.hpp"
#include "ap.hpp"
#include "options.hpp"
#include <boost/heap/pairing_heap.hpp>

#define FMT_HEADER_ONLY
#include "fmt/core.h"

#define time_double(x) (double(x)/precision)

using namespace boost;
using namespace std;

typedef pair<unsigned long,int> FState; // (workers,last done task)

struct OState{
  unsigned long workers;
  int lastTask;
  unsigned long c;
};

namespace std {
  inline bool operator<(const FState& a, const FState& b) {
    return a.second<b.second;
  }
}

struct compare_pq {
  bool operator()(const OState& a, const OState& b) {
    if(a.lastTask<b.lastTask) return true;
    if((a.lastTask==b.lastTask)&&(a.c<b.c)) return true;
    return false;
  }
};


namespace cbfsF {

typedef heap::pairing_heap<FState> PQ;

pair<Solution,bool> cbfsFixedC(const Instance& I, Solution &S, long unsigned int currentC) {
  int t;
  //memory to store solutions
  map<unsigned long,int> memory;
  //priority queues with solution
  vector<PQ> q(I.m);
  //current solution
  unsigned long current;
  int performedTasks;
  bool update,nodeSelected;

  //initialize queue and map
  q[0].push({0,-1}); 
  memory[0]= -1;
  do { //loop until queues are empty
    if(::S.elapsed()>opt.tlim) return {S,false};
    update=false;
    for(auto depth=0u;depth<I.m;depth++) { //loop through levels
      do {
        nodeSelected=false;
        current=LONG_MAX;
        if(q[depth].size()>0) {
          auto [c, p] = q[depth].top();
          current=c;
          performedTasks=p;
          q[depth].pop();
          nodeSelected=true;
        }else break;
        if(memory[current]<performedTasks) {
          current=LONG_MAX;
          nodeSelected=false;
        }
        else update=true;
      }while(nodeSelected==false);
      if(current<LONG_MAX) { //there is a state
        assert(memory[current]==performedTasks);
        for(auto w=0u;w<I.m;w++) { //pass through workers
          if (checkbit(current,w)==0) { // unassigned
            //obtain largest number of tasks to add
            long int maxValue=I.d[performedTasks][w]+currentC;
            for(t=performedTasks+1;t<int(I.n)&&(I.d[t][w]<maxValue);t++) { }
            t=t-1;
            if(t==(I.n-1)) { //complete solution
              auto currentOriginal=current;
              vprint(2,"solution found.\n");
              //rebuild solution
              vector<int> cumTasks(I.m,0); //cumulative number of tasks up to workstation
              vector<unsigned int> workersInUse(I.m,I.m); //which worker works at workstation 
              workersInUse[depth]=w;
              cumTasks[depth]=performedTasks;
              for(auto level=depth-1;level>0;level--) { //trace solution
                for(auto ww=0u;ww<I.m;ww++) {
                  if(checkbit(current,ww)) { //check each worker still available
                    unsetbit(current,ww);
                    auto it=memory.find(current);
                    if((it!=memory.end())&&((I.d[cumTasks[level+1]][ww]-I.d[it->second][ww])<currentC)) {
                      workersInUse[level]=ww;
                      cumTasks[level]=it->second;
                      break;
                    }
                    setbit(current,ww);
                  }
                }
                assert(workersInUse[level]<I.m);
              }
              for(auto ww=0u;ww<I.m;ww++) { //first worker
                if(checkbit(current,ww)>0) {
                  workersInUse[0]=ww;
                  break;
                }
              }
              auto currentWorkstation=depth;
              for(auto ww=0u;ww<I.m;ww++) {
                if(checkbit(currentOriginal,ww)==0) {
                  workersInUse[currentWorkstation]=ww;
                  currentWorkstation++;
                }
              }
              Solution S2=optW(I,workersInUse);
              vprint(1,"{}\n",S2.to_string());
              return {S2,true};
            } else {
              if((depth+1)<I.m) { //otherwise it's a no-end
                setbit(current,w);
                //check lower bound
                if(sadpO::checkLowerBoundSADP(I,current,t)<currentC) { 
                  auto it=memory.find(current);
                  if((it==memory.end())||(it->second<t)) { //store
                    memory[current]=t;
                    if(memory.size()>=opt.cbfsmemory) return {S,false};
                    q[depth+1].push({current,t});
                  }
                }
                unsetbit(current,w);
              }
            }
          }
        }
      }
    }
  }while(update==true);
  return {S,true};
}

pair<Solution,bool> cbfsFeasibility(const Instance& I, Solution &S) { 
  bool isOptimal;

  tie(S,isOptimal)=sadpO::sadpOptimization(I, S);
  if(isOptimal==false) {
    long unsigned int currentC=S.value*precision; //current cycle time (bad practice).
    vprint(2,"\nIncumbent objective {}\n",S.value);
    //start search for improving solutions
    while (true) {
      vprint(2,"currentBest: {}\n",S.value);
      tie(S,isOptimal)=cbfsFixedC(I,S,currentC);
      if(isOptimal==false) break;
      if (S.value*precision<currentC) {
        currentC=S.value*precision;
      } else {
        break;
      }
    }
    sadpO::freeMemoryCBFSO();
  }
  return {S,isOptimal};
}

} //end cbfsF

namespace cbfsO {

struct compare_map_cbfsOpt {
  bool operator()(const pair<unsigned long,int>& a, const pair<unsigned long,int>& b) const {
    if(a.first<b.first) return true;
    if((a.first==b.first)&&(a.second<b.second)) return true;
    return false;
  }
};

pair<Solution,bool> cbfsOpt(const Instance& I,Solution &S,unsigned int currentC) {
  int t,tMin,tMax;
  map<pair<unsigned long,int>,unsigned long int,compare_map_cbfsOpt> memory; //key <workers,lastTask> value c (either a bound or current c)
  vector< priority_queue<OState,std::vector<OState>,compare_pq> > q(I.m+1);  //I.m just for ease of reconstruction
  bool update,valid,nodeSelected;
  OState st;

  vprint(1,"In cbfsOpt currentC {}\n",currentC);

  //initialize queue and map
  auto bound=sadpO::checkLowerBoundSADP(I,0,-1);
  q[0].push({0,-1,bound});
  memory[{0,-1}]=bound;
  do{ //loop until queues are empty
    if(::S.elapsed()>opt.tlim) return {S,false};
    update=false;
    valid=false;
    for(auto depth=0;depth<I.m;depth++) {
      do {
        nodeSelected=false;
        if(q[depth].size()>0) {
          st=q[depth].top();
          q[depth].pop();
          valid=true;
          nodeSelected=true;
        }else break;
        if((currentC<=st.c)||(memory[make_pair(st.workers,st.lastTask)]<st.c)) nodeSelected=false;
        else update=true;
      }while(nodeSelected==false);
      if(valid==true) { //there is a state
        vprint(3,"depth: {} state {} performedTasks {} cycle: {}\n",depth,st.workers,st.lastTask,st.c);
        for(auto w=0u;w<I.m;w++) { //pass through workers
          if(checkbit(st.workers,w)==0) {
            setbit(st.workers,w);
            long int minValue=I.d[st.lastTask][w]+st.c;
            for(tMin=st.lastTask+1;tMin<int(I.n)&&(I.d[tMin][w]<minValue);tMin++) { }
            tMin=tMin-2;
            if(tMin<st.lastTask) tMin++;
            long int maxValue=I.d[st.lastTask][w]+currentC;
            for(tMax=st.lastTask+1;tMax<int(I.n)&&(I.d[tMax][w]<maxValue);tMax++) { }
            if(tMax==(I.n)) { //complete solution
              //reconstruction
              auto cycle=max(st.c,I.d[tMax-1][w]-I.d[st.lastTask][w]);
              vprint(2,"complete solution {} vs {}\n",cycle,currentC);
              if(cycle<currentC) { //it improves
                currentC=cycle;
                auto current=st.workers;
                unsetbit(current,w);
                //rebuild solution              
                vector<int> cumTasks(I.m,0); //cumulative number of tasks up to workstation
                vector<unsigned int> workersInUse(I.m,I.m); //which worker works at workstation 
                workersInUse[depth]=w;
                cumTasks[depth]=st.lastTask;
                for(auto level=depth-1;level>0;level--) { //trace solution
                  auto found=false;
                  for(auto ww=0u;ww<I.m;ww++) {
                    if(checkbit(current,ww)) { //check each worker still available
                      unsetbit(current,ww);
                      for(auto ct=cumTasks[level+1];(ct>0)&&((I.d[cumTasks[level+1]][ww]-I.d[ct][ww])<=currentC);ct--) {
                        pair<unsigned long,int> key={current,ct};
                        auto it=memory.find(key);
                        if((it!=memory.end())&&(it->second<=currentC)) {
                          found=true;
                          workersInUse[level]=ww;
                          cumTasks[level]=ct;
                          break;
                        }
                      }
                      if(found==false) setbit(current,ww);
                    }
                    if(found==true) break;
                  }
                  assert(workersInUse[level]<I.m);
                }
                for(auto ww=0u;ww<I.m;ww++) { //first worker
                  if(checkbit(current,ww)>0) {
                    workersInUse[0]=ww;
                    break;
                  }
                }
                auto currentWorkstation=depth+1;
                for(auto ww=0u;ww<I.m;ww++) { //remaining workers
                  if(checkbit(st.workers,ww)==0) {
                    workersInUse[currentWorkstation]=ww;
                    currentWorkstation++;
                  }
                }
                if(currentWorkstation!=I.m) {
                  vprint(1,"error {} vs {}\n",currentWorkstation,I.m);
                }
                Solution S2=optW(I,workersInUse);
                vprint(2,"{}\n",S2.to_string());
                currentC=S2.value*precision;
                //check if S2 is actually correct
                S=S2;
              } //end it improves
            } else {
              if((depth+1)<I.m) { //otherwise it's a no-end
                for(t=tMin;t<tMax;t++) { //pass through different number of tasks
                  //check lower bound
                  pair<unsigned long,int> key={st.workers,t};
                  bound=max(st.c,sadpO::checkLowerBoundSADP(I,st.workers,t));
                  bound=max(bound,I.d[t][w]-I.d[st.lastTask][w]);
                  if(bound<currentC) {
                    auto it=memory.find(key);
                    //if(it!=memory.end()) { vprint(1,"{} ** {} {}\n",key,it->first,it->second);}
                    if((it==memory.end())||(it->second>bound)) { //store (not yet ready)
                      memory[key]=bound;
                      if(memory.size()>=opt.cbfsmemory) return {S,false};
                      q[depth+1].push({st.workers,t,bound});
                    }
                  }
                }
              }
            }
            unsetbit(st.workers,w);
          }
        }
      }
    }
  }while(update==true);
  return {S,true};
}

pair<Solution,bool> cbfsOptimization(const Instance& I, Solution &S){
  bool isOptimal;
  tie(S,isOptimal)=sadpO::sadpOptimization(I, S);
  if(isOptimal==true) return {S,true};
  else {
    unsigned int currentC=S.value*precision;
    vprint(2,"\nInitial objective cbfsOptimization {}\n",S.value);
    tie(S,isOptimal)=cbfsOpt(I,S,currentC);
    sadpO::freeMemoryCBFSO();
    vprint(1,"\nFinal objective cbfsOptimization {}\n",S.value);
    return {S,isOptimal};
  }
}

} //end cbfsO
