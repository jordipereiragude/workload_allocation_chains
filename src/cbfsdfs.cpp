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
#include "cbfsdfs.hpp"
#include <boost/heap/pairing_heap.hpp>

#define FMT_HEADER_ONLY
#include "fmt/core.h"

#define time_double(x) (double(x)/precision)

using namespace boost;
using namespace std;

namespace cbfsDfsF {

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
  vector<unsigned int>workersVector;
  bool inCbfs=true;
  bool nodeSelected;
  int depth=(-1);

  //initialize queue and map
  q[0].push({0,-1,workersVector}); 
  memory[0]= -1;
  for(;;) {
    if(::S.elapsed()>opt.tlim) return {S,false};
    if(memory.size()>=opt.cbfsmemory) inCbfs=false;
    if(inCbfs==true) {
      auto initialDepth=depth;
      depth++;
      if(depth==I.m) depth=0;
      current=LONG_MAX;
      nodeSelected=false;
      for(;depth<I.m;depth++) {
        do {
          if(q[depth].size()>0) {
            auto [c, p, vect] = q[depth].top();
            current=c;
            performedTasks=p;
            workersVector=vect;
            q[depth].pop();
            if(memory[current]<performedTasks) {
              current=LONG_MAX;
            } else { 
              nodeSelected=true;
            }
          } else break; //empty queue
        }while(nodeSelected==false);
        if(current<LONG_MAX) break;
      }
      if(current==LONG_MAX) depth=0;
      for(;depth<=initialDepth;depth++) {
        do{
          if(q[depth].size()>0) {
            auto [c, p, vect] = q[depth].top();
            current=c;
            performedTasks=p;
            workersVector=vect;
            q[depth].pop();
            if(memory[current]<performedTasks) {
              current=LONG_MAX;
            } else {
              nodeSelected=true;
            }
          } else break; //empty queue
        }while(nodeSelected==false);
        if(nodeSelected==true) break;
        if(current<LONG_MAX) break;
      }
    } 
    else { //dfs
      current=LONG_MAX;
      for(depth=I.m-1;depth>=0;depth--) {
        if(q[depth].size()>0) {
          auto [c, p, vect] = q[depth].top();
          current=c;
          performedTasks=p;
          workersVector=vect;
          q[depth].pop();
          //while we do not store in memory anymore, it is still worth checking it
          auto it=memory.find(current);
          if((it!=memory.end())&&(it->second>performedTasks)) { //do not use
            depth++; //skip this partial solution but look again in the same queue
            current=LONG_MAX;
          } else { //use
            break;
          }
        }
      }
    }
    if(current==LONG_MAX) return {S,true}; //no state to explore
    //branch
    for(auto w=0u;w<I.m;w++) { //pass through workers
      if (checkbit(current,w)==0) { // unassigned
        //obtain largest number of tasks to add
        long int maxValue=I.d[performedTasks][w]+currentC;
        for(t=performedTasks+1;t<int(I.n)&&(I.d[t][w]<maxValue);t++) { }
        t=t-1;
        if(t==(I.n-1)) { //complete solution
          vprint(3,"solution found.\n");
          workersVector.push_back(w);
          if(workersVector.size()<I.m) {
            setbit(current,w);
            for(auto ww=0u;ww<I.m;ww++) {
              if(checkbit(current,ww)==0){
                workersVector.push_back(ww);
              }
            }
          }
          Solution S2=optW(I,workersVector);
          vprint(2,"{}\n",S2.to_string());
          return {S2,true};
        } else {
          if((depth+1)<I.m) { //otherwise it's a no-end
            setbit(current,w);
            workersVector.push_back(w);
            //check lower bound
            if(sadpO::checkLowerBoundSADP(I,current,t)<currentC) { 
              auto it=memory.find(current);
              if((it==memory.end())||(it->second<t)) { //store
                q[depth+1].push({current,t,workersVector});
                if(inCbfs==true) {
                  memory[current]=t;
                }
              }
            }
            unsetbit(current,w);
            workersVector.pop_back();
          }
        }
      }
    }
  }
  return {S,true};
}

pair<Solution,bool> cbfsDfsFeasibility(const Instance& I, Solution &S) { 
  bool isOptimal;

  tie(S,isOptimal)=sadpO::sadpOptimization(I, S);
  if(isOptimal==false) {
    long unsigned int currentC=S.value*precision; //current cycle time (bad practice).
    currentC=currentC+2;
    vprint(2,"Incumbent objective {}\n",S.value);
    //start search for improving solutions
    while (true) {
      vprint(2,"currentBest: {}\n",S.value);
      tie(S,isOptimal)=cbfsFixedC(I,S,currentC);
      if(isOptimal==false) break;
      if (S.value*precision<currentC)
        currentC=S.value*precision;
      else {
        break;
      }
    }
    sadpO::freeMemoryCBFSO();
  }
  return {S,isOptimal};
}

} //end cbfsAltF

namespace cbfsDfsO {

struct compare_map_cbfsOpt {
  bool operator()(const pair<unsigned long,int>& a, const pair<unsigned long,int>& b) const {
    if(a.first<b.first) return true;
    if((a.first==b.first)&&(a.second<b.second)) return true;
    return false;
  }
};

pair<Solution,bool> cbfsOpt(const Instance& I,Solution &S,unsigned int currentC)
{
  int t,tMin,tMax;
  //key <workers,lastTask> value c (either a bound or current c)
  map<pair<unsigned long,int>,unsigned long int,compare_map_cbfsOpt> memory; 
  vector< priority_queue<OState,std::vector<OState>> > q(I.m);  //I.m just for ease of reconstruction
  OState st;
  vector<unsigned int>workersVector;
  int depth=(-1);
  bool nodeSelected;
  bool inCbfs=true;

  //initialize queue and map
  auto bound=sadpO::checkLowerBoundSADP(I,0,-1);
  q[0].push({0,-1,bound,workersVector});
  memory[{0,-1}]=bound;
  for(;;) {
    if(::S.elapsed()>opt.tlim) return {S,false};
    if((inCbfs==true)&&(memory.size()>=opt.cbfsmemory)) {
      inCbfs=false;
      vprint(1,"moving to dfs\n");
    }
    if(inCbfs==true) {
      auto initialDepth=depth;
      depth++;
      if(depth==I.m) depth=0;
      nodeSelected=false;
      for(;depth<I.m;depth++) {
        do {
          if(q[depth].size()>0) {
            st = q[depth].top();
            q[depth].pop();
            if((st.c<currentC)&&(memory[make_pair(st.workers,st.lastTask)]==st.c)) {
              nodeSelected=true;
            }
          } else break;
        }while(nodeSelected==false);
        if(nodeSelected==true) break;
      }
      if(nodeSelected==false) {
        depth=0;
        for(;depth<=initialDepth;depth++) {
          do {
            if(q[depth].size()>0) {
              st=q[depth].top();
              q[depth].pop();
              if((st.c<currentC)&&(memory[make_pair(st.workers,st.lastTask)]==st.c)) {
                nodeSelected=true;
              }
            }else break;
          }while(nodeSelected==false);
          if(nodeSelected==true) break;
        }
      }
    } 
    else { //dfs 
      nodeSelected=false; 
      for(depth=I.m-1;depth>=0;depth--) { 
        do{ 
          if(q[depth].size()>0) { 
            st=q[depth].top();
            q[depth].pop(); 
            nodeSelected=true;
            if(st.c>=currentC) nodeSelected=false;
            auto it=memory.find(make_pair(st.workers,st.lastTask));
            if((it!=memory.end())&&(it->second<st.c)) nodeSelected=false;
          } else break;
        }while(nodeSelected==false);
        if(nodeSelected) break;
      }
    }
    if(nodeSelected==false) return {S,true}; //no state to explore (optimal solution)

    vprint(2,"depth: {} state {} performedTasks {} cycle: {} s: {}\n",depth,st.workers,st.lastTask,st.c,st.workersInUse);
    for(auto w=0u;w<I.m;w++) { //pass through workers
      if(checkbit(st.workers,w)==0) {
        setbit(st.workers,w);
        st.workersInUse.push_back(w);
        long int minValue=I.d[st.lastTask][w]+st.c;
        for(tMin=st.lastTask+1;tMin<int(I.n)&&(I.d[tMin][w]<minValue);tMin++) { }
        tMin=tMin-2;
        if(tMin<st.lastTask) tMin++;
        long int maxValue=I.d[st.lastTask][w]+currentC;
        for(tMax=st.lastTask+1;tMax<int(I.n)&&(I.d[tMax][w]<maxValue);tMax++) { }
        if(tMax==(I.n)) { //complete solution
          //reconstruction
          auto cycle=max(st.c,I.d[tMax-1][w]-I.d[st.lastTask][w]);
          vprint(1,"solution found with cycle {} vs {}\n",cycle,currentC);
          if(cycle<currentC) { //it improves
            if(st.workersInUse.size()<I.m) {
              for(auto ww=0u;ww<I.m;ww++) {
                if(checkbit(st.workers,ww)==0) {
                  st.workersInUse.push_back(ww);
                }
              }
            }
            Solution S2=optW(I,st.workersInUse);
            vprint(1,"{}\n",S2.to_string());
            currentC=S2.value*precision;
            S=S2;
          }
        } else {
          if((depth+1)<I.m) { //otherwise it's a no-end
            for(t=tMin;t<tMax;t++) { //pass through different number of tasks
              pair<unsigned long,int> key={st.workers,t};
              bound=max(st.c,sadpO::checkLowerBoundSADP(I,st.workers,t));
              bound=max(bound,I.d[t][w]-I.d[st.lastTask][w]);
              if(bound<currentC) {
                auto it=memory.find(key);
                if((it==memory.end())||(it->second>bound)) {
                  q[depth+1].push({st.workers,t,bound,st.workersInUse});
                  if(inCbfs==true) {
                    memory[key]=bound;
                  }
                }
              }
            }
          }
        }
        unsetbit(st.workers,w);
        st.workersInUse.pop_back();
      }
    }
  }
}

pair<Solution,bool> cbfsDfsOptimization(const Instance& I, Solution &S){
  bool isOptimal;

  tie(S,isOptimal)=sadpO::sadpOptimization(I, S);
  if(isOptimal==true) return {S,true};
  else {
    unsigned int currentC=S.value*precision;
    vprint(2,"\nInitial objective cbfsOptimization {}\n",S.value);
    tie(S,isOptimal)=cbfsOpt(I,S,currentC);
    vprint(1,"\nFinal objective cbfsOptimization {}\n",S.value);
    sadpO::freeMemoryCBFSO();
    return {S,isOptimal};
  }
}

} //end cbfsDfsO
