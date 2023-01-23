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
#include <boost/heap/pairing_heap.hpp>

#include <queue>
#include <map>
#include <assert.h>

#include "beamSearch.hpp"
#include "util.hpp"
#include "hs.hpp"
#include "options.hpp"

using namespace std;
using namespace boost;

//macros for generating key for the map d is the data, n is the bit number
#define setbit(d,n)   ((d) |=  (1<<(n)))
#define unsetbit(d,n) ((d) &= ~(1<<(n)))
#define checkbit(d,n) ((d) &   (1<<(n)))


// shared map to use in pqOrder
extern map<unsigned long,int> avoidRepetition;  //key is an integer bitset to consider only one copy of a DP states in the exploration.

Solution tryDP(const Instance &I, Solution &S, unsigned long currentC)
{
  //memory structures for beam search
  map<unsigned long,int>::iterator it;  //key is an integer bitset to consider only one copy of a DP state in the exploration.
  priority_queue<unsigned long,std::vector<unsigned long>, std::function<bool(unsigned long,unsigned long)>> oldStates(pqOrder); //must order from worst to best element according to value recorded in map
  vector<unsigned long> newStates;
  unsigned long current=0;	// current state (allocated workers)
  unsigned int maxValue, t;
  int performedTasks;		// number of last performed task

  //empty for repeated use
  avoidRepetition.clear();
  //completely unneeded assertions
  assert(oldStates.empty());
  assert(avoidRepetition.empty());
  assert(newStates.empty());

  vprint(2,"BS testing currentC: {}\n",currentC);

  oldStates.push(current);
  avoidRepetition[current]=(-1);
  for(unsigned depth=0;depth<I.m;depth++)
  {
    assert(newStates.empty());
    vprint(2,"setting up worker at position {} states: {}\n",depth+1,oldStates.size());
    while(!(oldStates.empty()))
    {
      current=oldStates.top();
      performedTasks=avoidRepetition[current];
      oldStates.pop();
      for(unsigned w=0;w<I.m;w++) //check all workers
      {
        if(checkbit(current,w)==0) //unassigned
        {
          setbit(current,w);
          //we must look at the number of tasks that the worker can perform for a given c
          maxValue=I.d[performedTasks][w]+currentC;
          for(t=performedTasks+1;(t<I.n)&&(I.d[t][w]<maxValue);t++) {  }
          t=t-1;
          vprint(3,"worker: {} Code: {} maxValue: {} PerformedTasks: {}\ttasks to perform: {}\n",w,current,maxValue,performedTasks,t);
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
          unsetbit(current,w);
        }
      }
    }
    vprint(2,"finished with size: {}\n",newStates.size());
    //now we have to move newStates into oldStates
    for (unsigned i=0; i<newStates.size(); i++) //for each generated stage
    {
      oldStates.push(newStates[i]);
    }
    newStates.clear();
  }
  assert(oldStates.size()==1); //check all workers have been assigned
  current=oldStates.top();
  t=avoidRepetition.find(current)->second; //check all tasks have been assigned (otherwise there is no improvement)
  assert((t+1)<=I.n);
  if((t+1)!=I.n) return S; //otherwise no improvement, hence return
#ifndef NDEBUG
  current=0;
  for(unsigned i=0;i<I.m;i++) setbit(current,i);
  assert(oldStates.top()==current);
#endif
  
  current=oldStates.top();
  S=decode(I,current,currentC);
  //sequence is ready. Now I want to double-check cycle time
  assert(currentC>(S.value*precision));
  return S;


}

//simple iterative DP (needs initial solution. No time limit)
Solution dynamicProgramming(const Instance& I, Solution &S)
{
  unsigned int currentC=S.value*precision;

  vprint(2,"Starting dynamic programming currentC {} numTasks: {} numWorkers: {}\n",currentC,I.n,I.m);

  while(1) //keep trying while solution improves
  {
    vprint(1,"{} ",time_double(currentC)); fflush(stdout);
    S=tryDP(I,S,currentC);
    if((S.value*precision)<currentC) currentC=S.value*precision; 
    else break;
  }
  vprint(1,"\n");
  
  return S;
}


