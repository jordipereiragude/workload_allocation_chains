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
#pragma once

#include <iostream>
#include <vector>

#include <boost/multi_array.hpp>

#include "dimensions.hpp"
#define time_double(x) (double(x)/precision)

const long precision = 1000000;

struct Instance {
  unsigned n;			// number of tasks
  unsigned m;			// number of workers
  unsigned c;			// cycle time
  unsigned B;			// batch size

  std::vector<Time> s,s0;	// standard times (expanded, regular)
  boost::multi_array<Time,2> t;	// execution times; t[i,w] for task i=[0,n-1], worker w=[0,m-1]
  int ninf;		        // number of infeasible task-worker pairs

  boost::multi_array<Time,2> T;	// prefix sums of execution times; T[i,w]=sum(t[0:i-1][w]) for (excl) task i=[0,n], worker w=[0,m-1]
  std::vector<Time> T_; // ditto, for smallest times
  
  boost::multi_array<unsigned long,2> d; // prefix sums of task times, with fixed precision, d[i][j]=sum_{k=0}^i t[k][j], for i\in [-1,n), j\in[0,m)
  boost::multi_array<unsigned long,2> db; // prefix sums of task times, with fixed precision (reverse)
  void computeD();

  void readSTD(std::istream&);	  // read from stream in STD format
  void readALWABP(std::istream&); // read from stream in ALWABP format
  double normtime(double);

  int getB();
  double getMinTime();
  double segCost(int i, int j, int w) const; // cost of segment [i,j] for worker w
  double segCost(int i, int j) const;        // cost of segment [i,j], cheapest over all workers
  Time sumStdTimes() const;		     // sum of operation times

  // helpers
  std::vector<Time> getTask(std::istream&); // process a line of time values
};
