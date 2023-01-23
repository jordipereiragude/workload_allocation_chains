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

#include <boost/program_options.hpp>
namespace po = boost::program_options;

// global options
struct Options {
  std::string instance;	  // input file name
  std::string format;	  // input format

  unsigned tlim;	  // time limit, seconds
  unsigned ilim;	  // iteration limit
  unsigned B;             // batch size
  unsigned sadpmemory;    // max memory for sadp (in table entries)
  unsigned cbfsmemory;    // max memory for cbfs (in table entries)
  unsigned sadpFmemory;   // max memory for sadpF (in MB)

  std::string initial, heuristic, exact; // initial and full heuristic algorithm, exact algorithm
  bool onlylb;			// to compute only lower bounds

  bool m2i, m2ibc, m3i, m3icg, m3if, corf; // models to solve
  bool m2iup;				   // use "up" constraints
  unsigned mcuts;      // maximum cuts in cut loops (m2ibc currently)
  unsigned ncols;      // number of columns per iteration (m3i via cg)
  unsigned nthreads;   // number of solver threads

  unsigned lb3w;		// number of workers for lb3
  bool showlb3;			// show values of lb3

  unsigned bsmin, bsmax;	// minimum and maximum beam size
  unsigned maxpasses; // maximum number of passes for cyclic best-first search

  unsigned nsamples;	      // number of samples in sampling methods
  unsigned pstrength;	      // perturbation strength in ILS

  double opt;			// optimal solution value (for tests with segments)
};

extern Options opt;

// process commandline options, return false, if some options are missing
bool process_options(int, char *[], Options&, po::variables_map&, bool=false);
