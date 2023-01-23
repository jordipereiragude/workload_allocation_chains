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

/*
 * Generate instances according to Battarra et al. (2020).
 * 1. Sample average operation times from N(0.4,0.2)
 * 2. For average performance (AP): add N(0,0.5) limit to at most 1.0
 * 3. For mixed performance (MP):
 *    a. Select a skill level l\in{-1,0,1}
 *    b. add N(l/10,0.3+0.2[l==0]), limit to at most 1.0
 */

#include <iostream>
#include <cassert>
#include <cstdlib>
using namespace std;

#include <sys/ioctl.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#define FMT_HEADER_ONLY
#include "fmt/core.h"

#include "util.hpp"
#include "random.hpp"
using namespace wap;

struct Options {
  unsigned n, m; // number of operations and workers
  string type;	 // instance type: AP or MP
};

// mean and variance of basic ("estimated") times
const double est_avg = 0.4;
const double est_var = 0.2;
normal_distribution<double> est(est_avg,est_var);

// upper bound on times
const double max_time = 0.1;

// variation of AP (average performance) instances
normal_distribution<double> ap(0.0,0.5);

// variation for MP (mixed performance) instances
vector<normal_distribution<double>> mp = {
	   normal_distribution<double>(-0.1,0.3), // level 1 in paper
	   normal_distribution<double>( 0.0,0.5), // level 0 in paper
	   normal_distribution<double>( 0.1,0.3)  // level 2 in paper
};

bool process_options(int argc, char *argv[], Options& opt, po::variables_map& vm) {
  struct winsize wsize;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &wsize);

  po::options_description desc("Options",wsize.ws_col);
  desc.add_options()
    ("help",                                                              "Show help.")
    ("version",                                                           "Show version.")
    ("verbose,v",      po::value(&verbosec)->zero_tokens(),               "Verbosity. If present, output is sent to screen. If -v is repeated, more output is given.")
    ("type",           po::value<string>(&opt.type)->default_value("ap"), "Instance type, either average performance (ap), or mixed performance (mp).")
    ("seed",           po::value<unsigned>()->default_value(1),           "Random seed (0: acquire from random source).")
    ("n",              po::value<unsigned>(&opt.n),                       "Number of operations.")
    ("m",              po::value<unsigned>(&opt.m),                       "Number of workers.")
    ;

  po::options_description all("Options",wsize.ws_col);
  all.add(desc);

  po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    cout << all << endl;
    return false;
  }
  return true;
}

Options opt;

int main(int argc, char *argv[]) {
  // (0) process commandline
  po::variables_map vm;
  if (!process_options(argc,argv,opt,vm))
    return 1;

  // (1) random seed
  [[maybe_unused]] unsigned seed = setupRandom(vm["seed"].as<unsigned>());

  // (2) generate instance
  fmt::print("{} Number of operations\n",opt.n);
  fmt::print("{} Number of workers\n",opt.m);
  fmt::print("Standard production times\n");

  // (2.1) produce estimated times
  vector<double> tm(opt.n,0.0);
  generate(tm.begin(),tm.end(),[&]() { return max(est(rng),0.1); });

  // (2.2) define worker skill levels (only used for MP instances)
  vector<unsigned> ws(opt.m,0);
  generate(ws.begin(),ws.end(),[]() { return getRandom(0,2); });

  // (3) output instance
  // (3.1) output basic times
  for(auto i=0ul, ie=tm.size(); i!=ie; ++i)
    fmt::print("{:.6f} ",tm[i]);
  fmt::print("\n");
  // (3.2) output worker-task time matrix
  fmt::print("Production times: workers row and ops columns\n");
  for(unsigned w=0; w!=opt.m; ++w) {
    for(unsigned i=0; i!=opt.n; ++i) {
      double t = tm[i];
      t += opt.type=="ap" ? ap(rng) : mp[ws[w]](rng);
      t = max(t,max_time);
      fmt::print("{:.6f} ",t);
    }
    fmt::print("\n");
  }
}
