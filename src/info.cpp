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

#define FMT_HEADER_ONLY
#include "fmt/core.h"

#include "util.hpp"
#include "options.hpp"
#include "random.hpp"
#include "instance.hpp"
#include "ip.hpp"
#include "hs.hpp"
#include "beamSearch.hpp"
#include "reducedIP.hpp"
#include "lb.hpp"
#include "dp.hpp"
#include "sadp.hpp"

using namespace wap;

//// globals
string iname;

void all_times(const Instance& I) {
  set<Time> all_times;
  for(auto k=0u; k!=I.m; ++k) {
    set<Time> w_times;
    for(auto i=0u; i!=I.n; ++i)
      for(auto j=i; j!=I.n; ++j)
	w_times.insert(I.segCost(i,j,k));
    fmt::print("Worker {}: {} times (out of {}).\n", k, w_times.size(), I.n*I.n);
    all_times.insert(w_times.begin(),w_times.end());
  }
  fmt::print("Overall: {} times (out of {}).\n", all_times.size(), I.n*I.n*I.m);
}

int main(int argc, char *argv[]) {
  // (0) process commandline
  S.start = Clock::now();
  S.iter = 0;

  po::variables_map vm;
  if (!process_options(argc,argv,opt,vm,true))
    return 1;
  
  // (1) random seed
  [[maybe_unused]] unsigned seed = setupRandom(vm["seed"].as<unsigned>());
  // (2) read the instance
  iname = string(fs::path(opt.instance).parent_path().filename())+"/"+string(fs::path(opt.instance).stem());

  vtprint(1,fg(fmt::color::green),"Read input {}.\n",iname);
  
  Instance I;
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
    return 1;
  } catch (const string s) {
    cerr << "Failed to read instance: " << s << endl;
    return 1;
  } catch (...) {
    cerr << "Failed to read instance." << endl;
    return 1;
  }
  in.close();

  vtprint(1,fg(fmt::color::green),"Instance read.\n",iname);
  vprint(1,"Instance has {} tasks, and {} workers. Batch size is {}.\n",I.n,I.m,I.B);

  fmt::print("{} {} {} {} {}\n", iname, I.n, I.m, I.B, I.sumStdTimes());
  //all_times(I);
}
