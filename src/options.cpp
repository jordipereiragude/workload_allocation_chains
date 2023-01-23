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
#include "options.hpp"

#include <iostream>
#include <cassert>
#include <cstdlib>
using namespace std;

#include <sys/ioctl.h>

#include "util.hpp"
#include "version.hpp"
#include "dimensions.hpp"

Options opt;

bool process_options(int argc, char *argv[], Options& opt, po::variables_map& vm, bool only_standard) {
  struct winsize wsize;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &wsize);

  po::options_description desc("General options",wsize.ws_col);
  desc.add_options()
    ("help",                                                                   "Show help.")
    ("version",                                                                "Show version.")
    ("verbose,v",      po::value(&verbosec)->zero_tokens(),                    "Verbosity. If present, output is sent to screen. If -v is repeated, more output is given.")
    ("instance",       po::value<string>(&opt.instance),                       "Instance.")
    ("format",         po::value<string>(&opt.format)->default_value("guess"), "Input format, either Standard (std), ALWABP (alwabp) or guess based on the extension.")
    ("seed",           po::value<unsigned>()->default_value(1),                "Random seed (0: get from random source).")
    ("batch",          po::value<unsigned>(&opt.B)->default_value(0),          "Batch size (0: default).")
    ("timelimit",      po::value<unsigned>(&opt.tlim)->default_value(600),     "Time limit (s).")
    ;

  po::options_description all("Options",wsize.ws_col);
  all.add(desc);

  if (!only_standard) {
    po::options_description algo("Algorithmic options",wsize.ws_col);
    algo.add_options()
      ("initial",   po::value<string>(&opt.initial)->default_value("seg"),         "Initial heuristic, rnd (random), sample, seg (segments), ils (iterated local search), or test.")
      ("heuristic", po::value<string>(&opt.heuristic)->default_value("bs"),        "Heuristic algorithm, either cbfs (CBFS), bs (beam search), rnd (random), or cmb (combined).")
      ("exact",     po::value<string>(&opt.exact)->default_value("ip"),            "Exact algorithm, either ip (integer programming), dp (dynamic programming), sadpF (successive approximation, dp, feasibility), sadpFS (ditto, with tails), sadpO (ditto, optimization version), cbfsF (exact version of cbfs, solves feasibility problems), cbfsO (exact version of cbfs, solves optimization problem), cbfsDfsF (exact feasibily cbfs/dfs hybrid), cbfsDfsO (exact optimization cbfs/dfs hybrid).")
      ("sadpMem",   po::value<unsigned>(&opt.sadpmemory)->default_value(UINT_MAX), "Maximum memory sadp (in table entries).")
      ("sadpFMem",  po::value<unsigned>(&opt.sadpFmemory)->default_value(2048),    "Maximum memory sadp (in MB).")
      ("cbfsMem",   po::value<unsigned>(&opt.cbfsmemory)->default_value(UINT_MAX), "Maximum memory cbfs (in table entries).")
      ("onlylb",    po::bool_switch(&opt.onlylb)->default_value(false),            "Compute only lower bounds.")
      ("lb3w",      po::value<unsigned>(&opt.lb3w)->default_value(0),              "Number of workers for LB3.")
      ("showlb3",   po::bool_switch(&opt.showlb3)->default_value(false),           "Show details of LB3.")
      ("opt",       po::value<double>(&opt.opt)->default_value(0),                 "Optimal solution value.")
      ;
  
    po::options_description model("Model options",wsize.ws_col);
    model.add_options()
      ("m2i",    po::bool_switch(&opt.m2i)->default_value(false),      "Run model M2I.")
      ("m2ibc",  po::bool_switch(&opt.m2ibc)->default_value(false),    "Run model M2I via Branch&Cut.")
      ("m3i",    po::bool_switch(&opt.m3i)->default_value(false),      "Run model M3I.")
      ("m3icg",  po::bool_switch(&opt.m3icg)->default_value(false),    "Run model M3I via Column generation.")
      ("m3if",   po::bool_switch(&opt.m3if)->default_value(false),     "Run model M3I feasibility model.")
      ("corf",   po::bool_switch(&opt.corf)->default_value(false),     "Run model CORF.")
      ("m2iup",  po::bool_switch(&opt.m2iup)->default_value(false),    "Use up constraints in M2I.")
      ("mcuts",  po::value<unsigned>(&opt.mcuts)->default_value(5),    "Maximum number of cuts in M2I cut loop.")
      ("ncols",  po::value<unsigned>(&opt.ncols)->default_value(200),  "Number of columns to add per iteration in M3I via CG.")
      ("nthread",po::value<unsigned>(&opt.nthreads)->default_value(1), "Number of solver threads.")
      ;
  
    po::options_description heur("Heuristic options",wsize.ws_col);
    heur.add_options()
      ("nsamples",  po::value<unsigned>(&opt.nsamples)->default_value(1),  "Number of samples in sample methods.")
      ("pstrength", po::value<unsigned>(&opt.pstrength)->default_value(2), "Perturbation strength in ILS.")
      ;
    
    po::options_description beam("Beam search options",wsize.ws_col);
    beam.add_options()
      ("bsmin",     po::value<unsigned>(&opt.bsmin)->default_value(25),       "Minimum beam size.")
      ("bsmax",     po::value<unsigned>(&opt.bsmax)->default_value(500),      "Maximum beam size.")
      ("maxpasses", po::value<unsigned>(&opt.maxpasses)->default_value(400),  "Maximum number of passes for cyclic best-first search.")
      ;
  
    all.add(algo).add(model).add(heur).add(beam);
  }
  
  po::positional_options_description pod;
  pod.add("instance", 1);

  po::store(po::command_line_parser(argc, argv).options(all).positional(pod).run(), vm);
  po::notify(vm);
  
  if (vm.count("version")) {
    cout << version << endl;
    return 0;
  }
  
  if (vm.count("help") || !vm.count("instance")) {
    if (!vm.count("instance"))
      cout << "Instance is missing." << endl;
    cout << all << endl;
    return false;
  }

  return true;
}
