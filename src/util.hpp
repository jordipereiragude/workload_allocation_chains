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

#include <chrono>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include <boost/any.hpp>

#define FMT_HEADER_ONLY
#include "fmt/format.h"
#include "fmt/ostream.h"
#include "fmt/ranges.h"
#include "fmt/color.h"

//macros for generating key for the map d is the data, n is the bit number
#define setbit(d,n)   ((d) |=  (1<<(n)))
#define unsetbit(d,n) ((d) &= ~(1<<(n)))
#define checkbit(d,n) ((d) &   (1<<(n)))

typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::duration<double> Duration;
typedef Clock::time_point Timepoint;

struct timer {
  using clock = std::chrono::high_resolution_clock;
  using time_point = clock::time_point;
  time_point tp_start;
  double tm_lim;
  timer(double timeLimitSecs, const timer& parent) {
    reset(std::min(timeLimitSecs, parent.secs_left()));
  }
  timer(double time_lim_secs = std::numeric_limits<double>::max()) {
    reset(time_lim_secs);
  }
  void reset(double time_lim_secs = std::numeric_limits<double>::max()) {
    tm_lim = time_lim_secs, tp_start = clock::now();
  }
  double elapsed_secs() const {
    return std::chrono::duration_cast<std::chrono::duration<double>>(clock::now()-tp_start).count();
  }
  double secs_left() const { return tm_lim - elapsed_secs(); }
  bool timed_out() const { return elapsed_secs() >= tm_lim; }
};

// global state
struct State {
  Timepoint start; // starting time
  unsigned iter;   // total number of iterations
  double elapsed() {
    return Duration(Clock::now()-start).count();
  }
};

extern State S;

struct options_counter { int count = 0; };
inline void validate(boost::any& v, std::vector<std::string> const& xs, options_counter*, long) {
  if (v.empty()) v = options_counter{1};
  else ++boost::any_cast<options_counter&>(v).count;
}
extern options_counter verbosec;

// helpers for verbose and timestamped output
template <typename S, typename... Args>
void vtprint(int level, std::ostream& os, const S& format_str, Args&&... args) {
  if (verbosec.count>=level) {
    fmt::print(os,"({:7.2f}) ",::S.elapsed());
    fmt::print(os,format_str,args...);
  }
}

template <typename S, typename... Args>
void vtprint(int level, const fmt::text_style& ts, const S& format_str, Args&&... args) {
  if (verbosec.count>=level) {
    fmt::print(ts,"({:7.2f}) ",::S.elapsed());
    fmt::print(ts,format_str,args...);
  }
}

template <typename S, typename... Args>
void vtprint(int level, const S& format_str, Args&&... args) {
  vtprint(level,std::cout,format_str,std::forward<Args>(args)...);
}

template <typename S, typename... Args>
void vprint(int level, std::ostream& os, const S& format_str, Args&&... args) {
  if (verbosec.count>=level)
    fmt::print(os,format_str,std::forward<Args>(args)...);
}

template <typename S, typename... Args>
void vprint(int level, const fmt::text_style& ts, const S& format_str, Args&&... args) {
  if (verbosec.count>=level)
    fmt::print(ts,format_str,std::forward<Args>(args)...);
}

template <typename S, typename... Args>
void vprint(int level, const S& format_str, Args&&... args) {
  if (verbosec.count>=level)
    fmt::print(format_str,std::forward<Args>(args)...);
}

template <typename S, typename... Args>
void tprint(std::ostream& os, const S& format_str, Args&&... args) {
  fmt::print(os,"({:7.2f}) ",::S.elapsed());
  fmt::print(os,format_str,args...);
}

template <typename S, typename... Args>
void tprint(const fmt::text_style& ts, const S& format_str, Args&&... args) {
  fmt::print(ts,"({:7.2f}) ",::S.elapsed());
  fmt::print(ts,format_str,args...);
}

template <typename S, typename... Args>
void tprint(const S& format_str, Args&&... args) {
  tprint(std::cout,format_str,std::forward<Args>(args)...);
}

std::string get_date_string(Timepoint t);
