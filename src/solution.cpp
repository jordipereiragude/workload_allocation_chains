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
#include "util.hpp"
#include "ff.hpp"
#include "solution.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
using namespace std;

string Solution::to_string() const {
  string s = "\n";

  if (!isValid())
    return s + "Infeasible  cycle Inf\n";

  const int wd = 5;
  for(int k=0,ke=w.size(); k!=ke; k++)
    s += fmt::format("{:{}} ", k+1, wd);
  s += "\n";
  for(int k=0,ke=w.size(); k!=ke; k++)
    s += fmt::format("{}", string(wd+1,'-'));
  s += "\n";
  for(int k=0,ke=w.size(); k!=ke; k++)
    s += fmt::format("{:{}} ", w[k]+1, wd);
  s += "\n";
  int nt = 0;
  for(int k=0,ke=w.size(); k!=ke; k++) {
    s += fmt::format("{:{}} ", l[k], wd);
    nt+=l[k];
  }
  s += "\n";
  nt = 0;
  for(int k=0,ke=w.size(); k!=ke; k++) {
    Time sc = I->segCost(nt,nt+l[k]-1,w[k]);
    if (ff(sc)==ff(value))
      s += fmt::format("{:{}.2f} ", fmt::styled(sc, fmt::fg(fmt::color::red)), wd);
    else
      s += fmt::format("{:{}.2f} ", sc, wd);
    nt+=l[k];
  }
  s += fmt::format("  cycle {:.2f}\n", value);
  return s;
}

ostream& operator<<(ostream &o, const Solution& s) {
  return o << s.to_string();
}
