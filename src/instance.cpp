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
#include "instance.hpp"

using namespace std;
using namespace boost;

#include "util.hpp"
#include "options.hpp"

double Instance::segCost(int i, int j, int w) const {
  return T[j+1][w]-T[i][w];
}

double Instance::segCost(int i, int j) const {
  return T_[j+1]-T_[i];
}

void Instance::readSTD(std::istream& in) {
  string line;
  in >> n; getline(in,line);
  in >> m; getline(in,line);
  getline(in,line);

  s0.resize(n);
  for(unsigned i=0; i!=n; i++) {
    in >> s0[i];
    vprint(3,"{} ",s0[i]);
  }
  vprint(3,"\n");
  B = opt.B!=0 ? opt.B : getB();
  vprint(3,"Batch size {}.\n",B);

  s.resize(n*B);
  for(unsigned i=0; i!=n; i++)
    fill(s.begin()+i*B,s.begin()+i*B+B,s0[i]);

  getline(in,line);
  getline(in,line);
  t.resize(extents[n*B][m]);
  for(unsigned w=0; w!=m; w++) {
    for(unsigned i=0; i!=n; i++) {
      Time ct;
      in >> ct;
      for(unsigned j=0; j!=B; j++)
	t[i*B+j][w]=ct;
      vprint(3," {}",t[i*B][w]);
    }
    vprint(3,"\n");
  }
  n *= B;
  computeD();
  //Malloc (maybe substitute by multi_array)
}

vector<Time> Instance::getTask(istream& in) {
  vector<Time> result;

  string s;
  getline(in,s);
  istringstream line(s);

  while (!line.eof()) {
    string token;
    line >> token;
    if (token.size()==0) // for DOS mode, since last element is CR
      continue;
    if (token[0]=='I' || token[0]=='i') {
      result.push_back(inf_time);
      ninf++;
    } else {
      assert(stod(token) < inf_time);
      if (stod(token)<1)
	cerr << "WARNING: Instance contains a 0-time task. Local search may fail." << endl;
      result.push_back(stod(token));
    }
  }
  return result;
}

void Instance::computeD() {
  T.resize(extents[boost::multi_array_types::extent_range(0,n+1)][m]);
  T_.resize(n+1);
  d.resize(extents[boost::multi_array_types::extent_range(-1,n)][m]);
  db.resize(extents[boost::multi_array_types::extent_range(-1,n)][m]);
  for(unsigned j=0;j!=m;j++) { //for each worker
    T[0][j]=0.0;
    for(unsigned i=1;i<=n;i++) // accumulate over tasks
      T[i][j]=T[i-1][j]+t[i-1][j];
  }
  T_[0]=0.0;
  for(unsigned i=1;i<=n;++i) {
    double min_time = t[i-1][0];
    for(unsigned j=1;j!=m;++j)
      min_time=min(min_time,t[i-1][j]);
    T_[i]=T_[i-1]+min_time;
  }
  for(unsigned j=0;j!=m;j++) { //for each worker
    d[-1][j]=0;
    for(int i=0;i<int(n);i++) { // accumulate over tasks
      d[i][j]=d[i-1][j]+t[i][j]*precision;
      assert(d[i][j]>d[i-1][j]);
    }
  }
  for(unsigned j=0;j!=m;j++) { //for each worker
    db[-1][j]=0;
    for(int i=(n-1);i>=0;i--)
      db[int(n)-1-i][j]=db[int(n)-2-i][j]+t[i][j]*precision;
  }
}

double Instance::normtime(double t) {
  double r = round(1.0e6*t)/1.0e6;
  return r<0.01?0.01:r;
}

void Instance::readALWABP(std::istream& in) {
  in >> n;
  string dummy;
  getline(in,dummy);

  // cache ALWABP times
  multi_array<Time,2> u;
  for(unsigned i=0; i!=n; i++) {
    auto t_ = getTask(in);
    if (i==0) { // first task defines number of workers
      m = t_.size();
      u.resize(extents[n][m]);
    }
    else if (t_.size()!=size_t(m))
      throw fmt::format("Number of workers {} in line {}  inconsistent with first line {}",t_.size(),i+1,m);
    for(unsigned w=0; w!=m; w++)
      u[i][w]=t_[w];
  }

  // dependencies are ignored

  // first worker defines standard times
  s.resize(n);

  for(unsigned i=0; i!=n; i++)
    s[i]=u[i][0];

  double scale = 0.4*double(m)/accumulate(s.begin(),s.end(),0.0);
  vprint(1,"ALWABP scaling factor {}.\n",scale);
  for(unsigned i=0; i!=n; i++) {
    s[i]=normtime(s[i]*scale);
    vprint(3,"{} ",s[i]);
  }
  s0=s;

  B = opt.B!=0 ? opt.B : getB();
  vprint(3,"Batch size {}.\n",B);

  vector<Time> s_(n*B);
  for(unsigned i=0; i!=n; i++)
    fill(s_.begin()+i*B,s_.begin()+i*B+B,s[i]);
  s.swap(s_);

  // assign times
  t.resize(extents[n*B][m]);
  for(unsigned w=0; w!=m; w++) {
    for(unsigned i=0; i!=n; i++) {
      for(unsigned j=0; j!=B; j++) {
	if (u[i][w]!=inf_time)
	  t[i*B+j][w]=normtime(u[i][w]*scale);
	else
	  t[i*B+j][w]=normtime(10*u[i][0]*scale);
      }
    }
  }
  n *= B;

  computeD();
}

int Instance::getB() {
  Time sum = sumStdTimes();
  return floor(double(60*m)/sum);
}

double Instance::getMinTime() {
  Time res = inf_time;
  for(unsigned w=0; w!=m; w++)
    for(unsigned i=0; i!=n; i++)
      res = min(res,t[i][w]);
  return res;
}

Time Instance::sumStdTimes() const {
  Time sum = accumulate(s0.begin(),s0.end(),0.0);
  return sum;
}
