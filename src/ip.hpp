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

#include <tuple>
#include <map>

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include "instance.hpp"
#include "solution.hpp"

void setupSolver(IloCplex solver);
IloCplex solveRelaxed(IloModel model, IloNumVarArray x);

std::string to_string(IloAlgorithm::Status status);

using Interval = std::tuple<unsigned,unsigned,unsigned>; // start, end, worker

// two-index model, with continuity
struct M2I {
  const Instance& I;

  IloEnv env;
  IloModel model;
  IloObjective obj;
  IloNumVarArray x,u;

  M2I(const Instance& I) : I(I), model(env), obj(IloAdd(model, IloMinimize(env))), x(env), u(env) {}

  ~M2I() {
    env.end();
  }
  void addVars(bool);		// variables
  void addCycle();		// cycle time constraints
  void addAssignment();		// assignment constraints
  void addContinuity();		// continuity constraints
  void addContinuity2();	// continuity constraints

  void buildBase(bool);		// omit full continuity constraints
  void buildFull(bool);

  void setSolution(IloCplex solver, const Solution& S);
  Solution getSolution(IloCplex solver);
  void showSolution(IloCplex solver);

  double solveRelaxed();	       // relax full model, solve it
  std::pair<Solution,IloAlgorithm::Status> solve(const Solution& S);   // solve full model
  std::pair<Solution,IloAlgorithm::Status> solveBC(const Solution& S); // solve via B&C
};

// three-index segment model
struct M3I {
  const Instance& I;

  IloEnv env;
  IloModel model;
  IloObjective obj;
  IloNumVarArray y;
  IloRangeArray cycle, assign, worker;

  int idx_c;
  std::map<std::tuple<unsigned,unsigned,unsigned>,int> idx;

  M3I(const Instance& I) : I(I), model(env), obj(IloAdd(model, IloMinimize(env))), y(env), cycle(env), assign(env), worker(env) {}

  ~M3I() {
    env.end();
  }

  double addVars(double,bool=true);  // variables, optional: disable cycle time
  void   addCycle();		   // cycle time constraints
  void   addAssignment();		   // assignment constraints
  void   addWorker();		   // worker constraints

  void   buildCore(const Solution& S); // subset of columns in `S`
  void   buildFull(const Solution& S, bool=false); // using upper bound `S`, optional: feasibility problem
  double buildFull(double ub, bool=false); // using upper bound `S`, optional: feasibility problem

  void setSolution(IloCplex solver, const Solution& S);
  Solution getSolution(IloCplex solver);
  void showSolution(IloCplex solver);

  double solveRelaxed();		// relax full model, solve it
  bool solveRelaxedFeasibility();	// relax feasibility model, solve it
  double solveRelaxedCG(double,double=0.0); // solve relaxation via CG with upper bound, and optional lower bound
  void addColumn(unsigned i, unsigned j, unsigned w);
  void addColumn(Interval iv);

  std::pair<Solution,IloAlgorithm::Status> solvedTenseCG();	     // solve IP with columns from CG
  std::pair<Solution,IloAlgorithm::Status> solve(const Solution& S); // solve full IP model, with MIP start `S`
};

// CORF
struct Corf {
  const Instance& I;

  Corf(const Instance& I) : I(I), model(env), obj(IloAdd(model, IloMinimize(env))), p(env), q(env), r(env), c(env), cycle(env), assign(env), continuity(env), packing(env), linking(env), fringe(env) {}

  ~Corf() {
    env.end();
  }

  IloEnv env;
  IloModel model;
  IloObjective obj;
  IloNumVarArray p,q,r;
  IloNumVar c;
  IloRangeArray cycle, assign, continuity, packing, linking, fringe;

  void addVars();
  void addCycle();
  void addAssignment();
  void addContinuity();
  void addPacking();
  void addLinking();
  void addFringe();
  
  void build();					  // build model
  std::pair<double,IloAlgorithm::Status> solve(); // solve it
};
