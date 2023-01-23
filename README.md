# Exact and heuristic methods for a workload allocation problem with chained precedence constraints

This repository contains the instance, the source code, and additional information for the paper [Exact and heuristic methods for a workload allocation problem with chained precedence constraints](https://doi.org/10.1016/j.ejor.2022.12.035).

## Instances

The instances can be found in folder `dat`.

## Source code

The source can be found in folder `src`. It is written in C++ and requires the [Boost libraries](https://www.boost.org/), [CMake](https://cmake.org), CPLEX 20.1 and [libfmt](https://github.com/fmtlib/fmt) (libfmt is a submodule and will be cloned automatically). On some platforms, other libraries may be required (i.e., in Ubuntu 20.04, installing the zlib library was needed). 

### How to clone

```bash
git clone https://github.com/jordipereiragude/workload_allocation_chains.git
cd workerallocation-paper
git submodule update --init --recursive
```

### How to compile

```sh
mkdir build
cmake ../src
make -j
```

*Note*: if cmake can't find Boost, set BOOST_ROOT to the root of your local installation.

### How to run

From the build directory run, for example

```sh
./ipsolve --heuristic cbfs --maxpasses 20000 --exact sadpF --initial seg --lb3w 4 --showlb3 ../dat/Battarra,etal/AP001.txt
```

To see all options run
```sh
./ipsolve --help
```

### How to reproduce the computational results

To reproduce the tables run the script in `src/scripts/tables.R` in GNU R. This script processes log files from the `test` folder and produces figures in tables in `doc/figures` and `doc/tables`.

To reproduce the runs, have a look at the options in `test/testNNN/config.dat` and run the script `src/scripts/test.sh` with these options.

## How to cite

```bibtex
@Article{Pereira.Ritt/2022,
  author = 	 {Jorge Pereira and Marcus Ritt},
  title = 	 {Exact and heuristic methods for a workload allocation problem with chain precedence constraints},
  journal =  {Eur. J. Oper. Res},
  year = 	 {2023},
  doi = 	 {10.1016/j.ejor.2022.12.035}
}
```
