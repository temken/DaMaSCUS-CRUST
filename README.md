[![Build Status](https://travis-ci.org/temken/DaMaSCUS-CRUST.svg?branch=master)](https://travis-ci.org/temken/DaMaSCUS-CRUST)
[![codecov](https://codecov.io/gh/temken/DaMaSCUS-CRUST/branch/master/graph/badge.svg)](https://codecov.io/gh/temken/DaMaSCUS-CRUST)


# DaMaSCUS - CRUST

Dark Matter Simulation Code for Underground Scatterings - Crust Edition

DaMaSCUS-CRUST Version 1.0 14/02/2018

<img src="https://user-images.githubusercontent.com/29034913/36166492-fa29a4fe-10f2-11e8-90c1-4b9fc3f96f2e.png" width="425">

## GENERAL NOTES

- DaMaSCUS-CRUST is a MC simulator of dark matter (DM) particle trajectories through the Earth crust, atmosphere and/or additional shielding layers undergoing elastic scatterings on nuclei. 
- It is a tool to systematically determine the critical cross-section of strongly interacting DM, above which a given experiment loses detection sensitivity. As a final result it returns a precise exclusion band.
- Users can easily change the layer structure in the config file and define layers by their density, thickness and composition.
- DaMaSCUS-CRUST is written in C++ and fully parallelized (openMPI).

For the underlying physics we refer to the corresponding paper [[arXiv:1802.04764]](https://arxiv.org/abs/1802.04764)

## CONTENT

The included folders are:

- bin: This folder contains the configuration files and the executables after compilation.
- build: This folder contains all object files.
- detectors: Detector specific data, such as tabulated efficiencies or recoil energy data are stored here.
- include: All header files can be found here. Necessary 3rd party libraries can also be placed here.
- results: The resulting limits are saved in this folder.
- src: Here you find the source code of DaMaSCUS-CRUST.


## INSTALLATION AND USAGE

### Dependencies:

Before the code can be compiled, three dependencies need to be taken care of:

- **libconfig v1.6**: For the input configuration files we use the *libconfig* library. Note that we need at least version 1.6. [(Link)](https://hyperrealm.github.io/libconfig/)
- **Eigen**: A linear algebra library for C++. [(Link)](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- **MPI**: For parallelization we use the Message Passing Interface [(Link)](https://www.open-mpi.org)

### Installation:

The code can be compiled using the makefile. It might be necessary to adjust the compiler lines and the path to the libraries

```
#Compiler and compiler flags
CXX := mpic++
CXXFLAGS := -Wall -std=c++11 
LIB := -lconfig++
INC := -I include
(...)
```

After that simply run
```
make
```
from the root directory in the terminal to compile DaMaSCUS-CRUST.

Running
```
make clean
```
deletes all objective files and executables.

### Usage:

After successful compilation there is only one final step before running the program: Adjusting the configuration file in */bin/*. This file specifies all the input parameter, the mass scan range, the experiment of interest, the shielding layer structure,... . This should be rather self-explanatory. However note that *libconfig* is not very forgiving with regards to the input parameter type. For example it might complain if an input parameter of type double is given as “1” instead of “1.0”.

Afterwards the program starts by

```
mpirun -n N DaMaSCUS-CRUST config.cfg
```
where N is the number of MPI processes.

After a successful run, the resulting constraints can be plotted with the included *Mathematica* notebook */Plot.nb*.

## CITING DaMaSCUS

If you decide to use this code, please cite

<!--- >Emken, T. & Kouvaris, C., 2017, DaMaSCUS, Astrophysics Source Code Library, record [ascl:xx](link)--->

<!---as well as the original publication, --->

>Emken, T. & Kouvaris, C., How blind are underground and surface detectors to strongly interacting Dark Matter?, (2018), [[arXiv:1802.04764]](https://arxiv.org/abs/1802.04764).

## AUTHORS & CONTACT

The authors of DaMaSCUS-CRUST are Timon Emken and Chris Kouvaris.

For questions, bug reports or other suggestions please contact Timon Emken (emken@cp3.sdu.dk).


## LICENCE

This project is licensed under the MIT License - see the LICENSE file.
