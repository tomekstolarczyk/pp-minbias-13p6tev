# pp-minbias-13p6tev

A complete pipeline for minimum-bias proton-proton collisions at √s = 13.6 TeV:
event generation with Pythia 8, analysis of final-state particles in C++/ROOT,
and visualization of the resulting distributions.

Author: Tomasz Stolarczyk

## Requirements

- ROOT 6, built with C++17
- Pythia 8
- C++17 compiler
- CMake 3.18 or newer

Tested on Ubuntu 22.04 LTS with ROOT 6.40.04, Pythia 8.318, CMake 4.4.5 and GCC 12.3.

*Identical histograms are expected only with the same Pythia version and settings.*

## Installing ROOT and Pythia

**ROOT:** follow [root.cern/install](https://root.cern/install/). Check the C++
standard of your build with `root-config --cflags`; it should contain `-std=c++17`.

**Pythia:** download the source from [pythia.org](https://pythia.org/) and build it
with C++17 (its default is C++11):

```bash
./configure --prefix=$HOME/pythia8 \
    --cxx-common='-O2 -std=c++17 -pedantic -W -Wall -Wshadow -fPIC -pthread'
make -j4 && make install
```

## Environment Setup

In every new shell, before building or running:

```bash
source <root-install-dir>/bin/thisroot.sh   # not needed if ROOT comes from a package manager
export PATH=$HOME/pythia8/bin:$PATH
```

Both `root-config --version` and `pythia8-config --version` should then work.

## Build

```bash
cmake -S . -B build
cmake --build build
```

## Event Generation

```bash
mkdir -p data results
./build/generate config/minbias.cmnd data/events.root | tee results/generate.log
```

All settings, including the event count and the random seed, are in
[config/minbias.cmnd](config/minbias.cmnd). 


