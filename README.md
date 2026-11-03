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

## Full Pipeline

```bash
scripts/run_all.sh
```

The script builds the programs and runs every step below in order
(about 2 minutes on the test machine). It overwrites `data/`, `results/` and `plots/`.

## Event Generation

```bash
mkdir -p data results
./build/generate config/minbias.cmnd data/events.root | tee results/generate.log
```

All settings, including the event count and the random seed, are in
[config/minbias.cmnd](config/minbias.cmnd).

The generator writes only final-state particles to `data/events.root`: a tree
`events` with one entry per event and the branches `pdg`, `charge`, `pt` (GeV),
`eta` and `phi` (rad), each holding one value per particle.

## Analysis

```bash
./build/analysis data/events.root results/analysis.root | tee results/analysis.log
```

The analysis writes 32 histograms to `results/analysis.root`: `hN_<cat>`,
`hPt_<cat>`, `hEta_<cat>` and `hPhi_<cat>` for eight categories of final-state particles.

| Category | Selection |
| --- | --- |
| `all` | every final-state particle, including neutral ones |
| `ch` | charge ≠ 0 |
| `pip`, `pim` | π⁺, π⁻ (PDG ±211) |
| `Kp`, `Km` | K⁺, K⁻ (PDG ±321) |
| `p`, `pbar` | p, p̄ (PDG ±2212) |

| Histogram | Bins | Range |
| --- | --- | --- |
| `hN` | 400 | [0, 400] |
| `hPt` | 100 | [0, 10] GeV |
| `hEta` | 100 | [−5, 5] |
| `hPhi` | 64 | [0, π] |


Values outside a range are kept in the underflow and overflow bins; `results/analysis.log` lists them for every histogram.

Two points of the task left room for interpretation:

- **φ range.** Pythia gives φ in (−π, π], the specified range is [0, π]. The histograms
  are filled with |φ|, which is equivalent because the distribution does not depend on φ.
- **`all` and `ch`.** The task names both categories, so `all` includes neutral particles
  and `ch` only charged ones. As a result 12.7% of the events have N ≥ 400 and fall into
  the overflow of `hN_all`, compared with 0.08% for `hN_ch`.




