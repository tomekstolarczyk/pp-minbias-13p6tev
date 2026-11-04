#!/usr/bin/env bash
# Runs the whole pipeline: build, event generation, analysis and plots.
# ROOT and Pythia must be set up first (see README).

# stop at the first failing command
set -euo pipefail

# run from the repository root, wherever the script is called from
cd "$(dirname "$0")/.."

# the user sets up the tools (see README); the script only checks that they are available
for tool in cmake root pythia8-config; do
    if ! command -v "$tool" > /dev/null; then
        echo "Error: $tool not found in PATH, set up ROOT and Pythia first (see README)" >&2
        exit 1
    fi
done

echo "=== Build"
cmake -S . -B build
cmake --build build

mkdir -p data results

echo "=== Generation"
./build/generate config/minbias.cmnd data/events.root | tee results/generate.log

echo "=== Analysis"
./build/analysis data/events.root results/analysis.root | tee results/analysis.log

echo "=== Plots"
root -l -b -q 'macros/draw_plots.C("results/analysis.root", "plots")'

echo "=== Check against the task"
root -l -b -q tests/check_task.C

echo "=== Extra: comparison with ALICE (not part of the task)"
./build/generate config/minbias_alice_primaries.cmnd data/events_alice_primaries.root > results/generate_alice_primaries.log
root -l -b -q 'macros/compare_alice.C("data/events.root", "data/events_alice_primaries.root", "plots/alice_comparison")' | tee results/compare_alice.log

echo "=== Done"
