#!/bin/bash -ex
perf script | ../FlameGraph/stackcollapse-perf.pl > out.perf-folded
../FlameGraph/flamegraph.pl --title=flamegraph out.perf-folded > flamegraph.svg
../FlameGraph/flamegraph.pl --title=flamegraph out.perf-folded --reverse > flamegraph_rev.svg
