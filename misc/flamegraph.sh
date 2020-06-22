#!/bin/bash -ex
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
perf script | $DIR/stackcollapse-perf.pl > out.perf-folded
$DIR/flamegraph.pl --title=flamegraph out.perf-folded > flamegraph.svg
$DIR/flamegraph.pl --title=flamegraph out.perf-folded --reverse > flamegraph_rev.svg
