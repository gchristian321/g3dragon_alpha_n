#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: $0 <dsinit.sh file> [<numthreads=4>]"
    exit 1
fi

if [ -z "$2" ]; then
    python docker_run_dsbatch_parallel.py $1
else
    python docker_run_dsbatch_parallel.py $1 $2
fi
