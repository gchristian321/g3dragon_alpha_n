#!/bin/bash


# Check if first argument ($1) is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <dsinit.sh file>"
    exit 1
fi

# If provided, continue with script
echo "Proceeding with dsinit file: $1"


docker run --rm -u $(id -u):$(id -g) \
       -v /data/experiments/g3dragon:/data/experiments/g3dragon \
       -w "$PWD" \
       ghcr.io/gchristian321/docker_g3dragon:v2 \
       bash -c "source ${1} && ./bin/dsbatch"
