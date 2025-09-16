#!/bin/bash


docker run --rm -u $(id -u):$(id -g) \
       -v /data/experiments/g3dragon:/data/experiments/g3dragon \
       -w "$PWD" \
       ghcr.io/gchristian321/docker_g3dragon:v2 \
       bash -c "make dsbatch"
