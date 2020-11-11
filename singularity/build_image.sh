#!/usr/bin/env bash

cd $(mktemp -d)
sudo singularity build reat.img REAT-Singularity.def