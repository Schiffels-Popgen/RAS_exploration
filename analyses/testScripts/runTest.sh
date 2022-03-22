#!/usr/bin/env bash

REPO=~/dev/poseidon-framework/published_data

# xerxes ras --maxSnps 100000 --noMaxAlleleCount -d $REPO --popConfigFile popConfig2.yml -f testOut
xerxes fstats --maxSnps 100000 -d $REPO --statFile fstats.txt
