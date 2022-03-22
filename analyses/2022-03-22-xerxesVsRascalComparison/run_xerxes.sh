#!/usr/bin/env bash

# xerxes ras -d all_vars_chr1 --noMinFreq --maxAC 5 --popConfigFile popConfig.yml -f xerxes_ras_indsVsPops.txt
# xerxes ras -d all_vars_chr1 --noMinFreq --maxAC 5 --popConfigFile popConfig2.yml -f xerxes_ras_popsVsPops.txt
xerxes ras -d all_vars_chr1 -j CHR --minAC 2 --maxAC 5 --popConfigFile popConfigSimple.yml -f xerxes_ras_simple.txt