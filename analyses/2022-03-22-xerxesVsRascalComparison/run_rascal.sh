#!/usr/bin/env bash

cat all_vars_m1.0_chr1.freqsum | RASCalculator.py -M 5 -o Ref -a ind0,ind1,ind2,ind3 -c ind0,ind1,ind2,ind3 -L ind0,ind1,ind2,ind3 --skipJackknife