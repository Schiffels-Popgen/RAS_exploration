#!/usr/bin/env bash

DATA_DIR=/mnt/archgen/MICROSCOPE/rarevar_sim_study/real_data
POS25_ANCIENT=$DATA_DIR/ancient/ancientBritish/hgdp/poseidon2.5
POS25_MODERN=$DATA_DIR/modern/HGDP/plink/poseidon/all_vars_2.5

qsub -b y -cwd xerxes ras -d $POS25_ANCIENT -d $POS25_MODERN --popConfigFile pop_config.yml -f AncientBritish_HGDP_ras.table.txt

