#!/usr/bin/env bash

# DATA_DIR=/mnt/archgen/MICROSCOPE/rarevar_sim_study/real_data
# POS25_ANCIENT=$DATA_DIR/ancient/ancientBritish/hgdp/poseidon
# POS25_MODERN=$DATA_DIR/modern/HGDP/plink/poseidon/all_vars

# qsub -b y -cwd -l h_vmem=16G xerxes ras -d $POS25_ANCIENT -d $POS25_MODERN -j 100000 --popConfigFile pop_config.yml -f AncientBritish_HGDP_ras.table.txt

DATA_DIR=/mnt/archgen/MICROSCOPE/rarevar_sim_study/real_data
POS_ANCIENT=$DATA_DIR/ancient/ancientBritish/1000g/poseidon
POS_MODERN=$DATA_DIR/modern/1000G/plink/poseidon/all_vars

qsub -b y -cwd -l h_vmem=16G xerxes ras -d $POS_ANCIENT -d $POS_MODERN -j 100000 --popConfigFile pop_config_TGP.yml -f AncientBritish_1000G_ras.table.txt
