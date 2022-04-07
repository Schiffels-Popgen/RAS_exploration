#!/usr/bin/env bash

DATA_DIR=/mnt/archgen/MICROSCOPE/rarevar_sim_study/real_data
POS_ANCIENT=$DATA_DIR/ancient/ancientBritish/hgdp/poseidon
POS_MODERN=$DATA_DIR/modern/HGDP/plink/poseidon/all_vars

# Rightpops ~200 individuals
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $POS_ANCIENT -d $POS_MODERN \
  -j 100000 --noMinFreq --maxFreq 0.02 --popConfigFile pop_config_HGDP.yml \
  -f ../Data/AncientBritish_HGDP_ras02.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $POS_ANCIENT -d $POS_MODERN \
  -j 100000 --noMinFreq --maxFreq 0.05 --popConfigFile pop_config_HGDP.yml \
  -f ../Data/AncientBritish_HGDP_ras05.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $POS_ANCIENT -d $POS_MODERN \
  -j 100000 --noMinFreq --maxFreq 0.1 --popConfigFile pop_config_HGDP.yml \
  -f ../Data/AncientBritish_HGDP_ras10.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $POS_ANCIENT -d $POS_MODERN \
  -j 100000 --noMinFreq --maxFreq 0.2 --popConfigFile pop_config_HGDP.yml \
  -f ../Data/AncientBritish_HGDP_ras20.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $POS_ANCIENT -d $POS_MODERN \
  -j 100000 --noMinFreq --noMaxFreq --popConfigFile pop_config_HGDP.yml \
  -f ../Data/AncientBritish_HGDP_rasAll.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $POS_ANCIENT -d $POS_MODERN \
  -j 100000 --minFreq 0.05 --maxFreq 0.95 --popConfigFile pop_config_HGDP.yml \
  -f ../Data/AncientBritish_HGDP_rasCommon.table.txt


DATA_DIR=/mnt/archgen/MICROSCOPE/rarevar_sim_study/real_data
POS_ANCIENT=$DATA_DIR/ancient/ancientBritish/1000g/poseidon
POS_MODERN=$DATA_DIR/modern/1000G/plink/poseidon/all_vars

# Rightpops ~800 individuals
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $POS_ANCIENT -d $POS_MODERN \
  -j 100000 --noMinFreq --maxFreq 0.01 --popConfigFile pop_config_TGP.yml \
  -f ../Data/AncientBritish_1000G_ras01.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $POS_ANCIENT -d $POS_MODERN \
  -j 100000 --noMinFreq --maxFreq 0.02 --popConfigFile pop_config_TGP.yml \
  -f ../Data/AncientBritish_1000G_ras02.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $POS_ANCIENT -d $POS_MODERN \
  -j 100000 --noMinFreq --maxFreq 0.05 --popConfigFile pop_config_TGP.yml \
  -f ../Data/AncientBritish_1000G_ras05.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $POS_ANCIENT -d $POS_MODERN \
  -j 100000 --noMinFreq --maxFreq 0.1 --popConfigFile pop_config_TGP.yml \
  -f ../Data/AncientBritish_1000G_ras10.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $POS_ANCIENT -d $POS_MODERN \
  -j 100000 --noMinFreq --maxFreq 0.2 --popConfigFile pop_config_TGP.yml \
  -f ../Data/AncientBritish_1000G_ras20.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $POS_ANCIENT -d $POS_MODERN \
  -j 100000 --noMinFreq --noMaxFreq --popConfigFile pop_config_TGP.yml \
  -f ../Data/AncientBritish_1000G_rasAll.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $POS_ANCIENT -d $POS_MODERN \
  -j 100000 --minFreq 0.05 --maxFreq 0.95 --popConfigFile pop_config_TGP.yml \
  -f ../Data/AncientBritish_1000G_rasCommon.table.txt

