#!/usr/bin/env bash

REPO=/mnt/archgen/poseidon/published_data

# 1000G
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $REPO \
  --noMinFreq --maxFreq 0.01 --popConfigFile pop_config_TGP_1240K.yml \
  -f ../Data/AncientBritish_1000G_1240K_ras01.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $REPO \
  --noMinFreq --maxFreq 0.02 --popConfigFile pop_config_TGP_1240K.yml \
  -f ../Data/AncientBritish_1000G_1240K_ras02.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $REPO \
  --noMinFreq --maxFreq 0.05 --popConfigFile pop_config_TGP_1240K.yml \
  -f ../Data/AncientBritish_1000G_1240K_ras05.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $REPO \
  --noMinFreq --maxFreq 0.1 --popConfigFile pop_config_TGP_1240K.yml \
  -f ../Data/AncientBritish_1000G_1240K_ras10.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $REPO \
  --noMinFreq --maxFreq 0.2 --popConfigFile pop_config_TGP_1240K.yml \
  -f ../Data/AncientBritish_1000G_1240K_ras20.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $REPO \
  --noMinFreq --noMaxFreq --popConfigFile pop_config_TGP_1240K.yml \
  -f ../Data/AncientBritish_1000G_1240K_rasAll.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $REPO \
  --minFreq 0.05 --maxFreq 0.95 --popConfigFile pop_config_TGP_1240K.yml \
  -f ../Data/AncientBritish_1000G_1240K_rasCommon.table.txt

# HGDP
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $REPO \
  --noMinFreq --maxFreq 0.02 --popConfigFile pop_config_HGDP_1240K.yml \
  -f ../Data/AncientBritish_HGDP_1240K_ras02.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $REPO \
  --noMinFreq --maxFreq 0.05 --popConfigFile pop_config_HGDP_1240K.yml \
  -f ../Data/AncientBritish_HGDP_1240K_ras05.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $REPO \
  --noMinFreq --maxFreq 0.1 --popConfigFile pop_config_HGDP_1240K.yml \
  -f ../Data/AncientBritish_HGDP_1240K_ras10.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $REPO \
  --noMinFreq --maxFreq 0.2 --popConfigFile pop_config_HGDP_1240K.yml \
  -f ../Data/AncientBritish_HGDP_1240K_ras20.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $REPO \
  --noMinFreq --noMaxFreq --popConfigFile pop_config_HGDP_1240K.yml \
  -f ../Data/AncientBritish_HGDP_1240K_rasAll.table.txt
qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d $REPO \
  --minFreq 0.05 --maxFreq 0.95 --popConfigFile pop_config_HGDP_1240K.yml \
  -f ../Data/AncientBritish_HGDP_1240K_rasCommon.table.txt

