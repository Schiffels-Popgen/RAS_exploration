#!/usr/bin/env bash

DATA_DIR=/mnt/archgen/MICROSCOPE/rarevar_sim_study/real_data
POS_ANCIENT=$DATA_DIR/ancient/ancientBritish/hgdp/poseidon #63341718 SNPs
POS_MODERN=$DATA_DIR/modern/HGDP/plink/poseidon/all_vars #63341718 SNPs

# TEMPORARY: Update to Poseidon 2.5

POS25_ANCIENT=$DATA_DIR/ancient/ancientBritish/hgdp/poseidon2.5
# cp -r $POS_ANCIENT $POS25_ANCIENT
# cat $POS_ANCIENT/Hinxton.hgdp.janno | cut -f1,19,20 | awk -v OFS="\t" '{if($1=="Individual_ID")$1="Poseidon_ID"}{print}' > $POS25_ANCIENT/Hinxton.hgdp.janno
# trident update -d $POS25_ANCIENT --logText "updated to Poseidon v2.5" --poseidonVersion "2.5.0"

POS25_MODERN=$DATA_DIR/modern/HGDP/plink/poseidon/all_vars_2.5
# cp -r $POS_MODERN $POS25_MODERN
# cat $POS_MODERN/hgdp_hg19.all_vars.janno | cut -f1,19,20 | awk -v OFS="\t" '{if($1=="Individual_ID")$1="Poseidon_ID"}{print}' > $POS25_MODERN/hgdp_hg19.all_vars.janno
# trident update -d $POS25_MODERN --logText "updated to Poseidon v2.5" --poseidonVersion "2.5.0"


