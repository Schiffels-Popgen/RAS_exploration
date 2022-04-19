#!/usr/bin/env stack
--stack script --resolver=lts-19.0 --package turtle
{-# LANGUAGE OverloadedStrings #-}

import Control.Monad (forM_)
import Turtle
import Prelude hiding (FilePath)

data_dir :: FilePath
data_dir = "/mnt/archgen/MICROSCOPE/rarevar_sim_study/real_data"

pos_ancient :: FilePath
pos_ancient = data_dir </> "ancient/ancientBritish/1000g/poseidon"

pos_modern :: FilePath
pos_modern = data_dir </> "modern/1000G/plink/poseidon/all_vars"

bed_file :: FilePath
bed_file = "/mnt/archgen/users/schiffels/hs37m_filt35_99.bed"

repo_dir :: FilePath
repo_dir = "/mnt/archgen/poseidon/published_data"

pos_ancient_hgdp :: FilePath
pos_ancient_hgdp = data_dir </> "ancient/ancientBritish/hgdp/poseidon"

pos_modern_hgdp :: FilePath
pos_modern_hgdp = data_dir </> "modern/HGDP/plink/poseidon/all_vars"

main :: IO ()
main = do
    let process cmd = do
            print cmd
            stdout (inshell cmd empty)
    sh $ do
        afStr <- select ["01", "02", "05", "10", "20", "Common", "All"]
        mapMasked <- select [False, True]
        dataset <- select ["1000G", "HGDP"]
        tf <- select [False, True]
        tvOnly <- select [False, True]
        let afCond = case afStr of
                "01"     -> "--noMinFreq    --maxFreq 0.01"
                "02"     -> "--noMinFreq    --maxFreq 0.02"
                "05"     -> "--noMinFreq    --maxFreq 0.05"
                "10"     -> "--noMinFreq    --maxFreq 0.1"
                "20"     -> "--noMinFreq    --maxFreq 0.2"
                "Common" -> "--minFreq 0.05 --maxFreq 0.95"
                "All"    -> "--noMinFreq    --noMaxFreq"
        let bedOpt = if mapMasked then
                format ("--bedFile "%fp%" ") bed_file else " "
        let base_dirs = case (dataset, tf) of
                ("1000G", False) -> format ("-d "%fp%" -d "%fp) pos_ancient pos_modern
                ("HGDP", False) -> format ("-d "%fp%" -d "%fp) pos_ancient_hgdp pos_modern_hgdp
                (_, True) -> format ("-d "%fp) repo_dir
        let config_file = case (dataset, tf) of
                ("1000G", False) -> "pop_config_TGP.yml"
                ("HGDP", False) -> "pop_config_HGDP.yml"
                ("1000G", True) -> "pop_config_TGP_1240K.yml"
                ("HGDP", True) -> "pop_config_HGDP_1240K.yml"
        let tvCmd = if tvOnly then "--noTransitions" else ""
        let bedStr = if mapMasked then "_mapMasked" else ""
        let tfStr = if tf then "_1240K" else ""
        let tvStr = if tvOnly then "_TVonly" else ""
        let outFN = fromText $ format("AncientBritish_"%s%s%"_ras"%s%s%s%".table.txt") dataset tfStr afStr tvStr bedStr
        e <- testfile ("../Data" </> outFN)
        when (not e) . liftIO . process $ 
            format ("qsub -V -b y -cwd -l h_vmem=16G xerxes ras "%s%
                " -j 100000 "%s%" "%s%" --popConfigFile "%s%" "%s%
                "-f ../Data/"%fp) base_dirs tvCmd afCond config_file bedOpt outFN

  