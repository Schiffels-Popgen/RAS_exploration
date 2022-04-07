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
        -- afStr <- select ["01", "02", "05", "10", "20", "Common", "All"]
        afStr <- select ["02", "05", "10", "20", "Common", "All"]
        mapMasked <- select [False, True]
        let dataset = "HGDP"
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
        let bedStr = if mapMasked then "_mapMasked" else ""
        let outFN =
                if dataset == "1000G" then
                        format("AncientBritish_1000G_ras"%s%"_TVonly"%s%".table.txt") afStr bedStr
                else
                        format("AncientBritish_HGDP_ras"%s%"_TVonly"%s%".table.txt") afStr bedStr
        let base_dirs =
                if dataset == "1000G" then
                        format ("-d "%fp%" -d "%fp) pos_ancient pos_modern
                else
                        format ("-d "%fp%" -d "%fp) pos_ancient_hgdp pos_modern_hgdp
        let config_file = if dataset == "1000G" then "pop_config_TGP.yml" else "pop_config_HGDP.yml"
        liftIO . process $ format ("qsub -V -b y -cwd -l h_vmem=16G xerxes ras "%s%
            " -j 100000 --noTransitions "%s%" --popConfigFile "%s%" "%s%
            "-f ../Data/"%s) base_dirs afCond config_file bedOpt outFN

  