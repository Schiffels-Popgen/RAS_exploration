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

main :: IO ()
main = do
    let process cmd = do
            print cmd
            stdout (inshell cmd empty)
    sh $ do
        afStr <- select ["01", "02", "05", "10", "20", "Common", "All"]
        let mapMasked = True
        -- mapMasked <- select [False, True]
        let afCond = case afStr of
                "01"     -> "--noMinFreq    --maxFreq 0.01"
                "02"     -> "--noMinFreq    --maxFreq 0.02"
                "05"     -> "--noMinFreq    --maxFreq 0.05"
                "10"     -> "--noMinFreq    --maxFreq 0.1"
                "20"     -> "--noMinFreq    --maxFreq 0.2"
                "Common" -> "--minFreq 0.95 --maxFreq 0.05"
                "All"    -> "--noMinFreq    --noMaxFreq"
        let bedOpt = if mapMasked then
                format ("--bedFile "%fp%" ") bed_file else " "
        let bedStr = if mapMasked then "_mapMasked" else ""
        liftIO . process $ format ("qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d "%fp%" -d "%fp%
            " -j 100000 --noTransitions "%s%" --popConfigFile pop_config_TGP.yml "%s%
            "-f AncientBritish_1000G_ras"%s%"_TVonly"%s%".table.txt") pos_ancient pos_modern afCond bedOpt afStr bedStr

  