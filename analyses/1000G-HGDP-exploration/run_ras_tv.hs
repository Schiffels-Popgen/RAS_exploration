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

main :: IO ()
main = do
    let process cmd = do
            print cmd
            stdout (inshell cmd empty)
    sh $ do
        (afStr, af) <- select [("01", 0.01), ("02", 0.02), ("05", 0.05), ("10", 0.1), ("20", 0.2)]
        liftIO . process $ format ("qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d "%fp%" -d "%fp%
            " -j 100000 --noTransitions --noMinFreq --maxFreq "%f%" --popConfigFile pop_config_TGP.yml "%
            "-f AncientBritish_1000G_ras"%s%"_TVonly.table.txt") pos_ancient pos_modern af afStr
    process $ format ("qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d "%fp%" -d "%fp%
        " -j 100000 --noTransitions --noMinFreq --noMaxFreq --popConfigFile pop_config_TGP.yml "%
        "-f AncientBritish_1000G_rasAll_TVonly.table.txt") pos_ancient pos_modern
    process $ format ("qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d "%fp%" -d "%fp%
        " -j 100000 --noTransitions --minFreq 0.05 --maxFreq 0.95 --popConfigFile pop_config_TGP.yml "%
        "-f AncientBritish_1000G_rasCommon_TVonly.table.txt") pos_ancient pos_modern
  