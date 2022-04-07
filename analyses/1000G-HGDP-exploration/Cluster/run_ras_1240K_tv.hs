#!/usr/bin/env stack
--stack script --resolver=lts-19.0 --package turtle
{-# LANGUAGE OverloadedStrings #-}

import Control.Monad (forM_)
import Turtle
import Prelude hiding (FilePath)


repo_dir :: FilePath
repo_dir = "/mnt/archgen/poseidon/published_data"

bed_file :: FilePath
bed_file = "/mnt/archgen/users/schiffels/hs37m_filt35_99.bed"

main :: IO ()
main = do
    let process cmd = do
            print cmd
            stdout (inshell cmd empty)
    sh $ do
        afStr <- select ["01", "02", "05", "10", "20", "Common", "All"]
        mapMasked <- select [False, True]
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
        liftIO . process $ format ("qsub -V -b y -cwd -l h_vmem=16G xerxes ras -d "%fp%
            " -j 100000 --noTransitions "%s%" --popConfigFile pop_config_TGP_1240K.yml "%s%
            "-f ../Data/AncientBritish_1000G_1240K_ras"%s%"_TVonly"%s%".table.txt") repo_dir afCond bedOpt afStr bedStr

  