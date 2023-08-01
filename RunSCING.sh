#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=2G,h_rt=6:00:00
#$ -pe shared 12
#$ -V
#$ -m bea

python3 RunSCING.py $1