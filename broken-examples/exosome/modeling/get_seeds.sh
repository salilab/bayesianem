#!/bin/bash

source ~/modules.sh
cd /baycells/scratch/shanot/EM_benchmark/exosome/modeling
~/imp-fast/setup_environment.sh python get_seeds.py PREFIX 1000 10.0 STEP NUM_REPLICAS
