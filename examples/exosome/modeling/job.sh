#!/bin/bash

source ~/modules.sh
cd /baycells/scratch/shanot/EM_benchmark/exosome/modeling
mpirun -n NUM_REPLICAS ~/imp-fast/setup_environment.sh python model.py PREFIX STEP
