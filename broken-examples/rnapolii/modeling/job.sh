#!/bin/bash

source ~/modules.sh
cd /baycells/scratch/shanot/EM_benchmark/imp_tutorial/rnapolii/modeling_em/
mpirun -n NUM_REPLICAS ~/imp-fast/setup_environment.sh python modeling.py PREFIX STEP
