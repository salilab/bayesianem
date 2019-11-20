import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.bayesianem
import IMP.bayesianem.restraint
import tempfile,os
import sys
import glob
import itertools
import math
import numpy as np
import os
from collections import defaultdict

import IMP.pmi.output


identifier=sys.argv[1]
data_file=sys.argv[2]
num_models=int(sys.argv[3])
sts=sys.argv[4:]

mdl = IMP.Model()

are=IMP.pmi.macros.AnalysisReplicaExchange(mdl,sts)
are.stath0.data.sort(key=lambda d:d.score)
are.stath1.data.sort(key=lambda d:d.score)
are.stath0.data=are.stath0.data[:num_models]
are.stath1.data=are.stath1.data[:num_models]
are.save_data(data_file)
