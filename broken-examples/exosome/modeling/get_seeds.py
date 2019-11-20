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


class MyAre(IMP.pmi.macros.AnalysisReplicaExchange):
    def find_different_structures(self, rmsd_cutoff, N, included_structs=set([])):
        ret=included_structs
        idxs = set(range(len(self.stath0)))-ret
        while len(ret)<N and len(idxs)>1:
            n1=idxs.pop()
            d1=self.stath1[n1]
            for n0 in ret:
                d0=self.stath0[n0]
                rmsd, _ = self.rmsd()
                if rmsd<rmsd_cutoff:
                    break
            else:
                print '%d is different from everything in '%(n1), ret
                ret.add(n1)
        return ret

prefix=sys.argv[1]
num_models=int(sys.argv[2])
rmsd_cutoff=float(sys.argv[3])
step=int(sys.argv[4])
num_replicas=int(sys.argv[5])

print sys.argv

if num_models<0:
    num_models=None

output_dir='./'
pattern='%s/%s_output_%d/rmfs/*rmf3'%(output_dir,prefix,step)
sts = sorted(glob.glob(pattern))
print pattern,sts

prefix_dir=os.path.dirname(os.path.dirname(sts[0]))

mdl = IMP.Model()



datafile='%s/%sdata_%d.pkl'%(output_dir, prefix, step)
 
if os.path.isfile(datafile):
    are=MyAre(mdl,datafile,num_models)
else:
    are=MyAre(mdl,sts,num_models)
    are.save_data(datafile)

are.alignment=False
are.stath0.data.sort(key=lambda d:d.score)
are.stath1.data.sort(key=lambda d:d.score)

different_structures=set([])
c=rmsd_cutoff
while len(different_structures)<num_replicas:
    different_structures=are.find_different_structures(c, num_replicas, different_structures)
    c/=2
    print c


rmf_name = '%sseed_%s'%(prefix,step) + ".rmf3"
o=IMP.pmi.output.Output()
o.init_rmf(rmf_name, [are.stath1])

frames_saved=0
index=0
if len(different_structures)<num_replicas:
    print 'not enough frames to feed replicas'
    exit()
for n in different_structures:
    are.stath1[n]
    o.write_rmf(rmf_name)

o.close_rmf(rmf_name)
