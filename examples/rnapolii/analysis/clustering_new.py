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

data_file=sys.argv[1]
cluster_file=sys.argv[2]
num_models=int(sys.argv[3])

pdb_id=''

def get_native_score():
    with open("ref.rmf.score") as f:
        for line in f:
            tag,val = line.rstrip().split('=')
            if tag == "Total_Score":
                return float(val)


def get_selections_dictionary(input_objects):
    moldict=IMP.pmi.tools.get_molecules_dictionary(input_objects)
    seldict=defaultdict(list)
    for name, mols in moldict.items():
        for m in mols:
            seldict[name].append(IMP.atom.Selection(m))
    return seldict

class MyAre(IMP.pmi.macros.AnalysisReplicaExchange):
    def __init__(self,
                 model,
                 stat_files,
                 best_models=None,
                 alignment=True):
        IMP.pmi.macros.AnalysisReplicaExchange.__init__(self, model, stat_files, best_models, alignment)
        self.rbs1, self.beads1 = IMP.pmi.tools.get_rbs_and_beads(IMP.pmi.tools.select_at_all_resolutions(self.stath1))
        self.rbs0, self.beads0 = IMP.pmi.tools.get_rbs_and_beads(IMP.pmi.tools.select_at_all_resolutions(self.stath0))
        self.set_rmsd_selection(representation_type=IMP.atom.BALLS, state_index=0)
        self.set_alignment_selection(representation_type=IMP.atom.BALLS, state_index=0)
        self.molcopydict0=IMP.pmi.tools.get_molecules_dictionary_by_copy(IMP.atom.Selection(self.stath0, representation_type=IMP.atom.BALLS, state_index=0).get_selected_particles())
        self.molcopydict1=IMP.pmi.tools.get_molecules_dictionary_by_copy(IMP.atom.Selection(self.stath1, representation_type=IMP.atom.BALLS, state_index=0).get_selected_particles())
        self.update_seldicts()
        self.coords_cache={}

    def rmsd_helper(self, sels0, sels1, metric):
        '''
        a function that returns the permutation best_sel of sels0 that minimizes metric
        '''
        best_rmsd2 = sys.float_info.max
        best_sel = None
        for sels in itertools.permutations(sels0):
            rmsd2=0.0
            for sel0, sel1 in itertools.takewhile(lambda x: rmsd2<best_rmsd2, zip(sels, sels1)):
                if type(sel0) is list:
                    if len(sel0)!=len(sel1):
                        print 'rmsd_helper len error 1', sel0, sel1
                        exit()
                elif len(sel0.get_selected_particles())!=len(sel1.get_selected_particles()):
                    print 'rmsd_helper len error 2', sel0, sel1, len(sel0.get_selected_particles()), len(sel1.get_selected_particles())
                    exit()

                r=metric(sel0, sel1)
                rmsd2+=r*r
            if rmsd2 < best_rmsd2:
                best_rmsd2 = rmsd2
                best_sel = sels
        return  best_sel, best_rmsd2

 
    def rmsd(self,metric=IMP.atom.get_rmsd):
        '''
        Computes the RMSD. Resolves ambiguous pairs assignments
        '''
        # here we memoize the rmsd and molecular assignment so that it's not done multiple times
        n0=self.stath0.current_index
        n1=self.stath1.current_index
        if ((n0,n1) in self.pairwise_rmsd) and ((n0,n1) in self.pairwise_molecular_assignment):
            return self.pairwise_rmsd[(n0,n1)], self.pairwise_molecular_assignment[(n0,n1)]

        if self.alignment:
            self.align()
        #if it's not yet memoized
        total_rmsd=0.0
        total_N=0
        # this is a dictionary which keys are the molecule names, and values are the list of IMP.atom.Selection for all molecules that share the molecule name
        molecular_assignment={}
        for molname, sels0 in self.seldict0.items():
            sels_best_order, best_rmsd2 = self.rmsd_helper(sels0, self.seldict1[molname], metric)

            Ncoords=0
            if type(sels_best_order[0]) is list:
                Ncoords = len(sels_best_order[0])
            else:
                Ncoords = len(sels_best_order[0].get_selected_particles())
            Ncopies = len(sels_best_order)
            total_rmsd += Ncoords*best_rmsd2
            total_N += Ncoords*Ncopies

            if len(sels_best_order)!=len(self.seldict1[molname]):
                print 'len error 1', len(sels_best_order), len(self.seldict1[molname])
            for sel0, sel1 in zip(sels_best_order, self.seldict1[molname]):
                if type(sel0) is list:
                    p0 = sel0[0]
                    p1 = sel1[0]
                else:
                    p0 = sel0.get_selected_particles()[0]
                    p1 = sel1.get_selected_particles()[0]
                m0 = IMP.pmi.tools.get_molecules([p0])[0]
                m1 = IMP.pmi.tools.get_molecules([p1])[0]
                c0 = IMP.atom.Copy(m0).get_copy_index()
                c1 = IMP.atom.Copy(m1).get_copy_index()
                molecular_assignment[(molname,c0)]=(molname,c1)

        total_rmsd = math.sqrt(total_rmsd/total_N)

        self.pairwise_rmsd[(n0,n1)]=total_rmsd
        self.pairwise_molecular_assignment[(n0,n1)]=molecular_assignment
        self.pairwise_rmsd[(n1,n0)]=total_rmsd
        self.pairwise_molecular_assignment[(n1,n0)]=molecular_assignment
        return total_rmsd, molecular_assignment

    def update_seldicts(self):
        """
        Update the seldicts
        """
        self.seldict0=get_selections_dictionary(self.sel0_rmsd.get_selected_particles())
        self.seldict1=get_selections_dictionary(self.sel1_rmsd.get_selected_particles())

    def exclude_from_selection(self,**kwargs):
        s0 = IMP.atom.Selection(self.stath0, **kwargs)
        s1 = IMP.atom.Selection(self.stath1, **kwargs)
        self.sel0_rmsd-=s0
        self.sel1_rmsd-=s1
        self.update_seldicts()


root_dir='/baycells/scratch/shanot/EM_benchmark/'
output_dir='./'

mdl = IMP.Model()

are=None
are=MyAre(mdl, data_file, num_models)
are.set_rmsd_selection(representation_type=IMP.atom.BALLS, state_index=0, resolution=1)
are.exclude_from_selection(molecule='Rpb8')
print are.sel0_rmsd.get_selected_particles()

print 'clustering...'
are.alignment=False
if os.path.isfile(cluster_file):
    are.load_clusters(cluster_file)
else:
    rmsd_cutoff=30.0
    are.cluster(rmsd_cutoff)
    are.save_clusters(cluster_file)

print(are)
