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
    if pdb_id=='1suv':
        seldict['1SUV:C']+=seldict['1SUV:D']
        seldict['1SUV:E']+=seldict['1SUV:F']
        del seldict['1SUV:D']
        del seldict['1SUV:F']
    return seldict

#fixit
def get_residue_indexes(sel):
    return [IMP.pmi.tools.get_residue_indexes(p) for p in sel.get_selected_particles()]
def get_min_max_residue_indexes(sel):
    idxs = get_residue_indexes(sel)
    return min(idxs)[0], max(idxs)[0]


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
        self.fix_seldicts()
        self.coords_cache={}

    def fix_seldicts(self):
        for name, sels0 in self.seldict0.items():
            sels1=self.seldict1[name]
            if any([len(x.get_selected_particles())!=len(sels0[0].get_selected_particles()) for x in sels0+sels1]):
                imin=1
                imax=sys.maxint
                mols0=[]
                mols1=[]

                for s0,s1 in zip(sels0, sels1):
                    mols0.append(IMP.pmi.tools.get_molecules(s0.get_selected_particles()[0]))
                    mols1.append(IMP.pmi.tools.get_molecules(s1.get_selected_particles()[0]))
                    for s in s0,s1:
                        m,M = get_min_max_residue_indexes(s)
                        imin=max(m,imin)
                        imax=min(M,imax)

                xyzs0 = [[] for s in sels0]
                prevs = [None for s in sels0]
                xyzs1 = [[] for s in sels1]

                for i in range(imin, imax+1):
                    sels0 = [IMP.atom.Selection(m, residue_index=i) for m in mols0]
                    idxs = [s.get_selected_particle_indexes()[0].get_index() for s in sels0]
                    #if any(idx!=prev for idx,prev in zip(idxs, prevs)):
                    sels1 = [IMP.atom.Selection(m, residue_index=i) for m in mols1]
                    for j, (s0, s1) in enumerate(zip(sels0, sels1)):
                        p0=s0.get_selected_particles()[0]
                        p1=s1.get_selected_particles()[0]
                        #if (not IMP.core.NonRigidMember.get_is_setup(p0)) and (not IMP.core.NonRigidMember.get_is_setup(p1)):
                        xyzs0[j].append(IMP.core.XYZ(p0))
                        xyzs1[j].append(IMP.core.XYZ(p1))
                self.seldict0[name]=xyzs0
                self.seldict1[name]=xyzs1

    def exclude_from_selection(self,**kwargs):
        s0 = IMP.atom.Selection(self.stath0, **kwargs)
        s1 = IMP.atom.Selection(self.stath1, **kwargs)
        self.sel0_rmsd-=s0
        self.sel1_rmsd-=s1
        self.update_seldicts()
        self.fix_seldicts()

    def get_cached_xyz(self, stath, sel):
        i = (stath.current_index, sel)
        if i not in self.coords_cache:
            self.coords_cache[i] = np.array([IMP.core.XYZ(p).get_coordinates() for p in sel.get_selected_particles()])
        return self.coords_cache[i]

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

                #xyz0 = self.get_cached_xyz(self.stath0, sel0)
                #xyz1 = self.get_cached_xyz(self.stath1, sel1)
                #delta=xyz1-xyz0
                #r2 = np.mean(np.square(delta).sum(axis=1))
                #r2 = np.square(delta).sum(axis=1).mean()
                #rmsd2+=r2
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
                m0,m1 = IMP.pmi.tools.get_molecules([p0,p1])
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
        if pdb_id=='1suv':
            self.seldict0['1SUV:C']+=self.seldict0['1SUV:D']
            self.seldict1['1SUV:C']+=self.seldict1['1SUV:D']
            self.seldict0['1SUV:E']+=self.seldict0['1SUV:F']
            self.seldict1['1SUV:E']+=self.seldict1['1SUV:F']
            del self.seldict0['1SUV:D']
            del self.seldict1['1SUV:D']
            del self.seldict0['1SUV:F']
            del self.seldict1['1SUV:F']


root_dir='/baycells/scratch/shanot/EM_benchmark/'
output_dir='./'

mdl = IMP.Model()

#print(output_dir)
#matchdir=['%s/*output_*_%s/'%(output_dir, suffix) for suffix in ["1","hw_all_hev"]]
#print matchdir
#dirs = []
#for d in matchdir:
#    dirs+=glob.glob(d)
#for d in dirs:
#    print d
#sts=[]
##print 'sts', sts
#for d in dirs: 
#    sts += ['%s/stat.%d.out'%(d, n) for n in range(48)]
#for s in sts:
#    print s

are=None
are=MyAre(mdl, data_file, num_models)

are.set_rmsd_selection(state_index=0)
are.exclude_from_selection(molecule='Rrp46_gfp', residue_indexes=range(247, 1000))

print 'clustering...'
are.alignment=False
if os.path.isfile(cluster_file):
    are.load_clusters(cluster_file)
else:
    rmsd_cutoff=20.0
    are.cluster(rmsd_cutoff)
    are.save_clusters(cluster_file)

print(are)
