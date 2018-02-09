#!/usr/bin/env python
# -*- coding: utf-8 -*-

import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.output
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.bayesianem
import IMP.bayesianem.restraint
import os
import sys
import glob
import itertools
import math
import numpy as np
import scipy
import os
from collections import defaultdict
import operator
try:
    import cPickle as pickle
except ImportError:
    import pickle

ref_file=sys.argv[1]
data_file=sys.argv[2]
cluster_file=sys.argv[3]
num_models=int(sys.argv[4])
max_clusters=int(sys.argv[5])

pdb_id=''
print '#', sys.argv

def get_native_score():
    try:
        with open("ref.rmf.score") as f:
            for line in f:
                tag,val = line.rstrip().split('=')
                if tag == "Total_Score":
                    return float(val)
    except:
        return 1.0


def get_selections_dictionary(input_objects):
    moldict=IMP.pmi.tools.get_molecules_dictionary(input_objects)
    seldict=defaultdict(list)
    for name, mols in moldict.items():
        for m in mols:
            seldict[name].append(IMP.atom.Selection(m))
    return seldict

#fixit
def get_residue_indexes(sel):
    ret = []
    for p in sel.get_selected_particles():
        ret+=IMP.pmi.tools.get_residue_indexes(p)
    return ret

def get_min_max_residue_indexes(sel):
    idxs = [IMP.pmi.tools.get_residue_indexes(p)[0] for p in sel.get_selected_particles()]
    return min(idxs), max(idxs)


class MyAre(IMP.pmi.macros.AnalysisReplicaExchange):
    def __init__(self,
                 model,
                 stat_files,
                 best_models=None,
                 alignment=True,
                 ref_rmf=None):
        IMP.pmi.macros.AnalysisReplicaExchange.__init__(self, model, stat_files, best_models, alignment)
        self.stath1 = IMP.pmi.output.RMFHierarchyHandler(mdl, ref_file)
        d0 = self.stath1[0]
        self.stath1.current_index=0
        model.update()
        self.rbs1, self.beads1 = IMP.pmi.tools.get_rbs_and_beads(IMP.pmi.tools.select_at_all_resolutions(self.stath1))
        self.rbs0, self.beads0 = IMP.pmi.tools.get_rbs_and_beads(IMP.pmi.tools.select_at_all_resolutions(self.stath0))
        self.set_rmsd_selection(representation_type=IMP.atom.BALLS, state_index=0)
        self.update_seldicts()
        self.set_alignment_selection(representation_type=IMP.atom.BALLS, state_index=0)
        self.molcopydict0=IMP.pmi.tools.get_molecules_dictionary_by_copy(IMP.atom.Selection(self.stath0, representation_type=IMP.atom.BALLS, state_index=0).get_selected_particles())
        self.molcopydict1=IMP.pmi.tools.get_molecules_dictionary_by_copy(IMP.atom.Selection(self.stath1, representation_type=IMP.atom.BALLS, state_index=0).get_selected_particles())

    def exclude_from_selection(self,**kwargs):
        s0 = IMP.atom.Selection(self.stath0, **kwargs)
        s1 = IMP.atom.Selection(self.stath1, **kwargs)
        self.sel0_rmsd-=s0
        self.sel1_rmsd-=s1
        self.update_seldicts()

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
        #if ((n0,n1) in self.pairwise_rmsd) and ((n0,n1) in self.pairwise_molecular_assignment):
        #    return self.pairwise_rmsd[(n0,n1)], self.pairwise_molecular_assignment[(n0,n1)]

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
        return total_rmsd, molecular_assignment

    def update_seldicts(self):
        """
        Update the seldicts
        """
        self.seldict0=get_selections_dictionary(self.sel0_rmsd.get_selected_particles())
        self.seldict1=get_selections_dictionary(self.sel1_rmsd.get_selected_particles())

root_dir='/baycells/scratch/shanot/EM_benchmark/'
output_dir='./'

print '#', root_dir
print '#', output_dir

mdl = IMP.Model()

are=None
are=MyAre(mdl, data_file, num_models, alignment=False, ref_rmf=ref_file)

if os.path.isfile(cluster_file):
    are.load_clusters(cluster_file)

native_score=get_native_score()


ref_gaussians=[]
target_gmm_file='../../../../divide_conquer_benchmark/1883/%d_imp.gmm'%2
IMP.isd.gmm_tools.decorate_gmm_from_text(target_gmm_file, ref_gaussians, mdl)

rmf_gmm_sel=IMP.atom.Selection(are.stath0, representation_type=IMP.atom.DENSITIES, state_index=0)
rmf_gaussians=rmf_gmm_sel.get_selected_particles()

tmp=[]
for p in rmf_gaussians:
    if not IMP.core.NonRigidMember.get_is_setup(p):
        tmp.append(p)

rmf_gaussians=tmp

num_subunits=len(are.stath1.get_children()[0].get_children())
num_clusters=len([c.members[0] for c in are.clusters if c.members[0]<len(are.stath0)])
num_gaussian=len(ref_gaussians)

try:
    time_per_frame=np.mean([float(d.features['Stopwatch_None_delta_seconds']) for d in are.stath0.data])
except:
    time_per_frame=-1

headers=['pdb', '#subunits', '#clusters', 'cluster dispersion', 'time per frame']
formats=['%s', '%d', '%d', '%.2f', '%.2f']
values=['.', num_subunits, num_clusters, are.get_shannon_entropy(), time_per_frame]

sep=' & '


scores=are.stath0.get_scores()
scores_sorted=sorted(scores)

best_scoring_index = scores.index(scores_sorted[0])


rmsd_file = '%s_rmsd.pkl'%(data_file)
def compute_rmsds():
    for i, d in enumerate(are.stath0):
        rmsd, _ = are.rmsd()
        rmsds.append(rmsd)
        f=open(rmsd_file,'wb')
        pickle.dump(rmsds,f)

rmsds=[]
if os.path.isfile(rmsd_file):
    f=open(rmsd_file, 'rb')
    rmsds=pickle.load(f)
    if len(rmsds)<num_models:
        compute_rmsds()
else:
    compute_rmsds()

best_rmsd_index, best_rmsd = min(enumerate(rmsds), key=operator.itemgetter(1))


def get_xyzs(mol0, mol1):
    
    s0=IMP.atom.Selection(mol0)
    s1=IMP.atom.Selection(mol1)

    idxs0 = set(get_residue_indexes(s0))
    idxs1 = set(get_residue_indexes(s1))
    idxs = list(idxs0 & idxs1)

    xyzs0 = []
    xyzs1 = []

    for i in idxs:
        s0 = IMP.atom.Selection(mol0, residue_index=i)
        s1 = IMP.atom.Selection(mol1, residue_index=i)
        p0=s0.get_selected_particles()[0]
        p1=s1.get_selected_particles()[0]
        if ((not IMP.core.NonRigidMember.get_is_setup(p0)) and (not IMP.core.NonRigidMember.get_is_setup(p1))) and (IMP.core.RigidMember.get_is_setup(p0) and IMP.core.RigidMember.get_is_setup(p1)):
            xyzs0.append(IMP.core.XYZ(p0))
            xyzs1.append(IMP.core.XYZ(p1))
    return xyzs0, xyzs1

def get_placement_score(mol0, mol1):
        xyz0, xyz1 = get_xyzs(mol0, mol1)

        xyz0 = np.array([p.get_coordinates() for p in xyz0])
        xyz1 = np.array([p.get_coordinates() for p in xyz1])
        centroid0 = np.mean(xyz0,axis=0)
        centroid1 = np.mean(xyz1,axis=0)

        dist = np.linalg.norm(centroid0-centroid1)

        xyz0-=centroid0
        xyz1-=centroid1

        A = np.dot(xyz0.transpose(), xyz1)
        U = scipy.linalg.sqrtm(np.dot(A.transpose(), A))*np.linalg.inv(A)
        [[a,b,c],[d,e,f],[g,h,i]] = U
        u=[h-f, c-g, d-b]
        arg = np.linalg.norm(u)/2
        arg2 = (a + e + i - 1)/2
        if arg<-1:
            arg=-1
        elif arg>1:
            arg=1
        anglesin = np.arcsin(arg)
        if arg2 < 0:
            angle=np.pi-anglesin
        else:
            angle=anglesin

        #dist, angle = IMP.atom.get_placement_score(xyz0, xyz1)
        return [dist, angle]
    
indexes_labels=['best rsmd'] + ['C%d'%c.cluster_id for c in are.clusters][:max_clusters]
indexes = [best_rmsd_index] + [c.members[0] for c in are.clusters][:max_clusters]


print '#>>>> ', indexes_labels
print '#>>>> ', indexes

subunit_strings=[]
APS_strings=[]
APS_header=[]

def get_p10():
    molecular_assignment = are.pairwise_molecular_assignment[(are.stath0.current_index,are.stath1.current_index)]
    N=0.0
    N_tot=0.0
    for (m0, c0), (m1,c1) in molecular_assignment.items():
        mol0 = are.molcopydict0[m0][c0]
        mol1 = are.molcopydict1[m1][c1]
        xyz0, xyz1 = get_xyzs(mol0,mol1)
        xyz0 = [x.get_coordinates() for x in xyz0]
        xyz1 = [x.get_coordinates() for x in xyz1]
        d = np.array(xyz0)-np.array(xyz1)
        d = np.sqrt(np.sum(d*d,axis=1))
        N+=np.sum(d<20.0)
        N_tot+=len(d)
    return N/N_tot


chains = {
"1CS4:A":"A",
"1CS4:B":"B",
"1CS4:C":"C",
"1GPQ:A":"A",
"1GPQ:A.1":"B",
"1GPQ:C":"C",
"1GPQ:C.1":"D",
"1MDA:A":"A",
"1MDA:A.1":"B",
"1MDA:H":"H",
"1MDA:H.1":"J",
"1MDA:L":"L",
"1MDA:L.1":"M",
"1SUV:A":"A",
"1SUV:A.1":"B",
"1SUV:C":"C",
"1SUV:D":"D",
"1SUV:E":"E",
"1SUV:F":"F",
"1TYQ:A":"A",
"1TYQ:B":"B",
"1TYQ:C":"C",
"1TYQ:D":"D",
"1TYQ:E":"E",
"1TYQ:F":"F",
"1TYQ:G":"G",
"1VCB:G":"A",
"1VCB:H":"B",
"1VCB:F":"C",
"1Z5S:A":"A",
"1Z5S:B":"B",
"1Z5S:C":"C",
"1Z5S:D":"D",
"2BBK:H":"H",
"2BBK:H.1":"J",
"2BBK:L":"L",
"2BBK:L.1":"M",
"2BO9:A":"A",
"2BO9:B":"B",
"2BO9:A.1":"C",
"2BO9:B.1":"D",
"2DQJ:H":"H",
"2DQJ:L":"L",
"2DQJ:Y":"Y",
"2GC7:I":"A",
"2GC7:J":"B",
"2GC7:G":"C",
"2GC7:H":"D",
"2UZX:A":"A",
"2UZX:B":"B",
"2WVY:A":"A",
"2WVY:A.1":"B",
"2WVY:A.2":"C",
"2Y7H:A":"A",
"2Y7H:B.1":"C",
"2Y7H:B":"B",
"2Y7H:D":"D",
"2Y7H:E":"E",
"3LU0:A":"A",
"3LU0:A.1":"B",
"3LU0:C":"C",
"3LU0:C":"C",
"3LU0:C":"C",
"3LU0:D":"D",
"3LU0:D":"D",
"3LU0:D":"D",
"3LU0:E":"E",
"3NVQ:A":"A",
"3NVQ:A.1":"E",
"3NVQ:B":"B",
"3NVQ:B.1":"F",
"3PDU:A":"A",
"3PDU:A.1":"B",
"3PDU:A.2":"C",
"3PDU:A.3":"D",
"3PUV:A":"A",
"3PUV:A.1":"B",
"3PUV:E":"E",
"3PUV:G":"G",
"3PUV:F":"F",
"3PUV:F":"F",
"3PUV:F":"F",
"3PUV:F":"F",
"3PUV:F":"F",
"3R5D:A":"A",
"3R5D:A.1":"B",
"3R5D:A.2":"C",
"3SFD:A":"A",
"3SFD:B":"B",
"3SFD:C":"C",
"3SFD:C":"C",
"3SFD:C":"C",
"3SFD:D":"D",
"3V6D:A":"A",
"3V6D:B":"B",
"3V6D:E":"T",
"3V6D:F":"P"}

for i,label in zip(indexes, indexes_labels):
    if i>len(are.stath0):
        print 'i>len(are.stath0) (i=%d)'%i
        continue

    di=are.stath0[i]
    mdl.update()

    pse = (native_score-di.score)/native_score
    are.pairwise_rmsd={}
    are.pairwise_molecular_assignment={}
    rmsd, molecular_assignment = are.rmsd()
    cc = IMP.bayesianem.gem_score_cc(ref_gaussians, rmf_gaussians)
    PS=0.0
    AS=0.0
    N=0
    _subunit_string = [' ']
    _APS_header = [' ']
    _APS_strings = [label]
    print '#>>', label, i
    print '#>>', len(sorted(molecular_assignment.keys(), key=lambda m:str(m[0])+str(m[1]))), molecular_assignment
    for m0,c0 in sorted(molecular_assignment.keys(), key=lambda m:str(m[0])+str(m[1])):
        m1, c1 = molecular_assignment[(m0,c0)]
        mol0 = are.molcopydict0[m0][c0]
        mol1 = are.molcopydict1[m1][c1]
        xyz0, xyz1 = get_xyzs(mol0, mol1)

        molrmsd = IMP.atom.get_rmsd(xyz0, xyz1)
        if 'Rpb1' in str(m0):
            domains = [ [1,1140], [1141,1274],[1275,1733] ]
            for a,b in domains:
                x0 = [x for x in xyz0 if (IMP.atom.Residue(x.get_particle()).get_index()<=b and a<=IMP.atom.Residue(x.get_particle()).get_index())]
                x1 = [x for x in xyz1 if (IMP.atom.Residue(x.get_particle()).get_index()<=b and a<=IMP.atom.Residue(x.get_particle()).get_index())]
                domrmsd = IMP.atom.get_rmsd(x0, x1)
                print '# >>>> domain', label, m0, [a,b], domrmsd


        xyz0 = np.array([p.get_coordinates() for p in xyz0])
        xyz1 = np.array([p.get_coordinates() for p in xyz1])


        #sel0 = IMP.atom.Selection(mol0)-IMP.atom.Selection(mol0, resolution=20)
        #sel1 = IMP.atom.Selection(mol1)-IMP.atom.Selection(mol1, resolution=20)
        #molrmsd = IMP.atom.get_rmsd(sel0, sel1)
        #print '# >>> sel0', sel0.get_selected_particles()
        #print '# >>> sel1', sel1.get_selected_particles()

        n = len(xyz0)
        dist, angle = get_placement_score(mol0, mol1)
        angle*=180.0/math.pi

        if ':' in str(m0):
            mname, subname = str(m0).split(":")
            mnamecopy = '%s.%d'%(str(m0),c0)
            if mnamecopy in chains:
                _subunit_string += ['%s:%s (%d)'%(mname,chains[mnamecopy], n),  ' ']
            else:
                _subunit_string += ['%s:%s (%d)'%(mname,chr(ord(subname)+c0), n),  ' ']
        else:
            _subunit_string += ['%s (%d)'%(str(m0), n),  ' ']

                
        _APS_header.append('rmsd [Å]')
        _APS_header.append('APS [Å,°]')
        _APS_strings.append("%.1f"%np.around(molrmsd, decimals=1))
        _APS_strings.append("%.1f, %.1f"%(np.around(dist, decimals=1), np.around(angle, decimals=1)))

        N+=n
        PS+=n*dist
        AS+=n*angle

    subunit_strings.append(_subunit_string)
    APS_strings.append(_APS_strings)
    APS_header.append(_APS_header)

    print '#>>> get_p10'
    p10 = get_p10()
    
    headers+=['%s rank'%label,  'rmsd [Å]', 'p(10)', 'CC', 'APS [Å,°]', '']
    formats+=['%d', '%.1f', '%.2f', '%.2f', '%.1f', '%.1f']
    values+=[i, np.around(rmsd, decimals=1), p10, cc, np.around(PS/N, decimals=1), np.around(AS/N, decimals=1)]

h=sep.join(headers)
v=sep.join(formats)%tuple(values)
print h
print v

print '\n%s\n'%('='*max(len(h), len(v)))

for n,(sub, h, s) in enumerate(zip(subunit_strings, APS_header, APS_strings)):
    if n==0:
        print sep.join(sub)
        print sep.join(h)
    print sep.join(s)

h=sep.join(headers)
v=sep.join(formats)%tuple(values)
print '\n%s\n'%('='*max(len(h), len(v)))


for i, di in enumerate(are.stath0.data):
    if i>=len(rmsds):
        break
    print i, di.rmf_name, di.rmf_index, rmsds[i], di.score
