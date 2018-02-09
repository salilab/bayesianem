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
import IMP.pmi.restraints.crosslinking
import IMP.bayesianem
import IMP.bayesianem.restraint
import tempfile,os
import sys
import glob
import itertools
import math
import numpy as np
import scipy
import os
from collections import defaultdict
import IMP.pmi.output
import operator
try:
    import cPickle as pickle
except ImportError:
    import pickle

ref_file=sys.argv[1]
data_file=sys.argv[2]
target_gmm_file=sys.argv[3]
cluster_file=sys.argv[4]
num_models=int(sys.argv[5])

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
                 alignment=True,
                 ref_rmf=None):
        IMP.pmi.macros.AnalysisReplicaExchange.__init__(self, model, stat_files, best_models, alignment)
        self.stath1 = IMP.pmi.output.RMFHierarchyHandler(mdl, ref_rmf)
        self.stath1.current_index=0
        self.rbs1, self.beads1 = IMP.pmi.tools.get_rbs_and_beads(IMP.pmi.tools.select_at_all_resolutions(self.stath1))
        self.rbs0, self.beads0 = IMP.pmi.tools.get_rbs_and_beads(IMP.pmi.tools.select_at_all_resolutions(self.stath0))
        self.set_rmsd_selection(representation_type=IMP.atom.BALLS, state_index=0)
        self.set_alignment_selection(representation_type=IMP.atom.BALLS, state_index=0)
        self.molcopydict0=IMP.pmi.tools.get_molecules_dictionary_by_copy(IMP.atom.Selection(self.stath0, representation_type=IMP.atom.BALLS, state_index=0).get_selected_particles())
        self.molcopydict1=IMP.pmi.tools.get_molecules_dictionary_by_copy(IMP.atom.Selection(self.stath1, representation_type=IMP.atom.BALLS, state_index=0).get_selected_particles())
        self.fix_seldicts()

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

are=None
ref_rmf='ref.rmf'
are=MyAre(mdl, data_file, num_models, alignment=False, ref_rmf=ref_rmf)

are.exclude_from_selection(molecule='Rrp46_gfp', residue_indexes=range(247, 1000))


if os.path.isfile(cluster_file):
    are.load_clusters(cluster_file)

native_score=get_native_score()

ref_gaussians=[]
IMP.isd.gmm_tools.decorate_gmm_from_text(target_gmm_file, ref_gaussians, mdl)

rmf_gmm_sel=IMP.atom.Selection(are.stath0, representation_type=IMP.atom.DENSITIES, state_index=0)
rmf_gaussians=rmf_gmm_sel.get_selected_particles()

tmp=[]
for p in rmf_gaussians:
    if not IMP.core.NonRigidMember.get_is_setup(p):
        tmp.append(p)

rmf_gaussians=tmp

num_subunits=len(are.stath1.get_children()[0].get_children())
num_clusters=len(are)
num_gaussian=len(ref_gaussians)

from IMP.pmi.io.crosslink import FilterOperator as FO
cldbkc=IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
cldbkc.set_protein1_key("Protein 1")
cldbkc.set_protein2_key("Protein 2")
cldbkc.set_residue1_key("Residue1")
cldbkc.set_residue2_key("Residue2")
cldbkc.set_id_score_key("Score")
cldbkc.set_unique_id_key("Unique ID")
cldb=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
cldb.create_set_from_file("../data/exosome_XLMS_column07012014.csv")

cldb.offset_residue_index('GFP', 246)
cldb.rename_proteins({"GFP":"Rrp46_gfp","Rrp46":"Rrp46_gfp"})


xl = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=are.stath0,
                            CrossLinkDataBase=cldb,
                            length=21.0,
                            slope=0.02,
                            resolution=1.0,
                            label="XL")
xl.add_to_model()
xl.set_psi_is_sampled(False)
psi=xl.psi_dictionary["PSI"][0]
psi.set_scale(0.01)

try:
    time_per_frame=np.mean([float(d.features['Stopwatch_None_delta_seconds']) for d in are.stath0.data])
except:
    time_per_frame=-1

headers=['pdb', '#subunits', '#clusters', 'cluster dispersion', 'time per frame']
formats=['%s', '%d', '%d', '%.2f', '%.2f']
values=['.', num_subunits, len(are), are.get_shannon_entropy(), time_per_frame]

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
    imin=1
    imax=sys.maxint
    
    s0=IMP.atom.Selection(mol0)
    s1=IMP.atom.Selection(mol1)

    for s in s0,s1:
        m,M = get_min_max_residue_indexes(s)
        imin=max(m,imin)
        imax=min(M,imax)

    xyzs0 = []
    xyzs1 = []
    prevs = [None for s in s0,s1]

    for i in range(imin, imax+1):
        s0 = IMP.atom.Selection(mol0, residue_index=i)
        s1 = IMP.atom.Selection(mol1, residue_index=i)
        #idxs = [s.get_selected_particle_indexes()[0].get_index() for s in s0,s1]
        #if any(idx!=prev for idx,prev in zip(idxs, prevs)):
        #prevs[:]=idxs[:]
        p0=s0.get_selected_particles()[0]
        p1=s1.get_selected_particles()[0]
        if (not IMP.core.NonRigidMember.get_is_setup(p0)) and (not IMP.core.NonRigidMember.get_is_setup(p1)):
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
        if arg<-1:
            print '#changing arg from %e to %e'%(arg,-1)
            arg=-1
        elif arg>1:
            print '#changing arg from %e to %e'%(arg,1)
            arg=1
        angle = np.arcsin(arg)
        print '#checksum=',np.sum(xyz1-np.dot(xyz0, U))

        #t=IMP.algebra.get_transformation_aligning_first_to_second(xyz0,xyz1)

        return [dist, angle]
        #return IMP.atom.get_placement_score(xyz0, xyz1)#, np.mean([p.get_coordinates() for p in xyz0],axis=0), np.mean([p.get_coordinates() for p in xyz1], axis=0)

def get_rmsd(mol0, mol1):
    xyz0, xyz1 = get_xyzs(mol0, mol1)
    return IMP.atom.get_rmsd(xyz0, xyz1)

    
indexes_labels=['best rsmd'] + ['C%d'%c.cluster_id for c in are.clusters]
indexes = [best_rmsd_index] + [c.members[0] for c in are.clusters]

APS_strings=[]
APS_header=[]

def get_p10():
    rmsd, molecular_assignment = are.rmsd()
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
        N+=np.sum(d<10.0)
        N_tot+=len(d)
    return N/N_tot


print '\n%s\n'%('*'*50)

for i,label in zip(indexes, indexes_labels):
    if i>len(are.stath0):
        continue

    di=are.stath0[i]
    mdl.update()
    pse = (native_score-di.score)/native_score
    rmsd, molecular_assignment = are.rmsd()
    cc = IMP.bayesianem.gem_score_cc(ref_gaussians, rmf_gaussians)
    PS=0.0
    AS=0.0
    N=0
    _APS_header = ['Protein', 'rmsd', 'Pos', 'Ang']
    _APS_strings = []
    for p0,p1 in xl.get_particle_pairs():
        print p0, p1, IMP.core.get_distance(IMP.core.XYZ(p0), IMP.core.XYZ(p1))
        

    for (m0, c0), (m1,c1) in molecular_assignment.items():
        mol0 = are.molcopydict0[m0][c0]
        mol1 = are.molcopydict1[m1][c1]
        aps = get_placement_score(mol0, mol1)
        
        rmsd01 = get_rmsd(mol0, mol1)

        _APS_strings.append(["%s"%m0, "%.2f"%rmsd01, "%.2f"%aps[0], "%.2f"%(aps[1]*180.0/math.pi)])

        n = len(IMP.atom.Selection(mol0).get_selected_particles())
        N+=n
        PS+=n*aps[0]
        AS+=n*aps[1]*180.0/math.pi

    #_APS_header[0]='%s %s'%(label, _APS_header[0])
    APS_strings.append(_APS_strings)
    APS_header.append(_APS_header)

    p10 = get_p10()
    
    headers+=['%s rank'%label, '%score error', 'rmsd', 'p10', 'cc', 'Pos', 'Ang']
    formats+=['%d', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f']
    values+=[i, pse, rmsd, p10, cc, PS/N, AS/N]

h=sep.join(headers)
v=sep.join(formats)%tuple(values)
print h
print v

print '\n%s\n'%('='*max(len(h), len(v)))






#for h, s in zip(APS_header, APS_strings):
#    print sep.join(h)
#    print sep.join(s)
#
#h=sep.join(headers)
#v=sep.join(formats)%tuple(values)
#print '\n%s\n'%('='*max(len(h), len(v)))
print sep.join(_APS_header)
for aps in _APS_strings:
    print sep.join(aps)

print '\n\n%s\n\n'%('='*max(len(h), len(v)))


for i, di in enumerate(are.stath0.data):
    if i>=len(rmsds):
        break
    print i, di.rmf_name, di.rmf_index, rmsds[i], di.score
