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
import numpy as np


output_objects=[]

root_dir='/baycells/scratch/shanot/EM_benchmark/'

output_prefix=""
if len(sys.argv)>1:
    output_prefix=sys.argv[1]
EV_weight=1.0
if len(sys.argv)>2:
    EV_weight=float(sys.argv[2])

print output_prefix, EV_weight
output_dir='%s/modeling/'%(root_dir)
gmm_file='%s/data/em/512_imp.gmm'%(root_dir)


###################### SYSTEM SETUP #####################
# Read sequences etc
# The TopologyReader reads the text file, and the BuildSystem macro constructs it
mdl = IMP.Model()
reader = IMP.pmi.topology.TopologyReader('%s/modeling/topology.dat'%(root_dir),
                                         pdb_dir =   '%s/data/'%(root_dir),
                                         fasta_dir = '%s/data/'%(root_dir),
                                         gmm_dir =   '%s/gmm'%(output_dir))

bs = IMP.pmi.macros.BuildSystem(mdl)
bs.add_state(reader) # note you can call this multiple times to create a multi-state system
hier, dof = bs.execute_macro()


# Connectivity keeps things connected along the backbone (ignores if inside same rigid body)
crs = []
moldict = bs.get_molecules()[0]
mols = []
for molname in moldict:
    for mol in moldict[molname]:
        IMP.pmi.tools.display_bonds(mol)
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
        cr.add_to_model()
        cr.set_label(molname)
        output_objects.append(cr)
        crs.append(cr)
        mols.append(mol)

# Excluded volume - automatically more efficient due to rigid bodies
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects = mols)
evr.add_to_model()
evr.set_weight(EV_weight)
output_objects.append(evr)


# EM restraint
densities = IMP.atom.Selection(hier,representation_type=IMP.atom.DENSITIES).get_selected_particles()
# for EM maps simulated from pdbs with gaps, we select only the Gaussians that belong to that pdb
# when using real data, just pass densities to the IMP.bayesianem.restraint.GaussianEMRestraintWrapper below
pdb_densities = []
for p in densities:
    if IMP.core.RigidBody.get_is_setup(p) and not IMP.core.NonRigidMember.get_is_setup(p):
        pdb_densities.append(p)
   
gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(pdb_densities,target_fn=gmm_file,
                                            scale_target_to_mass=True,
                                            slope=0.01,
                                            target_radii_scale=3.0,
                                            target_is_rigid_body=False)
gem.set_label("EM")
gem.add_target_density_to_hierarchy(hier)
gem.add_to_model()
output_objects.append(gem)

###################### SAMPLING #####################
# First shuffle the system
# Mind that you cannod shuffle the root hier, as it have the gaussian beads attached to
# by-state shuffling is safer

for state in bs.system.get_states():
    IMP.pmi.tools.shuffle_configuration(state.get_hierarchy(),
                                        max_translation=100)

# Quickly move all flexible beads into place
dof.optimize_flexible_beads(10)

## Run replica exchange Monte Carlo sampling
rex_obj=None
weights = gem.weight * np.array([0.1, 1]*100)
i=0;
for w in weights:
    gem.set_weight(w);
    print rex_obj
    num_frames=100
    if w==1:
        num_frames=1000
    if rex_obj is None:
        rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                            root_hier=hier,                          # pass the root hierarchy
                                            monte_carlo_sample_objects=dof.get_movers(),  # pass MC movers
                                            global_output_directory='%s/%soutput_%d_%g/'%(output_dir,output_prefix, i,w),
                                            output_objects=output_objects,
                                            monte_carlo_steps=10,
                                            number_of_best_scoring_models=0,      # set >0 to store best PDB files (but this is slow to do online)
                                            number_of_frames=num_frames)                   # increase number of frames to get better results!
        rex.execute_macro()
        rex_obj=rex.replica_exchange_object
    else:
        rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                                root_hier=hier,                          # pass the root hierarchy
                                                monte_carlo_sample_objects=dof.get_movers(),  # pass MC movers
                                                global_output_directory='%s/%soutput_%d_%g/'%(output_dir,output_prefix,i,w),
                                                output_objects=output_objects,
                                                monte_carlo_steps=10,
                                                number_of_best_scoring_models=0,      # set >0 to store best PDB files (but this is slow to do online)
                                                number_of_frames=num_frames,                   # increase number of frames to get better results!
                                                replica_exchange_object=rex_obj)
        rex.execute_macro()
    i=i+1
