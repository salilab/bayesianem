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


output_objects=[]


###################### SYSTEM SETUP #####################
# Read sequences etc
topology='''
|molecule_name|color|fasta_fn|fasta_id|pdb_fn|chain|residue_range|pdb_offset|bead_size|em_residues_per_gaussian|rigid_body|super_rigid_body|chain_of_super_rigid_bodies|
|Dis3  |blue      |exosome.fasta|Dis3 |4IFD.pdb|J|1,237  |0 |5|20 |1 |1,3| |
|Dis3  |blue      |exosome.fasta|Dis3 |4IFD.pdb|J|238,471|0 |5|20 |2 |1,3| |
|Dis3  |blue      |exosome.fasta|Dis3 |4IFD.pdb|J|472,END |0 |5|20 |3 |1,3| |
|Rrp45  |green    |exosome.fasta|Rrp45|4IFD.pdb|A|1,END   |0 |5|20 |4 |2,3| |
|Rrp4  |orange    |exosome.fasta|Rrp4 |4IFD.pdb|H|1,102  |0 |5|20 |4 |2,3| |
|Rrp4  |yellow    |exosome.fasta|Rrp4 |4IFD.pdb|H|103,END |0 |5|20 |4 |2,3| |
|Csl4  |salmon    |exosome.fasta|Csl4 |4IFD.pdb|I|1,END   |0 |5|20 |4 |2,3| |
|Mtr3  |gold      |exosome.fasta|Mtr3 |4IFD.pdb|F|1,30   |0 |5|20 |4 |2,3| |
|Mtr3  |gold      |exosome.fasta|Mtr3 |4IFD.pdb|F|31,END  |0 |5|20 |4 |2,3| |
|Rrp40 |pink      |exosome.fasta|Rrp40|4IFD.pdb|G|1,60   |0 |5|20 |4|2,3| |
|Rrp40 |pink      |exosome.fasta|Rrp40|4IFD.pdb|G|61,END  |0 |5|20 |4|2,3| |
|Rrp42 |red       |exosome.fasta|Rrp42|4IFD.pdb|E|1,END   |0 |5|20 |4|2,3| |
|Ski6  |white     |exosome.fasta|Ski6 |4IFD.pdb|B|1,END   |0 |5|20 |4|2,3| |
|Rrp46 |purple    |exosome.fasta|Rrp46|4IFD.pdb|D|1,END   |0 |5|20 |4|2,3| |
|Rrp43 |gray      |exosome.fasta|Rrp43|4IFD.pdb|C|1,END   |0 |5|20 |4|2,3| |
'''


# Normally the topology table is kept in a text file but here we just write it to a temporary one
tf = tempfile.NamedTemporaryFile(delete=False)
tf.write(topology)
tf.close()

# The TopologyReader reads the text file, and the BuildSystem macro constructs it
mdl = IMP.Model()
reader = IMP.pmi.topology.TopologyReader(tf.name,
                                         pdb_dir = '../data/',
                                         fasta_dir = '../data/',
                                         gmm_dir = '../data/')
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
output_objects.append(evr)


# EM restraint
densities = IMP.atom.Selection(hier,representation_type=IMP.atom.DENSITIES).get_selected_particles()
gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(densities,target_fn="../data/32_imp.txt",
                                            slope=0.000001,
                                            target_radii_scale=3.0,
                                            target_is_rigid_body=False)
gem.set_label("EM")
gem.add_target_density_to_hierarchy(hier)
gem.write_target_gmm_to_mrc()
gem.add_to_model()


###################### SAMPLING #####################
# First shuffle the system
# Mind that you cannod shuffle the root hier, as it have the gaussian beads attached to
# by-state shuffling is safer

for state in bs.system.get_states():
    IMP.pmi.tools.shuffle_configuration(state.get_hierarchy(),
                                        max_translation=100)

# Quickly move all flexible beads into place
dof.optimize_flexible_beads(10)


# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                    root_hier=hier,                          # pass the root hierarchy
                                    #crosslink_restraints=xls,                     # will display like XLs
                                    monte_carlo_sample_objects=dof.get_movers(),  # pass MC movers
                                    global_output_directory='output/',
                                    output_objects=output_objects,
                                    monte_carlo_steps=10,
                                    number_of_best_scoring_models=0,      # set >0 to store best PDB files (but this is slow to do online)
                                    number_of_frames=10000)                   # increase number of frames to get better results!
rex.execute_macro()

