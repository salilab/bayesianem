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
import IMP.bayesianem.movers
import tempfile,os
import sys
import numpy as np
import time

try:
    import IMP.mpi
    print('ReplicaExchange: MPI was found. Using Parallel Replica Exchange')
    rex_obj = IMP.mpi.ReplicaExchange()
except ImportError:
    print('ReplicaExchange: Could not find MPI. Using Serial Replica Exchange')
    rex_obj = _SerialReplicaExchange()

replica_number = rex_obj.get_my_index()

#time.sleep(replica_number)


output_prefix=""
if len(sys.argv)>1:
    output_prefix=sys.argv[1]

step=int(sys.argv[2])

output_objects=[]


###################### SYSTEM SETUP #####################
# Read sequences etc
topology='''
|  molecule_name  |  color   |  fasta_fn       |  fasta_id   |  pdb_fn                      |  chain  |  residue_range  |  pdb_offset  |  bead_size  |  em_residues_per_gaussian  |  rigid_body  |  super_rigid_body  |  chain_of_super_rigid_bodies  |
|  Dis3           |  blue    |  exosome.fasta  |  Dis3       |  4IFD.pdb                    |  J      |  1,237          |  0           |  RESB       |  RESB                      |  1           |  0,1               |  |
|  Dis3           |  blue    |  exosome.fasta  |  Dis3       |  4IFD.pdb                    |  J      |  238,471        |  0           |  RESB       |  RESB                      |  2           |  0,1               |  |
|  Dis3           |  blue    |  exosome.fasta  |  Dis3       |  4IFD.pdb                    |  J      |  472,END        |  0           |  RESB       |  RESB                      |  3           |  0,1               |  |
|  Rrp45          |  green   |  exosome.fasta  |  Rrp45      |  4IFD.pdb                    |  A      |  1,END          |  0           |  RESB       |  RESB                      |  4           |  0,2               |  |
|  Rrp4           |  orange  |  exosome.fasta  |  Rrp4       |  4IFD.pdb                    |  H      |  1,102          |  0           |  RESB       |  RESB                      |  5           |  0,3               |  |
|  Rrp4           |  orange  |  exosome.fasta  |  Rrp4       |  4IFD.pdb                    |  H      |  103,END        |  0           |  RESB       |  RESB                      |  6           |  0,3               |  |
|  Csl4           |  salmon  |  exosome.fasta  |  Csl4       |  4IFD.pdb                    |  I      |  1,END          |  0           |  RESB       |  RESB                      |  7           |  0,4               |  |
|  Mtr3           |  gold    |  exosome.fasta  |  Mtr3       |  4IFD.pdb                    |  F      |  1,30           |  0           |  RESB       |  RESB                      |  8           |  0,5               |  |
|  Mtr3           |  gold    |  exosome.fasta  |  Mtr3       |  4IFD.pdb                    |  F      |  31,END         |  0           |  RESB       |  RESB                      |  9           |  0,5               |  |
|  Rrp40          |  pink    |  exosome.fasta  |  Rrp40      |  4IFD.pdb                    |  G      |  1,60           |  0           |  RESB       |  RESB                      |  10          |  0,6               |  |
|  Rrp40          |  pink    |  exosome.fasta  |  Rrp40      |  4IFD.pdb                    |  G      |  61,END         |  0           |  RESB       |  RESB                      |  11          |  0,6               |  |
|  Rrp42          |  red     |  exosome.fasta  |  Rrp42      |  4IFD.pdb                    |  E      |  1,END          |  0           |  RESB       |  RESB                      |  12          |  0,7               |  |
|  Ski6           |  white   |  exosome.fasta  |  Ski6       |  4IFD.pdb                    |  B      |  1,END          |  0           |  RESB       |  RESB                      |  13          |  0,8               |  |
|  Rrp46_gfp      |  purple  |  exosome.fasta  |  Rrp46_gfp  |  4IFD.pdb                    |  D      |  1,246          |  0           |  RESB       |  RESB                      |  14          |  0,9               |  |
|  Rrp46_gfp      |  purple  |  exosome.fasta  |  Rrp46_gfp  |  GFP_1GFL.pdb                |  A      |  1,229          |  246         |  RESB       |  RESB                      |  16          |  0,9               |  |
|  Rrp43          |  gray    |  exosome.fasta  |  Rrp43      |  4IFD.pdb                    |  C      |  1,END          |  0           |  RESB       |  RESB                      |  15          |  0,10              |  |
'''.replace("RESB",str(10))

# Normally the topology table is kept in a text file but here we just write it to a temporary one
tf = tempfile.NamedTemporaryFile(delete=False)
tf.write(topology)
tf.close()

# The TopologyReader reads the text file, and the BuildSystem macro constructs it
mdl = IMP.Model()
reader = IMP.pmi.topology.TopologyReader(tf.name,
                                         pdb_dir = '../data/',
                                         fasta_dir = '../data/',
                                         gmm_dir = 'gmm/')
bs = IMP.pmi.macros.BuildSystem(mdl)
bs.add_state(reader) # note you can call this multiple times to create a multi-state system
hier, dof = bs.execute_macro(max_rb_trans=4.0, max_rb_rot=0.1, max_bead_trans=4.0, max_srb_trans=4.0,max_srb_rot=0.1)

#pm=IMP.bayesianem.movers.PCAMover(mdl,hier,resolution=10)

# Connectivity keeps things connected along the backbone (ignores if inside same rigid body)
crs = []
moldict = bs.get_molecules()[0]
mols = []
print moldict.keys()
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
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols)
evr.add_to_model()
evr.set_label("complex-complex")
evr.set_weight(0.1)
output_objects.append(evr)

densities=[]
molecules_in_em = ['Dis3', 'Rrp45', 'Rrp4', 'Csl4', 'Mtr3', 'Rrp40', 'Rrp42', 'Ski6', 'Rrp43']

for molname in molecules_in_em:
    sel = IMP.atom.Selection(hierarchy=hier, molecule=molname, representation_type=IMP.atom.DENSITIES)
    densities+=sel.get_selected_particles()
sel=IMP.atom.Selection(hierarchy=hier, molecule='Rrp46_gfp', residue_indexes=range(1,246+1),representation_type=IMP.atom.DENSITIES)
densities+=sel.get_selected_particles()

gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(densities,
                                                 target_fn='../data/%d_imp.gmm'%step,
                                                 scale_target_to_mass=True,
                                                 slope=0.01,
                                                 target_radii_scale=3.0,
                                                 target_is_rigid_body=False)
gem.add_to_model()
gem.set_label("Total")
output_objects.append(gem)

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


xl = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=hier,
                            CrossLinkDataBase=cldb,
                            length=21.0,
                            slope=0.02,
                            resolution=1.0,
                            label="XL")
xl.add_to_model()
output_objects.append(xl)
xl.set_psi_is_sampled(False)
psi=xl.psi_dictionary["PSI"][0]
psi.set_scale(0.01)

###################### SAMPLING #####################
# First shuffle the system
# Mind that you cannod shuffle the root hier, as it have the gaussian beads attached to
# by-state shuffling is safer




if step==1:
    for state in bs.system.get_states():
        IMP.pmi.tools.shuffle_configuration(state.get_hierarchy(),
                                            max_translation=500)
else:
    rh_ref = RMF.open_rmf_file_read_only('%sseed_%d.rmf3'%(output_prefix,step-1))
    IMP.rmf.link_hierarchies(rh_ref, [hier])
    IMP.rmf.load_frame(rh_ref, RMF.FrameID(replica_number))




# Quickly move all flexible beads into place
dof.optimize_flexible_beads(10)

output_dir='.'
global_output_directory='%s/%s_output_%d/'%(output_dir,output_prefix, step)

nf=1000

rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                    root_hier=hier,                          # pass the root hierarchy
                                    monte_carlo_sample_objects=dof.get_movers(),  # pass MC movers
                                    global_output_directory=global_output_directory,
                                    crosslink_restraints=[xl],
                                    output_objects=output_objects,
                                    monte_carlo_steps=10,
                replica_exchange_minimum_temperature=1.0,
                                    replica_exchange_maximum_temperature=5.0,
                            replica_exchange_swap=True,
                save_coordinates_mode="25th_score",
                                    number_of_best_scoring_models=0,      # set >0 to store best PDB files (but this is slow to do online)
                                    number_of_frames=nf,                   # increase number of frames to get better results!
                                    replica_exchange_object=rex_obj)
rex.execute_macro()
