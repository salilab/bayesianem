import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.plotting
import IMP.pmi.plotting.topology
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

###################### SYSTEM SETUP #####################
# Read sequences etc
#2 rigid bodies
topology='''
|  molecule_name  |  color   |  fasta_fn       |  fasta_id   |  pdb_fn          |  chain  |  residue_range  |  pdb_offset  |  bead_size  |  em_residues_per_gaussian  |  rigid_body  |  super_rigid_body  |  chain_of_super_rigid_bodies  |
|  Dis3           |  blue    |  exosome.fasta  |  Dis3       |  5g06_align_2.pdb  |  J      |  1,END          |  0           |  RESB       |  RESB                      |  1           |  0                 |  |
|  Rrp45          |  green   |  exosome.fasta  |  Rrp45      |  5g06_align_2.pdb  |  A      |  1,END          |  0           |  RESB       |  RESB                      |  4           |  0                 |  |
|  Rrp4           |  orange  |  exosome.fasta  |  Rrp4       |  5g06_align_2.pdb  |  H      |  1,102          |  0           |  RESB       |  RESB                      |  4           |  0                 |  |
|  Rrp4           |  yellow  |  exosome.fasta  |  Rrp4       |  5g06_align_2.pdb  |  H      |  103,END        |  0           |  RESB       |  RESB                      |  4           |  0                 |  |
|  Csl4           |  salmon  |  exosome.fasta  |  Csl4       |  5g06_align_2.pdb  |  I      |  1,END          |  0           |  RESB       |  RESB                      |  4           |  0                 |  |
|  Mtr3           |  gold    |  exosome.fasta  |  Mtr3       |  5g06_align_2.pdb  |  F      |  1,30           |  0           |  RESB       |  RESB                      |  4           |  0                 |  |
|  Mtr3           |  gold    |  exosome.fasta  |  Mtr3       |  5g06_align_2.pdb  |  F      |  31,END         |  0           |  RESB       |  RESB                      |  4           |  0                 |  |
|  Rrp40          |  pink    |  exosome.fasta  |  Rrp40      |  5g06_align_2.pdb  |  G      |  1,60           |  0           |  RESB       |  RESB                      |  4           |  0                 |  |
|  Rrp40          |  pink    |  exosome.fasta  |  Rrp40      |  5g06_align_2.pdb  |  G      |  61,END         |  0           |  RESB       |  RESB                      |  4           |  0                 |  |
|  Rrp42          |  red     |  exosome.fasta  |  Rrp42      |  5g06_align_2.pdb  |  E      |  1,END          |  0           |  RESB       |  RESB                      |  4           |  0                 |  |
|  Ski6           |  white   |  exosome.fasta  |  Ski6       |  5g06_align_2.pdb  |  B      |  1,END          |  0           |  RESB       |  RESB                      |  4           |  0                 |  |
|  Rrp46_gfp      |  purple  |  exosome.fasta  |  Rrp46_gfp  |  5g06_align_2.pdb  |  D      |  1,246          |  0           |  RESB       |  RESB                      |  7           |  0,2               |  |
|  Rrp43          |  gray    |  exosome.fasta  |  Rrp43      |  5g06_align_2.pdb  |  C      |  1,END          |  0           |  RESB       |  RESB                      |  6           |  0                 |  |
'''.replace("RESB",str(20))


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
hier, dof = bs.execute_macro()
IMP.pmi.plotting.topology.draw_component_composition(dof)

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
                                                 target_fn='../data/emd_3367_256.gmm',
                                                 scale_target_to_mass=True,
                                                 slope=0.01,
                                                 target_radii_scale=3.0,
                                                 target_is_rigid_body=False)
gem.add_to_model()
gem.set_label("Total")
output_objects.append(gem)

#from IMP.pmi.io.crosslink import FilterOperator as FO
#cldbkc=IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
#cldbkc.set_protein1_key("Protein 1")
#cldbkc.set_protein2_key("Protein 2")
#cldbkc.set_residue1_key("Residue1")
#cldbkc.set_residue2_key("Residue2")
#cldbkc.set_id_score_key("Score")
#cldbkc.set_unique_id_key("Unique ID")
#cldb=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
#cldb.create_set_from_file("../data/exosome_XLMS_column07012014.csv")
#
#cldb.offset_residue_index('GFP', 246)
#cldb.rename_proteins({"GFP":"Rrp46_gfp","Rrp46":"Rrp46_gfp"})
#
#
#xl = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=hier,
#                            CrossLinkDataBase=cldb,
#                            length=21.0,
#                            slope=0.02,
#                            resolution=1.0,
#                            label="XL")
#xl.add_to_model()
#output_objects.append(xl)
#xl.set_psi_is_sampled(False)
#psi=xl.psi_dictionary["PSI"][0]
#psi.set_scale(0.01)



output_dir='.'
dof.optimize_flexible_beads(100)
score=IMP.pmi.tools.get_restraint_set(mdl).unprotected_evaluate(None)
em_score=gem.rs.unprotected_evaluate(None)
print("score=%e"%(score))
print("em score=%e"%(em_score))
with open("%s/ref.rmf.score"%(output_dir), "w") as f:
    f.write("Total_Score=%e\n"%(score))
    f.write("EM_Score=%e\n"%(em_score))
rh = RMF.create_rmf_file("%s/ref.rmf"%(output_dir))
IMP.rmf.add_hierarchies(rh, [hier])
IMP.rmf.save_frame(rh)
del rh
