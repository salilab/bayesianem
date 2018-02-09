import IMP
import IMP.isd
import IMP.isd.gmm_tools
import IMP.bayesianem
import IMP.bayesianem.restraint

import sys
import os
import glob

m = IMP.Model()
init=False
plotscores=False
ngauss=30


def generate_data(num_gaussian,pdbfile,chain_names,outputdir="./",voxel_size=2.0):

    fit_coords = []

    for chain in chain_names:
        sel = IMP.atom.AndPDBSelector(IMP.atom.ChainPDBSelector(chain),
                                    IMP.atom.NonWaterNonHydrogenPDBSelector())

        t = IMP.atom.read_pdb(pdbfile, m, sel)
        #IMP.atom.show_molecular_hierarchy(t)

        fit_coords += [IMP.core.XYZ(p).get_coordinates()
                      for p in IMP.atom.get_leaves(t)]

    density_ps = []
    (score,akaikescore)=IMP.isd.gmm_tools.fit_gmm_to_points(fit_coords,
                                        num_gaussian,
                                        m,
                                        density_ps,
                                        min_covar=2.0)
    IMP.isd.gmm_tools.write_gmm_to_text(
        density_ps, outputdir + str(num_gaussian) + "_" + str(voxel_size) + '.txt')
    IMP.isd.gmm_tools.write_gmm_to_map(
        density_ps, outputdir + str(num_gaussian) + "_" + str(voxel_size) + '.mrc', voxel_size)
    
    return(akaikescore)


if init:
    #for f in glob.glob('./input/bayesian_em_restraint/gmm_data/*'):
    #    os.remove(f)
    generate_data(ngauss,"./input/1WCM.pdb",["C","K","L","J"],outputdir='./input/bayesian_em_restraint/gmm_data/')
    exit()


import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.crosslinking_new
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros
import IMP.pmi.io


import os

# setting up parameters

rbmaxtrans = 4.00
fbmaxtrans = 4.00
rbmaxrot=0.015
outputobjects = []
sampleobjects = []

# setting up topology

m = IMP.Model()
simo1 = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=False)
fastadirectory="./input"
pdbdirectory="./input/"

compactrepresentation=False


       # compname  hier_name    color         fastafile              fastaid          pdbname      chain    resrange      read    "BEADS"ize rigid_body super_rigid_body emnum_components emtxtfilename  emmrcfilename chain of super rigid bodies
domains=   [
    ("chainC",  "chainC",    0.0,  fastadirectory+"/1WCM.fasta.txt",  "1WCM:C", pdbdirectory+"/1WCM.pdb" ,   "C",   (3,268,0),  True,        10,      0,         [1],     -10,   './input/bayesian_em_restraint/C.txt',  './input/bayesian_em_restraint/C.mrc',   None),
    ("chainK",  "chainK",    0.10,  fastadirectory+"/1WCM.fasta.txt",  "1WCM:K", pdbdirectory+"/1WCM.pdb" ,   "K",   (1,114,0), True,        10,      1,         [1],    -10,   './input/bayesian_em_restraint/K.txt',  './input/bayesian_em_restraint/K.mrc',   None),
    ("chainL",  "chainL",    0.20,  fastadirectory+"/1WCM.fasta.txt",  "1WCM:L", pdbdirectory+"/1WCM.pdb" ,   "L",   (25,70,0), True,        10,      2,         [1],    -10,   './input/bayesian_em_restraint/L.txt',  './input/bayesian_em_restraint/L.mrc',   None),
    ("chainJ",  "chainJ",    0.30,  fastadirectory+"/1WCM.fasta.txt",  "1WCM:J", pdbdirectory+"/1WCM.pdb" ,   "J",   (1,65,0),  True,        10,      3,         [1],    -10,   './input/bayesian_em_restraint/J.txt',  './input/bayesian_em_restraint/J.mrc',   None),
    ]


bm1=IMP.pmi.macros.BuildModel1(simo1)
bm1.build_model(domains,sequence_connectivity_scale=1.0,sequence_connectivity_resolution=1)


simo1.set_rigid_bodies_max_rot(rbmaxrot)
simo1.set_floppy_bodies_max_trans(fbmaxtrans)
simo1.set_rigid_bodies_max_trans(rbmaxtrans)
simo1.set_current_coordinates_as_reference_for_rmsd()

outputobjects.append(simo1)
sampleobjects.append(simo1)


ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(simo1,resolution=10)
ev.set_weight(0.5)
ev.add_to_model()
outputobjects.append(ev)


#em


resdensities=bm1.get_density_hierarchies([t[1] for t in domains])
gemh = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(resdensities,"input/bayesian_em_restraint/gmm_data/"+str(ngauss)+"_2.0.txt",
                                                slope=0.0005,
                                                target_radii_scale=3.0)

gemh.add_to_model()
outputobjects.append(gemh)
#sampleobjects.append(gemh)
"""
target_ps = []
IMP.isd.gmm_tools.decorate_gmm_from_text(
                "input/bayesian_em_restraint/gmm_data/100_2.0.txt",
                target_ps,
                m)


resdensities=bm1.get_density_hierarchies([t[1] for t in domains])
gemh = BayesianPointWiseGaussianRestraint(m,resdensities,target_ps)


print(gemh.unprotected_evaluate(False))
gemh.add_to_model()


outputobjects.append(gemh)
"""

simo1.shuffle_configuration(50)

mc1=IMP.pmi.macros.ReplicaExchange0(m,
                                    simo1,
                                    monte_carlo_sample_objects=sampleobjects,
                                    output_objects=outputobjects,
                                    monte_carlo_temperature=1.0,
                 		            simulated_annealing=False,
                                    simulated_annealing_minimum_temperature=1.0,
                                    simulated_annealing_maximum_temperature=2.5,
                                    simulated_annealing_minimum_temperature_nframes=100,
                                    simulated_annealing_maximum_temperature_nframes=10,
                                    replica_exchange_minimum_temperature=1.0,
                                    replica_exchange_maximum_temperature=2.5,
                                    number_of_best_scoring_models=10,
                                    monte_carlo_steps=10,
                                    number_of_frames=100000,
                                    write_initial_rmf=True,
                                    initial_rmf_name_suffix="initial",
                                    stat_file_name_suffix="stat",
                                    best_pdb_name_suffix="model",
                                    do_clean_first=True,
                                    do_create_directories=True,
                                    global_output_directory="./input/bayesian_em_restraint/output",
                                    rmf_dir="rmfs/",
                                    best_pdb_dir="pdbs/",
                                    replica_stat_file_suffix="stat_replica")
mc1.execute_macro()




    
