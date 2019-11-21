from __future__ import print_function
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.test
import IMP.pmi.macros

import IMP.bayesianem
import IMP.bayesianem.restraint

import math

def setup_gaussian(m):
    std=[1,1,1]
    center=[0,0,0]
    var=[s**2 for s in std]
    trans=IMP.algebra.Transformation3D(center)
    shape=IMP.algebra.Gaussian3D(IMP.algebra.ReferenceFrame3D(trans),var)
    p=IMP.Particle(m)
    IMP.core.Gaussian.setup_particle(p,shape)
    IMP.atom.Mass.setup_particle(p,1.0)
    IMP.core.XYZR.setup_particle(p)
    IMP.core.XYZR(p).set_radius(1.0)
    return p


def create_gem_ref(root_hier, target_ps, density_ps, label):
    gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(
        density_ps, target_ps=target_ps, slope=0.000001,
        target_radii_scale=3.0, target_is_rigid_body=True)
    gem.set_label(label)
    gem.add_target_density_to_hierarchy(root_hier)
    gem.add_to_model()
    """
    rb=gem.rb
    rbxyz = (rb.get_x(), rb.get_y(), rb.get_z())

    transformation = IMP.algebra.get_random_local_transformation(
        rbxyz,
        100,
        math.pi)

    IMP.core.transform(rb, transformation)
    """
    return gem

def create_gem_cross(target_ps,density_ps,label):
    gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(density_ps,target_ps=target_ps,
                                                slope=0.000001,
                                                target_radii_scale=3.0,
                                                target_is_rigid_body=False)
    gem.set_label(label)
    gem.add_to_model()
    return gem


class Tests(IMP.test.TestCase):
    def test_multiple(self):
        # setting up parameters

        rbmaxtrans = 4.00
        fbmaxtrans = 5.00
        rbmaxrot=0.05
        outputobjects = []
        sampleobjects = []

        # setting up topology

        m = IMP.Model()
        s = IMP.pmi.topology.System(m)
        root_hier = s.build()

        #em
        p11=setup_gaussian(m)
        p21=setup_gaussian(m)
        p31=setup_gaussian(m)
        p12=setup_gaussian(m)
        p22=setup_gaussian(m)
        p32=setup_gaussian(m)

        #root_hier.add_child(IMP.atom.Hierarchy.setup_particle(p1))
        #root_hier.add_child(IMP.atom.Hierarchy.setup_particle(p2))
        #root_hier.add_child(IMP.atom.Hierarchy.setup_particle(p3))

        IMP.atom.show_molecular_hierarchy(root_hier)


        gem=create_gem_ref(root_hier, [p11,p12],[p21,p22],'12')
        print(gem.get_output())
        outputobjects.append(gem)
        sampleobjects.append(gem)

        gem=create_gem_ref(root_hier, [p11,p12],[p31,p32],'13')
        print(gem.get_output())
        outputobjects.append(gem)
        sampleobjects.append(gem)

        gem=create_gem_cross([p21,p22],[p31,p32],'23')
        print(gem.get_output())
        outputobjects.append(gem)
        sampleobjects.append(gem)



        mc1=IMP.pmi.macros.ReplicaExchange0(
            m, root_hier=root_hier,
            monte_carlo_sample_objects=sampleobjects,
            output_objects=outputobjects,
            monte_carlo_temperature=1.0,
            simulated_annealing=None,
            simulated_annealing_minimum_temperature=1.0,
            simulated_annealing_maximum_temperature=20.0,
            simulated_annealing_minimum_temperature_nframes=100,
            simulated_annealing_maximum_temperature_nframes=100,
            replica_exchange_minimum_temperature=1.0,
            replica_exchange_maximum_temperature=100.0,
            number_of_best_scoring_models=0,
            monte_carlo_steps=10,
            number_of_frames=100000,
            write_initial_rmf=True,
            initial_rmf_name_suffix="initial",
            stat_file_name_suffix="stat",
            best_pdb_name_suffix="model",
            do_clean_first=True,
            do_create_directories=True,
            global_output_directory="output_shuffle_rex_xl",
            rmf_dir="rmfs/",
            best_pdb_dir="pdbs/",
            replica_stat_file_suffix="stat_replica")
        mc1.execute_macro()


if __name__ == '__main__':
    IMP.test.main()
