import IMP
import IMP.test
import IMP.pmi
import IMP.isd
import IMP.bayesianem
import IMP.bayesianem.restraint


class Tests(IMP.test.TestCase):
    def make_restraint(self):
        mdl = IMP.Model()
        s = IMP.pmi.topology.System(mdl)
        st1 = s.create_state()
        hier = s.build()

        fname = self.get_input_file_name('2A73.pdb50.txt')
        target_ps = []
        IMP.isd.gmm_tools.decorate_gmm_from_text(
            fname, target_ps, mdl, radius_scale=3.0, mass_scale=1.0)
        gemh = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(
            target_ps, fname, target_mass_scale=1.0, slope=0.000001,
            target_radii_scale=3.0)
        gemh.set_label("Mobile")
        gemh.add_target_density_to_hierarchy(st1)
        gemh.add_to_model()
        gemh.set_weight(100.0)
        return mdl, gemh

    def test_info(self):
        """Test key:value restraint info"""
        mdl, gemh = self.make_restraint()
        score = gemh.rs.evaluate(False)
        gem = gemh.gaussianEM_restraint
        s = gem.get_static_info()
        s.set_was_used(True)
        # No density_fn set
        self.assertEqual(s.get_number_of_filename(), 1)
        self.assertEqual(s.get_filename_key(0), "filename")
        self.assertEqual(s.get_filename_value(0),
            self.get_input_file_name('2A73.pdb50.txt'))
        self.assertEqual(s.get_number_of_string(), 1)
        self.assertEqual(s.get_string_key(0), "type")
        self.assertEqual(s.get_string_value(0),
            "IMP.bayesianem.GaussianEMRestraint")

        self.assertEqual(s.get_number_of_int(), 0)
        self.assertEqual(s.get_number_of_float(), 2)
        self.assertEqual(s.get_float_key(0), "model cutoff")
        self.assertAlmostEqual(s.get_float_value(0), 10.0, delta=1e-4)
        self.assertEqual(s.get_float_key(1), "density cutoff")
        self.assertAlmostEqual(s.get_float_value(1), 10.0, delta=1e-4)

        s = gem.get_dynamic_info()
        s.set_was_used(True)
        self.assertEqual(s.get_number_of_float(), 1)
        self.assertEqual(s.get_float_key(0), "global sigma")
        self.assertAlmostEqual(s.get_float_value(0), 1.0, delta=1e-4)


if __name__ == '__main__':
    IMP.test.main()
