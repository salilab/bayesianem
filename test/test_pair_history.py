import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.isd
import IMP.test
import IMP.bayesianem


class Tests(IMP.test.TestCase):

    def test_score_independent_of_history(self):
        """Test that the score does not depend on close pair list history"""
        m = IMP.Model()

        def gaussian(x, radius):
            p = IMP.Particle(m)
            IMP.core.Gaussian.setup_particle(p, IMP.algebra.Gaussian3D(
                IMP.algebra.ReferenceFrame3D(IMP.algebra.Transformation3D(
                    IMP.algebra.Vector3D(x, 0, 0))),
                IMP.algebra.Vector3D(9, 9, 9)))
            IMP.atom.Mass.setup_particle(p, 1000.)
            IMP.core.XYZR.setup_particle(p, radius)
            return p.get_index()

        def make_restraint():
            return IMP.bayesianem.GaussianEMRestraint(
                m, model, density, sigma, 10., 10., 0., False, False)
        model = [gaussian(0., 3.), gaussian(-36.5, 3.)]
        density = [gaussian(0., 9.), gaussian(-60., 9.)]
        sigma = IMP.isd.Scale.setup_particle(
            IMP.Particle(m), 1.).get_particle_index()
        r = make_restraint()
        r.evaluate(False)
        # Move the second model Gaussian by less than the close pair
        # container's slack (1A), so that its pair with the second density
        # Gaussian crosses the list edge without a list update
        IMP.core.XYZ(m, model[1]).set_x(-35.6)
        self.assertAlmostEqual(r.evaluate(False),
                               make_restraint().evaluate(False), delta=1e-8)


if __name__ == '__main__':
    IMP.test.main()
