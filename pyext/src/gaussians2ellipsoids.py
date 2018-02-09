import IMP
import IMP.isd
import IMP.isd.gmm_tools
import sys

gmmfile=sys.argv[1]
mdl=IMP.Model()
ps=[]
IMP.isd.gmm_tools.decorate_gmm_from_text(gmmfile,ps,mdl)


maxmass=max([IMP.atom.Mass(p).get_mass() for p in ps])

for p in ps:
    g=IMP.core.Gaussian(p)
    gg=g.get_gaussian()
    rf=gg.get_reference_frame()
    tr=rf.get_transformation_to()
    c=tr.get_translation()
    r=tr.get_rotation()
    q=r.get_quaternion()
    v=gg.get_variances()
    mass=IMP.atom.Mass(p).get_mass()
    nm=mass/maxmass
    cl=(nm,1-nm,0)
    transp=1.0
    s=('shape ellipsoid radius %f,%f,%f qrotation %f,%f,%f,%f center %f,%f,%f color %f,%f,%f,%f\n')%(v[0]/10,v[1]/10,v[2]/10,q[1],q[2],q[3],q[0],c[0],c[1],c[2],cl[0],cl[1],cl[2],1.0)
    print(s)
