import IMP
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.io
import os
import glob

class ParticlesTool(object):
    def __init__(self,model,rmf_file_name,selection_tuples,nframe=0):
        self.model=model
        self.rmf_name=rmf_file_name
        self.nframe=nframe
        hiers,self.rs=IMP.pmi.analysis.get_hiers_and_restraints_from_rmf(self.model,nframe,rmf_file_name)
        self.hier=hiers[0]
        self.ps=[]
        for st in selection_tuples:
           if type(st) is str:
              s=IMP.atom.Selection(self.hier,molecule=st)
           elif type(st) is tuple:
              s=IMP.atom.Selection(self.hier,molecule=st[2],
                     residue_indexes=range(st[0],st[1]+1))
              self.ps+=s.get_selected_particles()        

    def get_coordinates(self):
        xyz=[IMP.core.XYZ(p).get_coordinates() for p in self.ps]
        return xyz

    def transform(self,transformation):
        rbs = set()
        for p in IMP.atom.get_leaves(self.hier):
            if not IMP.core.XYZR.get_is_setup(p):
                IMP.core.XYZR.setup_particle(p)
                IMP.core.XYZR(p).set_radius(0.0001)
                IMP.core.XYZR(p).set_coordinates((0, 0, 0))

            if IMP.core.RigidBodyMember.get_is_setup(p):
                rb = IMP.core.RigidBodyMember(p).get_rigid_body()
                rbs.add(rb)
            else:
                IMP.core.transform(IMP.core.XYZ(p),
                                   transformation)
        for rb in rbs:
            IMP.core.transform(rb,transformation)

    def save_rmf(self,rmf_name):
         o=IMP.pmi.output.Output()
         out_rmf_fn=os.path.join(rmf_name)
         o.init_rmf(out_rmf_fn,[self.hier],self.rs)
         o.write_rmf(out_rmf_fn)
         o.close_rmf(out_rmf_fn)

    def save_pdb(self,pdb_name):
         o=IMP.pmi.output.Output()
         out_pdb_fn=os.path.join(pdb_name)
         o.init_pdb(out_pdb_fn,self.hier)
         o.write_pdb(out_pdb_fn)




class ClusterTool(object):
    def __init__(self,GetModelDensity_custom_ranges):
        self.clusters={}
        self.cluster_directory={}
        self.density_objects={}
        self.gmdcr=GetModelDensity_custom_ranges
    
    def get_cluster_number(self):
        return len(self.clusters.keys())

    def add_cluster(self,member_ParticlesTool):
        nclusters=self.get_cluster_number()
        directory="cluster."+str(nclusters)
        if not os.path.exists(directory):
           os.makedirs(directory)
        else:
           print("Directory already existing")
           exit()
        self.cluster_directory[nclusters]=directory
        self.clusters[nclusters]=[member_ParticlesTool]
        self.density_objects[nclusters]=IMP.pmi.analysis.GetModelDensity(self.gmdcr,voxel=2)
        self.density_objects[nclusters].add_subunits_density(member_ParticlesTool.hier)

    def add_member(self,cluster_number,member_ParticlesTool):
        self.clusters[cluster_number].append(member_ParticlesTool)
        self.density_objects[cluster_number].add_subunits_density(member_ParticlesTool.hier)      

    def get_members(self,cluster_number):
        return self.clusters[cluster_number]

    def get_cluster_directory_name(self,cluster_number):
        return self.cluster_directory[cluster_number]





model = IMP.Model()

rmf_names=glob.glob("../output/rmfs/*.rmf3")
selection_tuples=["chainJ","chainK","chainC","chainL"]
em_density_map={"chainJ":["chainJ"],
                "chainK":["chainK"],
                "chainC":["chainC"],
                "chainL":["chainL"]}
#selection_tuples=[(1,500,"Nup82")]
ct=ClusterTool(em_density_map)
pt_dict={}
dist_cutoff=10

for n,rmf_name in enumerate(rmf_names):
    print n,rmf_name
    nframe=0
    while 1:
        try:
           pt_dict[(rmf_name,nframe)]=ParticlesTool(model,rmf_name,selection_tuples,nframe)
           pt_st=pt_dict[(rmf_name,nframe)]
           nframe+=1
        except:
           break

        if (ct.get_cluster_number() == 0):
            ct.add_cluster(pt_st)
            cdn=ct.get_cluster_directory_name(0)
            pt_st.save_rmf(cdn+"/"+str(n)+"."+str(nframe)+".rmf3")
            pt_st.save_pdb(cdn+"/"+str(n)+"."+str(nframe)+".pdb")
            print(rmf_name,"created new cluster",cdn)    
        else:
            xyz_st=pt_st.get_coordinates()
            assigned=False
            for cn in ct.clusters.keys():    
                cc=ct.get_members(cn)[0].rmf_name
                ccf=ct.get_members(cn)[0].nframe
                pt_cc=pt_dict[(cc,ccf)]
                xyz_cc=pt_cc.get_coordinates()
                transformation=IMP.algebra.IMP.algebra.get_transformation_aligning_first_to_second(xyz_st,xyz_cc)
                pt_st.transform(transformation)
                dist=IMP.algebra.get_rmsd(xyz_st, xyz_cc)
                print(dist)
                if dist < dist_cutoff:
                   cdn=ct.get_cluster_directory_name(cn)
                   pt_st.save_rmf(cdn+"/"+rmf_name)
                   pt_st.save_pdb(cdn+"/"+rmf_name+".pdb")
                   ct.add_member(cn,pt_st)
                   print(rmf_name,"RMSD",dist,"assigned to cluster",cdn)
                   assigned=True
                   break
                pt_st.transform(transformation.get_inverse())            
            if not assigned:
                ct.add_cluster(pt_st)
                cdn=ct.get_cluster_directory_name(ct.get_cluster_number()-1)
                pt_st.save_rmf(cdn+"/"+str(n)+"."+str(nframe)+".rmf3") 
                pt_st.save_pdb(cdn+"/"+str(n)+"."+str(nframe)+".pdb")
                print(rmf_name,"created new cluster",cdn)           
               
               
for cn in range(ct.get_cluster_number()):
    cdn=ct.get_cluster_directory_name(cn)    
    ct.density_objects[cn].write_mrc(cdn)
        
        

    

"""
first="0.0.rmf3"
second="199.0.rmf3"


pt_first=ParticlesTool(model,first,tuple_selections)
pt_second=ParticlesTool(model,second,tuple_selections)

xyz_first=pt_first.get_coordinates()
xyz_second=pt_second.get_coordinates()

dist=IMP.algebra.get_rmsd(xyz_second, xyz_first)
print("RMSD",dist)

transformation=IMP.algebra.IMP.algebra.get_transformation_aligning_first_to_second(xyz_first,xyz_second)

pt_first.transform(transformation)

xyz_first_new=pt_first.get_coordinates()
dist=IMP.algebra.get_rmsd(xyz_second, xyz_first_new)
print("RMSD",dist)

pt_first.save_rmf("first.rmf3")
"""
