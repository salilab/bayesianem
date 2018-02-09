import IMP
import IMP.core
import IMP.algebra
import itertools
import math
import random
import numpy as np
import scipy.spatial
import networkx as nx
from IMP.pmi.tools import get_restraint_set


class MoverScheme(object):
    def __init__(self, length):
        self.length = length
        self.mover_scheme = None
        self.score_scheme = {}
        self.score_initial_weights = {}
        self.position = 0
        self.numbermc = 0

    def add_mover_scheme(self, mover_scheme):
        if len(mover_scheme) != self.length:
            raise ("Error")
        else:
            mc_scheme = []
            for mvs in mover_scheme:
                smv = IMP.core.SerialMover(mvs)
                mc = IMP.core.MonteCarlo(mvs[0].get_model())
                mc.set_scoring_function(get_restraint_set(mvs[0].get_model()))
                mc.set_return_best(False)
                mc.set_kt(1.0)
                mc.add_mover(smv)
                mc_scheme.append(mc)
                mc.set_name(str(self.numbermc))
                self.numbermc += 1
            self.mover_scheme = mc_scheme

    def add_montecarlo_temperatures(self, temps):
        self.temps = temps

    def add_montecarlo_steps(self, steps):
        self.steps = steps

    def add_score_scheme(self, score, weigth_scheme):
        if len(weigth_scheme) != self.length:
            raise ("Error")
        else:
            self.score_scheme[score] = weigth_scheme
        self.score_initial_weights[score] = score.get_weight()

    def get_length(self):
        return self.length

    def get_position(self, position):
        movers = self.mover_scheme[position]
        score_tmp_dict = {}
        for score in self.score_scheme:
            score_tmp_dict[score] = self.score_scheme[score][position]
        return movers, score_tmp_dict

    def get_and_update(self):
        mc, score_tmp_dict = self.get_position(self.position)
        temp = self.temps[self.position]
        mc.set_kt(temp)
        steps = self.steps[self.position]
        self.position += 1
        if self.position > self.length - 1: self.position = 0
        return mc, score_tmp_dict, steps


class HybridMonteCarloMover(IMP.core.MonteCarloMover):
    def __init__(self, m, hier, MoverScheme):
        IMP.core.MonteCarloMover.__init__(
            self,
            m,
            "HybridMonteCarloMover")

        self.mols = IMP.pmi.tools.get_molecules(hier)
        self.refframes = {}
        self.xyzs = {}
        self.ps = []
        self.m = m
        self.rejects = []
        self.lenrejects = 20

        self.rbs, self.beads = IMP.pmi.tools.get_rbs_and_beads([hier])
        for rb in self.rbs:
            self.ps.append(rb.get_particle())
        for bead in self.beads:
            self.ps.append(bead.get_particle())

        self.indexes = [p.get_index() for p in self.ps]

        self.do_store_current_coordinates()

        self.mover_scheme = MoverScheme

    def do_store_current_coordinates(self):
        for rb in self.rbs:
            self.refframes[rb] = rb.get_reference_frame()
        for bead in self.beads:
            self.xyzs[bead] = IMP.core.XYZ(bead).get_coordinates()

    def do_propose(self):

        self.rejects.append(0)
        if len(self.rejects) > self.lenrejects:
            self.rejects.pop(0)
        print
        self.rejects, float(sum(self.rejects)) / len(self.rejects)

        self.do_store_current_coordinates()

        rs = get_restraint_set(self.m)
        # print "Initial",rs.evaluate(False)

        for i in range(self.mover_scheme.get_length()):
            mc, score_tmp_dict, steps = self.mover_scheme.get_and_update()

            mc.reset_statistics()

            for score in score_tmp_dict:
                score.set_weight(score_tmp_dict[score])

            mc.optimize(steps)

            print
            mc.get_name(), float(mc.get_number_of_accepted_steps()) / mc.get_number_of_proposed_steps()
        # print i,steps,rs.evaluate(False)

        # print "Final",rs.evaluate(False)

        return IMP.core.MonteCarloMoverResult(self.indexes, 1.0)

    def do_reject(self):
        print
        "REJECT"
        self.rejects[-1] = 1
        for rb in self.refframes:
            rb.set_reference_frame(self.refframes[rb])
        for bead in self.xyzs:
            IMP.core.XYZ(bead).set_coordinates(self.xyzs[bead])

    def do_get_inputs(self):
        return self.ps


class PCAMover(IMP.core.MonteCarloMover):
    def __init__(self, m, hier, resolution=10, depth=4):

        IMP.core.MonteCarloMover.__init__(
            self,
            m,
            "PCAMover")
        self.depth = depth
        self.propose = 0
        self.reject = 0
        self.rejects = []
        self.lenrejects = 20
        self.mols = IMP.pmi.tools.get_molecules(hier)
        self.resolution = resolution
        self.rbsd = {}
        self.beadsd = {}
        self.refframes = {}
        self.xyzs = {}
        self.ps = []
        self.particle_blocks = []
        for mol in self.mols:
            rbs, beads = IMP.pmi.tools.get_rbs_and_beads([mol])
            self.rbsd[mol] = rbs
            self.beadsd[mol] = beads
            for rb in self.rbsd[mol]:
                self.ps.append(rb.get_particle())
            for bead in self.beadsd[mol]:
                self.ps.append(bead.get_particle())
            ps = IMP.atom.Selection(mol, resolution=self.resolution).get_selected_particles()
            self.particle_blocks.append(ps)

        self.do_store_current_coordinates()

    def do_store_current_coordinates(self):
        for mol in self.mols:
            for rb in self.rbsd[mol]:
                self.refframes[rb] = rb.get_reference_frame()
            for bead in self.beadsd[mol]:
                self.xyzs[bead] = IMP.core.XYZ(bead).get_coordinates()

    def get_full_graph(self):
        '''
        get the full graph of distances between every particle pair
        '''

        pdist_array = np.array(
            IMP.pmi.get_list_of_bipartite_minimum_sphere_distance(self.particle_blocks))
        pdist_mat = scipy.spatial.distance.squareform(pdist_array)
        pdist_mat[pdist_mat <= 10] = 1
        pdist_mat[pdist_mat > 10] = 0
        graph = nx.Graph(pdist_mat)
        return graph

    def get_list_of_interacting_molecules(self, succ_dict, pointer):
        if pointer in succ_dict:
            paths = [path for p in succ_dict[pointer] for path in self.get_list_of_interacting_molecules(succ_dict, p)]
            ret = [[pointer] + path for path in paths]
            # print pointer,paths,ret
            return ret
        else:
            return [[pointer]]

    def get_minimum_spanning_tree(self):
        """
        return the minimum spanning tree
        """
        graph = self.get_full_graph()
        graph = nx.minimum_spanning_tree(graph)
        return graph

    def get_sublists(self, lst, size):
        result = [sublist for sublist in
                  (lst[x:x + size] for x in range(len(lst) - size + 1))]
        return result

    def get_list_of_connected_mols(self, depth):
        gr = self.get_minimum_spanning_tree()
        all_groups = set()
        if len(gr.edges()) > 0:

            for node in gr.nodes():
                succ_dict = nx.dfs_successors(gr, node)
                paths = [path for p in succ_dict[node] for path in self.get_list_of_interacting_molecules(succ_dict, p)]
                ret = [[node] + path for path in paths]
                for path in ret:
                    for s in range(2, min([depth + 1, len(path) + 1, len(self.mols)])):
                        subls = self.get_sublists(path, s)
                        for subl in subls: all_groups.add(frozenset(subl))
        return all_groups
        # succ=nx.bfs_successors(gr,0)

    # succmol=dict([(self.mols[i].get_name()+str(IMP.atom.Copy(self.mols[i]).get_copy_index()),[self.mols[k].get_name()+str(IMP.atom.Copy(self.mols[i]).get_copy_index()) for k in ks]) for i,ks in succ.iteritems()])
    # print(succmol)

    def do_compute_pcas(self):
        self.singletons = {}
        for i, ps in enumerate(self.particle_blocks):
            mol = self.mols[i]
            pca = IMP.algebra.get_principal_components([IMP.core.XYZ(p).get_coordinates() for p in ps])
            ## MB normalize
            pcan = IMP.bayesianem.NormalizePCA(pca, ps)
            self.singletons[(mol,)] = pcan

        if self.depth > 1:
            all_groups = self.get_list_of_connected_mols(self.depth)

            for group in all_groups:
                mols = [self.mols[n] for n in group]
                ps = IMP.atom.Selection(mols, resolution=self.resolution).get_selected_particles()
                pca = IMP.algebra.get_principal_components([IMP.core.XYZ(p).get_coordinates() for p in ps])
                pcan = IMP.bayesianem.NormalizePCA(pca, ps)
                self.singletons[tuple(mols)] = pcan

    def do_get_pairs(self, temp=20):
        sk = self.singletons.keys()

        pairs = []
        for s1, s2 in itertools.combinations(sk, 2):
            rms = IMP.algebra.get_distance(self.singletons[s1].get_principal_values(),
                                           self.singletons[s2].get_principal_values())
            pairs.append(((s1, s2), math.exp(-rms / temp)))

        return pairs

    def do_compute_transforms(self, pair):
        amols = pair[0]
        bmols = pair[1]

        arbs = []
        abeads = []
        for mol in amols:
            arbs += self.rbsd[mol]
            abeads += self.beadsd[mol]
        apca = self.singletons[amols]

        brbs = []
        bbeads = []
        for mol in bmols:
            brbs += self.rbsd[mol]
            bbeads += self.beadsd[mol]
        bpca = self.singletons[bmols]

        # transatob=IMP.algebra.get_alignments_from_first_to_second(apca,bpca)
        # MB get single transformation
        transatob = IMP.bayesianem.PCAalign(apca, bpca)[0]
        transbtoa = transatob.get_inverse()

        return (transatob, transbtoa, arbs, abeads, brbs, bbeads)

    def weighted_choice(self, choices):
        total = sum(w for c, w in choices)
        r = random.uniform(0, total)
        upto = 0
        for c, w in choices:
            if upto + w >= r:
                return c
            upto += w
        assert False, "Shouldn't get here"

    def do_propose(self):
        self.rejects.append(0)
        if len(self.rejects) > self.lenrejects:
            self.rejects.pop(0)

        self.propose += 1
        self.do_store_current_coordinates()
        self.do_compute_pcas()
        pairs = self.do_get_pairs()
        pair = self.weighted_choice(pairs)
        # print pair,self.reject,self.propose,float(self.propose-self.reject)/self.propose,float(sum(self.rejects))/len(self.rejects)
        # print self.rejects
        t = self.do_compute_transforms(pair)
        # (ta2b,tb2a)=random.choice(zip(t[0],t[1]))
        # MB no need to randomly select
        (ta2b, tb2a) = (t[0], t[1])

        moved_indexes = []
        for arb in t[2]:
            IMP.core.transform(arb, ta2b)
            moved_indexes.append(arb.get_particle_index())
        for abead in t[3]:
            IMP.core.transform(IMP.core.XYZ(abead), ta2b)
            moved_indexes.append(arb.get_particle_index())
        for brb in t[4]:
            IMP.core.transform(brb, tb2a)
        for bbead in t[5]:
            IMP.core.transform(IMP.core.XYZ(bbead), tb2a)
        return IMP.core.MonteCarloMoverResult(moved_indexes, 1.0)

    def do_reject(self):
        self.rejects[-1] = 1

        self.reject += 1
        for rb in self.refframes:
            rb.set_reference_frame(self.refframes[rb])
        for bead in self.xyzs:
            IMP.core.XYZ(bead).set_coordinates(self.xyzs[bead])

    def do_get_inputs(self):
        return self.ps
