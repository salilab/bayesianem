import sys
import IMP
import IMP.pmi
import IMP.pmi.output
import numpy
from collections import defaultdict
import matplotlib.pyplot as plt


import pickle

clusterfile=sys.argv[1]
fl=open(clusterfile,"rb")
clusters=pickle.load(fl)

datas=[]

for nc,cluster in enumerate(clusters):
    distances=defaultdict(list)
    found=False
    for model in cluster:
        for k in model.features.keys():
            if "Distance" in k:
                distances[k].append(float(model.features[k]))
                found=True
        if found:
            break

    print '#cluster', nc
    for k,vs in distances.items():
        k=' '.join(k.split('|')[3:7])
        k=k.replace('_gfp','')
        for v in vs:
            print k,v

    distances_tuple=(list(),list(),list())
    [(distances_tuple[0].append(distances[k]),distances_tuple[1].append(k),distances_tuple[2].append(numpy.median(distances[k]))) for k in distances]

    distances=[]
    keys=[]
    medians=[]
    for n,dists in enumerate(distances_tuple[0]):
        if any([dists[0]!=x for x in dists]):
            distances.append(dists)
            key=distances_tuple[1][n]
            key=' '.join(key.split('|')[3:7])
            key=key.replace('_gfp','')
            keys.append(key)
            medians.append(distances_tuple[2][n])

    distances_sorted_by_median=[x for _,x in sorted(zip(medians,distances))]
    datas.append(distances_sorted_by_median)
    keys_sorted_by_median=[x for _,x in sorted(zip(medians,keys))]
#    IMP.pmi.output.plot_fields_box_plots("XL_Distances_"+str(cluster.cluster_id),
#                                         distances_sorted_by_median,
#                                         range(len(distances_sorted_by_median)),
#                                         xlabels=keys_sorted_by_median,scale_plot_length=0.2, valuename='crosslink length')
#

#f, ax = plt.subplots(figsize=(18, 7))
#v1 = ax.violinplot(datas[0], points=50, positions=np.arange(0, len(datas[0])), widths=0.85,
#               showmeans=False, showextrema=False, showmedians=False)
#for b in v1['bodies']:
#    m = np.mean(b.get_paths()[0].vertices[:, 0])
#    b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], -np.inf, m)
#    b.set_color('r')
#v2 = ax.violinplot(datas[1], points=50, positions=np.arange(0, len(datas[1])), widths=0.85,
#               showmeans=False, showextrema=False, showmedians=False)
#for b in v2['bodies']:
#    m = np.mean(b.get_paths()[0].vertices[:, 0])
#    b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], m, np.inf)
#    b.set_color('b')
#
#plt.show()
