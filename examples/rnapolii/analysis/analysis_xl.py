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

for nc,cluster in enumerate(clusters[:1]):
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
        k=' '.join(k.split('|'))
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
