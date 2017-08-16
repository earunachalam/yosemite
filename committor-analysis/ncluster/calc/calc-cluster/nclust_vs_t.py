#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import sys

# plot over time 
# 1. mean & max cluster size
# 2. number of icelike clusters
# 3. number of icelike mws (total in system across all clusters)

mean_size       = []
max_size        = []
num_clusters    = []
num_icelike     = []
frames          = range(0, 5000)

for f in frames:
    ifname = sys.argv[1] + '/f' + str(f) + '.dat'
    try:
        with open(ifname, 'r') as f:
            dat = [d.split() for d in f.readlines()]
    except:
        print("Syntax: nclust_vs_t.py cluster_stats_directory.")
        sys.exit()

    sizes = [float(len(member)) for member in dat]
    mean_size.append(np.mean(sizes))
    max_size.append(np.max(sizes))
    num_clusters.append(len(sizes))

    tot_icelike = 0
    for s in sizes: tot_icelike += s
    num_icelike.append(tot_icelike)

frames = [0.001*float(f) for f in frames]

plt.plot(frames,mean_size)
plt.xlabel(r'$t$/ns')
plt.ylabel(r'mean cluster size')
plt.savefig('meansize_v_t.pdf',bbox_inches='tight')
plt.close()

plt.plot(frames,max_size)
plt.xlabel(r'$t$/ns')
plt.ylabel(r'max cluster size')
plt.savefig('max_v_t.pdf',bbox_inches='tight')
plt.close()

plt.plot(frames,num_clusters)
plt.xlabel(r'$t$/ns')
plt.ylabel(r'number of icelike clusters')
plt.savefig('nClusters_v_t.pdf',bbox_inches='tight')
plt.close()

plt.plot(frames,num_icelike)
plt.xlabel(r'$t$/ns')
plt.ylabel(r'number of icelike mWs')
plt.savefig('nIceLike_v_t.pdf',bbox_inches='tight')
plt.close()
