#!/usr/bin/python

import numpy as np
import scipy
import kmapper as km
import sklearn
from sklearn import ensemble
from sklearn.datasets import load_digits
from sklearn.manifold import Isomap
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.decomposition import PCA
import matplotlib as mpl
import matplotlib.pyplot as plt
import sklearn
from sklearn import datasets
from kmapper.plotlyviz import plotlyviz
import math
import networkx as nx
import pygraphviz
import pdfkit
import re
import os
import subprocess
import shlex
import gudhi as gd
#from sklearn_tda import *
import sklearn_tda
import statmapper as stm
from sklearn.manifold import MDS
import mapper # Daniel Müller, for lenses.
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA

# Load summarized OTU table.
df1 = pd.read_csv("../../RNA_genus_species.csv", index_col="Genus")
df1.columns = df1.columns.str.replace('r.PS79.', '')
df1.columns = df1.columns.str.replace('m_p.fastq', '')
df1.columns = df1.columns.str.replace('.', '_')

names = [c for c in df1.columns]
runtype = "SVD_k0_RNA_10_04_1lens"
X = np.array(df1[names].fillna(0)) # doesen't work with NaNs

print(names)
print(df1)

## Load labels for later coloring. Matching key should be the station!
df2 = pd.read_csv("../../nonscale_imputed_meinhard_full_province_alphadiv_mld.txt")
feature_names = [c for c in df2.columns if c not in ["ID","Station","Province","Cruise","DNA_div", "DNA_ENS", "RNA_div", "RNA_ENS"]] #, "Lon",$
print(feature_names)

##df["diagnosis"] = df["diagnosis"].apply(lambda x: 1 if x == "M" else 0) ## good for binary categories
#Y0 = np.array(df2[feature_names].fillna(0)) # doesen't work with NaNs

#Turn categorical data to numbers.
plot_labels = np.array(df2["Province"]) #z
x = df2.iloc[:, 1:].values
labelencoder_x = LabelEncoder()
x[:, 2] = labelencoder_x.fit_transform(x[:, 2])
newval = list(x[:, 2])
df2.iloc[:, 3] =newval
filtered = df2[df2.ID.isin(df1.columns)]

# Only for coloring. Here NaNs will be zeros, and that will mark the missing values.
#labels = ['Province', 'Depth [m]'] #, 'Bot. Depth [m]', 'DNA_div','DNA_ENS', 'RNA_div', 'RNA_ENS']
labels = [c for c in df2.columns if c not in ["ID","Station","Cruise"]] #, 'Depth_m']] #, 'Bot_Depth..m']]
Y = np.array(filtered[labels].fillna(0))

# run Keplermapper with filters -  question is optimisation of parameters. Should be run with GUDHI beforehand?
kmmapper = km.KeplerMapper(verbose=2)

# choose 1 lens here

# Custom 1-D lens with Isolation Forest - detecting anomalies
##model = ensemble.IsolationForest(random_state=1729)
##model.fit(X)
##lens = model.decision_function(X).reshape((X.shape[0], 1))

# 1-D lens with Isomap - to disperse and reduce dimensions.
##embedding = Isomap(n_components=1) # not sure of n_components are good
##lens = embedding.fit_transform(X)

#Normalization, get some distance with l2norm
#lens = kmmapper.fit_transform(X, projection="l2norm")

# If I want to use an environmental parameter as a lens
##lens = mapper.fit_transform(X, projection=[5]) #choose from the environmental table

# Combine all lenses: (not for this setup)
##lens = np.c_[lens1,lens2]

#Multidimensional scaling
##embedding = MDS(n_components=1, random_state=1)
##lens = embedding.fit_transform(X)

# Row mean
##lens = kmmapper.fit_transform(X, projection="mean")

# Row median
##lens = kmmapper.fit_transform(X, projection="median")

# eccenticity
##lens=mapper.filters.eccentricity(X)
##lens = np.reshape(lens, (9370, 1))

# Distance to measure
##lens=mapper.filters.distance_to_measure(X,2)
##lens = np.reshape(lens, (9370, 1))

# L-inf centrality with Bray-Curtis
##pairwise_dist = squareform(pdist(X, "braycurtis"))
##lens=np.amax(pairwise_dist, axis = 1)
##lens=np.nan_to_num(lens)
##lens = np.reshape(lens, (9370, 1))

#PCA
##lens = kmmapper.fit_transform(X, projection=PCA(n_components=1, random_state=1))

# Laplacian - doesen't work, mode has weird format input
##lens=mapper.filters.graph_Laplacian(X, eps=30000, k=2)
##lens = np.reshape(lens, (9370, 1))

#Eigenvectors
lens=mapper.filters.dm_eigenvector(X,mean_center=True,k=0,metricpar={'metric': 'braycurtis'})
lens = np.reshape(lens, (9370, 1))

print(lens)

# Transforming mapper simplicial complex to simplex tree. Gain has to be set beforehand, between 0.1 and 0.4.
# Resolution however can be estimated, but it will give simplier graphs.


params = {"filters": lens, "filter_bnds": np.array([[np.nan,np.nan]]), "colors": X[:,2:3],\
 "resolutions": np.array([10]), "gains": np.array([0.4]), "inp": "point cloud",\
 "clustering": sklearn.cluster.KMeans(n_clusters=2,random_state=1618033)}

mymapper = sklearn_tda.clustering.MapperComplex(**params).fit(X)

mapper_simplex_tree = mymapper.mapper_

print("Number of simplices and vertices:")
print(mymapper.mapper_.num_simplices())
print(mymapper.mapper_.num_vertices())
print(mymapper.resolutions)

diag = mapper_simplex_tree.persistence(min_persistence=-1,persistence_dim_max=3) # show all elements to Betti number 2.
print(diag)

barcode_pdf = (('Meinhard_%s_barcode.pdf') % runtype)
#plt.figure()
gd.plot_persistence_barcode(diag)
plt.savefig(barcode_pdf)
plt.clf()

persdiag_pdf = (('Meinhard_%s_persistencediag.pdf') % runtype)
#plt.figure()
gd.plot_persistence_diagram(persistence=diag)
plt.savefig(persdiag_pdf)
plt.clf()

# here statmapper
statnum_pdf = (('Meinhard_%s_node_nums.pdf') % (runtype))
G = stm.mapper2networkx(mymapper)
nx.draw(G,with_labels=True)
#plt.show()
plt.savefig(statnum_pdf)
plt.clf()

function = lens # filter or sum? or color?
topo = ["connected_components", "loop", "downbranch", "upbranch"]
confidence = 0.95 # check if enough
bootstrap  = 100

for t in topo:
    dgm, bnd = stm.compute_topological_features(mymapper, function, "data", t)
    print(dgm, bnd)
    topo_pdf = (('Meinhard_%s_%s_topo_elements.pdf') % (runtype,t))
    plt.figure(figsize=(40,10))
    for idx, bd in enumerate(bnd):
        plt.subplot(1,len(bnd),idx+1)
        nx.draw(G, pos=nx.nx_agraph.graphviz_layout(G, prog='neato'), with_labels=True, node_color=[1 if node in bd else 0 for node in G.nodes()])
        plt.title(t)

    plt.savefig(topo_pdf)
    plt.clf()

    sdgm, sbnd = stm.evaluate_significance(dgm, bnd, X, mymapper, function, params, t, confidence, bootstrap)

    #features, pv = stm.compute_DE_features(X, mapper, sbnd[0]) this part doesn't work
    #print(features, pv)
    
    stat_pdf = (('Meinhard_%s_%s_topo_stat.pdf') % (runtype,t))
    plt.figure(figsize=(40,10))
    for idx, bd in enumerate(sbnd):
        plt.subplot(1,len(sbnd),idx+1)
        nx.draw(G, pos=nx.nx_agraph.graphviz_layout(G, prog='neato'), with_labels=True, node_color=[1 if node in bd else 0 for node in G.nodes()])
    plt.savefig(stat_pdf)
    plt.clf()

plt.close('all')

## stamapper and kepler-mapper to networkx solution - should be the same input! Later: try params?
keplermapper_mapper = kmmapper.map(lens,X,cover=km.Cover(n_cubes=10, perc_overlap=0.4),clusterer=sklearn.cluster.KMeans(n_clusters=2,random_state=1618033))
kepler_name = (('Meinhard_%s_keplermapper.html') % runtype)
keplerhtml = kmmapper.visualize(keplermapper_mapper, path_html=kepler_name)

kepler_pdf = (('Meinhard_%s_keplermapper.pdf') % runtype)
#plt.figure()
plt.title("kepler-mapper output")
keplergraph = km.adapter.to_nx(keplermapper_mapper)
keplerpos = nx.nx_agraph.graphviz_layout(keplergraph, prog='neato')
nodedict = keplermapper_mapper.get('nodes')

nodesize= []
color_map=[]

for k, v in nodedict.items():
    nodelist=v
    nodecolor=np.mean(nodelist)
    color_map.append(nodecolor)
    nodesize.append(len(nodelist))
print(nodesize)
node_normalized=(color_map-min(color_map))/float((max(color_map)-min(color_map)))
vmin = min(node_normalized)
vmax = max(node_normalized)
nx.draw_networkx(keplergraph, keplerpos, with_labels=False, node_size=nodesize, node_color=node_normalized, vmin=vmin, vmax=vmax, cmap=plt.cm.get_cmap('viridis'), width=0.25, edge_color="grey")
plt.savefig(kepler_pdf)
plt.clf()

# this is for the later loop
nx_graph = stm.mapper2networkx(mymapper)

#pos = nx.spring_layout(nx_graph, dim=2, k=None, pos=None, fixed=None, iterations=100, weight='weight', scale=1.0)
pos = nx.nx_agraph.graphviz_layout(nx_graph, prog='neato')

# betweenes centrality
cent_pdf = (('Meinhard_%s_betweeness.pdf') % runtype)
plt.figure()
plt.title("Betweeness")
#nx.draw_networkx(nx_graph, pos, with_labels=False, node_size=5,width=0.25, edge_color="grey")
degree_centrality = nx.degree_centrality(nx_graph)
#sorted(dict(nx_graph.degree()).items(), key=lambda x: x[1], reverse=True)[:10]

nx.set_node_attributes(nx_graph, name='degree centrality', values=degree_centrality)
betweeness_centrality = nx.betweenness_centrality(nx_graph)
vmin = min(betweeness_centrality.values())
vmax = max(betweeness_centrality.values())
nx.draw_networkx(nx_graph, pos, with_labels=False, node_size=5, node_color=list(betweeness_centrality.values()),nodelist=list(betweeness_centrality.keys()), vmin=vmin, vmax=vmax, cmap=plt.cm.get_cmap('viridis'), width=0.25, edge_color="grey")
nodes1 = nx.draw_networkx_nodes(nx_graph, pos, with_labels=False, node_size=5, node_color=list(betweeness_centrality.values()),nodelist=list(betweeness_centrality.keys()), vmin=vmin, vmax=vmax, cmap=plt.cm.get_cmap('viridis'), width=0.25, edge_color="grey")

plt.colorbar(nodes1)
plt.axis('off')
plt.savefig(cent_pdf)
plt.clf()

eig_pdf = (('Meinhard_%s_eigvec.pdf') % runtype)
plt.figure()
plt.title("Eigenvector centrality")
#eigenvector_centrality = nx.eigenvector_centrality(nx_graph) # or
eigenvector_centrality = nx.eigenvector_centrality_numpy(nx_graph) 
max_value = max(eigenvector_centrality.items(), key=lambda x: x[1])
ec_scaled = {}
for k in eigenvector_centrality.keys():
    ec_scaled[k] = eigenvector_centrality[k] / max_value[1]

vmin = min(ec_scaled.values())
vmax = max(ec_scaled.values())
nx.draw_networkx(nx_graph, pos, with_labels=False, node_size=5, node_color=list(ec_scaled.values()), vmin=vmin, vmax=vmax, cmap=plt.cm.get_cmap('viridis'), width=0.25, edge_color="grey")
nodes2 = nx.draw_networkx_nodes(nx_graph, pos, with_labels=False, node_size=5, node_color=list(ec_scaled.values()),nodelist=list(ec_scaled.keys()), vmin=vmin, vmax=vmax, cmap=plt.cm.get_cmap('viridis'), width=0.25, edge_color="grey")

plt.colorbar(nodes2)
plt.axis('off')
plt.savefig(eig_pdf)
plt.clf()

print(labels)
num = 0

for i in labels:
    print(i)
    extractkey = i
    num = num +1

    i = re.sub('[^A-Za-z0-9]+', '', i)
    out_pdf = ('Meinhard_tmp_%s_%s.pdf' % (runtype, i))

    # extracting values (mean) for colors

    if i == "Province":

        # pair provinces with correspontding value, define names for later barplot label
        pairs = sorted(set(zip(plot_labels, newval)))
        names =[]
        for m in range(0,len(pairs)):
            names.append(pairs[m][0])
        print(pairs)
        # set node colors and node size:
        cmap = mpl.cm.get_cmap('viridis', 8)
        colors=[]
        for l in range(cmap.N):
            rgb = cmap(l)[:3] # will return rgba, we take only first 3 so we get rg
            colors.append(mpl.colors.rgb2hex(rgb))

        color_map = []
        nodesize = []
        colordict={}
        newv = []

        for node in nx_graph:
            nodelist=mymapper.node_info_[node]['indices']
            ids_tmp = df1.iloc[nodelist,].sum(0) > 1
            ids=ids_tmp[ids_tmp].index.tolist() # full list of stations
            nodesize.append(len(ids))
            #nodecolor=np.mean(df2.loc[df2.ID.isin(ids),[extractkey]][extractkey].tolist()) # average value for a node
            nodecolor = list(scipy.stats.mode(df2.loc[df2.ID.isin(ids),[extractkey]][extractkey].tolist()))[0]
            color_map.append(nodecolor)
            nom_elements= df2.loc[df2.ID.isin(ids),[extractkey]][extractkey].tolist()
            print(nom_elements)
            print(scipy.stats.mode(nom_elements))

            newv= []
            for e in nom_elements:
                #newv.append(colors[e-1])
                newv.append(colors[e])
            colordict.update( {node :newv} )

        color_map = np.concatenate(color_map).ravel().tolist()
        
        vmin=min(color_map)
        vmax=max(color_map)
        modelist=[]
        print(color_map)
        for y in color_map:
            modelist.append(colors[y])

        plt.figure()
        plt.title(i)

        # First nodes, because otherwise the colorbar is not discrete.
        nodes = nx.draw_networkx_nodes(nx_graph, pos, with_labels=False, node_size=nodesize, node_color=color_map, vmin=vmin, vmax=vmax,cmap=cmap,width=0.25, edge_color="grey")
        nx.draw_networkx(nx_graph, pos, with_labels=False, node_size=nodesize, node_color=modelist, width=0.25, edge_color="grey")
        cb=plt.colorbar(nodes, cmap=cmap)
        tick_locs=np.linspace(cb.vmin+((cb.vmax-cb.vmin)/16),cb.vmax-((cb.vmax-cb.vmin)/16),8)
        cb.set_ticks(tick_locs)
        cb.set_ticklabels(names)
        cb.ax.set_ylabel('Modal value of the nodes')
        plt.axis('off')
        plt.savefig(out_pdf)
        plt.clf()

        # with pie charts:
        plt.figure()
        plt.title(i)
        nx.draw_networkx_edges(nx_graph, pos, with_labels=False, node_size=5, node_color=color_map, vmin=vmin, vmax=vmax,cmap=cmap,width=0.25, edge_color="grey")

        for node in nx_graph.nodes:
            #a = plt.pie(extractdict[node],center=pos[node],radius=25,colors=[cmap(a) for a in extractdict[node]])
            a = plt.pie(ids,center=pos[node],radius=25, colors=colordict[node])

        plt.axis('scaled')
        nodes = nx.draw_networkx_nodes(nx_graph, pos, with_labels=False, node_size=0.005, node_color=color_map, vmin=vmin, vmax=vmax,cmap=cmap, width=0.25, edge_color="grey")
        cb=plt.colorbar(nodes, cmap=cmap)
        tick_locs = np.linspace(cb.vmin+((cb.vmax-cb.vmin)/16),cb.vmax-((cb.vmax-cb.vmin)/16),8)
        cb.set_ticks(tick_locs)
        cb.set_ticklabels(names)
        plt.axis('off')
        plt.savefig(("Meinhard_%s_piechart.pdf") % runtype)

    else:

        color_map = []
        nodesize = []
        newv = []

        for node in nx_graph:
            nodelist=mymapper.node_info_[node]['indices']
            ids_tmp = df1.iloc[nodelist,].sum(0) > 1
            ids=ids_tmp[ids_tmp].index.tolist() # full list of stations
            nodesize.append(len(ids))

            # Non normalized values
            nodecolor=np.mean(df2.loc[df2.ID.isin(ids),[extractkey]][extractkey].tolist()) # average value for a node
            color_map.append(nodecolor)
            #num_elements= df2.loc[df2.ID.isin(ids),[extractkey]][extractkey].tolist()

        vmin = min(color_map)
        vmax = max(color_map)

        # If it should be normalized:
        numlist_normalized=[x-min(color_map)/(max(color_map)-float(min(color_map))) for x in color_map]
        numlist_normalized=[0 if math.isnan(x) else x for x in numlist_normalized]

        vmin = min(numlist_normalized)
        vmax = max(numlist_normalized)
        numlist=[0 if math.isnan(x) else x for x in color_map]
        nums= np.around(np.linspace(min(color_map), max(color_map), 6),decimals=3)
        plt.figure()
        plt.title(i)
        print(color_map)
        #color_map = numlist_normalized or color_map
        nx.draw_networkx(nx_graph, pos, with_labels=False, node_size=nodesize, node_color=numlist_normalized, vmin=vmin, vmax=vmax, cmap=plt.cm.get_cmap('viridis'), width=0.25, edge_color="grey")
        nodes =nx.draw_networkx_nodes(nx_graph, pos, with_labels=False, node_size=nodesize, node_color=color_map, vmin=vmin, vmax=vmax, cmap=plt.cm.get_cmap('viridis'), width=0.25, edge_color="grey")
        #sm = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap('viridis'), norm=plt.Normalize(vmin = vmin, vmax=vmax))
        cb=plt.colorbar(nodes)
        tick_locs = np.linspace(cb.vmin,cb.vmax,6)
        cb.set_ticks(tick_locs)
        cb.set_ticklabels(nums)
        cb.ax.set_ylabel('Average value of the nodes')
        plt.axis('off')
        plt.savefig(out_pdf)
        #nx.draw_networkx_labels(nx_graph, pos=pos, font_size=5)
        #plt.show()
        #nx.write_gml(nx_graph, 'ga_graph.gexf')

        # save values of nodes into a table
        #nodetable = df2.loc[df2.ID.isin(ids)]
        #nodetable.insert(0,"Node",node, True)
        #endtable = pd.concat([endtable, nodetable], ignore_index=True)

nodetable = None
endtable = None

for node in nx_graph:
    nodelist=mymapper.node_info_[node]['indices']
    ids_tmp = df1.iloc[nodelist,].sum(0) > 1
    ids=ids_tmp[ids_tmp].index.tolist()
    nodetable = df2.loc[df2.ID.isin(ids)]
    nodetable.insert(0,"Node",node, True)
    endtable = pd.concat([endtable, nodetable], ignore_index=True)

    spectable = df1.iloc[nodelist]
    spectable.to_csv(runtype + "_node_" + str(node) + "_spectable.txt", index=True, header= True)

endtable.to_csv(runtype + "_node_table.txt", index=False, header= True)

command = "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=Meinhard_" + runtype + ".pdf Meinhard_tmp_" + runtype + "_*.pdf"
os.system(command)

command = "rm  Meinhard_tmp_" + runtype + "_*.pdf"
os.system(command)

print(labels)
print(mymapper.resolutions)
