## Performing Topological Data Analysis, statistical evaulation of topological elements, visualization and data extraxtion. ##

# Details:
# title           :keplermapper_meinhard_networkx_biodiv_stations_both.py  
# author          :Emese Xochitl Szabo
# email:	  :emese.szabo@uni-oldenburg.de
# date            :25/02/2021
# version         :0.1
# license         :GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# usage           :python keplermapper_meinhard_networkx_biodiv_stations_both.py
# notes           :GUDHI environment has to be activated first
# python_version  :3.8.2

#!/usr/bin/python

## Import basic modules for data loading, preprocessing visualisation. ##
import numpy as np
import scipy
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import Counter
import math
import pdfkit
import re
import os
import subprocess
import shlex

## Import modules for mandatory TDA construction. ##
import kmapper as km
from kmapper.plotlyviz import plotlyviz
import networkx as nx
import pygraphviz

## Import optional modules, if the lens needs it. ##
import sklearn
from sklearn import ensemble
from sklearn.datasets import load_digits
from sklearn.manifold import Isomap
from sklearn.preprocessing import LabelEncoder
from sklearn.decomposition import PCA
from sklearn import datasets
from sklearn.manifold import MDS
import mapper # Daniel Müller, for lenses.
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA

## Import optional modules for TDA statistical evaulation and persistent homology. ##
import gudhi as gd
import sklearn_tda
import statmapper as stm

## Give the name of the run. ## 
runtype = "Linf_BC_10_04_new_both_species_scaled_filtered_1lens_stations"

# Load OTU table.
df1 = pd.read_csv("../../Full_scaled_filtered_both_Count.tsv", sep="\t")


# Modifications for Genus runs:
##df1 = df1.drop(df1.columns[0:6], axis=1)
##df1 = df1.drop(df1.columns[1], axis=1)
##df1 = df1.set_index('Genus')

# Modifications for Species run:
df1 = df1.drop(df1.columns[0:7], axis=1)
df1 = df1.set_index('Species')
print(df1)
#df1.columns = df1.columns.str.replace('d.PS79.', '')
#df1.columns = df1.columns.str.replace('m_p.fastq', '')
#df1.columns = df1.columns.str.replace('.', '_')
df1 = df1.select_dtypes([np.number])
df1 = df1.T
names = [c for c in df1.columns]
X = np.array(df1) #(df1[names]) # .fillna(0)) # safety step, because TDA doesen't work with NaNs

## Load labels for later coloring. Matching key should be the station! ##
df_ori = pd.read_csv("../../both_merged_scale_imputed_meinhard_ocean.csv", sep=",") # two tables were separately imputed.
feature_names = [c for c in df_ori.columns if c not in ["Cruise","Station","Date","Time"]] # don't forget to deselect unnecessary values
labeldf = df_ori.copy()
labeldf = labeldf[labeldf.Sample_ID.isin(df1.index)]
labeldf = labeldf[feature_names]
print("labels")
print(labeldf)

## Optional step for binary categories ##
#df["diagnosis"] = df["diagnosis"].apply(lambda x: 1 if x == "M" else 0)

##Turn categorical data to numbers. ##

df2 = df_ori.copy() # labels have to stay the same, copy provides separate table to manipulate.
for o in  0, 3: # choose columns manually with categorical data
    x = df2.iloc[:, 0:].values
    labelencoder_x = LabelEncoder()
    print(x)
    x[:, o] = labelencoder_x.fit_transform(x[:, o])
    newval = list(x[:, o])
    df2.iloc[:, o] =newval
filtered = df2[df2.Sample_ID.isin(df1.index)]
filtered = filtered[feature_names]
Y=np.array(filtered)
labels= feature_names

## initialize Keplermapper ##
kmmapper = km.KeplerMapper(verbose=2)

## choose 1 lens here: it is recommened to test out all of them and see which lens separates the patterns best. Advice: if it does not work, check the array shape! ##

## Custom 1-D lens with Isolation Forest - detecting anomalies ##
#model = ensemble.IsolationForest(random_state=1729)
#model.fit(X)
#lens = model.decision_function(X).reshape((X.shape[0], 1))

## 1-D lens with Isomap - to disperse and reduce dimensions. ##
#embedding = Isomap(n_components=1) # not sure of n_components are good 
#lens = embedding.fit_transform(X)

## Normalization, get some distance with l2norm ##
#lens = kmmapper.fit_transform(X, projection="l2norm")

## If I want to use an environmental parameter as a lens (make sure to choose the right column)##
#Z = filtered.fillna(0)
#Z = np.array(Z)
#lens = kmmapper.fit_transform(Z, projection=[4]) # It was Latitude for the subproject

## Combine all lenses: (not for this setup, if you want to do, choose two lenses and name them lens1 and lens2)
#lens = np.c_[lens1,lens2]

## Multidimensional scaling ##
#embedding = MDS(n_components=1, random_state=1)
#lens = embedding.fit_transform(X)

## Row mean ##
#lens = kmmapper.fit_transform(X, projection="mean")

## Row median ##
#lens = kmmapper.fit_transform(X, projection="median")

## Eccenticity ##
#lens=mapper.filters.eccentricity(X)
#lens = np.reshape(lens, (len(lens), 1))

## Distance to measure ##
#lens=mapper.filters.distance_to_measure(X,2)
#lens = np.reshape(lens, (len(lens), 1))

## L-inf centrality with Bray-Curtis ##
#pairwise_dist = squareform(pdist(X, "braycurtis"))
#lens=np.amax(pairwise_dist, axis = 1)
#lens=np.nan_to_num(lens)
#lens = np.reshape(lens, (len(lens), 1))

## PCA ##
#lens = kmmapper.fit_transform(X, projection=PCA(n_components=1, random_state=1))

## Laplacian ##
#lens=mapper.filters.graph_Laplacian(X, eps=30000, k=2)

## Eigenvectors - SVD ##
#lens=mapper.filters.dm_eigenvector(X,mean_center=True,k=0,metricpar={'metric': 'braycurtis'})
lens = np.reshape(lens, (len(lens), 1))

## Transforming mapper simplicial complex to simplex tree. Gain has to be set beforehand, between 0.1 and 0.4. Resolution however can be estimated, but it will give simplier graphs. ##

params = {"filters": lens, "filter_bnds": np.array([[np.nan,np.nan]]), "colors": X[:,2:3],\
 "resolutions": np.array([10]), "gains": np.array([0.4]), "inp": "point cloud",\
 "clustering": sklearn.cluster.KMeans(n_clusters=2,random_state=1618033)}

mymapper = sklearn_tda.clustering.MapperComplex(**params).fit(X)
mapper_simplex_tree = mymapper.mapper_

## Check structure of the simplicial complex ##
print("Number of simplices and vertices:")
print(mymapper.mapper_.num_simplices())
print(mymapper.mapper_.num_vertices())
print(mymapper.resolutions)

## Create persistence barcode and diagrams: visualize the levels of topological elements. ##
diag = mapper_simplex_tree.persistence(min_persistence=-1,persistence_dim_max=3) # show all elements to Betti number 2.
barcode_pdf = (('Meinhard_%s_barcode.pdf') % runtype)
gd.plot_persistence_barcode(diag)
plt.savefig(barcode_pdf)
plt.clf()

persdiag_pdf = (('Meinhard_%s_persistencediag.pdf') % runtype)
gd.plot_persistence_diagram(persistence=diag)
plt.savefig(persdiag_pdf)
plt.clf()

## Transform the simplicial complex to a NetworkX graph for reprodicuble representation. ##
nx_graph = stm.mapper2networkx(mymapper)
pos = nx.nx_agraph.graphviz_layout(nx_graph, prog='neato')
statnum_pdf = (('Meinhard_%s_node_nums.pdf') % (runtype))
nx.draw(nx_graph,with_labels=True)
plt.savefig(statnum_pdf)
plt.clf()

## stamapper and kepler-mapper to networkx solution - should be the same input! ##
keplermapper_mapper = kmmapper.map(lens,X,cover=km.Cover(n_cubes=10, perc_overlap=0.4),clusterer=sklearn.cluster.KMeans(n_clusters=2,random_state=1618033))
kepler_name = (('Meinhard_%s_keplermapper.html') % runtype)
keplerhtml = kmmapper.visualize(keplermapper_mapper, path_html=kepler_name)

kepler_pdf = (('Meinhard_%s_keplermapper.pdf') % runtype)
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

node_normalized=(color_map-min(color_map))/float((max(color_map)-min(color_map)))
vmin = min(node_normalized)
vmax = max(node_normalized)
nx.draw_networkx(keplergraph, keplerpos, with_labels=False, node_size=nodesize, node_color=node_normalized, vmin=vmin, vmax=vmax, cmap=plt.cm$
plt.savefig(kepler_pdf)
plt.clf()

## Statistical evaulation  of the simplicial complex. Beware, when it tries to evaulate the down/upbranch it analyzes the values IN the graph. ##
function = lens 
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
        nx.draw(nx_graph, pos=nx.nx_agraph.graphviz_layout(G, prog='neato'), with_labels=True, node_color=[1 if node in bd else 0 for node in nx_graph.nodes()])
        plt.title(t)

    plt.savefig(topo_pdf)
    plt.clf()

    sdgm, sbnd = stm.evaluate_significance(dgm, bnd, X, mymapper, function, params, t, confidence, bootstrap)
    
    stat_pdf = (('Meinhard_%s_%s_topo_stat.pdf') % (runtype,t))
    plt.figure(figsize=(40,10))
    for idx, bd in enumerate(sbnd):
        plt.subplot(1,len(sbnd),idx+1)
        nx.draw(nx_graph, pos=nx.nx_agraph.graphviz_layout(G, prog='neato'), with_labels=True, node_color=[1 if node in bd else 0 for node in nx_graph.nodes()])
    plt.savefig(stat_pdf)
    plt.clf()

plt.close('all')

## Add Node numbers and names ##
color_map =[]
nodelabel = {}
stationlabel = {}
newv= None

for node in nx_graph:
    #print(mymapper.node_info_[node])
    nodenum=np.mean(mymapper.node_info_[node]['indices'])
    stations = mymapper.node_info_[node]['indices']
    stations_names = df1.index[stations] #.Sample_ID.tolist()
    pairs = sorted(set(zip(stations, stations_names)))
    names =[]
    for m in range(0,len(pairs)):
        names.append(pairs[m][1])

    color_map.append(nodenum)
    newv = mymapper.node_info_[node]['indices'].tolist()
    if len(newv) > 7:
        reps = int(len(newv)/4)
        new = [newv[i: i+reps] for i in range(0, len(newv), reps)]
        new_names = [names[i: i+reps] for i in range(0, len(names), reps)]
        updated = []
        for j in new:
            updated.append(str(j))
            updated.append('\n')
        upd = ''.join(updated)
        nodelabel.update( {node : upd} )

        names_updated = []
        for k in new_names:
            names_updated.append(str(k))
            names_updated.append('\n')
        n_upd = ''.join(names_updated)
        stationlabel.update( {node : n_upd} )

    else:
        nodelabel.update( {node : newv} )
        stationlabel.update( {node : names } )

node_pdf = (('Meinhard_%s_nodes.pdf') % runtype)
nx.draw_networkx(nx_graph,pos,labels= nodelabel, node_color="grey", vmin=vmin, vmax=vmax, cmap=plt.cm.get_cmap('viridis'), width=0.25, edge_color="grey", font_size=6)
plt.axis('off')
plt.tight_layout()
plt.savefig(node_pdf)
plt.clf()

station_node_pdf = (('Meinhard_%s_nodes_stations.pdf') % runtype)
nx.draw_networkx(nx_graph,pos,labels= stationlabel, node_color="grey", vmin=vmin, vmax=vmax, cmap=plt.cm.get_cmap('viridis'), width=0.25, edge_color="grey", font_size=3)
plt.axis('off')
plt.tight_layout()
plt.savefig(station_node_pdf)
plt.clf()

## Betweenes centrality - ##
cent_pdf = (('Meinhard_%s_betweeness.pdf') % runtype)
plt.figure()
plt.title("Betweeness")
degree_centrality = nx.degree_centrality(nx_graph)

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

## Eigenvector centrality - ##
eig_pdf = (('Meinhard_%s_eigvec.pdf') % runtype)
plt.figure()
plt.title("Eigenvector centrality")
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

num = 0

for i in labels:
    extractkey = i
    num = num +1
    ##i = re.sub('[^A-Za-z0-9]+', '', i)
    out_pdf = ('Meinhard_tmp_%s_%s.pdf' % (runtype, i))
    out_png = ('Meinhard_tmp_%s_%s.png' % (runtype, i))
    print(i)
    # extracting values (mean) for colors

    if i == "Sample_ID": pass

    elif i == "Province" or i == "Ocean":
        plot_labels = labeldf[i]
        numvalcat = list(filtered[i])
        print(filtered[i])
        # pair provinces with correspontding value, define names for later barplot label
        pairs = sorted(set(zip(plot_labels, numvalcat)))

        names =[]
        for m in range(0,len(pairs)):
            names.append(pairs[m][0])

        # set node colors and node size:
        numcol = len(pairs)
        cmap = mpl.cm.get_cmap('viridis', numcol)
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
            ids = df1.iloc[nodelist,].index.tolist() # full list of stations
            nodesize.append(len(ids))
            nodesize_norm =  [math.sqrt(x)*25 for x in nodesize]
            #nodecolor=np.mean(df2.loc[df2.ID.isin(ids),[extractkey]][extractkey].tolist()) # average value for a node

            nodecolor = list(scipy.stats.mode(filtered.loc[filtered.Sample_ID.isin(ids),[extractkey]][extractkey].tolist()))[0]
            color_map.append(nodecolor)
            nom_elements= filtered.loc[filtered.Sample_ID.isin(ids),[extractkey]][extractkey].tolist()
            ##print(scipy.stats.mode(nom_elements))

            newv= []
            for e in nom_elements:
                newv.append(colors[e])
            newv.sort()
            colordict.update( {node :newv} )

        # mandatory modification of modelist, because it is an array.
        color_map = np.concatenate(color_map).ravel().tolist()
        
        vmin=min(color_map)
        vmax=max(color_map)
        modelist=[]

        for y in color_map:
            modelist.append(colors[y])

        plt.figure()
        plt.title(i)

        # First nodes, because otherwise the colorbar is not discrete.
        nodes = nx.draw_networkx_nodes(nx_graph, pos, with_labels=False, node_size=nodesize_norm, node_color=color_map, vmin=vmin, vmax=vmax,cmap=cmap,width=0.25, edge_color="grey")
        nx.draw_networkx(nx_graph, pos, with_labels=False, node_size=nodesize_norm, node_color=modelist, width=0.25, edge_color="grey")
        cb=plt.colorbar(nodes, cmap=cmap)
        tick_locs=np.linspace(cb.vmin+((cb.vmax-cb.vmin)/(2*numcol)),cb.vmax-((cb.vmax-cb.vmin)/(2*numcol)),numcol)
        cb.set_ticks(tick_locs)
        cb.set_ticklabels(names)
        cb.ax.set_ylabel('Modal value of the nodes')
        plt.axis('off')
        plt.savefig(out_pdf)
        plt.savefig(out_png)
        plt.clf()

        # with pie charts:
        plt.figure()
        plt.title(i)
        nx.draw_networkx_edges(nx_graph, pos, with_labels=False, node_size=5, node_color=color_map, vmin=vmin, vmax=vmax,cmap=cmap,width=0.25, edge_color="grey")
        
        for node in nx_graph.nodes:
            #a = plt.pie(extractdict[node],center=pos[node],radius=25,colors=[cmap(a) for a in extractdict[node]])
            c = Counter(colordict[node])
            print(c)
            ratios = [(c[y] / len(colordict[node]) * 100.0) for y in c]
            #a = plt.pie(ratios,center=pos[node],radius=25, colors=colordict[node], wedgeprops = {'linewidth': 0}) # sorting ids?
            a = plt.pie(ratios,center=pos[node],radius=25, colors=c.keys(), wedgeprops = {'linewidth': 0}) # sorting ids? -> no, because the Counter does it

        plt.axis('scaled')
        nodes = nx.draw_networkx_nodes(nx_graph, pos, with_labels=True, node_size=0.005, node_color=color_map, vmin=vmin, vmax=vmax,cmap=cmap, width=0.25, edge_color="grey")
        cb=plt.colorbar(nodes, cmap=cmap)
        tick_locs = np.linspace(cb.vmin+((cb.vmax-cb.vmin)/(2*numcol)),cb.vmax-((cb.vmax-cb.vmin)/(2*numcol)),numcol)
        cb.set_ticks(tick_locs)
        cb.set_ticklabels(names)
        plt.axis('off')
        plt.savefig(("Meinhard_%s_%s_piechart.pdf") % (runtype, i))

    else:

        color_map = []
        nodesize = []
        #newv = []

        for node in nx_graph:
            nodelist=mymapper.node_info_[node]['indices']
            ids = df1.iloc[nodelist,].index.tolist() # full list of stations
            nodesize.append(len(ids))
            nodesize_norm =  [math.sqrt(x)*25 for x in nodesize]
            # Non normalized values
            nodecolor=np.mean(filtered.loc[filtered.Sample_ID.isin(ids),[extractkey]][extractkey].tolist()) # average value for a node
            color_map.append(nodecolor)
            #num_elements= df2.loc[df2.ID.isin(ids),[extractkey]][extractkey].tolist()

        vmin = min(color_map)
        vmax = max(color_map)

        # If it should be normalized:
        ##numlist_normalized=[x-min(color_map)/(max(color_map)-float(min(color_map))) for x in color_map]
        ##numlist_normalized=[0 if math.isnan(x) else x for x in numlist_normalized]
        ##vmin = min(numlist_normalized)
        ##vmax = max(numlist_normalized)

        numlist=[0 if math.isnan(x) else x for x in color_map]
        nums= np.around(np.linspace(min(color_map), max(color_map), 6),decimals=3)
        nums=[np.format_float_scientific(x,precision=2) for x in nums] # pretty format
        plt.figure()
        plt.title(i)
        print(color_map)
        #color_map = numlist_normalized or color_map
        nx.draw_networkx(nx_graph, pos, with_labels=False, node_size=nodesize_norm, node_color=color_map, vmin=vmin, vmax=vmax, cmap=plt.cm.get_cmap('viridis'), width=0.25, edge_color="grey")
        nodes =nx.draw_networkx_nodes(nx_graph, pos, with_labels=False, node_size=nodesize_norm, node_color=color_map, vmin=vmin, vmax=vmax, cmap=plt.cm.get_cmap('viridis'), width=0.25, edge_color="grey")
        #sm = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap('viridis'), norm=plt.Normalize(vmin = vmin, vmax=vmax))
        cb=plt.colorbar(nodes,format='%.0e')
        tick_locs = np.linspace(cb.vmin,cb.vmax,6)
        cb.set_ticks(tick_locs)
        cb.set_ticklabels(nums)
        cb.ax.set_ylabel('Average value of the nodes')
        plt.axis('off')
        plt.savefig(out_pdf)
        plt.savefig(out_png)

        #nx.draw_networkx_labels(nx_graph, pos=pos, font_size=5)
        #plt.show()
        #nx.write_gml(nx_graph, 'ga_graph.gexf')

nodetable = None
endtable = None

for node2 in nx_graph:
    nodelist2=mymapper.node_info_[node2]['indices']
    stations_names = df1.index[nodelist2] #.Sample_ID.tolist()
    nodetable = labeldf.loc[labeldf.Sample_ID.isin(stations_names)]
    nodetable.insert(0,"Node",node2, True)
    endtable = pd.concat([endtable, nodetable], ignore_index=True)

    spectable = df1.iloc[nodelist2]
    spectable.to_csv(runtype + "_node_" + str(node2) + "_spectable.csv", index=True, header= True)

endtable.to_csv(runtype + "_node_table.csv", index=False, header= True)

command = "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=Meinhard_" + runtype + ".pdf Meinhard_tmp_" + runtype + "_*.pdf"
os.system(command)

command = "rm  Meinhard_tmp_" + runtype + "_*.pdf"
os.system(command)

command = "ffmpeg -pattern_type glob  -i 'Meinhard_tmp_" + runtype + "_*.png' -vf palettegen " + runtype + "_palette.png"
os.system(command)

#command = "ffmpeg -pattern_type glob -framerate 1 -i 'Meinhard_tmp_" + runtype + "_*.png' -i " + runtype + "_palette.png -lavfi paletteuse -loop 1 " + runtype + "_video.gif"
#os.system(command)

command = "ffmpeg -pattern_type glob -framerate 1 -i 'Meinhard_tmp_" + runtype + "_*.png' -i " + runtype + "_palette.png -lavfi paletteuse -loop 1 " + runtype + "_video.mp4"
os.system(command)

command = "rm  Meinhard_tmp_" + runtype + "_*.png"
os.system(command)

# color according to species:

specieslabels = df1.columns
num2 = 0

for p in specieslabels:
    extractkey = p
    p_name = re.sub('[^A-Za-z0-9]+', '', p)

    out_pdf = ('Meinhard_tmp_%s_%s%d.pdf' % (runtype, p_name, num2))
    out_png = ('Meinhard_tmp_%s_%s%d.png' % (runtype, p_name, num2))
    print(p)
    color_map = []
    nodesize = []
    for node in nx_graph:
        nodelist=mymapper.node_info_[node]['indices']
        ids = df1.iloc[nodelist,].index.tolist() # full list of stations
        values = df1.iloc[nodelist,num2]
        nodesize.append(np.mean(values))
    nodesize_norm =  [math.sqrt(x)*2500 for x in nodesize]
    # Non normalized values
    color_map = nodesize # modify zeros
    vmin = min(color_map)
    vmax = max(color_map)
    numlist=[0 if math.isnan(x) else x for x in color_map]
    nums= np.around(np.linspace(min(color_map), max(color_map), 6),decimals=3)
    nums=[np.format_float_scientific(x,precision=2) for x in nums] # pretty format
    plt.figure()
    plt.title(p)
    print(color_map)
    #color_map = numlist_normalized or color_map
    nx.draw_networkx(nx_graph, pos, with_labels=False, node_size=5,node_color="grey", width=0.25, edge_color="grey") # base
    nx.draw_networkx(nx_graph, pos, with_labels=False, node_size=nodesize_norm,\
     node_color=color_map, vmin=vmin, vmax=vmax, cmap=plt.cm.get_cmap('viridis'), width=0.25, edge_color="grey")
    nodes =nx.draw_networkx_nodes(nx_graph, pos, with_labels=False, node_size=nodesize_norm,\
     node_color=color_map, vmin=vmin, vmax=vmax, cmap=plt.cm.get_cmap('viridis'), width=0.25, edge_color="grey")
    #sm = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap('viridis'), norm=plt.Normalize(vmin = vmin, vmax=vmax))
    cb=plt.colorbar(nodes,format='%.0e')
    tick_locs = np.linspace(cb.vmin,cb.vmax,6)
    cb.set_ticks(tick_locs)
    cb.set_ticklabels(nums)
    cb.ax.set_ylabel('Average abundance')
    plt.axis('off')
    plt.savefig(out_pdf)
    plt.savefig(out_png)
    num2 = num2 +1

plt.close('all')

command = "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=Meinhard_" + runtype + "_species.pdf Meinhard_tmp_" + runtype + "*.pdf"
os.system(command)

command = "rm  Meinhard_tmp_" + runtype + "*.pdf"
os.system(command)

command = "ffmpeg -pattern_type glob  -i 'Meinhard_tmp_" + runtype + "*.png' -vf palettegen " + runtype + "_species_palette.png"
os.system(command)

command = "ffmpeg -pattern_type glob -framerate 1 -i 'Meinhard_tmp_" + runtype + "*.png' -i " + runtype + "_species_palette.png -lavfi paletteuse -loop 1 " + runtype + "_species_video.mp4"
os.system(command)

command = "rm  Meinhard_tmp_" + runtype + "*.png"
os.system(command)

