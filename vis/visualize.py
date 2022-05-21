import pandas as pd
import numpy as np
from tqdm import tqdm, trange
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
import matplotlib.pyplot as plt
from sklearn import manifold, datasets
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sb
import numpy as np

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE

from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.DataManip.Metric.rdMetricMatrixCalc import GetTanimotoSimMat

from math import pi

vischemmnx = pd.read_csv('mnx_chem_bioreachable.csv', sep=',')
# vischemmnx['intra'] = '1'
# vischemmnx['database'] = 'MetaNetX'
# #print(vischemmnx)
# vischemmnx.iloc[:,[8,9,10]].to_csv('vischemmnx.csv', index=False)
vischemmnx.iloc[:,8].to_csv('vischemmnx.csv', index=False)

vischemPub = pd.read_csv('mnx_chem_bioreachable_PubChemBioCyc.csv')
# vischemPub['intra'] = '1'
# vischemPub['database'] = 'PubChemBioCyc'
#print(vischemmnx)
# vischemPub.iloc[:,[8,9,10]].to_csv('vischemPub.csv', index=False)
vischemPub.iloc[:,8].to_csv('vischemPub.csv', index=False)

def ClusterFps(fps,cutoff=0.2):

    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    return cs

bioreachable =vischemmnx.dropna(subset=['SMILES'])
bioreachable.index = range(len(bioreachable))
# metabolites_list = []
# for i in trange(len(reac)):
#     metabolites_list += reac.loc[i, 'metabolites']
# metabolites_list = list(set(metabolites_list))
# bioreachable = np.array(bioreachable, dtype=object)
# bioreachable = bioreachable.tolist()

#
# bioreachable_mols = [Chem.MolFromSmiles(s) for s in bioreachable["SMILES"] if type(s) is not None]
# fps = [AllChem.GetMorganFingerprintAsBitVect(m, radius=4, nBits=2048) for m in bioreachable_mols]
# clusters = ClusterFps(fps,cutoff=0.5)

smi = list(bioreachable['SMILES'])
bioreachable_mols = [Chem.MolFromSmiles(s) for s in bioreachable["SMILES"] if type(s) is not None]
fps = [AllChem.GetMorganFingerprintAsBitVect(m, radius=4, nBits=2048) for m in bioreachable_mols]
#fps = [MACCSkeys.GenMACCSKeys(x) for x in smi] #will use MACCSKeys for this
tanimoto_sim_mat = GetTanimotoSimMat(fps) #computes a similartity matrix between all the molecules
n_mole = len(fps)
similarity_matrix = np.ones([n_mole,n_mole])
i_lower = np.tril_indices(n = n_mole, m = n_mole, k = -1)
i_upper = np.triu_indices(n = n_mole, m = n_mole, k = 1)
similarity_matrix[i_lower] = tanimoto_sim_mat
similarity_matrix[i_upper] = similarity_matrix.T[i_upper]
distance_matrix = np.subtract(1,similarity_matrix) #similarity matrix of all vs all molecules in our table

TSNE_sim = TSNE(n_components = 2, init = 'pca', random_state = 90, angle = 0.3,
                perplexity = 50).fit_transform(distance_matrix) #tune the parameters according to your dataset
tsne_result = pd.DataFrame(data = TSNE_sim , columns = ["TC1","TC2"]) #new table containing tSNE results
print("t-SNE analysis results\n", tsne_result.head(5))

print("Beginning K-means clustering")

range_n_clusters = [2, 3, 4, 5, 6, 7, 8, 9, 10] #explore a range of cluster sizes to best clasify our molecules
best_silhouette_score = -1
best_cluster_size = 0

for n_clusters in range_n_clusters:

    kmeans = KMeans(n_clusters = n_clusters, random_state =  10)
    cluster_labels = kmeans.fit_predict(tsne_result[['TC1','TC2']])
    silhouette_avg = silhouette_score(tsne_result[['TC1','TC1']], cluster_labels) ###check this

    if silhouette_avg > best_silhouette_score:
        best_silhouette_score = silhouette_avg
        best_cluster_size = n_clusters

    #print silhouette scores. scores between [-1,1] with 1 being the best, hence the better our data is
    #distributed inside the clusters
    print("For n_clusters =", n_clusters,
            "The average silhouette_score is :", silhouette_avg)


    kmeans = KMeans(n_clusters = best_cluster_size, random_state = 10) ##add automation for n-clusters for best performing cluster
    clusters = kmeans.fit(tsne_result[['TC1','TC2']])
    tsne_result['Cluster'] = pd.Series(clusters.labels_, index = tsne_result.index)
    print("t-SNE results with cluster number for each element\n", tsne_result.head(5)) #tSNE table now contains the number of clusters for each element

    print("Plotting t-SNE and K-means results")
    plt.rcParams['axes.linewidth'] = 1.5
    fig, ax = plt.subplots(figsize = (6, 6))
    ax = sb.scatterplot(x = 'TC1', y = 'TC2', data = tsne_result, hue = 'Cluster', s = 22, palette = sb.color_palette("Set2", 2),
                        linewidth = 0.2, alpha = 1)

    plt.xlabel('tSNE 1', fontsize = 24, fontweight = 'bold')
    plt.ylabel('tSNE 2', fontsize = 24, fontweight = 'bold')
    plt.tick_params('both', width = 2, labelsize = 18)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    handles, labels = ax.get_legend_handles_labels()

    ax.legend(handles = handles[1:], labels = labels[1:])
    plt.legend(loc = 'best', frameon = False, prop = {'size': 16}, ncol = 2)

    #plt.tight_layout()
    plt.show()

# def highest_tanimoto_precalc_fps(mol, fps):
#
#     if fps is None or len(fps) == 0:
#         return 0
#
#     fp1 = AllChem.GetMorganFingerprintAsBitVect(mol, 4, 2048)
#     sims = np.array(DataStructs.BulkTanimotoSimilarity(fp1, fps))
#
#     return sims.max()
#
#
#
#
# '''X是特征，不包含target;X_tsne是已经降维之后的特征'''
# tsne = manifold.TSNE(n_components=2, init='pca', random_state=501)
# X_tsne = tsne.fit_transform(X)
# print("Org data dimension is {}. Embedded data dimension is {}".format(X.shape[-1], X_tsne.shape[-1]))
#
# '''嵌入空间可视化'''
# x_min, x_max = X_tsne.min(0), X_tsne.max(0)
# X_norm = (X_tsne - x_min) / (x_max - x_min)  # 归一化
# plt.figure(figsize=(8, 8))
# for i in range(X_norm.shape[0]):
#     plt.text(X_norm[i, 0], X_norm[i, 1], str(y[i]), color=plt.cm.Set1(y[i]),
#              fontdict={'weight': 'bold', 'size': 9})
# plt.xticks([])
# plt.yticks([])
# plt.show()


df_a = pd.read_csv('vischemmnx.csv', sep=',')
df_b = pd.read_csv('vischemPub.csv', sep=',')
x_list = pd.concat([df_a,df_b])

#merge = x.drop_duplicates(subset=['SMILES'], keep=False)
x_list.to_csv('merge.csv', index=False)

# #jiaoji insec; chaji differa=mnx differb=pub;bingji merge mergex include yuan
# df_insec = a[a['SMILES'].isin(b['SMILES'])]
# df_insec['secdatabase'] = 'PubChemBioCyc'
# df_differa = a.append(df_insec).drop_duplicates(subset=['SMILES'], keep=False)
# df_differb = b.append(df_insec).drop_duplicates(subset=['SMILES'], keep=False)
# df_mergex = pd.concat([df_insec,df_differb,df_differa])
# #df_mergex = df_mergex.drop_duplicates(subset=['SMILES'], keep=False)
# df_differa.to_csv('differa.csv', index=False)
# df_differb.to_csv('differb.csv', index=False)
# df_insec.to_csv('insec.csv', index=False)
# df_mergex.to_csv('mergex.csv', index=False)


