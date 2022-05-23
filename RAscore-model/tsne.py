import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import datasets
from sklearn.manifold import TSNE
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit import DataStructs
import numpy as np

def transform_mol(mol, radius, nBits):
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)
    arr = np.zeros((nBits))
    DataStructs.ConvertToNumpyArray(fp, arr)
    return list(arr)


# 对样本进行预处理并画图
def plot_embedding(data, label, title):
    """
    :param data:数据集
    :param label:样本标签
    :param title:图像标题
    :return:图像
    """
    x_min, x_max = np.min(data, 0), np.max(data, 0)
    data = (data - x_min) / (x_max - x_min)     # 对数据进行归一化处理
    fig = plt.figure()      # 创建图形实例
    ax = plt.subplot(111)       # 创建子图
    # 遍历所有样本
    print(data.shape)
    for i in range(data.shape[0]):
        # 在图中为每个数据点画出标签
        print(data[i, 0], data[i, 1], str(label[i]))
        plt.text(data[i, 0], data[i, 1], str(label[i]), color=plt.cm.Set1(label[i] / 2),
                 fontdict={'weight' :  'bold' ,  'size' : 10})
    plt.xticks()        # 指定坐标的刻度
    plt.yticks()
    plt.title(title, fontsize=14)
    # 返回值
    return fig


# 主函数，执行t-SNE降维
if __name__ == '__main__':
    # 数据集需要换下
    df_bioreachable = pd.read_csv("bioreachable_molesules_fromql.csv", index_col=0)[:1500]
    df_inbioreachable = pd.read_csv("inbioreachable_molesules_fromql.csv", index_col=0)[:1500]
    bioreachable_mols = [Chem.MolFromSmiles(s) for s in df_bioreachable["SMILES"] if type(s) is not None]
    bioreachable_fps = np.array([transform_mol(m, radius=4, nBits=2048) for m in bioreachable_mols])  # 可保存加速检索
    inbioreachable_mols = [Chem.MolFromSmiles(s) for s in df_inbioreachable["SMILES"] if type(s) is not None]
    inbioreachable_fps = np.array([transform_mol(m, radius=4, nBits=2048) for m in inbioreachable_mols])  # 可保存加速检索
    all_fps = np.vstack([bioreachable_fps, inbioreachable_fps])
    bioreachable_label = np.array([1 for i in range(len(bioreachable_fps))])
    inbioreachable_label = np.array([0 for i in range(len(inbioreachable_fps))])
    all_label = np.hstack([bioreachable_label, inbioreachable_label])
    data, label = all_fps, all_label    # 调用函数，获取数据集信息
    print(data.shape)
    print(label)
    print( 'Starting compute t-SNE Embedding...' )
    ts = TSNE(n_components=2, init= 'pca' , random_state=0)
    # t-SNE降维
    reslut = ts.fit_transform(data)
    # 调用函数，绘制图像

    fig = plot_embedding(reslut, label,  't-SNE Embedding of ECFP' )
    # 显示图像
    plt.show()

