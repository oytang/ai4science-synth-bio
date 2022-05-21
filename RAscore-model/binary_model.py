from pip import main
from RAscore import RAscore_XGB #For XGB based models
xgb_scorer = RAscore_XGB.RAScorerXGB()

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect 
from rdkit import DataStructs
import numpy as np

import pandas as pd
df_bioreachable = pd.read_csv("bioreachable_molesules_fromql.csv",index_col=0) #可合成的分子库，文件可替换

bioreachable_mols = [Chem.MolFromSmiles(s) for s in df_bioreachable["SMILES"] if type(s) is not None]
fps = [AllChem.GetMorganFingerprintAsBitVect(m, radius=4, nBits=2048) for m in bioreachable_mols] #可保存加速检索


def highest_tanimoto_precalc_fps(mol, fps):

    if fps is None or len(fps) == 0:
        return 0

    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol, 4, 2048)
    sims = np.array(DataStructs.BulkTanimotoSimilarity(fp1, fps))

    return sims.max() 


def is_bioreachable(molecule):
    sim_score = highest_tanimoto_precalc_fps(mol=Chem.MolFromSmiles(molecule), fps=fps)
    if sim_score>0.5:
        ra_score = xgb_scorer.predict(molecule)
        if ra_score>0.5:
            return 1
        else:
            return 0
    else:
        return 0

if __name__ == '__main__':
    test_molecule = 'CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC'
    res = is_bioreachable(test_molecule)
    print(res)