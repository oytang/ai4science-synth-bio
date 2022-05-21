import logging
import os
import time
import pkg_resources
import shutil
import copy
import pandas as pd
import sys
import sklearn
import pickle
import warnings
import os
import glob
import numpy as np
import pandas as pd
from rdkit import RDLogger
import warnings , os
import pickle
from deeprfc import utils
from keras.models import model_from_json 
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit import Chem

import torch
from deeprfc.featurizer import OneHotFeaturizer
from deeprfc.models import MolecularVAE


 
# os.environ["CUDA_VISIBLE_DEVICES"] = "2"
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

def read_input_file(filename):
    input_info = {}
    with open(filename, 'r') as fp:
        fp.readline()
        for line in fp:
            sptlist = line.strip().split('\t')
            label = sptlist[0].strip()
            reactant = sptlist[1].strip()
            product = sptlist[2].strip()
            if len(reactant) <= 120:
                if len(product) <=120:
                    mol1 = Chem.MolFromSmiles(reactant)
                    mol2 = Chem.MolFromSmiles(product)
                    mol1 = Chem.RemoveHs(mol1)
                    mol2 = Chem.RemoveHs(mol2)
                    reactant = Chem.MolToSmiles(mol1)
                    product = Chem.MolToSmiles(mol2)
                    input_info[label] = [reactant, product]
    return input_info

def calc_z(model, smi):
    start = smi
    start = start.ljust(120)
    oh = OneHotFeaturizer()
    start_vec = torch.from_numpy(oh.featurize([start]).astype(np.float32)).to(device)
    z = model(start_vec)[1].cpu().detach().numpy()
    return z

def calculate_features(input_info, vae_model):
    r=2
    feature_info = {}
    for each_label in input_info:
        
        reactant_smi = input_info[each_label][0]
        product_smi = input_info[each_label][1]
        
        z1 = calc_z(vae_model, reactant_smi)
        z1 = z1[0]
        
        z2 = calc_z(vae_model, product_smi)
        z2 = z2[0]
        
        reactant_mol = Chem.MolFromSmiles(reactant_smi)
        fp = AllChem.GetMorganFingerprintAsBitVect(reactant_mol, radius=r, nBits=1024)
        bits = fp.ToBitString()
        reactant_feature = []
        for f in z1:
            reactant_feature.append(float(f))
        for f in bits:
            reactant_feature.append(int(f))
            
        product_mol = Chem.MolFromSmiles(product_smi)
        fp = AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=r, nBits=1024)
        bits = fp.ToBitString()
        product_feature = []
        for f in z2:
            product_feature.append(float(f))
        for f in bits:
            product_feature.append(int(f))
            
        reactant_feature = [reactant_feature]
        product_feature = [product_feature]
        reactant_feature = np.asarray(reactant_feature)
        product_feature = np.asarray(product_feature)
        
        feature_info[each_label]=[reactant_feature, product_feature]    
    return feature_info

def predict_reaction_feasibility(feature_info, model):    
    results = {}
    for each_label in feature_info:
        X1 = feature_info[each_label][0]
        X2 = feature_info[each_label][1]
        
        val_list = []
        for i in range(10):
            pred_result = model.predict([X1, X2])
            val = pred_result[0][0]
            val_list.append(val)
        mean_val = np.mean(val_list)
        std_val = np.std(val_list)
        final_val = mean_val - (std_val*0.5)
        if final_val >= 0.32:
            feasibility = 'feasible'
        else:
            feasibility = 'infeasible'

        results[each_label] = [mean_val, std_val,feasibility]
        
    return results

def main():   
    warnings.filterwarnings(action='ignore')
    
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    
    warnings.filterwarnings('ignore')
    start = time.time()
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
    
    weight_file = pkg_resources.resource_filename('deeprfc', 'data/final_model.h5') 
    model_file = pkg_resources.resource_filename('deeprfc', 'data/final_model.json')
    vae_model_file = pkg_resources.resource_filename('deeprfc', 'data/vae_model.pth')
    
    parser = utils.argument_parser()
    
    options = parser.parse_args()
    input_file = options.input_file
    output_dir = options.output_dir
        
    try:
        os.mkdir(output_dir)
    except:
        pass
    
    vae_model = MolecularVAE()
    vae_model.load_state_dict(torch.load(vae_model_file))
    vae_model.to(device)
    vae_model.eval()
    
    input_info = read_input_file(input_file)
    feature_info = calculate_features(input_info, vae_model)
    
    ## load model
    json_file = open(model_file, "r")
    loaded_model_json = json_file.read() 
    json_file.close()
    loaded_model = model_from_json(loaded_model_json)
    loaded_model.load_weights(weight_file)
    
    results = predict_reaction_feasibility(feature_info, loaded_model)
    
    with open(output_dir+'/result.txt', 'w') as fp:
        fp.write('ID\tPredictive_mean\tStd\tFeasibility\n')
        for label in results:
            fp.write('%s\t%s\t%s\t%s\n'%(label, results[label][0], results[label][1],results[label][2]))
    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))
    
    
if __name__ == '__main__':
    main()
