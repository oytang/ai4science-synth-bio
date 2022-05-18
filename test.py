import pandas as pd
import numpy as np
from tqdm import tqdm, trange
reac = pd.read_csv('reac_prop.tsv', sep='\t', header=351) # skip 351 lines of documentation

def preprocess(self):
    """
    rex: one reaction, such as '1 MNXM01@MNXD1 = 1 MNXM1@MNXD1'

    output:
        - reaction formula without compartment (MNXM01 instead of MNXM01@MNXD1)
        - list of metabolites involved in the reaction
        - list of substrates
        - list of products
    """
    # compartment can only be @MNXD1 or @MNXD2
    rex_clean = rex.replace('@MNXD1', '').replace('@MNXD2', '')

    metabolites = take_MNXM(rex_clean.split(' '))

    # drop duplicates
    metabolites = list(set(metabolites))

    substrates, products = rex_clean.split('=')
    substrates = take_MNXM(substrates.split(' '))
    products = take_MNXM(products.split(' '))

 #   return rex_clean, metabolites, substrates, products
    print(rex_clean,metabolites, substrates, products)
p= preprocess(self)
print(p)
def take_MNXM(str_list):
    """
    Helper function for preprocessing.
    """
    return [mol for mol in str_list if mol.startswith('MNXM')]