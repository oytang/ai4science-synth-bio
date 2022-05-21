import numpy as np
from sklearn import model_selection
import torch
import torch.utils.data
from torch import nn, optim
from torch.nn import functional as F
from deeprfc.featurizer import OneHotFeaturizer
from rdkit import Chem

class MolecularVAE(nn.Module):
    def __init__(self):
        super(MolecularVAE, self).__init__()

        self.conv1d1 = nn.Conv1d(120, 9, kernel_size=9)
        self.conv1d2 = nn.Conv1d(9, 9, kernel_size=9)
        self.conv1d3 = nn.Conv1d(9, 10, kernel_size=11)
        self.fc0 = nn.Linear(90, 435)
        self.fc11 = nn.Linear(435, 292)
        self.fc12 = nn.Linear(435, 292)

        self.fc2 = nn.Linear(292, 292)
        self.gru = nn.GRU(292, 501, 3, batch_first=True)
        self.fc3 = nn.Linear(501, 35)

    def encode(self, x):
        h = F.relu(self.conv1d1(x))
        h = F.relu(self.conv1d2(h))
        h = F.relu(self.conv1d3(h))
        h = h.view(h.size(0), -1)
        h = F.selu(self.fc0(h))
        return self.fc11(h), self.fc12(h)

    def reparametrize(self, mu, logvar):
        if self.training:
            std = torch.exp(0.5 * logvar)
            eps = 1e-2 * torch.randn_like(std)
            w = eps.mul(std).add_(mu)
            return w
        else:
            return mu

    def decode(self, z):
        z = F.selu(self.fc2(z))
        z = z.view(z.size(0), 1, z.size(-1)).repeat(1, 120, 1)
        out, h = self.gru(z)
        out_reshape = out.contiguous().view(-1, out.size(-1))
        y0 = F.softmax(self.fc3(out_reshape), dim=1)
        y = y0.contiguous().view(out.size(0), -1, y0.size(-1))
        return y

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparametrize(mu, logvar)
        return self.decode(z), mu, logvar
    
    def sample(self, x, decode_attempts=1000, noise_norm=1.0):
        target_atoms = ['C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'H']
        
        oh = OneHotFeaturizer()
        mu, logvar = self.encode(x)
        z = self.reparametrize(mu, logvar)
        z2 = z.cpu().detach().numpy()
        smiles_results = {}
        for i in range(decode_attempts):
            Z = self.perturb_z(z2, noise_norm, False)
            z3 = torch.from_numpy(Z.astype(np.float32)).to('cuda')
            X = self.decode(z3)
            X = X.cpu().detach().numpy()
            y = np.argmax(X, axis=2)
            smiles = oh.decode_smiles_from_index(y[0])
            mol = Chem.MolFromSmiles(smiles)
            if mol != None:
                smiles_results[smiles]=mol
                
        for i in range(decode_attempts):
            Z = self.perturb_z(z2, noise_norm, True)
            z3 = torch.from_numpy(Z.astype(np.float32)).to('cuda')
            X = self.decode(z3)
            X = X.cpu().detach().numpy()
            y = np.argmax(X, axis=2)
            smiles = oh.decode_smiles_from_index(y[0])
            mol = Chem.MolFromSmiles(smiles)
            if mol != None:
                smiles_results[smiles]=mol
                
        for smiles in smiles_results:
            flag = True
            each_tmp_mol = smiles_results[smiles]
            for each_atom in each_tmp_mol.GetAtoms():
                if each_atom.GetSymbol() not in target_atoms:
                    flag = False
            if flag:        
                print (smiles)
        return
    

    def perturb_z(self, z, noise_norm, constant_norm=False):
        if noise_norm > 0.0:
            noise_vec = np.random.normal(0, 1, size=z.shape)
            noise_vec = noise_vec / np.linalg.norm(noise_vec)
            if constant_norm:
                return z + (noise_norm * noise_vec)
            else:
                noise_amp = np.random.uniform(0, noise_norm, size=(z.shape[0], 1))
                return z + (noise_amp * noise_vec)
        else:
            return z