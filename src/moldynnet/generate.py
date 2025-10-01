import torch
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from .models.diffusion import TinyDiffusion, time_embedding
from .geometry import pdb_to_coords

def latent_to_smiles(z):
    # silly mapping for demo: hash latent to fragments
    ints = (z[0].detach().cpu().numpy() * 1000).astype(int).tolist()
    # keep this simple so it's reproducible in a portfolio
    seed = abs(sum(ints)) % 3
    return ["CCO", "c1ccccc1", "CCN(CC)CC"][seed]

def generate_main(protein: str, num: int, out: str):
    coords = pdb_to_coords(protein)  # not used deeply, but gives conditioning in real impl
    model = TinyDiffusion()
    with open(out, 'w') as f:
        for i in range(num):
            z = torch.randn(1, 128)
            t = torch.rand(1)
            te = time_embedding(t, 128)
            zhat = model(z, te)
            smi = latent_to_smiles(zhat)
            f.write(smi + "\n")
