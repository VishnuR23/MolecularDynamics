import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

def load_labels_csv(path: str):
    df = pd.read_csv(path)
    # expecting columns: protein_id, ligand_id, dG
    return df

def load_ligand_sdf(path: str):
    mol = Chem.SDMolSupplier(path, removeHs=False, sanitize=True)[0]
    if mol is None:
        raise ValueError(f"Failed to read SDF: {path}")
    AllChem.EmbedMolecule(mol, randomSeed=42)  # if no 3D
    return mol
