import json
import numpy as np
from .geometry import pdb_to_coords
from .utils import save_json
import moldyn_core as core

def predict_main(protein: str, ligand: str, out: str):
    coords = pdb_to_coords(protein)
    flat = coords.reshape(-1).astype(np.float64).tolist()
    feats = core.compute_features(flat, coords.shape[0])
    # toy ΔG prediction: inverse mean distance heuristic (placeholder for model call)
    mean_dist = sum(feats.distances)/max(len(feats.distances),1)
    dG = -4.0 - 3.0*(1.0/max(mean_dist,1e-6))
    curve = [dG + 0.2*i for i in range(-5,6)]
    save_json({"dG": dG, "curve": curve, "pocket_center": feats.pocket_center}, out)
