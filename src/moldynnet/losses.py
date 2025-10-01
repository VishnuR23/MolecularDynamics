import torch
import torch.nn.functional as F
from .physics import physics_regularizer

def supervised_energy_loss(pred_dG, true_dG, distances_tensor):
    rmse = torch.sqrt(F.mse_loss(pred_dG, true_dG))
    reg = physics_regularizer(distances_tensor)
    return rmse + 0.05 * reg  # small physics weight by default

def contrastive_bound_loss(z_bound, z_unbound, margin=1.0):
    # Pull bound close to anchor; push unbound apart.
    pos = torch.norm(z_bound[0] - z_bound[1], p=2)
    neg = torch.clamp(margin - torch.norm(z_bound[0] - z_unbound[0], p=2), min=0.0)
    return pos + neg
