import torch
import torch.nn as nn

# Tiny diffusion-like stub for ligand SMILES latent; good enough to demo API.
class TinyDiffusion(nn.Module):
    def __init__(self, latent_dim=128):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(latent_dim, 256),
            nn.SiLU(),
            nn.Linear(256, latent_dim)
        )
    def forward(self, z, t_embed):
        return self.net(z + t_embed)

def time_embedding(t, dim=128):
    # positional-style embedding for time step t (scalar 0..1)
    import math, torch
    half = dim // 2
    freqs = torch.exp(torch.linspace(0, math.log(10000), half))
    angles = t[:, None] * freqs[None, :]
    return torch.cat([torch.sin(angles), torch.cos(angles)], dim=-1)
