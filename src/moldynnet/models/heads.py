import torch.nn as nn

class BindingEnergyHead(nn.Module):
    def __init__(self, hidden=128):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(hidden, hidden),
            nn.SiLU(),
            nn.Linear(hidden, 1)
        )
    def forward(self, g):
        return self.mlp(g).squeeze(-1)  # ΔG (scalar)
