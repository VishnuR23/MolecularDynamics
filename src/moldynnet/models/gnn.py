import torch
import torch.nn as nn
import torch.nn.functional as F

# Note: keeping a compact GNN; enough for a portfolio but easy to read.
class SimpleGNN(nn.Module):
    def __init__(self, in_dim=64, hidden=128, layers=3):
        super().__init__()
        self.embed = nn.Linear(in_dim, hidden)
        self.layers = nn.ModuleList([nn.Linear(hidden, hidden) for _ in range(layers)])
        self.readout = nn.Linear(hidden, hidden)

    def forward(self, x, edge_index=None, edge_attr=None):
        h = F.silu(self.embed(x))
        for layer in self.layers:
            # simple residual MLP (replace with GAT/GIN/EGNN if desired)
            h = h + F.silu(layer(h))
        g = h.mean(dim=0, keepdim=True)  # graph-level
        return self.readout(g)  # [1, hidden]
