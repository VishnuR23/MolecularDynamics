import torch
from moldynnet.models.gnn import SimpleGNN
from moldynnet.models.heads import BindingEnergyHead

def test_forward_shapes():
    model = SimpleGNN()
    head = BindingEnergyHead()
    x = torch.randn(16, 64)
    g = model(x)
    y = head(g)
    assert y.shape == torch.Size([1])
