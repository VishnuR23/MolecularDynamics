import torch
from torch.optim import AdamW
from .utils import load_config
from .models.gnn import SimpleGNN
from .models.heads import BindingEnergyHead
from .losses import supervised_energy_loss

# note: I'm intentionally keeping this lightweight so it's easy to run on a laptop.
def train_main(config: str):
    cfg = load_config(config)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = SimpleGNN(in_dim=64, hidden=cfg.hidden).to(device)
    head = BindingEnergyHead(hidden=cfg.hidden).to(device)
    params = list(model.parameters()) + list(head.parameters())
    opt = AdamW(params, lr=cfg.lr)

    # toy data: random embeddings + synthetic distances/labels
    # (In practice, feed real protein-ligand graphs and distances from C++ backend)
    for epoch in range(cfg.epochs):
        x = torch.randn(32, 64, device=device)
        distances = torch.rand(32, device=device) * 5.0 + 1.0
        true_dG = torch.randn(1, device=device) * 2.0 - 8.0  # around -8 kcal/mol

        g = model(x)
        pred = head(g)
        loss = supervised_energy_loss(pred, true_dG, distances)
        opt.zero_grad()
        loss.backward()
        opt.step()

        if (epoch+1) % 2 == 0:
            print(f"[epoch {epoch+1}] loss={loss.item():.4f}")  # keeping it human-readable
