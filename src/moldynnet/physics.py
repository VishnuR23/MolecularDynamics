import torch

def lennard_jones(distance, epsilon=0.1, sigma=3.5):
    # I'm purposefully keeping this simple for stability in the loss.
    r = distance.clamp(min=1e-6)
    sr6 = (sigma / r) ** 6
    return 4 * epsilon * (sr6 * sr6 - sr6)

def coulomb(qi, qj, r, k=8.9875517923e9):
    r = r.clamp(min=1e-6)
    return k * qi * qj / (r * r)  # toy; not unit-consistent

def physics_regularizer(distances):
    # Encourage learned distances to avoid unphysical collapse/explosion.
    lj = lennard_jones(distances).mean()
    return torch.relu(lj)  # bias toward reasonable basin
