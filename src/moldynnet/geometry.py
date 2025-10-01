import numpy as np

def pdb_to_coords(pdb_path: str):
    # super-minimal PDB parser for demo; good enough for the toy file
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    arr = np.array(coords, dtype=np.float64)
    return arr
