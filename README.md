# MolDynamicsNet — Hybrid ML + Physics for Drug Binding Prediction

> Hybrid **physics + ML** pipeline for protein–ligand binding: C++ backend (pybind11) for fast feature extraction, **PyTorch** GNNs for ΔG prediction and binding curves, **self-supervised** pretraining, **physics-regularized** loss, and an optional **diffusion generator** for ligand proposals. Includes Streamlit dashboard, CLI, Docker/CUDA, and Linux-cluster friendly Makefiles.

## TL;DR
- **C++/pybind11** backend (`cpp/`) computes fast structural features (distance, HBonds, RMSD-ish proxies) and can be extended to mini-MD.
- **PyTorch** models (`src/moldynnet/`) with GNN encoder + binding-energy head + contrastive pretrain + physics regularizer.
- **Ligand diffusion generator** (`src/moldynnet/generative/`) conditioned on pocket embeddings.
- **Command-line** pipeline: `moldynnet` (train / predict / generate / explain).
- **Dashboard**: Streamlit app for energy landscape and attribution.
- **Docker** + **Make** workflows for local or GPU machines.
- **Tests** and small **toy dataset** to run smoke checks.

> ⚠️ Note: For a demo, the C++ backend computes **physically motivated features** without heavy external MD deps. If you have OpenMM/GROMACS, plug it into `cpp/backend.cpp` (hooks included).

---

## Quickstart

```bash
# 1) Python environment (CUDA optional but recommended)
python3 -m venv .venv && source .venv/bin/activate
pip install -U pip
pip install -r requirements.txt

# 2) Build C++ extension
make cpp

# 3) Run smoke tests
make test

# 4) Train (toy data)
moldynnet train --config configs/train_toy.yaml

# 5) Predict ΔG and curve for a protein-ligand pair
moldynnet predict --protein data/toy/protein.pdb --ligand data/toy/ligand.sdf --out runs/demo_pred.json

# 6) Generate ligands (diffusion; conditioned on pocket embedding)
moldynnet generate --protein data/toy/protein.pdb --num 5 --out runs/gen_ligands.smi

# 7) Explainability (attributions/IG on residues)
moldynnet explain --protein data/toy/protein.pdb --ligand data/toy/ligand.sdf --out runs/explain.html

# 8) Streamlit dashboard
streamlit run app/app.py
```

### Docker (CUDA 12.x)
```bash
# Build
docker build -t moldynnet:cuda -f docker/Dockerfile.cuda .

# Run with GPU
docker run --gpus all -it --rm -p 8501:8501 -v $PWD:/workspace moldynnet:cuda bash

# Inside container
make cpp && make test
moldynnet train --config configs/train_toy.yaml
```

---

## Project Layout
```
MolDynamicsNet/
├─ README.md
├─ LICENSE
├─ .gitignore
├─ requirements.txt
├─ setup.py
├─ pyproject.toml
├─ Makefile
├─ docker/
│  ├─ Dockerfile.cuda
│  └─ Dockerfile.cpu
├─ cpp/
│  ├─ CMakeLists.txt
│  ├─ bindings.cpp
│  ├─ backend.hpp
│  ├─ backend.cpp
│  └─ utils.hpp
├─ src/moldynnet/
│  ├─ __init__.py
│  ├─ cli.py
│  ├─ data.py
│  ├─ geometry.py
│  ├─ physics.py
│  ├─ models/
│  │  ├─ gnn.py
│  │  ├─ heads.py
│  │  └─ diffusion.py
│  ├─ losses.py
│  ├─ train.py
│  ├─ predict.py
│  ├─ generate.py
│  ├─ explain.py
│  └─ utils.py
├─ app/
│  └─ app.py
├─ configs/
│  ├─ train_toy.yaml
│  └─ model.yaml
├─ tests/
│  ├─ test_cpp_backend.py
│  ├─ test_models.py
│  └─ test_cli.py
└─ data/
   └─ toy/
      ├─ protein.pdb
      ├─ ligand.sdf
      └─ labels.csv
```

---

## Highlights
- **Physics-regularized loss** encourages the model to respect coarse Lennard–Jones/Coulombic patterns using learned embeddings.
- **Contrastive pretraining**: bound vs. unbound snapshots; masked-atom-style denoising.
- **Attribution**: Integrated Gradients & residue-level saliency with HTML report.
- **Linux-first**: Make targets (`make cpp`, `make test`, `make train`) and zero-UI CLI for cluster use.

---

## Extending to Real MD
- Implement short restrained MD via OpenMM in `cpp/backend.cpp` (hooks and stubs included).
- Stream XYZ or DCD frames; pass into Python via pybind11 NumPy buffers.
- Aggregate time-series features into `timepool` encodings (mean/max/attention).

---

## License
MIT (see `LICENSE`). If you adapt for proprietary use, please comply with third-party deps (pybind11, PyTorch).
