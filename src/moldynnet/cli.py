import typer
from rich import print
from .train import train_main
from .predict import predict_main
from .generate import generate_main
from .explain import explain_main

app = typer.Typer(help="MolDynamicsNet CLI.")

@app.command()
def train(config: str = "configs/train_toy.yaml"):
    """Train models (supervised + pretraining)."""
    # personal note: keeping the interface minimal; don't overfit flags for v0.
    train_main(config)

@app.command()
def predict(protein: str, ligand: str, out: str = "runs/pred.json"):
    """Predict ΔG and (toy) binding curve."""
    predict_main(protein, ligand, out)

@app.command()
def generate(protein: str, num: int = 4, out: str = "runs/gen.smi"):
    """Generate ligands conditioned on pocket embedding (diffusion stub)."""
    generate_main(protein, num, out)

@app.command()
def explain(protein: str, ligand: str, out: str = "runs/explain.html"):
    """Attributions/IG; residue-level HTML report."""
    explain_main(protein, ligand, out)

if __name__ == "__main__":
    app()
