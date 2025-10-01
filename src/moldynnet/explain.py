from .geometry import pdb_to_coords
from .utils import save_json

def explain_main(protein: str, ligand: str, out: str):
    # Minimal HTML report for the portfolio; shows where you'd dump IG/SHAP.
    coords = pdb_to_coords(protein)
    html = f"""
    <html><body>
    <h2>MolDynamicsNet Explainability Report</h2>
    <p>Protein atoms: {coords.shape[0]}</p>
    <p>(Demo) Residue saliency would be visualized here with color mapping.</p>
    </body></html>
    """
    import os
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f:
        f.write(html)
