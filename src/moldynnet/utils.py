import json, os
from dataclasses import dataclass
import yaml

@dataclass
class Config:
    batch_size: int = 8
    lr: float = 1e-3
    epochs: int = 10
    hidden: int = 128

def load_config(path: str) -> Config:
    with open(path) as f:
        raw = yaml.safe_load(f)
    return Config(**raw)

def save_json(obj, path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        json.dump(obj, f, indent=2)
