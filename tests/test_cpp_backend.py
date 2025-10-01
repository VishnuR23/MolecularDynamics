import numpy as np
import moldyn_core as core

def test_compute_features_shapes():
    coords = np.array([[0,0,0],[1,0,0],[0,1,0]], dtype=np.float64)
    feats = core.compute_features(coords.reshape(-1).tolist(), 3)
    assert len(feats.pocket_center) == 3
    assert len(feats.distances) == 3  # 3 pairs
