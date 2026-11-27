"""Two idealized acceptor directions about a neutral sp3 alcohol oxygen."""

import numpy as np


def _unit(x):
    n = np.linalg.norm(x, axis=-1, keepdims=True)
    if np.any(n < 1.0e-10):
        raise ValueError("Undefined alcohol acceptor geometry")
    return x / n


def lp_gate(oxygen, carbon, hydrogens, incoming_h):
    # Both lone pairs point away from the C/O/H bond bisector. Their
    # dot products with each ideal O–C and O–H arm equal cos(109.47°).
    oc = _unit(carbon - oxygen)
    oh = _unit(hydrogens - oxygen)
    away = -_unit(oc + oh) / np.sqrt(3)
    normal = np.sqrt(2 / 3) * _unit(np.cross(oc, oh))
    toward = _unit(incoming_h - oxygen)
    return (
        np.maximum(0, np.maximum((away + normal) @ toward, (away - normal) @ toward))
        ** 2
    )
