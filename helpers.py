import numpy as np


def step(x, x_ref):
    return np.clip(x-x_ref, 0, np.inf)