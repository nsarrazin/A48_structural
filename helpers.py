import numpy as np


def step(x, x_ref, pow=1):
    val = np.clip(x-x_ref, 0, np.inf)
    if val == 0 and pow == 0:
        return 0
    return val**pow