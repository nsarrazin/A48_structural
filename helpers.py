import numpy as np


def step(x, x_ref, power=1):
    val = np.clip(x-x_ref, 0, np.inf)
    if val == 0 and power == 0:
        return 0
    return val**power