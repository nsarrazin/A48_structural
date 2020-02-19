import numpy as np
from geometry import Geometry
from matrix import LoadCase


class Simulation:
    def __init__(self, **kwargs):
        self.geo = Geometry(**kwargs)
        self.case = LoadCase(parent, **kwargs)
        
    def run():
        pass