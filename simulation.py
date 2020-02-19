import numpy as np
from geometry import Geometry
from matrix import LoadCase
from interpolation import Interpolation

class Simulation:
    def __init__(self, **kwargs):
        self.geo = Geometry(**kwargs)
        self.case = LoadCase(parent, **kwargs)
        self.interp = Interpolation()
        
    def run():
        pass