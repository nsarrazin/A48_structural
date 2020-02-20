import numpy as np
from geometry import Geometry
from load import LoadCase
from interpolation import Interpolation
from data.consts import parameters_case, parameters_geometry
class Simulation:
    def __init__(self, **kwargs):
        self.geo = Geometry(**kwargs)
        self.interp = Interpolation()
        self.case = LoadCase(self, **kwargs)
        
    def run():
        pass

if __name__ == "__main__":
    # MMoI, SC, J
    params = {}
    params.update(parameters_case)
    params.update(parameters_geometry)
    # dic = parameters_geometry.update(parameters_case)
    # print(dic)
    sim = Simulation(**params)
