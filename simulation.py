import numpy as np
from geometry import Geometry
from load import LoadCase
from interpolation import Interpolation
from solution import Solution
from data.consts import parameters_case, parameters_geometry
from helpers import step
import matplotlib.pyplot as plt

class Simulation:
    def __init__(self, **kwargs):
        self.geo = Geometry(**kwargs)
        self.interp = Interpolation()
        self.case = LoadCase(self, **kwargs)
        self.solution = Solution(self)
        self.x = np.zeros(13)

    def run(self):
        self.x = np.linalg.solve(sim.case.A, sim.case.B)
        return self.x
    
    @property
    def BCs(self):
        BCs = ["Vy(la)", "Vz(la)", "My(la)", "Mz(la)", "T(la)","vy'(x1)", "vz'(x1)", "vy'(x2)","vz'(x2)", "vy'(x3)","vz'(x3)","vz(act)"]

        A = list(sim.case.A)

        dic = {}
        for BC, row in zip(BCs, A):
            dic[BC] = row
        
        return dic

if __name__ == "__main__":
    params = {}
    params.update(parameters_case)
    params.update(parameters_geometry)

    sim = Simulation(**params)
    sim.run()

    # print(sim.x)
    # print(sim.BCs)
    # print(sim.solution.sol)

    # sols = sim.x
    for key in sim.BCs.keys():
        print(f"{key} - {sim.BCs[key]}")
    print(sim.BCs)
    
    # print(f"Sum of fo/r v_z:{-1/(sim.case.E*sim.geo.MMoI[1])}")
    # print(sim.case.B)

    sim.solution.plot()


    # [print(n[-1]) for n in sim.BCs.values()]