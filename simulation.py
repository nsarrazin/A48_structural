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

    print(sim.x)
    print(sim.BCs)
    print(sim.solution.sol)

    sols = sim.x

    print(f"Sum of forces in the y-direction : {sols[0]+sols[2]+sols[4]+sols[6]*sim.case.a_y+sim.case.P*sim.case.a_y:.2f}N")
    print(f"Sum of forces in the z-direction : {sols[1]+sols[3]+sols[5]+sols[6]*sim.case.a_z+sim.case.P*sim.case.a_z+sim.interp.integrate_q(sim.geo.l_a,ord=1)[-1]:.2f}N")
    print(f"factor for v_y:{-1/(sim.case.E*sim.geo.MMoI[0])}")
    print(f"factor for v_z:{-1/(sim.case.E*sim.geo.MMoI[1])}")
    print(sim.case.B)

    def v_y(x, sol):
        return -1/(sim.case.E*sim.geo.MMoI[0])*(sol["Fa"]*sim.case.a_z/6*step(x,sim.geo.x_2-sim.geo.x_a/2)**3\
            +sim.case.P*sim.case.a_z*step(x,sim.geo.x_2+sim.geo.x_a/2)**3\
                -sol["Fz_1"]/6*step(x,sim.geo.x_1)**3-sol["Fz_2"]/6*step(x,sim.geo.x_2)**3-sol["Fz_3"]/6*step(x,sim.geo.x_3)**3\
                    +sol["C3"]*x+sol["C4"]*0)        
    xs = np.linspace(0, sim.geo.l_a, 100)
    ys = [v_y(i, sol) for i in xs]
    
    offset = v_y(sim.geo.x_2, sol)*0

    plt.plot(xs, ys-offset)
    
    for n,x in enumerate([sim.geo.x_1, sim.geo.x_2, sim.geo.x_3]):
        plt.axvline(x=x, linestyle="dashed", linewidth=1.2, color=f"C{n}")

    for n,y in enumerate([sim.case.d_1, 0, sim.case.d_3]):
        plt.axhline(y=y, linestyle="dashed", linewidth=1.2, color=f"C{n}")

    plt.show()

    def v_z(x,sol):
        return -1/(sim.case.E*sim.geo.MMoI[1])* (sol["Fy_1"]*-1/6*step(x,sim.geo.x_1)**3-sol["Fy_2"]*1/6*step(x,sim.geo.x_2)**3\
                  -sol["Fy_3"]*1/6*step(x,sim.geo.x_3)**3+sol["Fa"]*sim.case.a_y/6*step(x,sim.geo.x_2-sim.geo.x_a/2)**3\
                      +sol["C1"]*x+sol["C2"]*0\
                          +sim.interp.integrate_q(x, ord=4)[-1] + (sim.case.P*sim.case.a_y)/6*step(x,sim.geo.x_2+sim.geo.x_a/2)**3)
    def theta(x,sol):
        return 1/(sim.case.G*sim.geo.J)*(sol["Fy_1"]*-sim.case.z_sc*step(x,sim.geo.x_1)**1-sol["Fy_2"]*sim/case.z_sc*step(x,sim.geo.x_2)**1\
                  -sol["Fy_3"]*sim.case.z_sc*step(x,sim.geo.x_3)**1\
                      +sol["Fa"]*sim.case.a_y*sim.case.z_sc*step(x,sim.geo.x_2-sim.geo.x_a/2)**1-sol["Fa"]*sim.case.a_m*step(x,sim.geo.x_2-sim.geo.x_a/2)**1\
                          +sol["C5"]\
                              +sim.interp.integrate_tau(x, z_sc=sim.geo.h/2+sim.case.z_sc, ord=2)[-1] + sim.case.P*sim.case.a_y*sim.case.z_sc*step(x,sim.geo.x_2+sim.geo.x_a/2)**1-sim.case.P*sim.case.a_m*step(x,sim.geo.x_2+sim.geo.x_a/2)**1)
    
    xs = np.linspace(0, sim.geo.l_a, 100)
    ys = [v_z(i, sol) for i in xs]

    plt.plot(xs,ys)
    plt.show()
