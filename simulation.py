import numpy as np
from geometry import Geometry
from load import LoadCase
from interpolation import Interpolation
from data.consts import parameters_case, parameters_geometry
from helpers import step
import matplotlib.pyplot as plt

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
    
    # print(sim.case.T(sim.geo.l_a))
    # print(sim.case.B.shape)
    BCs = ["Vy(la)", "Vz(la)", "My(la)", "Mz(la)", "T(la)","vy'(x1)", "vz'(x1)", "vy'(x2)","vz'(x2)", "vy'(x3)","vz'(x3)","vz(act)"]

    A = list(sim.case.A)

    for BC, row in zip(BCs, A):
        print(BC, row)
    
    sols = np.linalg.solve(sim.case.A, sim.case.B)

    

    # names = ["Fy_1", "Fz_1", "Fx_2", "Fy_2", "Fz_2", "Fy_3", "Fz_3", "Fa", "C1", "C2", "C3", "C4", "C5"]
    names = ["Fy_1","Fz_1","Fy_2","Fz_2","Fy_3","Fz_3","Fa","C1","C2","C3","C4","C5","CST"]
    sol = {}
    for n, (val, name) in enumerate(zip(sols, names)):
        sol[name] = val
    
    print(sol)

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

        