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
    
    # print(sim.case.T(sim.geo.l_a))
    # print(sim.case.B.shape)
    sols = np.linalg.solve(sim.case.A, sim.case.B)

    

    # names = ["Fy_1", "Fz_1", "Fx_2", "Fy_2", "Fz_2", "Fy_3", "Fz_3", "Fa", "C1", "C2", "C3", "C4", "C5"]
    names = ["Fy_1","Fz_1","Fy_2","Fz_2","Fy_3","Fz_3","Fa","C1","C2","C3","C4","C5","CST"]
    for n, (sol, name) in enumerate(zip(sols, names)):
        print(str(n) + " - "+name + " = " + str(sol))

    print(f"Sum of forces in the y-direction : {sols[0]+sols[2]+sols[4]+sols[6]*sim.case.a_y+sim.case.P*sim.case.a_y:.2f}N")
    print(f"Sum of forces in the z-direction : {sols[1]+sols[3]+sols[5]+sols[6]*sim.case.a_z+sim.case.P*sim.case.a_z+sim.interp.integrate_q(sim.geo.l_a,ord=1)[-1]:.2f}N")