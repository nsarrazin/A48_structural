import numpy as np
import matplotlib.pyplot as plt
from helpers import step


class Solution:
    def __init__(self, parent):
        self.parent = parent
        self.case = parent.case
        self.geo = parent.geo
        self._sol = {}

    @property
    def sol(self):
        if self._sol == {}:
            names = ["Fy_1","Fz_1","Fy_2","Fz_2","Fy_3","Fz_3","Fa","C1","C2","C3","C4","C5","CST"]
            self._sol = {}
            for n, (val, name) in enumerate(zip(self.parent.x, names)):
                self._sol[name] = val
    
        return self._sol


    def v_y(self, x):
        sol = self.sol
        return -1/(self.case.E*self.geo.MMoI[0])*(sol["Fa"]*self.case.a_z/6*step(x,self.geo.x_2-self.geo.x_a/2)**3\
            +self.case.P*self.case.a_z*step(x,self.geo.x_2+self.geo.x_a/2)**3\
                -sol["Fz_1"]/6*step(x,self.geo.x_1)**3-sol["Fz_2"]/6*step(x,self.geo.x_2)**3-sol["Fz_3"]/6*step(x,self.geo.x_3)**3\
                    +sol["C3"]*x+sol["C4"]*0)       


    def plot(self):

        xs = np.linspace(0, self.geo.l_a, 100)
        ys = [self.v_y(i) for i in xs]

        offset = self.v_y(self.geo.x_2)*0

        plt.plot(xs, ys-offset)
        
        for n,x in enumerate([self.geo.x_1, self.geo.x_2, self.geo.x_3]):
            plt.axvline(x=x, linestyle="dashed", linewidth=1.2, color=f"C{n}")

        for n,y in enumerate([self.case.d_1, 0, self.case.d_3]):
            plt.axhline(y=y, linestyle="dashed", linewidth=1.2, color=f"C{n}")

        plt.show()

            