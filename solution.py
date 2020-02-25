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


    def v_y_prime(self, x):
        sol = self.sol
        return (-sol["Fy_1"]/6*step(x,self.geo.x_1)**3-sol["Fy_2"]/6*step(x,self.geo.x_2)**3-sol["Fy_3"]/6*step(x,self.geo.x_3)**3\
            -(sol["Fa"]*self.case.a_y)/6*step(x,self.case.x_I)**3\
                -(self.case.P*self.case.a_y)/6*step(x,self.case.x_II)**3-self.case.interp.integrate_q(x,ord=4)[-1])*-1/(self.case.E*self.geo.MMoI[1])\
                    +sol["C3"]*x+sol["C4"]        

    def v_z_prime(self, x):
        sol = self.sol
        return (-sol["Fz_1"]/6*step(x,self.geo.x_1)**3-sol["Fz_2"]/6*step(x,self.geo.x_2)**3-sol["Fz_3"]/6*step(x,self.geo.x_3)**3\
            -(sol["Fa"]*self.case.a_z)/6*step(x,self.case.x_I)**3\
                -(self.case.P*self.case.a_z)/6*step(x,self.case.x_II))*-1/(self.case.E*self.geo.MMoI[0])\
                    +sol["C1"]*x+sol["C2"]


    def theta(self, x):
        sol = self.sol
        return sol["Fy_1"]*self.case.z_sc*step(x,self.geo.x_1)**1+sol["Fy_2"]*self.case.z_sc*step(x,self.geo.x_2)**1+sol["Fy_3"]*self.case.z_sc*step(x,self.geo.x_3)**1\
            +sol["Fa"]*self.case.a_y*self.case.z_sc*step(x,self.case.x_I)**1+sol["Fa"]*self.case.a_m*step(x,self.case.x_I)**1\
                +self.case.P*self.case.a_y*self.case.z_sc*step(x,self.case.x_II)**1+self.case.P*self.case.a_m*step(x,self.case.x_II)**1\
                    +self.case.interp.integrate_tau(x,z_sc=self.case.z_sc,ord=2)[-1])*1/(self.case.G*self.geo.J)\
                        +sol["C5"]
    
    def v_y(self,x):
        return self.v_y_prime(x)*np.cos(self.case.defl)-self.v_z_prime(x)*np.sin(self.case.defl)
    
    def v_z(self,x): 
        return +self.v_y(x)*np.sin(self.case.defl)+self.v_z_prime(x)*np.cos(self.case.defl)
    
    def plot(self):

        xs = np.linspace(0, self.geo.l_a, 100)
        ys = [self.v_y(i)+self.theta(i)*self.case.z_sc for i in xs]

        print(self.v_y(self.geo.x_1)+self.theta(self.geo.x_1)*self.case.z_sc)

        offset = self.v_y(self.geo.x_1)+self.theta(self.geo.x_1)*self.case.z_sc*0

        plt.plot(xs, ys-offset)
        
        for n,x in enumerate([self.geo.x_1, self.geo.x_2, self.geo.x_3]):
            plt.axvline(x=x, linestyle="dashed", linewidth=1.2, color=f"C{n}")

        for n,y in enumerate([self.case.d_1, 0, self.case.d_3]):
            plt.axhline(y=y, linestyle="dashed", linewidth=1.2, color=f"C{n}")

        plt.show()