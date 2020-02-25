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

    def V_y_prime(self,x):
        sol = self.sol
        return sol["Fy_1"]*step(x,self.geo.x_1)**0+sol["Fy_2"]*step(x,self.geo.x_2)**0+sol["Fy_3"]*step(x,self.geo.x_3)**0\
            +sol["Fa"]*self.case.a_y*step(x,self.case.x_I)**0\
                +self.case.P*self.case.a_y*step(x,self.case.x_II)**0+self.case.interp.integrate_q(x,ord=1)[-1]

    def V_z_prime(self,x):
        sol = self.sol
        return sol["Fz_1"]*step(x,self.geo.x_1)**0+sol["Fz_2"]*step(x,self.geo.x_2)**0+sol["Fz_3"]*step(x,self.geo.x_3)**0\
            +sol["Fa"]*self.case.a_z*step(x,self.case.x_I)**0+self.case.P*self.case.a_z*step(x,self.case.x_II)**0

    def M_y_prime(self,x):
        sol = self.sol
        return -sol["Fz_1"]*step(x,self.geo.x_1)**1-sol["Fz_2"]*step(x,self.geo.x_2)**1-sol["Fz_3"]*step(x,self.geo.x_3)**1\
            -sol["Fa"]*self.case.a_z*step(x,self.case.x_I)**1-self.case.P*self.case.a_z*step(x,self.case.x_II)**1

    def M_z_prime(self,x):
        sol = self.sol
        return -sol["Fy_1"]*step(x,self.geo.x_1)**1-sol["Fy_2"]*step(x,self.geo.x_2)**1-sol["Fy_3"]*step(x,self.geo.x_3)**1\
            -sol["Fa"]*self.case.a_y*step(x,self.case.x_I)**1\
                -self.case.P*self.case.a_y*step(x,self.case.x_II)**1+self.case.interp.integrate_q(x,ord=2)[-1]

    def T(self,x):
        sol = self.sol
        return sol["Fy_1"]*self.case.z_sc*step(x,self.geo.x_1)**0+sol["Fy_2"]*self.case.z_sc*step(x,self.geo.x_2)**0+sol["Fy_3"]*self.case.z_sc*step(x,self.geo.x_3)**0\
            +sol["Fa"]*self.case.a_y*self.case.z_sc*step(x,self.case.x_I)**0+sol["Fa"]*self.case.a_m*step(x,self.case.x_I)**0\
                +self.case.P*self.case.a_y*self.case.z_sc*step(x,self.case.x_II)**0+self.case.P*self.case.a_m*step(x,self.case.x_II)**0\
                    +self.case.interp.integrate_tau(x,self.case.z_sc,ord=1)[-1]

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
        return (sol["Fy_1"]*self.case.z_sc*step(x,self.geo.x_1)**1+sol["Fy_2"]*self.case.z_sc*step(x,self.geo.x_2)**1+sol["Fy_3"]*self.case.z_sc*step(x,self.geo.x_3)**1\
            +sol["Fa"]*self.case.a_y*self.case.z_sc*step(x,self.case.x_I)**1+sol["Fa"]*self.case.a_m*step(x,self.case.x_I)**1\
                +self.case.P*self.case.a_y*self.case.z_sc*step(x,self.case.x_II)**1+self.case.P*self.case.a_m*step(x,self.case.x_II)**1\
                    +self.case.interp.integrate_tau(x,z_sc=self.case.z_sc,ord=2)[-1])*1/(self.case.G*self.geo.J)\
                        +sol["C5"]
    
    def v_y(self,x):
        return self.v_y_prime(x)*np.cos(self.case.defl)-self.v_z_prime(x)*np.sin(self.case.defl)
    
    def v_z(self,x): 
        return +self.v_y(x)*np.sin(self.case.defl)+self.v_z_prime(x)*np.cos(self.case.defl)
    
    def plot_1(self):
        xs = np.linspace(0, self.geo.l_a, 100)

        ys = []

        for x in xs:
            ys.append(self.v_y(x)+self.theta(x)*self.case.z_sc)

        plt.plot(xs, ys)
        
        for n,x in enumerate([self.geo.x_1, self.geo.x_2, self.geo.x_3]):
            plt.axvline(x=x, linestyle="dashed", linewidth=1.2, color=f"C{n}")

        for n,y in enumerate([self.case.d_1, 0, self.case.d_3]):
            plt.axhline(y=y, linestyle="dashed", linewidth=1.2, color=f"C{n}")

        plt.show()


    def plot_2(self):
        xs = np.linspace(0, self.geo.l_a, 100)

        ys_1 = []
        ys_2 = []

        for x in xs:
            ys_1.append(self.V_y_prime(x))
            ys_2.append(self.V_z_prime(x))

        plt.plot(xs, ys_1)
        plt.plot(xs, ys_2)

        for n,x in enumerate([self.geo.x_1, self.geo.x_2, self.geo.x_3, self.case.x_I, self.case.x_II]):
            plt.axvline(x=x, linestyle="dashed", linewidth=1.2, color=f"C{n}")

        # for n,y in enumerate([self.case.d_1, 0, self.case.d_3]):
            # plt.axhline(y=y, linestyle="dashed", linewidth=1.2, color=f"C{n}")

        plt.show()
