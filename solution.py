import numpy as np
import matplotlib.pyplot as plt
from helpers import step

class Solution:
    def __init__(self, parent):
        self.parent = parent
        self.case = parent.case
        self.geo = parent.geo
        self._sol = {}
        self.labels=["Hinge 1", "Hinge 2" , "Hinge 3", "Actuator I", "Actuator II"]

    @property
    def sol(self):
        if self._sol == {}:
            names = ["Fy'_1","Fz'_1","Fy'_2","Fz'_2","Fy'_3","Fz'_3","Fa","C1","C2","C3","C4","C5","CST"]
            self._sol = {}
            for n, (val, name) in enumerate(zip(self.parent.x, names)):
                self._sol[name] = val
    
        return self._sol

    def V_y_prime(self,x):
        sol = self.sol
        return -1*(sol["Fy'_1"]*step(x,self.geo.x_1,power=0)\
                  +sol["Fy'_2"]*step(x,self.geo.x_2,power=0)\
                  +sol["Fy'_3"]*step(x,self.geo.x_3,power=0)\
                  +sol["Fa"]*self.case.a_y*step(x,self.case.x_I,power=0)\
                  -self.case.P*self.case.a_y*step(x,self.case.x_II,power=0)\
                  +self.case.interp.integrate_q(x,ord=1)[-1])

    def V_z_prime(self,x):
        sol = self.sol
        return -1*(sol["Fz'_1"]*step(x,self.geo.x_1,power=0)\
                  +sol["Fz'_2"]*step(x,self.geo.x_2,power=0)\
                  +sol["Fz'_3"]*step(x,self.geo.x_3,power=0)\
                  +sol["Fa"]*self.case.a_z*step(x,self.case.x_I,power=0)\
                  -self.case.P*self.case.a_z*step(x,self.case.x_II,power=0))

    def M_y_prime(self,x):
        sol = self.sol
        return -sol["Fz'_1"]*step(x,self.geo.x_1,power=1)\
               -sol["Fz'_2"]*step(x,self.geo.x_2,power=1)\
               -sol["Fz'_3"]*step(x,self.geo.x_3,power=1)\
               -sol["Fa"]*self.case.a_z*step(x,self.case.x_I,power=1)\
               +self.case.P*self.case.a_z*step(x,self.case.x_II,power=1)

    def M_z_prime(self,x):
        sol = self.sol
        return -sol["Fy'_1"]*step(x,self.geo.x_1,power=1)\
               -sol["Fy'_2"]*step(x,self.geo.x_2,power=1)\
               -sol["Fy'_3"]*step(x,self.geo.x_3,power=1)\
               -sol["Fa"]*self.case.a_y*step(x,self.case.x_I,power=1)\
               +self.case.P*self.case.a_y*step(x,self.case.x_II,power=1)\
               -self.case.interp.integrate_q(x,ord=2)[-1]

    def T(self,x):
        sol = self.sol
        return (sol["Fy'_1"]*self.case.z_sc*step(x,self.geo.x_1,power=0)\
                  +sol["Fy'_2"]*self.case.z_sc*step(x,self.geo.x_2,power=0)\
                  +sol["Fy'_3"]*self.case.z_sc*step(x,self.geo.x_3,power=0)\
                  +sol["Fa"]*self.case.a_y*self.case.z_sc*step(x,self.case.x_I,power=0)\
                  +sol["Fa"]*self.case.a_m*step(x,self.case.x_I,power=0)\
                  -self.case.P*self.case.a_y*self.case.z_sc*step(x,self.case.x_II,power=0)\
                  -self.case.P*self.case.a_m*step(x,self.case.x_II,power=0)\
                  -self.case.interp.integrate_tau(x,self.case.z_sc,ord=1)[-1])

    def v_y_prime(self, x):
        sol = self.sol
        return (-sol["Fy'_1"]/6*step(x,self.geo.x_1,power=3)-sol["Fy'_2"]/6*step(x,self.geo.x_2,power=3)-sol["Fy'_3"]/6*step(x,self.geo.x_3,power=3)\
                -(sol["Fa"]*self.case.a_y)/6*step(x,self.case.x_I,power=3)\
                +(self.case.P*self.case.a_y)/6*step(x,self.case.x_II,power=3)\
                -self.case.interp.integrate_q(x,ord=4)[-1])\
                *-1/(self.case.E*self.geo.MMoI[1])\
                +sol["C3"]*x+sol["C4"]       

    def v_z_prime(self, x):
        sol = self.sol
        return (-sol["Fz'_1"]/6*step(x,self.geo.x_1,power=3)-sol["Fz'_2"]/6*step(x,self.geo.x_2,power=3)-sol["Fz'_3"]/6*step(x,self.geo.x_3,power=3)\
                -(sol["Fa"]*self.case.a_z)/6*step(x,self.case.x_I,power=3)\
                +(self.case.P*self.case.a_z)/6*step(x,self.case.x_II,power=3))*-1/(self.case.E*self.geo.MMoI[0])\
                +sol["C1"]*x+sol["C2"]


    def theta(self, x):
        sol = self.sol
        return ((sol["Fy'_1"]*self.case.z_sc*step(x,self.geo.x_1,power=1)\
                   +sol["Fy'_2"]*self.case.z_sc*step(x,self.geo.x_2,power=1)\
                   +sol["Fy'_3"]*self.case.z_sc*step(x,self.geo.x_3,power=1)\
                   +sol["Fa"]*self.case.a_y*self.case.z_sc*step(x,self.case.x_I,power=1)\
                   +sol["Fa"]*self.case.a_m*step(x,self.case.x_I,power=1)\
                   -self.case.P*self.case.a_y*self.case.z_sc*step(x,self.case.x_II,power=1)\
                   -self.case.P*self.case.a_m*step(x,self.case.x_II,power=1)\
                   -self.case.interp.integrate_tau(x,z_sc=self.case.z_sc,ord=2)[-1])\
                   *1/(self.case.G*self.geo.J)\
                   +sol["C5"])
    
    def v_y(self,x):
        return self.v_y_prime(x)*np.cos(self.case.defl)-self.v_z_prime(x)*np.sin(self.case.defl)
    
    def v_z(self,x): 
        return self.v_y_prime(x)*np.sin(self.case.defl)+self.v_z_prime(x)*np.cos(self.case.defl)

    def tau(self,x,z_sc):
        return self.case.interp.tau(x,self.case.z_sc)

    def slope_y_prime(self,x):
        sol = self.sol
        return (-sol["Fy'_1"]/2*step(x,self.geo.x_1,power=2)\
                -sol["Fy'_2"]/2*step(x,self.geo.x_2,power=2)\
                -sol["Fy'_3"]/2*step(x,self.geo.x_3,power=2)\
                -(sol["Fa"]*self.case.a_y)/2*step(x,self.case.x_I,power=2)\
                +(self.case.P*self.case.a_y)/2*step(x,self.case.x_II,power=2)\
                -self.case.interp.integrate_q(x,ord=3)[-1])\
                *-1/(self.case.E*self.geo.MMoI[1])\
                +sol["C3"]
    
    def slope_z_prime(self,x):
        sol = self.sol
        return (-sol["Fz'_1"]/2*step(x,self.geo.x_1,power=2)\
                -sol["Fz'_2"]/2*step(x,self.geo.x_2,power=2)\
                -sol["Fz'_3"]/2*step(x,self.geo.x_3,power=2)\
                -(sol["Fa"]*self.case.a_z)/2*step(x,self.case.x_I,power=2)\
                +(self.case.P*self.case.a_z)/2*step(x,self.case.x_II,power=2))\
                *-1/(self.case.E*self.geo.MMoI[0])\
                +sol["C1"]


    def rotate(self, x, y):
        x_prime = x * np.cos(self.case.defl) - y * np.sin(self.case.defl)
        y_prime = x * np.sin(self.case.defl) + y * np.cos(self.case.defl)

        return x_prime, y_prime
    
    def plot_defl(self):
        xs = np.linspace(0, self.geo.l_a, 100)

        ys_1 = []
        ys_2 = []

        for x in xs:
            ys_1.append(self.v_y_prime(x))
            ys_2.append(self.v_z_prime(x))

        plt.plot(xs, ys_1, label="v_y'")
        plt.plot(xs, ys_2, label="v_z'")
        
        for n,x in enumerate([self.geo.x_1, self.geo.x_2, self.geo.x_3, self.case.x_I, self.case.x_II]):
            plt.axvline(x=x, linestyle="dashed", linewidth=1.2, color=f"C{n}", label=self.labels[n])

        plt.xlabel("Span-wise length [m]")
        plt.ylabel("Deflection [m]")
        plt.legend()
        plt.tight_layout()
        plt.show()

    def plot_twist(self):
        xs = np.linspace(0, self.geo.l_a, 100)

        ys = []

        for x in xs:
            ys.append(self.theta(x))

        plt.plot(xs, ys, label="Twist")

        for n,x in enumerate([self.geo.x_1, self.geo.x_2, self.geo.x_3, self.case.x_I, self.case.x_II]):
            plt.axvline(x=x, linestyle="dashed", linewidth=1.2, color=f"C{n}", label=self.labels[n])
        
        plt.xlabel("Span-wise length [m]")
        plt.ylabel("Twist [rad]")
        plt.legend()

        plt.tight_layout()
        plt.show()


    def plot_shear(self):
        xs = np.linspace(0, self.geo.l_a, 100)

        ys_1 = []
        ys_2 = []

        for x in xs:
            ys_1.append(self.V_y_prime(x))
            ys_2.append(self.V_z_prime(x))

        plt.plot(xs, ys_1, label="Vy'")
        plt.plot(xs, ys_2, label="Vz'")
        
        for n,x in enumerate([self.geo.x_1, self.geo.x_2, self.geo.x_3, self.case.x_I, self.case.x_II]):
            plt.axvline(x=x, linestyle="dashed", linewidth=1.2, color=f"C{n}", label=self.labels[n])

        plt.xlabel("Span-wise length [m]")
        plt.ylabel("Shear [N]")
        plt.legend()
        plt.tight_layout()
        plt.show()

    def plot_moment(self):
        xs = np.linspace(0, self.geo.l_a, 100)

        ys_1 = []
        ys_2 = []

        for x in xs:
            ys_1.append(self.M_y_prime(x))
            ys_2.append(self.M_z_prime(x))

        plt.plot(xs, ys_1, label="My'")
        plt.plot(xs, ys_2, label="Mz'")

        for n,x in enumerate([self.geo.x_1, self.geo.x_2, self.geo.x_3, self.case.x_I, self.case.x_II]):
            plt.axvline(x=x, linestyle="dashed", linewidth=1.2, color=f"C{n}", label=self.labels[n])

        plt.xlabel("Span-wise length [m]")
        plt.ylabel("Bending moment [N.m]")
        plt.legend()
        plt.tight_layout()
        plt.show()

    def plot_torque(self):
        xs = np.linspace(0, self.geo.l_a, 100)

        ys_1 = []

        for x in xs:
            ys_1.append(self.T(x))

        plt.plot(xs, ys_1, label="Torque")
        
        for n,x in enumerate([self.geo.x_1, self.geo.x_2, self.geo.x_3, self.case.x_I, self.case.x_II]):
            plt.axvline(x=x, linestyle="dashed", linewidth=1.2, color=f"C{n}", label=self.labels[n])

        plt.xlabel("Span-wise length [m]")
        plt.ylabel("Torque [N.m]")
        plt.legend()
        plt.tight_layout()
        plt.show()

    def plot_tau(self):
        xs = np.linspace(0, self.geo.l_a,100)
        ys = []
        for x in xs:
            ys.append(self.tau(x,self.case.z_sc))
        plt.plot(xs,ys)

        for n,x in enumerate([self.geo.x_1, self.geo.x_2, self.geo.x_3, self.case.x_I, self.case.x_II]):
            plt.axvline(x=x, linestyle="dashed", linewidth=1.2, color=f"C{n}")

        plt.tight_layout()
        plt.show()

    def plot_slope(self):
        xs = np.linspace(0, self.geo.l_a, 100)

        ys_1 = []
        ys_2 = []

        for x in xs:
            ys_1.append(self.slope_y_prime(x))
            ys_2.append(self.slope_z_prime(x))

        plt.plot(xs, ys_1, label=r"$\frac{d_{v_{y'}}}{d_x}$")
        plt.plot(xs, ys_2, label=r"$\frac{d_{v_{z'}}}{d_x}$")

        for n,x in enumerate([self.geo.x_1, self.geo.x_2, self.geo.x_3, self.case.x_I, self.case.x_II]):
            plt.axvline(x=x, linestyle="dashed", linewidth=1.2, color=f"C{n}", label=self.labels[n])

        plt.xlabel("Span-wise length [m]")
        plt.ylabel("Deflection [rad]")
        plt.legend()
        plt.tight_layout()
        plt.show()
    
        
    def plot_solution(self):
        xs = np.linspace(0, self.geo.l_a, 100)

        ys_1 = []

        for x in xs:
            ys_1.append(self.theta(x)*self.case.z_sc+self.v_y(x))

        plt.plot(xs, ys_1, label="v_y'(x)+theta(x)*z_sc")

        bc = [(self.geo.x_1, self.case.d_1),
              (self.geo.x_2, 0),
              (self.geo.x_3, self.case.d_3)]
        
        for n,x in enumerate([self.geo.x_1, self.geo.x_2, self.geo.x_3, self.case.x_I, self.case.x_II]):
            plt.axvline(x=x, linestyle="dashed", linewidth=1.2, color=f"C{n}", label=self.labels[n])

        for n,point in enumerate(bc):
            plt.scatter(point[0], point[1], c=f"C{n}",s=40, label=f"BC @ Hinge {n}")
        
        plt.xlabel("Span-wise length [m]")
        plt.ylabel("Deflection [m]")
        plt.legend()
        plt.tight_layout()
        plt.show()
