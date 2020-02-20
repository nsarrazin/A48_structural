import numpy as np
from helpers import step

class LoadCase:
    def __init__(self, parent, **kwargs):
        self.geo = parent.geo
        self.interp = parent.interp

        self.d_1 = kwargs.get("d_1")
        self.d_3 = kwargs.get("d_3")
        self.theta = kwargs.get("defl")
        self.P = kwargs.get("load")
        self.E = kwargs.get("e_mod")
        self.G = kwargs.get("g_mod")

        self.a_y = np.sin(self.theta)
        self.a_z = np.cos(self.theta)
        self.a_m = self.geo.h/2 * (np.sin(theta)- np.cos(theta))

        # self.y_sc, self.z_sc = self.geo.shearcenter
        self.y_sc, self.z_sc = 0, 0
    def V_y(self, x):
        np.array([-step(x, self.geo.x_1)**0,
                  0,
                  -step(x, self.geo.x_2)**0,
                  0, 
                  -step(x, self.geo.x_3)**0,
                  0,
                  self.a_y*step(x, self.geo.x_2-self.geo.x_a/2)**0
                  0,
                  0,
                  0,
                  0,
                  0,
                  # int(q(x)) + self.P*self.a_y*step(x,self.geo.x_2-self.geo.x_a/2)**0
                  ])
        pass

    def V_z(self,x):
        np.array([0,
                 -step(x,self.geo.x_1)**0,
                 0,
                 -step(x,self.geo.x_2)**0,
                 0,
                 -step(x,self.geo.x_3)**0,
                 self.a_z*step(x,self.geo.x_2-self.geo.x_a/2)**0,
                 0,
                 0,
                 0,
                 0,
                 0,
                 self.P*self.a_z*step(x,self.geo.x_2+self.geo.x_a/2)**0])
        pass

    def M_y(self, x):
        np.array([0,
                  -step(x,self.geo.x_1)**1,
                  0,
                  -step(x,self.geo.x_2)**1,
                  0,
                  -step(x,self.geo.x_3)**1,
                  self.a_z*step(x,sleg.geo.x_2-self.geo.x_a/2)**1,
                  0,
                  0,
                  0,
                  0,
                  0,
                  self.P*self.a_z*step(x,self.geo.x_2+self.geo.x_a/2)**1])
        pass

    def M_z(self, x):
        np.array([-step(x,self.geo.x_1)**1,
                  0,
                  -step(x,self.geo.x_2)**1,
                  0,
                  -step(x,self.geo.x_3)**1,
                  0,
                  self.a_y*step(x,self.geo.x_2-self.geo.x_a/2)**1,
                  0,
                  0,
                  0,
                  0,
                  0,
                  #int(int(q(x))) + self.P*self.a_y*step(x,self.geo.x_2+self.geo.x_a/2)**1,
                  ])
        pass

    def T(self, x):
        np.array([-self.z_sc*step(x,self.geo.x_1)**0,
                  0,
                  -self.z_sc*step(x,self.geo.x_2)**0,
                  0,
                  -self.z_sc*step(x,self.geo.x_3)**0,
                  0,
                  self.a_y*self.z_sc*step(x,self.geo.x_2-self.ge0.x_a/2)**0-self.a_m*step(x,self.geo.x_2-self.geo.x_a/2)**0,
                  0,
                  0,
                  0,
                  0,
                  0,
                  #int(tau(x))+self.P*self.a_y*self.z_sc*step(x,self.geo.x_2+self.geo.x_a/2)**0-self.P*self.a_m*step(x,self.geo.x_2+self.geo.x_a/2)**0,
                  ])
        pass

    def v_y(self, x):
        np.array([0,
                  -1/6*step(x,self.geo.x_1)**3,
                  0,
                  -1/6*step(x,self.geo.x_2)**3,
                  0,
                  -1/6*step(x,self.geo.x_3)**3,
                  self.a_z/6*step(x,self.geo.x_2-self.geo.x_a/2)**3,
                  0,
                  0,
                  x,
                  1,
                  0,
                  self.P*self.a_z/6*step(x,self.geo.x_2+self.geo.x_a/2)**3])*-1/(self.E*self.geo.MMoI[0])
        pass

    def v_z(self, x):
        np.array([-1/6*step(x,self.geo.x_1)**3,
            	  0,
                  -1/6*step(x,self.geo.x_2)**3,
                  0,
                  -1/6*step(x,self.geo.x_3)**3,
                  0,
                  self.a_y/6*step(x,self.geo.x_2-self.geo.x_a/2)**3,
                  x,
                  1,
                  0,
                  0,
                  0,
                  #int(int(int(int(q(x)))))+self.P*self.a_y/6*step(x,self.geo.x_2+self.geo.x_a/2)**3,
                  ])*-1/(self.E*self.geo.MMoI[1])
        pass

    def theta(self, x):
        np.array([-self.z_sc*step(x,self.geo.x_1)**1,
                  0,
                  -self.z_sc*step(x,self.geo.x_2)**1,
                  0,
                  -self.z_sc*step(x,self.geo.x_3)**1,
                  0,
                  self.a_y*self.z_sc*step(x,self.geo.x_2-self.geo.x_a/2)**1-self.a_m*step(x,self.geo.x_2-self.geo.x_a/2)**1,
                  0,
                  0,
                  0,
                  0,
                  1,
                  #int(int(tau(x)))+self.P*self.a_y*self.z_sc*step(x,self.geo.x_2+self.geo.x_a/2)**1-self.P*self.a_m*step(x,self.geo.x_2+self.geo.x_a/2)**1,
                  ])*1/(self.G*self.J)
        pass