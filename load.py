import numpy as np
from helpers import step

class LoadCase:
    def __init__(self, parent, **kwargs):
        self.geo = parent.geo

        self.d_1 = kwargs.get("d_1")
        self.d_3 = kwargs.get("d_3")
        self.theta = kwargs.get("defl")
        self.P = kwargs.get("load")

        self.a_y = np.sin(self.theta)
        self.a_z = np.cos(self.theta)
        self.a_m = self.geo.h/2 * (np.sin(theta)- np.cos(theta))

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
                  # wip
                  ])
        pass

    def V_z(self, x):
        pass

    def M_y(self, x):
        pass

    def M_z(self, x):
        pass

    def T(self, x):
        pass

    def v_y(self, x):
        pass

    def v_z(self, x):
        pass

    def theta(self, x):
        pass