from interpolation import interpolation
from geometry import Geometry
import numpy as np
interp = interpolation()
geo = Geometry()

class Shearcenter:
    def __init__(self):
        Iyy2, Izz2 = geo.MMoI()
        self.Iyy = Iyy2
        self.Izz = Izz2
        self.h = geo.h/2
        self.t = geo.t_sk
        self.y_stringer = geo.crosssection()


    def q1(self, V):
        """
        :return: array; q at each value of theta
        """
        dtheta = (np.pi/2)/100
        theta = np.arange(0, ((np.pi/2) + dtheta), dtheta )
        y = self.t * self.h**2 * np.sin(theta)
        q = interpolation.trapezoidalrule(y, theta)
        q = q + (geo.Ast() * self.y_stringer[1][1])
        q = q * (-V/self.Izz)
        return q


    def q6(self, theta):
        return q


## Test implementation
test = Shearcenter()
print(test.q1(1))