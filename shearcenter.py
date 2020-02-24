from interpolation import interpolation
import numpy as np

class Shearcenter:
    def __init__(self, parent):
        self.geometry = parent
        self.interp = interpolation()

        I = self.geometry.MMoI
        self.Iyy = I[0]
        self.Izz = I[1]
        self.h = self.geometry.h/2
        self.t = self.geometry.t_sk
        self.y_stringer = self.geometry.crosssection


    def q1(self, V):
        """
        :return: array; q at each value of theta
        """
        dtheta = (np.pi/2)/100
        theta = np.arange(0, ((np.pi/2) + dtheta), dtheta)
        y = self.t * self.h**2 * np.sin(theta)
        q = self.interp.trapezoidalrule(y, theta)
        q = q + (self.geometry.Ast * self.y_stringer[1][1])
        q = q * (-V/self.Izz)
        return q

    #
    # def q6(self, theta):
    #     return q


## Test implementation
# test = Shearcenter()
# print(test.q1(1))