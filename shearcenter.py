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
        self.t_skin = self.geometry.t_sk
        self.t_spar = self.geometry.t_sp
        self.y_stringer = self.geometry.crosssection
        self.lsk = self.geometry.l1
        self.n_steps = 100 # Number of steps for the trapezoidal rule

    def q1(self, V):
        """
        :param V: scalar; virtual shear-force value
        :return: array; q1 at each value of theta
        """
        dtheta = (np.pi/2)/self.n_steps
        theta = np.arange(0, ((np.pi/2) + dtheta), dtheta)
        f = self.t_skin * self.h**2 * np.sin(theta)
        q = self.interp.trapezoidalrule(f, theta)
        q = q + (self.geometry.Ast * (self.h - self.y_stringer[1][1]))
        q = q * (-V/self.Izz)
        return q

    def q2(self, V):
        """
        :param V: scalar; virtual shear-force value
        :return: array; q2 at each value of y
        """
        dy = self.h/self.n_steps
        y = np.arange(0, (self.h+dy), dy)
        f = self.t_spar*y
        q = self.interp.trapezoidalrule(f, y)
        q = q * (-V/self.Izz)
        return q

    def q3(self, V):
        """
        :param V: scalar; virtual shear-force value
        :return: array; q3 at each value of s
        """
        ds = self.lsk/self.n_steps
        s = np.arange(0, (self.lsk + ds), ds)
        f = self.t_skin*(self.h - (self.h/self.lsk)*s)
        q = self.interp.trapezoidalrule(f, s)
        for i in range(2, 6):
            q = q + self.geometry.Ast*(self.h - self.y_stringer[1][i])
        q = q * (-V/self.Izz)
        q = q + self.q1(V)[-1] + self.q2(V)[-1]
        return q

    def q4(self, V):
        """
        :param V: scalar; virtual shear-force value
        :return: array; q4 at each value of s
        """
        ds = self.lsk/self.n_steps
        s = np.arange(0, (self.lsk + ds), ds)
        f = self.t_skin * ((-self.h/self.lsk)*s)
        q = self.interp.trapezoidalrule(f, s)
        for i in range(6, 10):
            q = q + self.geometry.Ast*(self.h - self.y_stringer[1][i])
        q = q * (-V/self.Izz)
        q = q + self.q3(V)[-1]
        return q

    def q5(self, V):
        """
        :param V: scalar; virtual shear-force value
        :return: array; q5 at each value of y
        """
        dy = -self.h/self.n_steps
        y = np.arange(0, (-self.h + dy), dy)
        f = self.t_spar*y
        q = self.interp.trapezoidalrule(f, y)
        q = q * (-V/self.Izz)
        return q

    def q6(self, V):
        """
        :param V: scalar; virtual shear-force value
        :return: array; q2 at each value of theta
        """
        dtheta = (np.pi/2)/self.n_steps
        theta = np.arange((-np.pi/2), (0+dtheta), dtheta)
        f = self.t_skin * self.h**2 * np.sin(theta)
        q = self.interp.trapezoidalrule(f, theta)
        q = q + (self.geometry.Ast * (self.h - self.y_stringer[1][-1]))
        q = q * (-V/self.Izz)
        q = q + self.q4(1)[-1] - self.q5(1)[-1]
        return q

    def q0_redundant(self, V):
        ds3 = self.lsk/self.n_steps
        s3 = np.arange(0, (self.lsk + ds3), ds3)

        dy = self.h/self.n_steps
        y2 = np.arange(0, (self.h+dy), dy)
        y5 = np.arange(0, (-self.h - dy), -dy)

        dtheta = (np.pi/2)/self.n_steps
        theta6 = np.arange((-np.pi/2), (0+dtheta), dtheta)
        theta1 = np.arange(0, ((np.pi / 2) + dtheta), dtheta)

        Rest_I = (self.h/self.t_skin)*self.interp.trapezoidalrule(self.q1(1), theta1)[-1] + (self.h/self.t_skin)*self.interp.trapezoidalrule(self.q6(1), theta6)[-1] - (1/self.t_spar)*self.interp.trapezoidalrule(self.q2(1), y2)[-1] - (1/self.t_spar)*self.interp.trapezoidalrule(self.q5(1), y5)[-1]
        Rest_II = (1/self.t_skin)*self.interp.trapezoidalrule(self.q3(1), s3)[-1] + (1/self.t_skin)*self.interp.trapezoidalrule(self.q4(1), s3)[-1] - (1/self.t_spar)*self.interp.trapezoidalrule(self.q5(1), y5)[-1] + (1/self.t_spar)*self.interp.trapezoidalrule(self.q2(1), y2)[-1]

        A = (np.pi*self.h)/(self.t_skin) + (2*self.h)/(self.t_spar)
        B = (-self.h)/(self.t_spar) * 2
        C = (2*self.lsk)/(self.t_skin) + (2*self.h)/(self.t_spar)

        q01 = ((B/C)*Rest_II - Rest_I)/(A - (B**2/C))
        q02 = (-Rest_II - B*q01)/C
        return q01, q02

    def shearcenter(self, V):
        """
        :param V: scalar; virtual shear-force value
        :return: scalar; distance of the shear center from the spar (positive means to the right)
        """
        T_I = np.pi * self.h**2 * self.q0_redundant(1)[0]
        T_II = (2 * self.h) * (self.geometry.c_a - self.h) * self.q0_redundant(1)[1]

        beta = np.arctan2(self.h, (self.geometry.c_a - self.h))

        ds3 = self.lsk/self.n_steps
        s3 = np.arange(0, (self.lsk + ds3), ds3)

        dtheta = (np.pi/2)/self.n_steps
        theta6 = np.arange((-np.pi/2), (0+dtheta), dtheta)
        theta1 = np.arange(0, ((np.pi / 2) + dtheta), dtheta)

        dz = (1/V) * (T_I + T_II + (self.h * np.cos(beta) * self.interp.trapezoidalrule(self.q3(V), s3)[-1]) + (self.h * np.cos(beta) * self.interp.trapezoidalrule(self.q4(V), s3)[-1]) + (self.h**2 * self.interp.trapezoidalrule(self.q6(1), theta6)[-1]) + (self.h**2 * self.interp.trapezoidalrule(self.q1(1), theta1)[-1]))
        return dz