from interpolation import interpolation
import matplotlib.pyplot as plt
import numpy as np


class Shearcenter:
    def __init__(self, parent):
        self.geometry = parent
        self.interp = interpolation()

        I = self.geometry.MMoI
        self.y_bar, self.z_bar = self.geometry.centroid
        self.z_bar = -1*self.z_bar
        self.Iyy = I[0]
        self.Izz = I[1]
        self.h = self.geometry.h/2
        self.t_skin = self.geometry.t_sk
        self.t_spar = self.geometry.t_sp
        self.y_stringer = self.geometry.crosssection
        self.lsk = self.geometry.l1
        self.n_steps = 100 # Number of steps for the trapezoidal rule
        self.ca = self.geometry.c_a

    def q1(self, V):
        """
        :param V: scalar; virtual shear-force value
        :return: array; q1 at each value of theta
        """
        dtheta = (np.pi/2)/self.n_steps
        theta = np.arange(0, ((np.pi/2) + dtheta), dtheta)
        f = self.t_skin * self.h**2 * np.sin(theta)
        q = self.interp.trapezoidalrule(f, theta)

        theta_compare = np.arctan2((self.h - self.y_stringer[1][1]),(self.h - self.y_stringer[0][1]))
        for i, theta_i in enumerate(theta):
            if theta_i > theta_compare:
                q[i] = q[i] + (self.geometry.Ast * (self.h - self.y_stringer[1][1]))

        q = q * (-V/self.Izz)
        return q

    def q1_z(self, V):
        """
        :param V: scalar; virtual shear-force value
        :return: array; q1 at each value of theta due to the Z component of V
        """
        dtheta = (np.pi / 2) / self.n_steps
        theta = np.arange(0, ((np.pi / 2) + dtheta), dtheta)
        f = self.t_skin * (-self.z_bar - (self.h - self.h*np.cos(theta))) * self.h
        q = self.interp.trapezoidalrule(f, theta)

        theta_compare = np.arctan2((self.h - self.y_stringer[1][1]), (self.h - self.y_stringer[0][1]))
        for i, theta_i in enumerate(theta):
            if theta_i > theta_compare:
                q[i] = q[i] + (self.geometry.Ast * (-self.z_bar - self.y_stringer[0][1]))

        q = q * (-V/self.Iyy)
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

    def q2_z(self, V):
        """
        :param V: scalar; virtual shear-force value
        :return: array; q2 at each value of y due to the z component of V
        """
        dy = self.h/self.n_steps
        y = np.arange(0, (self.h+dy), dy)
        f = self.t_spar * (-self.h - self.z_bar) * np.ones(len(y))
        q = self.interp.trapezoidalrule(f, y)
        q = q * (-V/self.Iyy)
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

        index = 2
        s_boom = (self.y_stringer[1][index])**2 + (self.h - self.y_stringer[0][index])**2

        for i, si in enumerate(s):
            if index == 6:
                break
            if si > s_boom:
                q[i:] = q[i:] + self.geometry.Ast*(self.h - self.y_stringer[1][index])
                index += 1
                s_boom = (self.y_stringer[1][index]) ** 2 + (self.h - self.y_stringer[0][index])**2

        q = q * (-V/self.Izz)
        q = q + self.q1(V)[-1] + self.q2(V)[-1]
        return q

    def q3_z(self, V):
        """
        :param V: scalar; virtual shear-force value
        :return: array; q3 at each value of s due to z component of V
        """
        ds = self.lsk/self.n_steps
        s = np.arange(0, (self.lsk + ds), ds)
        f = self.t_skin * (-self.h - self.z_bar - (self.ca-self.h)/self.lsk * s)
        q = self.interp.trapezoidalrule(f, s)

        index = 2
        s_boom = np.sqrt((self.y_stringer[1][index])**2 + (self.h - self.y_stringer[0][index])**2)

        for i, si in enumerate(s):
            if index == 6:
                break
            if si > s_boom:
                q[i:] = q[i:] + self.geometry.Ast*(-self.z_bar - self.y_stringer[0][index])
                index += 1
                s_boom = np.sqrt((self.y_stringer[1][index]) ** 2 + (self.h - self.y_stringer[0][index])**2)

        for i in range(2, 6):
            q = q + self.geometry.Ast*(-self.z_bar - self.y_stringer[0][i])
        q = q * (-V/self.Iyy)
        q = q + self.q1_z(V)[-1] + self.q2_z(V)[-1]
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

        index = 6
        s_boom = np.sqrt((self.ca - self.y_stringer[0][index])**2 + (self.y_stringer[1][index] - self.h)**2)

        for i, si in enumerate(s):
            if index == 10:
                break
            if si > s_boom:
                q[i:] = q[i:] + self.geometry.Ast*(self.h - self.y_stringer[1][index])
                index += 1
                s_boom = np.sqrt((self.ca - self.y_stringer[0][index])**2 + (self.y_stringer[1][index] - self.h)**2)

        q = q * (-V/self.Izz)
        q = q + self.q3(V)[-1]
        return q

    def q4_z(self, V):
        """
        :param V: scalar; virtual shear-force value
        :return: array; q4 at each value of s due to the z component of V
        """
        ds = self.lsk/self.n_steps
        s = np.arange(0, (self.lsk + ds), ds)
        f = self.t_skin * ((-self.ca - self.z_bar) + ((self.ca - self.h)/self.lsk * s))
        q = self.interp.trapezoidalrule(f, s)

        index = 6
        s_boom = np.sqrt((self.ca - self.y_stringer[0][index])**2 + (self.y_stringer[1][index] - self.h)**2)

        for i, si in enumerate(s):
            if index == 10:
                break
            if si > s_boom:
                q[i:] = q[i:] + self.geometry.Ast*(-self.z_bar - self.y_stringer[0][index])
                index += 1
                s_boom = np.sqrt((self.ca - self.y_stringer[0][index])**2 + (self.y_stringer[1][index] - self.h)**2)

        q = q * (-V/self.Iyy)
        q = q + self.q3_z(V)[-1]
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

    def q5_z(self, V):
        """
        :param V: scalar; virtual shear-force value
        :return: array; q5 at each value of y due to z component of V
        """
        dy = -self.h/self.n_steps
        y = np.arange(0, (-self.h + dy), dy)
        f = self.t_spar * (-self.h - self.z_bar) * np.ones(len(y))
        q = self.interp.trapezoidalrule(f, y)
        q = q * (-V/self.Iyy)
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

        theta_compare = -1*np.arctan2((self.h - self.y_stringer[1][1]),(self.h - self.y_stringer[0][1]))
        for i, theta_i in enumerate(theta):
            if theta_i > theta_compare:
                q[i] = q[i] + (self.geometry.Ast * (self.h - self.y_stringer[1][-1]))

        q = q * (-V/self.Izz)
        q = q + self.q4(V)[-1] - self.q5(V)[-1]
        return q

    def q6_z(self, V):
        """
        :param V: scalar; virtual shear-force value
        :return: array; q2 at each value of theta due to z component of V
        """
        dtheta = (np.pi/2)/self.n_steps
        theta = np.arange((-np.pi/2), (0+dtheta), dtheta)
        f = self.t_skin * (-1*(1-np.cos(theta))*self.h - self.z_bar) * self.h
        q = self.interp.trapezoidalrule(f, theta)

        theta_compare = -1*np.arctan2((self.h - self.y_stringer[1][1]),(self.h - self.y_stringer[0][1]))
        for i, theta_i in enumerate(theta):
            if theta_i > theta_compare:
                q[i] = q[i] + (self.geometry.Ast * (-self.z_bar - self.y_stringer[0][-1]))

        q = q * (-V/self.Iyy)
        q = q + self.q4_z(V)[-1] - self.q5_z(V)[-1]
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

        Rest_I = (self.h/self.t_skin)*self.interp.trapezoidalrule(self.q1(V), theta1)[-1] + (self.h/self.t_skin)*self.interp.trapezoidalrule(self.q6(V), theta6)[-1] - (1/self.t_spar)*self.interp.trapezoidalrule(self.q2(V), y2)[-1] - (1/self.t_spar)*self.interp.trapezoidalrule(self.q5(V), y5)[-1]
        Rest_II = (1/self.t_skin)*self.interp.trapezoidalrule(self.q3(V), s3)[-1] + (1/self.t_skin)*self.interp.trapezoidalrule(self.q4(V), s3)[-1] - (1/self.t_spar)*self.interp.trapezoidalrule(self.q5(V), y5)[-1] + (1/self.t_spar)*self.interp.trapezoidalrule(self.q2(V), y2)[-1]

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

        dz = (1/V) * (T_I + T_II + (self.h * np.cos(beta) * self.interp.trapezoidalrule(self.q3(V), s3)[-1]) + (self.h * np.cos(beta) * self.interp.trapezoidalrule(self.q4(V), s3)[-1]) + (self.h**2 * self.interp.trapezoidalrule(self.q6(V), theta6)[-1]) + (self.h**2 * self.interp.trapezoidalrule(self.q1(V), theta1)[-1]))
        return dz

    def shearstress(self, Vz, Vy):
        """
        :param Vz: z component of shearforce
        :param Vy: y component of shearforce
        :return: shearstress curve at each region --> for regions see verification model figure 3.1
        """
        tau1 = (self.q1(Vy) + self.q1_z(Vz) + qs1)/self.t_skin
        tau2 = (self.q2(Vy) + self.q2_z(Vz) + qs1 + qs2)/self.t_spar
        tau3 = (self.q3(Vy) + self.q3_z(Vz) + qs2)/self.t_skin
        tau4 = (self.q4(Vy) + self.q4_z(Vz) + qs2)/self.t_skin
        tau5 = (self.q5(Vy) + self.q5_z(Vz) + qs1 + qs2)/self.t_spar
        tau6 = (self.q6(Vy) + self.q6_z(Vz) + qs1)/self.t_skin
        return tau1, tau2, tau3, tau4, tau5, tau6


    def shearflow_plot(self, Vz, Vy):
        qs1, qs2 = self.q0_redundant(Vy)

        q1 = self.q1(Vy) + self.q1_z(Vz) + qs1
        q2 = self.q2(Vy) + self.q2_z(Vz) - qs1 + qs2
        q3 = self.q3(Vy) + self.q3_z(Vz) + qs2
        q4 = self.q4(Vy) + self.q4_z(Vz) + qs2
        q5 = self.q5(Vy) + self.q5_z(Vz) - qs1 + qs2
        q6 = self.q6(Vy) + self.q6_z(Vz) + qs1

        ds = self.lsk/self.n_steps
        s = np.arange(0, (self.lsk + ds), ds)

        dtheta = (np.pi/2)/self.n_steps
        theta6 = np.arange((-np.pi/2), (0+dtheta), dtheta)
        theta1 = np.arange(0, ((np.pi / 2) + dtheta), dtheta)

        dy = self.h/self.n_steps
        y2 = np.arange(0, (self.h+dy), dy)
        y5 = np.arange(0, (-self.h - dy), -dy)

        plt.plot(s, q3, label = "region 3")
        plt.plot(s, q4, label = "region 4")
        #plt.plot(y2, q2, label = "region 2")
        #plt.plot(y5, q5, label = "region 5")
        #plt.plot(theta1, q1, label = "region 1")
        #plt.plot(theta6, q6, label = "region 6")
        plt.xlabel("s")
        plt.ylabel("q")
        plt.legend()
        plt.show()
        return 0
