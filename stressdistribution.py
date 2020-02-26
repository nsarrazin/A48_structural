import numpy as np
import matplotlib.pyplot as plt
from data.consts import parameters_geometry
import geometry

class stressdistribution:
    def __init__(self, **kwargs):
        self.c_a = kwargs.get("c_a")
        self.l_a = kwargs.get("l_a")

        self.x_1 = kwargs.get("x_1")
        self.x_2 = kwargs.get("x_2")
        self.x_3 = kwargs.get("x_3")

        self.x_a = kwargs.get("x_a")

        self.h = kwargs.get("h")
        self.t_sk = kwargs.get("t_sk")
        self.t_sp = kwargs.get("t_sp")

        self.t_st = kwargs.get("t_st")
        self.h_st = kwargs.get("h_st")
        self.w_st = kwargs.get("w_st")
        self.n_st = kwargs.get("n_st")
        pass

    def sigmacalcduetoyforce(self,F,xF):
        c = self.c_a
        b = self.l_a
        Iyy = 4.594369333645184e-05
        Izz = 4.753851442684437e-06

        #direct stress due to force
        sigmay1 = F/(c*b)

        #
        x = b*np.linspace(0,1,41)
        sigmay2 = (F*(xF - x)*x)/Izz

        #
        sigmay = []
        for a in range(len(sigmay2)):
            sigmay.append(sigmay2[a] + sigmay1)



        y = c*np.linspace(0,1,81)
        sigmax = []
        sigmaxmax = []
        for b in range(len(x)):
            sigmax.append((F*(xF - x[b])*y)/Izz)
            sigmaxmax.append(max(abs(sigmax[b])))

        q = max(sigmaxmax)

        
        return sigmaxmax.index(q), x, y, sigmax, sigmay


test = stressdistribution(**parameters_geometry)

x_hinge1 = 0.125
x_hinge2 = 0.498
x_hinge3 = 1.494
x_actuator1 = 0.3755
x_actuator2 = 0.6205

Py = 7589.171334469135
Ry1 = -7589
Ry2 = 0
Ry3 = 0
Rya = 0

index1, x1, y1, sigmax1, sigmay1 = test.sigmacalcduetoyforce(Py,x_actuator2)
index2, x2, y2, sigmax2, sigmay2 = test.sigmacalcduetoyforce(Ry1,x_hinge1)
index3, x3, y3, sigmax3, sigmay3 = test.sigmacalcduetoyforce(Ry2,x_hinge2)
index4, x4, y4, sigmax4, sigmay4 = test.sigmacalcduetoyforce(Ry3,x_hinge3)
index5, x5, y5, sigmax5, sigmay5 = test.sigmacalcduetoyforce(Rya,x_actuator1)

sigmaytotal = []
for x in range(len(x1)):
    sigmaytotal.append(sigmay1[x] + sigmay2[x] + sigmay3[x] + sigmay4[x] + sigmay5[x])

plt.plot(x1,sigmaytotal)
plt.show()

#plt.plot(y,sigmax[index])
#plt.show()