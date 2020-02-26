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

        self.defl = kwargs.get('max_defl')
        pass

    def sigmaduetoyforce(self,F,xF): 
        c = self.c_a
        b = self.l_a
        Iyy = 4.594369333645184e-05
        Izz = 4.753851442684437e-06
        
        #spanwise locations
        x = b*np.linspace(0,1,10)

        #moment z varying spanwise
        Mz = F*(xF - x)

        #direct stress in y direction due to force (constant throughout span and chord)
        sigmay1 = np.matrix(len(x)*[F/(c*b)])

        #bending stress in y direction due to force (varying spanwise)
        sigmay2 = np.matrix(Mz/Izz*x)

        sigmaytotal = sigmay1 + sigmay2

        #heightwise locations
        y = c*np.linspace(0,1,10)

        #bending stress in x direction due to force (varying in height and spanwise)
        sigmax = []
        for a in range(len(x)):
            sigmax.append(Mz[a]/Izz*y)

        sigmax = np.matrix(sigmax)

        return x,y,sigmaytotal, sigmax
    
    def sigmaduetoxmoment(self, Mx):
        c = self.c_a
        h = self.h
        Ixx = 1./12.*self.l_a*self.c_a**3

        #chord locations
        z = c*np.linspace(0,1,10)

        #sigmay due to moment in x varies chordwise
        sigmay = Mx/Ixx*z

        y = h*np.linspace(0,1,10)
        #sigmaz due to moment in x varies heightwise
        sigmaz = Mx/Ixx*y

        return z, y, sigmay, sigmaz

test = stressdistribution(**parameters_geometry)

x_hinge1 = 0.125
x_hinge2 = 0.498
x_hinge3 = 1.494
x_actuator1 = 0.3755
x_actuator2 = 0.6205

Py = 7589.171334469135
Ry1 = 5 * 10**4
Ry2 = -7.8 * 10**4
Ry3 = 2.7 * 10**4
Rya = -3.6 * 10**4
Mxa = Rya*test.h/2.*(np.sin(np.radians(30))- np.cos(np.radians(30)))



x1, y1, sigmaytotal1, sigmax1 = test.sigmaduetoyforce(Py, x_actuator2)
x2, y2, sigmaytotal2, sigmax2 = test.sigmaduetoyforce(Ry1, x_hinge1)
x3, y3, sigmaytotal3, sigmax3 = test.sigmaduetoyforce(Ry2, x_hinge2)
x4, y4, sigmaytotal4, sigmax4 = test.sigmaduetoyforce(Ry3, x_hinge3)
x5, y5, sigmaytotal5, sigmax5 = test.sigmaduetoyforce(Rya, x_actuator1)

z1, y1, sigmaytotal6, sigmax6 = test.sigmaduetoxmoment(Mxa)

sigmaytotalarray = np.array(sigmaytotal1 + sigmaytotal2 + sigmaytotal3 + sigmaytotal4 + sigmaytotal5)
sigmaytotal = sigmaytotalarray[0]

sigmaxtotal = np.array(sigmax1 + sigmax2 + sigmax3 + sigmax4 + sigmax5)



fig = plt.figure()
cs = plt.imshow(sigmaxtotal, extent = (0., 1.611, 0., 0.161), cmap='Blues')
cbar= plt.colorbar(cs, shrink =0.5)
cbar.ax.set_ylabel("Bending in x direction [Pa]")
plt.xlabel("Span-wise direction [m]")
plt.ylabel("Height direction [m]")
plt.title("Bending stress over the aileron")
plt.tight_layout()
plt.show()