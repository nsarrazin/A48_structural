import numpy as np
import matplotlib.pyplot as plt
from data.consts import parameters_geometry
import geometry
import interpolation as inter

class stressdistribution2:
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

test = stressdistribution2(**parameters_geometry)

x_hinge1 = 0.125
x_hinge2 = 0.498
x_hinge3 = 1.494
x_actuator1 = 0.3755
x_actuator2 = 0.6205

Iyy = 4.594369333645184e-05
Izz = 4.753851442684437e-06
Ixx = 1./12.*test.l_a*test.c_a**3 - 1./12.*(test.l_a-2*test.t_sk)*(test.c_a-2*test.t_sk)**3

P = -49200
Py = np.sin(np.radians(30))*P
Pz = np.cos(np.radians(30))*P
Mxp = P*test.h/2.*(np.sin(np.radians(30))- np.cos(np.radians(30)))

Rz1 = -240804.32138928352
Rz2 = 312802.0499279723 
Rz3 = -82448.55420075239  

Ry1 = 31493.90656864783
Ry2 = -49484.785056661785
Ry3 = 18142.991298283927

Ra = 61267.574018492516
Rya = Ra*np.sin(np.radians(30))
Rza = Ra*np.cos(np.radians(30))
Rxa = 0

cz = 0.203625910851571

Mz = Py*abs(test.l_a/2. - x_actuator2)
Mx = -Py*abs(cz - test.h/2.)

x = 0.5
y = test.h*np.linspace(-0.5,0.5,100)
z = np.linspace(-cz,test.c_a-cz,100)

#due to Py

sigmax = Mz/Izz*y
sigmay = Mx/Ixx*z + Mz/Izz*x
sigmaz = Mx/Ixx*y

#due to Pz
My = Pz*abs(test.l_a/2. - x_actuator2)

sigmaz2 = My/Iyy*x
sigmax2 = My/Iyy*z

#due Mxp
sigmay3 = Mxp/Ixx*z
sigmaz3 = Mxp/Ixx*y

#total due to P
sigmaxP = []
for h in range(len(y)):
    sigmaxP.append(sigmax[h] + sigmax2)

sigmayP = sigmay + sigmay3
sigmazP = sigmaz + sigmaz2 + sigmaz3


interp = inter.interpolation()
p = 50

l = test.l_a/p
A = (test.l_a*test.c_a)/(p*p) 

sigmaxq = 0
sigmayq = 0
sigmazq = 0

for zi in z:
    qk, a = interp.bilinear_interpolation(x,zi)
    Fk = qk*A

    Mx = Fk*z
    Mz = Fk*x  

    sigmaz = Mx/Ixx*y
    sigmay = Mx/Ixx*z + Mz/Izz*x
    sigmax = Mz/Izz*y

    sigmazq += sigmaz
    sigmayq += sigmay
    sigmaxq += sigmax

#sigmaxtotal = []
#for h in range(len(y)):
#    sigmaxtotal.append(sigmaxP[h] + sigmaxq[h])


#sigmaytotal = sigmayP + sigmayq
#sigmaztotal = sigmazP + sigmazq



######### due to R1 #################
MzR1 = Ry1*abs(test.l_a/2. - x_hinge1)
MxR1 = -Ry1*abs(cz - test.h/2.)

#due to Ry1

sigmaxR1 = MzR1/Izz*y
sigmayR1 = MxR1/Ixx*z + MzR1/Izz*x
sigmazR1 = MxR1/Ixx*y

#due to Rz1
MyR1 = Rz1*abs(test.l_a/2. - x_hinge1)

sigmaz2R1 = MyR1/Iyy*x
sigmax2R1 = MyR1/Iyy*z

#total due to R1
sigmaxR1total = []
for h in range(len(y)):
    sigmaxR1total.append(sigmaxR1[h] + sigmax2R1)

sigmayR1total = sigmayR1
sigmazR1total = sigmazR1 + sigmaz2R1

######### due to Ra #################
MzRa = Rya*abs(test.l_a/2. - x_actuator1)
MxRa = -Rya*abs(cz - test.h/2.)

#due to Rya

sigmaxRa = MzRa/Izz*y
sigmayRa = MxRa/Ixx*z + MzRa/Izz*x
sigmazRa = MxRa/Ixx*y

#due to Rza
MyRa = Rza*abs(test.l_a/2. - x_actuator1)

sigmaz2Ra = MyRa/Iyy*x
sigmax2Ra = MyRa/Iyy*z

sigmaxRatotal = []
for h in range(len(y)):
    sigmaxRatotal.append(sigmaxRa[h] + sigmax2Ra)

sigmayRatotal = sigmayRa
sigmazRatotal = sigmazRa + sigmaz2Ra

######### due to R2 #################
MzR2 = Ry2*abs(test.l_a/2. - x_hinge2)
MxR2 = -Ry2*abs(cz - test.h/2.)

#due to Ry2

sigmaxR2 = MzR2/Izz*y
sigmayR2 = MxR2/Ixx*z + MzR2/Izz*x
sigmazR2 = MxR2/Ixx*y

#due to Rz2
MyR2 = Rz2*abs(test.l_a/2. - x_hinge2)

sigmaz2R2 = MyR2/Iyy*x
sigmax2R2 = MyR2/Iyy*z

sigmaxR2total = []
for h in range(len(y)):
    sigmaxR2total.append(sigmaxR2[h] + sigmax2R2)

sigmayR2total = sigmayR2
sigmazR2total = sigmazR2 + sigmaz2R2

######### due to R3 #################
MzR3 = Ry3*abs(test.l_a/2. - x_hinge3)
MxR3 = -Ry3*abs(cz - test.h/2.)

#due to Ry3

sigmaxR3 = MzR3/Izz*y
sigmayR3 = MxR3/Ixx*z + MzR3/Izz*x
sigmazR3 = MxR3/Ixx*y

#due to Rz3
MyR3 = Rz3*abs(test.l_a/2. - x_hinge3)

sigmaz2R3 = MyR3/Iyy*x
sigmax2R3 = MyR3/Iyy*z

sigmaxR3total = []
for h in range(len(y)):
    sigmaxR3total.append(sigmaxR3[h] + sigmax2R3)

sigmayR3total = sigmayR3
sigmazR3total = sigmazR3 + sigmaz2R3


####total due to reaction forces
sigmaxtotalRP = sigmaxR3total + sigmaxR2total + sigmaxR1total + sigmaxRatotal + sigmaxP
sigmaytotal = sigmayR3total + sigmayR2total + sigmayR1total + sigmayRatotal + sigmayP + sigmayq
sigmaztotal = sigmazR3total + sigmazR2total + sigmazR1total + sigmazRatotal + sigmazP + sigmazq

sigmaxtotal = []
for h in range(len(y)):
    sigmaxtotal.append(sigmaxtotalRP[h] + sigmaxq[h])

fig = plt.figure()
cs = plt.imshow(sigmaxtotal, extent = (z[0], z[-1], y[0], y[-1]), cmap='RdYlBu')
cbar= plt.colorbar(cs, shrink =0.5)
cbar.ax.set_ylabel("Sigmaxx in crosssection [Pa]")
plt.xlabel("Chord-wise direction [m]")
plt.ylabel("Height-wise direction [m]")
plt.title("Bending stress in x direction (sigma_x)")
plt.tight_layout()
plt.show()

plt.plot(z, sigmaytotal)
plt.show()

plt.plot(y, sigmaztotal)
plt.show()
