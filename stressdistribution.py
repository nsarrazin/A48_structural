import numpy as np
import matplotlib.pyplot as plt
from data.consts import parameters_geometry
import geometry
import interpolation as inter

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
        h = self.h
        Iyy = 4.594369333645184 * 10**(-5)
        Izz = 4.753851442684437 * 10**(-6)
        
        #spanwise locations
        x = b*np.linspace(0,1,10)

        #moment z varying spanwise
        Mz = F*(xF - x)
        print(Mz)

        #direct stress in y direction due to force (constant throughout span and chord)
        sigmay1 = np.matrix(len(x)*[F/(c*b)])

        #bending stress in y direction due to force (varying spanwise)
        sigmay2 = np.matrix(Mz/Izz*x)

        sigmay = sigmay1 + sigmay2

        #heightwise locations
        y = 0.5*h*np.linspace(-1,1,10)

        #bending stress in x direction due to force (varying in height and spanwise)
        sigmax = []
        for a in range(len(x)):
            sigmax.append(Mz[a]/Izz*y)

        sigmax = np.matrix(sigmax)

        return x,y,sigmay, sigmax
    
    def sigmaduetoxmoment(self, Mx):
        c = self.c_a
        h = self.h
        Ixx = 1./12.*self.l_a*self.c_a**3 - 1./12.*(self.l_a-2*self.t_sk)*(self.c_a-2*self.t_sk)**3

        #chord locations
        z = np.linspace(-self.h/2.,self.c_a-self.h/2.,10)

        #sigmay due to moment in x varies chordwise
        sigmay = np.matrix(Mx/Ixx*z)

        #height locations
        y = 0.5*h*np.linspace(-1,1,10)

        #sigmaz due to moment in x varies heightwise
        sigmaz = np.matrix(Mx/Ixx*y)

        return z, y, sigmay, sigmaz

    def sigmaduetoq(self,F,xF,zF): 
        c = self.c_a
        b = self.l_a
        h = self.h
        Iyy = 4.594369333645184e-05
        Izz = 4.753851442684437e-06
        Ixx = 1./12.*self.l_a*self.c_a**3 - 1./12.*(self.l_a-2*self.t_sk)*(self.c_a-2*self.t_sk)**3

        print(Ixx)


        #spanwise locations
        x = b*np.linspace(0,1,10)

        #chordwise locations
        z = np.linspace(-self.h/2.,self.c_a-self.h/2.,10)

        #moment z varying spanwise
        Mz = F*(xF - x)

        Mx = F*(zF -z)

        #direct stress in y direction due to force (constant throughout span and chord)
        sigmay1 = np.matrix(len(x)*[F/(c*b)])

        #bending stress in y direction due to force (varying spanwise)
        sigmay2 = np.matrix(Mz/Izz*x)

        #bending stress in y direction due to force (varying chordwise)
        sigmay3 = np.matrix(Mx/Ixx*z)

        sigmay = sigmay1 + sigmay2 + sigmay3

        #heightwise locations
        y = 0.5*h*np.linspace(-1,1,10)

        #bending stress in x direction due to force (varying in height and spanwise)
        sigmax = []
        for a in range(len(x)):
            sigmax.append(Mz[a]/Izz*y)


        sigmax = np.matrix(sigmax)


        #bending stress in z direction due to force(varying in height and chord)
        sigmaz = []
        for b in range(len(z)):
            sigmaz.append(Mx[b]/Ixx*y)
        
        sigmaz = np.matrix(sigmaz)

        return x,y,z, sigmax, sigmay, sigmaz

    def sigmaduetozforce(self,F,xF):
        c = self.c_a
        b = self.l_a
        h = self.h
        Iyy = 4.594369333645184e-05
        Izz = 4.753851442684437e-06

        #span locations
        x = b*np.linspace(0,1,10)

        #direct stress in z direction due to force
        sigmaz1 = np.matrix(len(x)*[F/(b*h)])

        #moment in y (My) varying over span
        My = Pz*(xF - x)

        #bending stress in z (sigmaz) varying over span
        sigmaz2 = np.matrix(My/Iyy*x)

        sigmaz = sigmaz1 + sigmaz2

        #chord locations
        z = np.linspace(-self.h/2.,self.c_a-self.h/2.,10)

        #bending stress in x direction (sigmax) varying over span and chord
        sigmax = []
        for a in range(len(x)):
            sigmax.append(My[a]/Izz*z)

        sigmax = np.matrix(sigmax)

        return x, z, sigmaz, sigmax

    def sigmaduetoxforce(self, F):
        h = self.h
        c = self.c_a
        sigmax = F/(c*h)

        return sigmax



test = stressdistribution(**parameters_geometry)

x_hinge1 = 0.125
x_hinge2 = 0.498
x_hinge3 = 1.494
x_actuator1 = 0.3755
x_actuator2 = 0.6205

P = -49200
Py = np.sin(np.radians(30))*P
Pz = np.cos(np.radians(30))*P

Rz1 = 0 #-145131.4967561881
Rz2 = 0 #219109.49055686826
Rz3 = 0 #-62191.58125357197

Ry1 = 0 #84700.72749871935
Ry2 = 0 #-99334.78822558284
Ry3 = 0 #15253.049364200167

Ra = 0 #-62809.776913705806
Rya = 0 #Ra*np.sin(np.radians(30))
Rza = 0 #Ra*np.cos(np.radians(30))
Rxa = 0

Mxa = 0 #Ra*test.h/2.*(np.sin(np.radians(30))- np.cos(np.radians(30)))
Mxp = P*test.h/2.*(np.sin(np.radians(30))- np.cos(np.radians(30)))


x1, y1, sigmay1, sigmax1 = test.sigmaduetoyforce(Py, x_actuator2)
x1, y1, sigmay2, sigmax2 = test.sigmaduetoyforce(Ry1, x_hinge1)
x1, y1, sigmay3, sigmax3 = test.sigmaduetoyforce(Ry2, x_hinge2)
x1, y1, sigmay4, sigmax4 = test.sigmaduetoyforce(Ry3, x_hinge3)
x1, y1, sigmay5, sigmax5 = test.sigmaduetoyforce(Rya, x_actuator1)

z1, y2, sigmay6, sigmaz1 = test.sigmaduetoxmoment(Mxa)
z1, y2, sigmay7, sigmaz2 = test.sigmaduetoxmoment(Mxp)

x2, z2, sigmaz3, sigmax6 = test.sigmaduetozforce(Pz, x_actuator2)
x2, z2, sigmaz4, sigmax7 = test.sigmaduetozforce(Rz1, x_hinge1)
x2, z2, sigmaz5, sigmax8 = test.sigmaduetozforce(Rz2, x_hinge2)
x2, z2, sigmaz6, sigmax9 = test.sigmaduetozforce(Rz3, x_hinge3)
x2, z2, sigmaz7, sigmax10 = test.sigmaduetozforce(Rza, x_actuator1)


#########################   sigma y #############################################################

#these vary over span
sigmaytotal1array = np.array(sigmay1 + sigmay2 + sigmay3 + sigmay4 + sigmay5)
sigmaytotal1 = sigmaytotal1array[0]

plt.plot(x1,sigmaytotal1)
plt.show()

#vary chordwise
sigmaytotal2array = np.array(sigmay6 + sigmay7)
sigmaytotal2 = sigmaytotal2array[0]

plt.plot(z1,sigmaytotal2)
plt.show()


###################### sigma z #############################################################
#vary heightwise
sigmaztotal1array = np.array(sigmaz1 + sigmaz2)
sigmaztotal1 = sigmaztotal1array[0]

plt.plot(y1,sigmaztotal1)
plt.show()

#these vary over span
sigmaztotal2array = np.array(sigmaz3 + sigmaz4 + sigmaz5 + sigmaz6 + sigmaz7)
sigmaztotal2 = sigmaztotal2array[0]

plt.plot(x1,sigmaztotal2)
plt.show()


####################################sigma x#########################################
#direct force
sigmax10 = test.sigmaduetoxforce(Rxa)


#these vary over span and height
sigmaxtotal1 = sigmax1 + sigmax2 + sigmax3 + sigmax4 + sigmax5

#take critical spanwise cut
sigmaxtotal1criticallst = np.array(sigmaxtotal1[-1])
sigmaxtotal1critical = sigmaxtotal1criticallst[0]




#take critical crosssectional cut
#sigmatotal1crosssection = sigmaxtotal1[:][-1]

#print(sigmaxtotal1)

#fig = plt.figure()
#cs = plt.imshow(sigmaxtotal1, extent = (0., 1.611, -0.0805, 0.0805), cmap='Blues')
#cbar= plt.colorbar(cs, shrink =0.5)
#cbar.ax.set_ylabel("Bending in x direction [Pa]")
#plt.xlabel("Span-wise direction [m]")
#plt.ylabel("Height direction [m]")
#plt.title("Bending stress in x direction (sigma_x)")
#plt.tight_layout()
#plt.show()

#these vary over span and chord
sigmaxtotal2 = np.array(sigmax6 + sigmax7 + sigmax8 + sigmax9 + sigmax10)

maximumtotal2 = sigmaxtotal2[-1][:]


fig = plt.figure()
cs = plt.imshow(sigmaxtotal2, extent = (0., 1.611, 0, 0.505), cmap='Blues')
cbar= plt.colorbar(cs, shrink =0.5)
cbar.ax.set_ylabel("Bending in x direction [Pa]")
plt.xlabel("Span-wise direction [m]")
plt.ylabel("Chord-wise direction [m]")
plt.title("Bending stress in x direction (sigma_x)")
plt.tight_layout()
plt.show()



#fig = plt.figure()
#cs = plt.imshow(criticalcrosssection, extent = (0., 0.505, -0.0805, 0.0805), cmap='Blues')
#cbar= plt.colorbar(cs, shrink =0.5)
#cbar.ax.set_ylabel("Sigmaxx in crosssection [Pa]")
#plt.xlabel("Chord-wise direction [m]")
#plt.ylabel("Height-wise direction [m]")
#plt.title("Bending stress in x direction (sigma_x)")
#plt.tight_layout()
#plt.show()


#################### aeroload #######################
interp = inter.interpolation()
p = 10
xlst = test.l_a*np.linspace(0,1,p)
zlst = np.linspace(-test.h/2.,test.c_a-test.h/2.,p)

l = test.l_a/p
A = (test.l_a*test.c_a)/(p*p) 

sigmaxtotal = 0
sigmaytotal = 0
sigmaztotal = 0


for x in xlst:
    for z in zlst:
        qk, a = interp.bilinear_interpolation(x,z)
        Fk = qk*l
        xn, yn, zn, sigmax, sigmay, sigmaz = test.sigmaduetoq(Fk,x,z)

        sigmaxtotal += sigmax
        sigmaytotal += sigmay
        sigmaztotal += sigmaz

maximumq = sigmaxtotal[-1][:]

#fig = plt.figure()
#cs = plt.imshow(sigmaztotal, extent = (0., 1.611, 0, 0.505), cmap='Blues')
#cbar= plt.colorbar(cs, shrink =0.5)
#cbar.ax.set_ylabel("Sigmaxx in crosssection [Pa]")
#plt.xlabel("Chord-wise direction [m]")
#plt.ylabel("Height-wise direction [m]")
#plt.title("Bending stress in x direction (sigma_x)")
#plt.tight_layout()
#plt.show()


############## sigma x total ########################


####### crosssection #######

#varies over chord
newtotallst = np.array(maximumq + maximumtotal2)
newtotal = newtotallst[0]

print(maximumq)

#varies over height
#sigmaxtotal1critical 

sigmafinalx = []

for c in range(10):
    sigmafinalx.append([])
    for d in range(10):
        sigmafinalx[c].append(newtotal[c] + sigmaxtotal1critical[d] + sigmax10)


fig = plt.figure()
cs = plt.imshow(sigmafinalx, extent = (0., 0.505, -0.0805, 0.0805 ), cmap='Blues')
cbar= plt.colorbar(cs, shrink =0.5)
cbar.ax.set_ylabel("Sigmaxx in crosssection [Pa]")
plt.xlabel("Chord-wise direction [m]")
plt.ylabel("Height-wise direction [m]")
plt.title("Bending stress in x direction (sigma_x)")
plt.tight_layout()
plt.show()

