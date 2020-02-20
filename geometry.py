import numpy as np
import matplotlib.pyplot as plt
from data.consts import parameters_geometry

class Geometry:
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
    

    @property

    def beta(self):
        beta = np.arctan2(self.h/2,(self.c_a-self.h/2))
        #beta = np.arctan(self.h/2/(self.c_a-self.h/2.))
    
        return beta

    @property

    def l1(self):
        l1 = np.sqrt((self.c_a-self.h/2.)**2 + (self.h/2.)**2)

        return l1

    @property
    def ltotal(self):
        ltotal = np.pi*self.h/2 + 2*self.l1

        return ltotal
    
    @property
    def strint(self):
        strint = self.ltotal/11.

        return strint

    @property
    def ltop(self):
        ltop = 1./2.*np.pi*self.h/2. + self.l1

        return ltop
    
    @property

    def strint2(self):
        strint2 = self.ltop/5.5

        return strint2
    
    @property
    def theta2(self):
        theta2 = self.strint/(self.h/2.)
        return theta2

    @property
    def Ast(self):
        Ast = (self.w_st+self.h_st)*self.t_st
        return Ast

    @property
    def crosssection(self):
        beta = self.beta
        theta2 = self.theta2

        zst1 = 0
        zst2 = self.h/2. - self.h/2.*np.cos(theta2)
        zst3 = (2*self.strint-np.pi*self.h/4.)*np.cos(beta) + self.h/2.
        zst4 = (3*self.strint-np.pi*self.h/4.)*np.cos(beta) + self.h/2.
        zst5 = (4*self.strint-np.pi*self.h/4.)*np.cos(beta) + self.h/2.
        zst6 = (5*self.strint-np.pi*self.h/4.)*np.cos(beta) + self.h/2.
        zst7 = (5*self.strint-np.pi*self.h/4.)*np.cos(beta) + self.h/2.
        zst8 = (4*self.strint-np.pi*self.h/4.)*np.cos(beta) + self.h/2.
        zst9 = (3*self.strint-np.pi*self.h/4.)*np.cos(beta) + self.h/2.
        zst10 = (2*self.strint-np.pi*self.h/4.)*np.cos(beta) + self.h/2.
        zst11 = self.h/2. - self.h/2.*np.cos(theta2)

        yst1 = self.h/2.
        yst2 = self.h/2. - self.h/2.*np.sin(theta2)
        yst3 = (2*self.strint-np.pi*self.h/4.)*np.sin(beta)
        yst4 = (3*self.strint-np.pi*self.h/4.)*np.sin(beta)
        yst5 = (4*self.strint-np.pi*self.h/4.)*np.sin(beta)
        yst6 = (5*self.strint-np.pi*self.h/4.)*np.sin(beta)
        yst7 = (6*self.strint-np.pi*self.h/4.)*np.sin(beta)
        yst8 = (7*self.strint-np.pi*self.h/4.)*np.sin(beta)
        yst9 = (8*self.strint-np.pi*self.h/4.)*np.sin(beta)
        yst10 = (9*self.strint-np.pi*self.h/4.)*np.sin(beta)
        yst11 = self.h/2. + self.h/2.*np.sin(theta2)

        return ([zst1,zst2,zst3,zst4,zst5,zst6,zst7,zst8,zst9,zst10,zst11],[yst1,yst2,yst3,yst4,yst5,yst6,yst7,yst8,yst9,yst10,yst11])

    @property
    #the y coordinate of the centroid is given w.r.t to the top of airfoil (but the value would be the same if the bottom is taken as reference line)
    #the z coordinate of the centroid is given w.r.t to the left side of the airfoil
    #1 is top skin
    #2 is bottom skin
    #3 is the spar
    #4 is the circular skin
    #st is stringer
    #first stringer located on leading edge

    def centroid(self):
        z1 = (self.c_a-self.h/2.)/2. + self.h/2.
        z2 = (self.c_a-self.h/2.)/2. + self.h/2.
        z3 = self.h/2.
        z4 = self.h/2. - 2*(self.h/2.)/(np.pi)

        zst1 = self.crosssection[0][0]
        zst2 = self.crosssection[0][1]
        zst3 = self.crosssection[0][2]
        zst4 = self.crosssection[0][3]
        zst5 = self.crosssection[0][4]
        zst6 = self.crosssection[0][5]
        zst7 = self.crosssection[0][6]
        zst8 = self.crosssection[0][7]
        zst9 = self.crosssection[0][8]
        zst10 = self.crosssection[0][9]
        zst11 = self.crosssection[0][10]

        y1 = self.h/4.
        y2 = self.h/4.*3.
        y3 = self.h/2.
        y4 = self.h/2.

        yst1 = self.crosssection[1][0]
        yst2 = self.crosssection[1][1]
        yst3 = self.crosssection[1][2]
        yst4 = self.crosssection[1][3]
        yst5 = self.crosssection[1][4]
        yst6 = self.crosssection[1][5]
        yst7 = self.crosssection[1][6]
        yst8 = self.crosssection[1][7]
        yst9 = self.crosssection[1][8]
        yst10 = self.crosssection[1][9]
        yst11 = self.crosssection[1][10]

        A1 = self.l1*self.t_sk
        A2 = self.l1*self.t_sk
        A3 = self.h*self.t_sp
        A4 =  np.pi*self.h/2.*self.t_sk
        Ast = self.Ast

        c_y = (y1*A1+y2*A2+y3*A3+y4*A4+(yst1+yst2+yst3+yst4+yst5+yst6+yst7+yst8+yst9+yst10+yst11)*Ast)/(A1+A2+A3+A4+11*Ast)
        c_z = (z1*A1+z2*A2+z3*A3+z4*A4+(zst1+zst2+zst3+zst4+zst5+zst6+zst7+zst8+zst9+zst10+zst11)*Ast)/(A1+A2+A3+A4+11*Ast)
        return c_y, c_z

    @property
    # stringers only have a steiner term (assumption!)
    #1 is top skin
    #2 is bottom skin
    #3 is the spar
    #4 is the circular skin
    #st is stringer


    def MMoI(self):
        beta = self.beta
        l1 = self.l1

        Iyy1 = self.t_sk*(l1)**3*(np.cos(beta))**2/12. + l1*self.t_sk* (self.centroid[1]-(self.h/2.+(self.c_a-self.h/2.)/2.))**2
        Iyy2 = self.t_sk*(l1)**3*(np.cos(beta))**2/12. + l1*self.t_sk* (self.centroid[1]-(self.h/2.+(self.c_a-self.h/2.)/2.))**2
        Iyy3 = self.h*(self.t_sp**3)/12. + self.h*self.t_sp*(self.centroid[1]-self.h/2.)**2
        Iyy4 = 1./2.*(self.h/2.)**3*np.pi*self.t_sk - (np.pi*self.h/2.*self.t_sk)*(self.h/np.pi)**2 + np.pi*self.t_sk*self.h/2.*(self.centroid[1]-(self.h/2. - 2*(self.h/2.)/(np.pi)))**2
        Iyyst1 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-self.crosssection[0][0])**2
        Iyyst2 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-self.crosssection[0][1])**2
        Iyyst3 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-self.crosssection[0][2])**2
        Iyyst4 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-self.crosssection[0][3])**2
        Iyyst5 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-self.crosssection[0][4])**2
        Iyyst6 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-self.crosssection[0][5])**2
        Iyyst7 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-self.crosssection[0][6])**2
        Iyyst8 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-self.crosssection[0][7])**2
        Iyyst9 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-self.crosssection[0][8])**2
        Iyyst10 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-self.crosssection[0][9])**2
        Iyyst11 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-self.crosssection[0][10])**2
        

        Izz1 = self.t_sk*(l1)**3*(np.sin(beta))**2/12. + l1*self.t_sk*(self.centroid[0]-self.h/4.)**2
        Izz2 = self.t_sk*(l1)**3*(np.sin(beta))**2/12. + l1*self.t_sk*(self.centroid[0]-self.h/4.*3.)**2
        Izz3 = self.t_sp*(self.h**3)/12. + self.h*self.t_sp*(self.centroid[0]-self.h/2.)**2
        Izz4 = 1./2.*(self.h/2.)**3*np.pi*self.t_sk + np.pi*self.t_sk*self.h/2.*(self.centroid[0]-self.h/2.)**2
        Izzst1 = (self.w_st+self.h_st)*self.t_st* (self.centroid[0]-self.crosssection[1][0])**2
        Izzst2 = (self.w_st+self.h_st)*self.t_st *(self.centroid[0]-self.crosssection[1][1])**2
        Izzst3 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-self.crosssection[1][2])**2
        Izzst4 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-self.crosssection[1][3])**2
        Izzst5 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-self.crosssection[1][4])**2
        Izzst6 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-self.crosssection[1][5])**2
        Izzst7 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-self.crosssection[1][6])**2
        Izzst8 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-self.crosssection[1][7])**2
        Izzst9 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-self.crosssection[1][8])**2
        Izzst10 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-self.crosssection[1][9])**2
        Izzst11 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-self.crosssection[1][10])**2

        Iyy = Iyy1 + Iyy2 + Iyy3 + Iyy4  + Iyyst1 + Iyyst2 + Iyyst3 + Iyyst4 + Iyyst5 + Iyyst6 + Iyyst7 + Iyyst8 + Iyyst9 + Iyyst10 + Iyyst11
        Izz = Izz1 + Izz2 + Izz3 + Izz4  + Izzst1 + Izzst2 + Izzst3 + Izzst4 + Izzst5 + Izzst6 + Izzst7 + Izzst8 + Izzst9 + Izzst10 + Izzst11

        return Iyy, Izz


    @property
    def shearcenter(self):
        #Iyy_verification = 4.5943507864451845e-05, Izz_verification = 4.753851442684436e-06

        l = self.c_a - self.h/2.
        alpha = np.arctan2(self.h/2.,l)
        d = l/np.cos(alpha)
        #shear center distance calculated from the leading edge
        dz = l-(self.t_sk*self.h*self.h*d*l)/(12.*self.MMoI[1])+(d*d*self.h*self.h*self.t_sk*np.cos(alpha))/(4.*self.MMoI[1])+(d*self.h*self.h*self.h*self.t_sk*np.cos(alpha))/(12.*self.MMoI[1])-(self.h*self.h*self.h*self.h*self.t_sk*l)/(8.*d*self.MMoI[1])
        dz = -dz
        return dz

        #y_sc = 0
        #z_sc = -0.22554053758032344

    @property
    def torsionalstiffness(self):
        



if __name__ == "__main__": # is called when you run the script
    # call an instance of the class
    geo = Geometry(**parameters_geometry) 

    print(geo.MMoI)
   