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
    def crosssection(self):
        beta = self.beta
        theta2 = self.strint/(self.h/2.)

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
        beta = self.beta
        #l1 = self.l1

        z1 = (self.c_a-self.h/2.)/2. + self.h/2.
        z2 = (self.c_a-self.h/2.)/2. + self.h/2.
        z3 = self.h/2.
        z4 = self.h/2. - 2*(self.h/2.)/(np.pi)

        theta2 = self.strint/(self.h/2.)

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

        y1 = self.h/4.
        y2 = self.h/4.*3.
        y3 = self.h/2.
        y4 = self.h/2.

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

        A1 = self.l1*self.t_sk
        A2 = self.l1*self.t_sk
        A3 = self.h*self.t_sp
        A4 =  np.pi*self.h/2.*self.t_sk
        Ast = (self.w_st-self.t_st+self.h_st)*self.t_st

        c_y = (y1*A1+y2*A2+y3*A3+y4*A4+(yst1+yst2+yst3+yst4+yst5+yst6+yst7+yst8+yst9+yst10+yst11)*Ast)/(A1+A2+A3+A4+11*Ast)
        c_z = (z1*A1+z2*A2+z3*A3+z4*A4+(zst1+zst2+zst3+zst4+zst5+zst6+zst7+zst8+zst9+zst10+zst11)*Ast)/(A1+A2+A3+A4+11*Ast)
        return c_y, c_z, Ast

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
        Iyy4 = 1./2.*(self.h/2.)**3*np.pi*self.t_sk # + np.pi*self.t_sk*self.h/2.*(self.centroid[1]-self.h/2.)**2
        Iyyst1 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-self.h/2.)**2
        Iyyst2 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-(self.h/2.+l1/4.*np.cos(beta)))**2
        Iyyst3 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-(self.h/2.+l1/2.*np.cos(beta)))**2
        Iyyst4 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-(self.h/2.+l1*3/4.*np.cos(beta)))**2
        Iyyst5 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-(self.h/2.+l1*3/4.*np.cos(beta)))**2
        Iyyst6 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-(self.h/2.+l1/2.*np.cos(beta)))**2
        Iyyst7 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-(self.h/2.+l1/4.*np.cos(beta)))**2
        Iyyst8 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-self.h/2.)**2
        Iyyst9 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-self.h/4.)**2
        Iyyst10 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1])**2
        Iyyst11 = (self.w_st+self.h_st)*self.t_st * (self.centroid[1]-self.h/4.)**2
        

        Izz1 = self.t_sk*(l1)**3*(np.sin(beta))**2/12. + l1*self.t_sk*(self.centroid[0]-self.h/4.)**2
        Izz2 = self.t_sk*(l1)**3*(np.sin(beta))**2/12. + l1*self.t_sk*(self.centroid[0]-self.h/4.*3.)**2
        Izz3 = self.t_sp*(self.h**3)/12. + self.h*self.t_sp*(self.centroid[0]-self.h/2.)**2
        Izz4 = 1./2.*(self.h/2.)**3*np.pi*self.t_sk # + np.pi*self.t_sk*self.h/2.*(self.centroid[0]-self.h/2.)**2
        Izzst1 = (self.w_st+self.h_st)*self.t_st* (self.centroid[0])**2
        Izzst2 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-(l1/4.*np.sin(beta)))**2
        Izzst3 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-(l1/2.*np.sin(beta)))**2
        Izzst4 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-(l1*3./4.*np.sin(beta)))**2
        Izzst5 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-(l1*5./4.*np.sin(beta)))**2
        Izzst6 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-(l1*3./2.*np.sin(beta)))**2
        Izzst7 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-(l1*7./4.*np.sin(beta)))**2
        Izzst8 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-self.h)**2
        Izzst9 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-self.h*3./4.)**2
        Izzst10 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-self.h/2.)**2
        Izzst11 = (self.w_st+self.h_st)*self.t_st * (self.centroid[0]-self.h/4.)**2

        Iyy = Iyy1 + Iyy2 + Iyy3 + Iyy4  + Iyyst1 + Iyyst2 + Iyyst3 + Iyyst4 + Iyyst5 + Iyyst6 + Iyyst7 + Iyyst8 + Iyyst9 + Iyyst10 + Iyyst11
        Izz = Izz1 + Izz2 + Izz3 + Izz4  + Izzst1 + Izzst2 + Izzst3 + Izzst4 + Izzst5 + Izzst6 + Izzst7 + Izzst8 + Izzst9 + Izzst10 + Izzst11

        return Iyy4, Izz4


    @property
    def shearcenter(self):
        #Iyy_verification = 4.5943507864451845e-05, Izz_verification = 4.753851442684436e-06
        l = self.c_a - self.h/2.
        alpha = np.arctan2(self.h/2.,l)
        d = l/np.cos(alpha)
        #shear center distance calculated from the leading edge
        d_z =((self.t_sk*self.h*self.h*d*d*(1./3.-np.cos(alpha)/l-(self.h*np.cos(alpha))/(3.*l))+self.h*self.h*((self.t_sk*self.h*self.h)/4.-2*self.t_sk*self.h*self.h+(4.*d*np.cos(alpha))/(self.h*self.h)+(self.t_sk*self.h*self.h)/4.))/(4.*4.753851442684436e-06))+self.h/2.
        return d_z

        #y_sc = 0
        #z_sc = -0.22554053758032344


if __name__ == "__main__": # is called when you run the script
    # call an instance of the class
    geo = Geometry(**parameters_geometry) 

    print(geo.MMoI)
   