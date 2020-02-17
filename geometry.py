import numpy as np
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
    # stringers only have a steiner term (assumption!)

    def MMoI(self):
        beta = np.arctan2(self.h/2,(self.c_a-self.h/2))
        l1 = np.sqrt((self.c_a-self.h/2.)**2 + (self.h/2.)**2)

        # Iyy1 =  
        # Iyy2 =
        Iyy3 = self.h*(self.t_sp**3)/12.
        # Iyy4 = 



        return 
    @property
    #the y coordinate of the centroid is given w.r.t to the top of airfoil (but the value would be the same if the bottom is taken as reference line)
    #the z coordinate of the centorid is given w.r.t to the left side of the airfoil
    #1 is top skin
    #2 is bottom skin
    #3 is the spar
    #4 is the circular skin
    #st is stringer

    def centroid(self):
        beta = np.arctan2(self.h/2,(self.c_a-self.h/2))

        z1 = (self.c_a-self.h/2.)/2. + self.h/2
        z2 = (self.c_a-self.h/2.)/2. + self.h/2
        z3 = self.h/2 - (4./3.)*self.h/(2*np.pi)
        z4 = self.h/2

        zst1 = self.h/2
        zst2 = self.h/2 + np.cos(beta)/4.*l1
        zst3 = self.h/2 + np.cos(beta)/2.*l1
        zst4 = self.h/2 + np.cos(beta)/4.*l1*3.
        zst5 = self.h/2 + np.cos(beta)/4.*l1*3.
        zst6 = self.h/2 + np.cos(beta)/2.*l1
        zst7 = self.h/2 + np.cos(beta)/4.*l1
        zst8 = self.h/2.
        zst9 = self.h/4.
        zst10 = 0
        zst11 = self.h/4.

        y1 = self.h/4.
        y2 = self.h/4.*3.
        y3 = self.h/2.
        y4 = self.h/2.

        yst1 = 0
        yst2 = np.sin(beta)/4.*l1
        yst3 = np.sin(beta)/2.*l1
        yst4 = np.sin(beta)/4.*l1*3.
        yst5 = np.sin(beta)/4.*l1*5
        yst6 = np.sin(beta)/4.*l1*6
        yst7 = np.sin(beta)/4.*l1*7
        yst8 = self.h
        yst9 = 3./4.*self.h
        yst10 = 1./2.*self.h
        yst11 = 1./4.*self.h

        A1 = np.sqrt((self.c_a-self.h/2.)**2 + (self.h/2.)**2)*self.t_sk
        A2 = np.sqrt((self.c_a-self.h/2.)**2 + (self.h/2.)**2)*self.t_sk
        A3 = np.pi*self.h*self.t_sk
        A4 = self.h*self.t_sp
        Ast = (self.w_st-self.t_st+self.h_st)*self.t_st

        c_y = (y1*A1+y2*A2+y3*A3+y4*A4+(yst1+yst2+yst3+yst4+yst5+yst6+yst7+yst8+yst9+yst10+yst11)*Ast)/(A1+A2+A3+A4+11*Ast)
        c_z = (z1*A1+z2*A2+z3*A3+z4*A4+(zst1+zst2+zst3+zst4+zst5+zst6+zst7+zst8+zst9+zst10+zst11)*Ast)/(A1+A2+A3+A4+11*Ast)
        return c_y, c_z


    @property
    def shearcenter(self):
        return


if __name__ == "__main__": # is called when you run the script
    # call an instance of the class
    geo = Geometry(**parameters_geometry) 

    print(geo.centroid)

