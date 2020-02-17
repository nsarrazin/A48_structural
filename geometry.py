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
    def MMoI(self):
        print(self.c_a) # we can access things now
        print(self.x_a)

        return self.c_a + self.x_a

    @property
    #the y coordinate of the centroid is given w.r.t to the top of airfoil (but the value would be the same if the bottom is taken as reference line)
    #the z coordinate of the centorid is given w.r.t to the left side of the airfoil
    
    def centroid(self):
        z1 = (self.c_a-self.h/2.)/2. + self.h/2
        z2 = (self.c_a-self.h/2.)/2. + self.h/2
        z3 = self.h/2 - (4./3.)*self.h/(2*np.pi)
        z4 = self.h/2

        y1 = self.h/4.
        y2 = self.h/4.*3.
        y3 = self.h/2.
        y4 = self.h/2.

        A1 = np.sqrt((self.c_a-self.h/2.)**2 + (self.h/2.)**2)*self.t_sk
        A2 = np.sqrt((self.c_a-self.h/2.)**2 + (self.h/2.)**2)*self.t_sk
        A3 = np.pi*self.h*self.t_sk
        A4 = self.h*self.t_sp

        c_y = (y1*A1+y2*A2+y3*A3+y4*A4)/(A1+A2+A3+A4)
        c_z = (z1*A1+z2*A2+z3*A3+z4*A4)/(A1+A2+A3+A4)
        return c_y, c_z


    @property
    def shearcenter(self):
        return


if __name__ == "__main__": # is called when you run the script
    # call an instance of the class
    geo = Geometry(**parameters_geometry) 

    #print(geo.MMoI)
    print(geo.centroid)
