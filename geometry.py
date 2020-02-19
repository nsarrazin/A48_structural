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
    #the y coordinate of the centroid is given w.r.t to the top of airfoil (but the value would be the same if the bottom is taken as reference line)
    #the z coordinate of the centroid is given w.r.t to the left side of the airfoil
    #1 is top skin
    #2 is bottom skin
    #3 is the spar
    #4 is the circular skin
    #st is stringer

    def centroid(self):
        beta = np.arctan2(self.h/2,(self.c_a-self.h/2))
        l1 = np.sqrt((self.c_a-self.h/2.)**2 + (self.h/2.)**2)

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
    # stringers only have a steiner term (assumption!)
    #1 is top skin
    #2 is bottom skin
    #3 is the spar
    #4 is the circular skin
    #st is stringer


    def MMoI(self):
        beta = np.arctan2(self.h/2,(self.c_a-self.h/2))
        l1 = np.sqrt((self.c_a-self.h/2.)**2 + (self.h/2.)**2)

        Iyy12 = self.t_sk*(l1*2)**3*(np.cos(beta))**2/12. + 2*l1*self.t_sk*(self.centroid[1]-(self.h/2+l1*np.cos(beta)))**2
        Iyy3 = self.h*(self.t_sp**3)/12. + self.h*self.t_sp*(self.centroid[1]-self.h/2.)**2
        Iyy4 = 1./2.*(self.h/2.)**3*np.pi*self.t_sk + np.pi*self.t_sk*self.h/2.*(self.centroid[1]-self.h/2.)**2
        Iyyst1 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[1]-self.h/2.)**2
        Iyyst2 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[1]-(self.h/2.+l1/4.*np.cos(beta))**2
        Iyyst3 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[1]-(self.h/2.+l1/2.*np.cos(beta))**2
        Iyyst4 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[1]-(self.h/2.+l1*3/4.*np.cos(beta))**2
        Iyyst5 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[1]-(self.h/2.+l1*3/4.*np.cos(beta))**2
        Iyyst6 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[1]-(self.h/2.+l1/2.*np.cos(beta))**2
        Iyyst7 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[1]-(self.h/2.+l1/4.*np.cos(beta))**2
        Iyyst8 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[1]-self.h/2.)**2
        Iyyst9 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[1]-self.h/4.)**2
        Iyyst10 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[1])**2
        Iyyst11 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[1]-self.h/4.)**2
        

        Izz12 = self.t_sk*(l1*2)**3*(np.sin(beta))**2/12.
        Izz3 = self.t_sp*(self.h**3)/12. + self.h*self.t_sp*(self.centroid[0]-self.h/2.)**2
        Izz4 = 1./2.*(self.h/2.)**3*np.pi*self.t_sk + np.pi*self.t_sk*self.h/2.*(self.centroid[0]-self.h/2.)**2
        Izzst1 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[0])**2
        Izzst2 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[0]-(l1/4.*np.sin(beta))**2
        Izzst3 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[0]-(l1/2.*np.sin(beta))**2
        Izzst4 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[0]-(l1*3./4.*np.sin(beta))**2
        Izzst5 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[0]-(l1*5./4.*np.sin(beta))**2
        Izzst6 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[0]-(l1*3./2.*np.sin(beta))**2
        Izzst7 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[0]-(l1*7./4.*np.sin(beta))**2
        Izzst8 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[0]-self.h)**2
        Izzst9 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[0]-self.h*3./4.)**2
        Izzst10 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[0]-self.h/2.)**2
        Izzst11 = self.w_st-self.t_st+self.h_st*self.t_st * (self.centroid[0]-self.h/4.)**2

        Iyy = Iyy12 + Iyy3 + Iyy4 + Iyyst1 + Iyyst2 + Iyyst3 + Iyyst4 + Iyyst5 + Iyyst6 + Iyyst7 + Iyyst8 + Iyyst9 + Iyyst10 + Iyyst11
        Izz = Izz12 + Izz3 + Izz4 + Izzst1 + Izzst2 + Izzst3 + Izzst4 + Izzst5 + Izzst6 + Izzst7 + Izzst8 + Izzst9 + Izzst10 + Izzst11

        return Iyy, Izz


    @property
    def shearcenter(self):
        return


if __name__ == "__main__": # is called when you run the script
    # call an instance of the class
    geo = Geometry(**parameters_geometry) 

    print(geo.MMoI)

