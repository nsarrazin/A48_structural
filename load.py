import numpy as np
from helpers import step

class LoadCase:
    def __init__(self, parent, **kwargs):
        self.geo = parent.geo
        self.interp = parent.interp

        self.d_1 = kwargs.get("d_1")
        self.d_3 = kwargs.get("d_3")
        self.defl = kwargs.get("defl") # theta
        self.P = kwargs.get("load")
        self.E = kwargs.get("E")
        self.G = kwargs.get("G")

        self.a_y = np.sin(self.defl)
        self.a_z = np.cos(self.defl)
        self.a_m = self.geo.h/2 * (np.sin(self.defl)- np.cos(self.defl))
        
        _ , self.z_sc = self.geo.shearcenter
        # _ , self.z_sc = 0, 0
        self.x_I = self.geo.x_2-self.geo.h/2
        self.x_II = self.geo.x_2+self.geo.h/2

    def V_y_prime(self, x):
        return np.array([step(x,self.geo.x_1)**0,     #Fy_1
                         0, #Fz_1
                        step(x,self.geo.x_2)**0,     #Fy_2
                        0,                          #Fz_2
                        step(x,self.geo.x_3),      #Fy_3
                        0,   #Fz_3
                        self.a_y*step(x,self.x_I)**0,    #Fa
                        0,  #C1
                        0,  #C2
                        0,  #C3
                        0,  #C4
                        0  #C5
                        ])

    def V_z_prime(self,x):
        return np.array([0,    #Fy_1
                         step(x,self.geo.x_1)**0,     #Fz_1
                         0, #Fy_2
                         step(x,self.geo.x_2)**0,     #Fz_2
                         0, #Fy_3
                         step(x,self.geo.x_3)**0,  #Fz_3
                         self.a_z*step(x,self.x_I)**0,   #Fa
                         0, #C1
                         0, #C2
                         0, #C3
                         0, #C4
                         0 #C5
                         ])

    def M_y_prime(self, x):
        return np.array([0,    #Fy_1
                         -step(x,self.geo.x_1)**1, #Fz_1
                         0, #Fy_2
                         -step(x,self.geo.x_2)**1,  #Fz_2
                         0,     #Fy_3
                         -step(x,self.geo.x_3)**1,   #Fz_3
                         -self.a_y*step(x,self.x_I)**1, #Fa
                         0, #C1
                         0, #C2
                         0, #C3
                         0, #C4
                         0 #C5
                         ])

    def M_z_prime(self, x):
        return np.array([-step(x,self.geo.x_1)**1,    #Fy_1
                         0,  #Fz_1
                         -step(x,self.geo.x_2)**1,   #Fy_2
                         0,    #Fz_2
                         -step(x,self.geo.x_3)**1,     #Fy_3
                         0,     #Fz_3
                         -self.a_y*step(x,self.x_I)**1, #Fa
                         0, #C1
                         0, #C2
                         0, #C3
                         0, #C4
                         0  #C5
                         ])

    def T(self, x):
        return np.array([self.z_sc*step(x,self.geo.x_1)**0,    #Fy_1
                         0, #Fz_1
                         self.z_sc*step(x,self.geo.x_2)**0,  #Fy_2
                         0,  #Fz_2
                         self.z_sc*step(x,self.geo.x_3)**0,  #Fy_3
                         0, #Fz_3
                         self.a_y*self.z_sc*step(x,self.x_I)**0 + self.a_m*step(x,self.x_I)**0,  #Fa
                         0, #C1
                         0, #C2
                         0, #C3
                         0, #C4
                         0 #C5
                         ])
        

    def v_y_prime(self, x):
        return np.array([-1/6*step(x,self.geo.x_1)**3,    #Fy_1
                         0, #Fz_1
                         -1/6*step(x,self.geo.x_2)**3,   #Fy_2
                         0,  #Fz_2
                         -1/6*step(x,self.geo.x_3)**3,   #Fy_3
                         0, #Fz_3
                         -self.a_y/6*step(x,self.x_I)**3, #Fa
                         0, #C1
                         0, #C2
                         -self.E*self.geo.MMoI[1]*x,    #C3
                         -self.E*self.geo.MMoI[1],  #C4
                         0 #C5
                         ])*-1/(self.E*self.geo.MMoI[1])
        

    def v_z_prime(self, x):
        return np.array([0,    #Fy_1
                         -1/6*step(x,self.geo.x_1)**3, #Fz_1
                         0, #Fy_2
                         -1/6*step(x,self.geo.x_2)**3,  #Fz_2
                         0, #Fy_3
                         -1/6*step(x,self.geo.x_3)**3,  #Fz_3
                         -self.a_z/6*step(x,self.x_I)**3,   #Fa
                         -self.E*self.geo.MMoI[0]*x,    #C1
                         -self.E*self.geo.MMoI[0],  #C2
                         0, #C3
                         0, #C4
                         0, #C5
                         ])*-self.E*self.geo.MMoI[0]
        

    def theta(self, x):
        return np.array([self.z_sc*step(x,self.geo.x_1)**1,    #Fy_1
                         0, #Fz_1
                         self.z_sc*step(x,self.geo.x_2)**1, #Fy_2
                         0, #Fz_2
                         self.z_sc*step(x,self.geo.x_3)**1, #Fy_3
                         0, #Fz_3
                         self.a_y*self.z_sc*step(x,self.x_I)**1 + self.a_m*step(x,self.x_I)**1, #Fa
                         0, #C1
                         0, #C2
                         0, #C3
                         0, #C4
                         self.G*self.geo.J #C5
                         ])*1/(self.G*self.geo.J)
        

    @property
    def A(self):
        rows = [self.V_y_prime(self.geo.l_a),
                self.V_z_prime(self.geo.l_a),
                self.M_y_prime(self.geo.l_a),
                self.M_z_prime(self.geo.l_a),
                self.T(self.geo.l_a),
                self.v_y_prime(self.geo.x_1)*np.cos(self.defl)\
                    -self.v_z_prime(self.geo.x_1)*np.sin(self.defl)\
                        +self.theta(self.geo.x_1)*self.z_sc, #vy(x1)+theta(x1)*z_sc
                self.v_z_prime(self.geo.x_1)*np.cos(self.defl)\
                    +self.v_y_prime(self.geo.x_1)*np.sin(self.defl),  #vz(x1)=0
                self.v_y_prime(self.geo.x_2)*np.cos(self.defl)\
                    -self.v_z_prime(self.geo.x_2)*np.sin(self.defl)\
                        +self.theta(self.geo.x_2)*self.z_sc, #vy(x2)+theta(x2)*z_sc
                self.v_z_prime(self.geo.x_2)*np.cos(self.defl)\
                    +self.v_y_prime(self.geo.x_2)*np.sin(self.defl), #vz(x2)=0
                self.v_y_prime(self.geo.x_3)*np.cos(self.defl)\
                    -self.v_z_prime(self.geo.x_3)*np.sin(self.defl)\
                        +self.theta(self.geo.x_3)*self.z_sc,    #vy(x3)+theta(x3)*z_sc
                self.v_z_prime(self.geo.x_3)*np.cos(self.defl)\
                    +self.v_y_prime(self.geo.x_3)*np.sin(self.defl), #vz(x3)=0
                self.v_z_prime(self.x_I)*np.cos(self.defl)\
                    +self.v_y_prime(self.x_I)*np.sin(self.defl) #vz(xI)=0
                    ]
        return np.vstack(rows)

    @property
    def B(self):
        return np.array([-self.P*self.a_y*step(self.geo.l_a,self.x_II)**0\
                            -self.interp.integrate_q(self.geo.l_a,ord=1)[-1], #V_y_prime
                         -self.P*self.a_z*step(self.geo.l_a,self.x_II)**0,  #V_z_prime
                         self.P*self.a_z*step(self.geo.l_a,self.x_II)**1, #M_y_prime
                         self.P*self.a_y*step(self.geo.l_a,self.x_II)**1\
                             +self.interp.integrate_q(self.geo.l_a,ord=2)[-1],     #M_z_prime 
                         -self.P*self.a_y*self.z_sc*step(self.geo.l_a,self.x_II)**0\
                             -self.P*self.a_m*step(self.geo.l_a,self.x_II)**0\
                                 -self.interp.integrate_tau(self.geo.l_a,self.z_sc,ord=1)[-1], #T
                         self.d_1\
                             +((self.P*self.a_y)/6*step(self.geo.x_1,self.x_II)**3+self.interp.integrate_q(self.geo.x_1,ord=4)[-1])*-np.cos(self.defl)/(self.E*self.geo.MMoI[0])\
                                +(self.P*self.a_z)/6*step(self.geo.x_1,self.x_II)**3*-np.sin(self.defl)/(self.E*self.geo.MMoI[1])\
                                    -(self.P*self.a_y*self.z_sc*step(self.geo.x_1,self.x_II)**1+self.P*self.a_m*step(self.geo.x_1,self.x_II)**1+self.interp.integrate_tau(self.geo.x_1,self.z_sc,ord=2)[-1])*self.z_sc/(self.G*self.geo.J), #vy(x1)+theta(x1)*z_sc
                         (self.P*self.a_z)/6*step(self.geo.x_1,self.x_II)**3*-np.cos(self.defl)/(self.E*self.geo.MMoI[0])\
                             +((self.P*self.a_y)/6*step(self.geo.x_1,self.x_II)**3+self.interp.integrate_q(self.geo.x_1,ord=4)[-1])*-np.sin(self.defl)/(self.E*self.geo.MMoI[1]), #vz(x1)=0
                         ((self.P*self.a_y)/6*step(self.geo.x_2,self.x_II)**3+self.interp.integrate_q(self.geo.x_2,ord=4)[-1])*-np.cos(self.defl)/(self.E*self.geo.MMoI[0])\
                                +(self.P*self.a_z)/6*step(self.geo.x_2,self.x_II)**3*-np.sin(self.defl)/(self.E*self.geo.MMoI[1])\
                                    -(self.P*self.a_y*self.z_sc*step(self.geo.x_2,self.x_II)**1+self.P*self.a_m*step(self.geo.x_2,self.x_II)**1+self.interp.integrate_tau(self.geo.x_2,self.z_sc,ord=2)[-1])*self.z_sc/(self.G*self.geo.J), #vy(x2)+theta(x2)*z_sc
                         (self.P*self.a_z)/6*step(self.geo.x_2,self.x_II)**3*-np.cos(self.defl)/(self.E*self.geo.MMoI[0])\
                             +((self.P*self.a_y)/6*step(self.geo.x_2,self.x_II)**3+self.interp.integrate_q(self.geo.x_2,ord=4)[-1])*-np.sin(self.defl)/(self.E*self.geo.MMoI[1]), #vz(x2)=0
                         self.d_3\
                             +((self.P*self.a_y)/6*step(self.geo.x_3,self.x_II)**3+self.interp.integrate_q(self.geo.x_3,ord=4)[-1])*-np.cos(self.defl)/(self.E*self.geo.MMoI[0])\
                                +(self.P*self.a_z)/6*step(self.geo.x_3,self.x_II)**3*-np.sin(self.defl)/(self.E*self.geo.MMoI[1])\
                                    -(self.P*self.a_y*self.z_sc*step(self.geo.x_3,self.x_II)**1+self.P*self.a_m*step(self.geo.x_3,self.x_II)**1+self.interp.integrate_tau(self.geo.x_3,self.z_sc,ord=2)[-1])*self.z_sc/(self.G*self.geo.J), #vy(x3)+theta(x3)*z_sc
                         (self.P*self.a_z)/6*step(self.geo.x_3,self.x_II)**3*-np.cos(self.defl)/(self.E*self.geo.MMoI[0])\
                             +((self.P*self.a_y)/6*step(self.geo.x_3,self.x_II)**3+self.interp.integrate_q(self.geo.x_3,ord=4)[-1])*-np.sin(self.defl)/(self.E*self.geo.MMoI[1]), #vz(x3)=0
                         (self.P*self.a_z)/6*step(self.x_I,self.x_II)**3*-np.cos(self.defl)/(self.E*self.geo.MMoI[0])\
                             +((self.P*self.a_y)/6*step(self.x_I,self.x_II)**3+self.interp.integrate_q(self.x_I,ord=4)[-1])*-np.sin(self.defl)/(self.E*self.geo.MMoI[1]), #vz(xI)=0
                         ], dtype=np.float64)  
    

if __name__ == "__main__":
    pass