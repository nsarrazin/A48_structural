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
    def V_y(self, x):
        return np.array([-step(x, self.geo.x_1)**0,
                  0,
                  -step(x, self.geo.x_2)**0,
                  0, 
                  -step(x, self.geo.x_3)**0,
                  0,
                  self.a_y*step(x, self.geo.x_2-self.geo.x_a/2)**0,
                  0,
                  0,
                  0,
                  0,
                  0,
                  self.interp.integrate_q(x, ord=1)[-1]+self.P*self.a_y*step(x,self.geo.x_2-self.geo.x_a/2)**0
                  ])

    def V_z(self,x):
        return np.array([0,
                 -step(x,self.geo.x_1)**0,
                 0,
                 -step(x,self.geo.x_2)**0,
                 0,
                 -step(x,self.geo.x_3)**0,
                 self.a_z*step(x,self.geo.x_2-self.geo.x_a/2)**0,
                 0,
                 0,
                 0,
                 0,
                 0,
                 self.P*self.a_z*step(x,self.geo.x_2+self.geo.x_a/2)**0])

    def M_y(self, x):
        return np.array([0,
                  -step(x,self.geo.x_1)**1,
                  0,
                  -step(x,self.geo.x_2)**1,
                  0,
                  -step(x,self.geo.x_3)**1,
                  self.a_z*step(x,self.geo.x_2-self.geo.x_a/2)**1,
                  0,
                  0,
                  0,
                  0,
                  0,
                  self.P*self.a_z*step(x,self.geo.x_2+self.geo.x_a/2)**1])

    def M_z(self, x):
        return np.array([-step(x,self.geo.x_1)**1,
                  0,
                  -step(x,self.geo.x_2)**1,
                  0,
                  -step(x,self.geo.x_3)**1,
                  0,
                  self.a_y*step(x,self.geo.x_2-self.geo.x_a/2)**1,
                  0,
                  0,
                  0,
                  0,
                  0,
                  self.interp.integrate_q(x, ord=2)[-1] + self.P*self.a_y*step(x,self.geo.x_2+self.geo.x_a/2)**1
                  ])

    def T(self, x):
        return np.array([-self.z_sc*step(x,self.geo.x_1)**0,
                  0,
                  -self.z_sc*step(x,self.geo.x_2)**0,
                  0,
                  -self.z_sc*step(x,self.geo.x_3)**0,
                  0,
                  self.a_y*self.z_sc*step(x,self.geo.x_2-self.geo.x_a/2)**0-self.a_m*step(x,self.geo.x_2-self.geo.x_a/2)**0,
                  0,
                  0,
                  0,
                  0,
                  0,
                  self.interp.integrate_tau(x, z_sc=self.geo.h/2+self.z_sc, ord=1)[-1] + self.P*self.a_y*self.z_sc*step(x,self.geo.x_2+self.geo.x_a/2)**0-self.P*self.a_m*step(x,self.geo.x_2+self.geo.x_a/2)**0,
                  ])
        

    def v_y(self, x):
        return np.array([0,                         #Fy_1
                  -1/6*step(x,self.geo.x_1)**3,     #Fz_1
                  0,                                #Fy_2
                  -1/6*step(x,self.geo.x_2)**3,     #Fz_2
                  0,                                #Fy_3
                  -1/6*step(x,self.geo.x_3)**3,     #Fz_3
                  self.a_z/6*step(x,self.geo.x_2-self.geo.x_a/2)**3,    #F_a
                  0,                                #C1
                  0,                                #C2
                  x,                                #C3
                  1,                                #C4
                  0,                                #C5
                  (self.P*self.a_z)/6*step(x,self.geo.x_2+self.geo.x_a/2)**3])*-1/(self.E*self.geo.MMoI[0]) #const*factor
        

    def v_z(self, x):
        return np.array([-1/6*step(x,self.geo.x_1)**3,   #Fy_1
            	  0,                                     #Fz_1
                  -1/6*step(x,self.geo.x_2)**3,          #Fy_2
                  0,                                     #Fz_2
                  -1/6*step(x,self.geo.x_3)**3,          #Fy_3   
                  0,                                     #Fz_3
                  self.a_y/6*step(x,self.geo.x_2-self.geo.x_a/2)**3,  #Fa
                  x,    #C1
                  1,    #C2 
                  0,    #C3
                  0,    #C4
                  0,    #C5 
                  self.interp.integrate_q(x, ord=4)[-1] + (self.P*self.a_y)/6*step(x,self.geo.x_2+self.geo.x_a/2)**3,
                  ])*-1/(self.E*self.geo.MMoI[1])
        

    def theta(self, x):
        return np.array([-self.z_sc*step(x,self.geo.x_1)**1,
                  0,
                  -self.z_sc*step(x,self.geo.x_2)**1,
                  0,
                  -self.z_sc*step(x,self.geo.x_3)**1,
                  0,
                  self.a_y*self.z_sc*step(x,self.geo.x_2-self.geo.x_a/2)**1-self.a_m*step(x,self.geo.x_2-self.geo.x_a/2)**1,
                  0,
                  0,
                  0,
                  0,
                  1,
                  self.interp.integrate_tau(x, z_sc=self.geo.h/2+self.z_sc, ord=2)[-1] + self.P*self.a_y*self.z_sc*step(x,self.geo.x_2+self.geo.x_a/2)**1-self.P*self.a_m*step(x,self.geo.x_2+self.geo.x_a/2)**1,
                  ])*1/(self.G*self.geo.J)
        

    @property
    def A(self):
        rows = [self.V_y(self.geo.l_a),
                self.V_z(self.geo.l_a),
                self.M_y(self.geo.l_a),
                self.M_z(self.geo.l_a),
                self.T(self.geo.l_a),

                self.v_y(self.geo.x_1)*np.cos(-self.defl) + self.v_z(self.geo.x_1)*np.sin(-self.defl) \
                                                        + self.theta(self.geo.x_1)*self.z_sc, # x1 v_y`
                -self.v_y(self.geo.x_1)*np.sin(-self.defl) + self.v_z(self.geo.x_1)*np.cos(-self.defl), # x1 v_z`


                self.v_y(self.geo.x_2)*np.cos(-self.defl) + self.v_z(self.geo.x_2)*np.sin(-self.defl) \
                                                        + self.theta(self.geo.x_2)*self.z_sc, # x2 v_y`
                -self.v_y(self.geo.x_2)*np.sin(-self.defl) + self.v_z(self.geo.x_2)*np.cos(-self.defl), # x2 v_z`


                self.v_y(self.geo.x_3)*np.cos(-self.defl) + self.v_z(self.geo.x_3)*np.sin(-self.defl) \
                                                        + self.theta(self.geo.x_3)*self.z_sc, # x3 v_y`
                -self.v_y(self.geo.x_3)*np.sin(-self.defl) + self.v_z(self.geo.x_3)*np.cos(-self.defl), # x3 v_z`

                -self.v_y(self.geo.x_2 - self.geo.x_a/2)*np.sin(-self.defl) + self.v_z(self.geo.x_2 - self.geo.x_a/2)*np.cos(-self.defl), # actuator I v_z`
                np.array([0,0,0,0,0,0,0,0,0,0,0,0,1], dtype=np.float64)
                ]
        return np.vstack(rows)

    @property
    def B(self):
        return np.array([0, #shear
                         0, #shear
                         0, #moment
                         0, #moment 
                         0, #torque
                         self.d_1, #x1 v_y`
                         0, #x1 v_z`
                         0, #x2 v_y`
                         0, # x2 v_z`
                         self.d_3, # x3 v_y`
                         0, # x3 v_z`
                         0, #actuator
                         1], dtype=np.float64)  #const
    

if __name__ == "__main__":
    pass