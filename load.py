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
        self.E = kwargs.get("e_mod")
        self.G = kwargs.get("g_mod")

        self.a_y = np.sin(self.defl)
        self.a_z = np.cos(self.defl)
        self.a_m = self.geo.h/2 * (np.sin(self.defl)- np.cos(self.defl))
        
        _ , self.z_sc = self.geo.shearcenter
        # _ , self.z_sc = 0, 0
    def V_y(self, x):
        np.array([-step(x, self.geo.x_1)**0,
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
        pass

    def V_z(self,x):
        np.array([0,
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
        pass

    def M_y(self, x):
        np.array([0,
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
        pass

    def M_z(self, x):
        np.array([-step(x,self.geo.x_1)**1,
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
        pass

    def T(self, x):
        np.array([-self.z_sc*step(x,self.geo.x_1)**0,
                  0,
                  -self.z_sc*step(x,self.geo.x_2)**0,
                  0,
                  -self.z_sc*step(x,self.geo.x_3)**0,
                  0,
                  self.a_y*self.z_sc*step(x,self.geo.x_2-self.ge0.x_a/2)**0-self.a_m*step(x,self.geo.x_2-self.geo.x_a/2)**0,
                  0,
                  0,
                  0,
                  0,
                  0, # XXX: add z_sc
                  self.interp.integrate_tau(x, z_sc=self.geo.h/2+self.z_sc, ord=1) + self.P*self.a_y*self.z_sc*step(x,self.geo.x_2+self.geo.x_a/2)**0-self.P*self.a_m*step(x,self.geo.x_2+self.geo.x_a/2)**0,
                  ])
        pass

    def v_y(self, x):
        np.array([0,
                  -1/6*step(x,self.geo.x_1)**3,
                  0,
                  -1/6*step(x,self.geo.x_2)**3,
                  0,
                  -1/6*step(x,self.geo.x_3)**3,
                  self.a_z/6*step(x,self.geo.x_2-self.geo.x_a/2)**3,
                  0,
                  0,
                  x,
                  1,
                  0,
                  self.P*self.a_z/6*step(x,self.geo.x_2+self.geo.x_a/2)**3])*-1/(self.E*self.geo.MMoI[0])
        pass

    def v_z(self, x):
        np.array([-1/6*step(x,self.geo.x_1)**3,
            	  0,
                  -1/6*step(x,self.geo.x_2)**3,
                  0,
                  -1/6*step(x,self.geo.x_3)**3,
                  0,
                  self.a_y/6*step(x,self.geo.x_2-self.geo.x_a/2)**3,
                  x,
                  1,
                  0,
                  0,
                  0,
                  self.interp.integrate_q(x, ord=4) + self.P*self.a_y/6*step(x,self.geo.x_2+self.geo.x_a/2)**3,
                  ])*-1/(self.E*self.geo.MMoI[1])
        pass

    def theta(self, x):
        np.array([-self.z_sc*step(x,self.geo.x_1)**1,
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
                  self.interp.integrate_tau(x, z_sc=self.geo.h/2+self.z_sc, ord=2) + self.P*self.a_y*self.z_sc*step(x,self.geo.x_2+self.geo.x_a/2)**1-self.P*self.a_m*step(x,self.geo.x_2+self.geo.x_a/2)**1,
                  ])*1/(self.G*self.geo.J)
        pass

    @property
    def A(self):
        rows = [V_y(self.geo.l_a),
                V_z(self.geo.l_a),
                M_y(self.geo.l_a),
                M_z(self.geo.l_a),
                T(self.geo.l_a),

                v_y(self.geo.x_1)*np.cos(-self.defl) + v_z(self.geo.x_1)*np.sin(-self.defl) \
                                                        + self.theta(self.geo.x_1)*self.z_sc, # x1 v_y`
                -v_y(self.geo.x_1)*np.sin(-self.defl) + v_z(self.geo.x_1)*np.cos(-self.defl), # x1 v_z`


                v_y(self.geo.x_2)*np.cos(-self.defl) + v_z(self.geo.x_2)*np.sin(-self.defl) \
                                                        + self.theta(self.geo.x_2)*self.z_sc, # x2 v_y`
                -v_y(self.geo.x_2)*np.sin(-self.defl) + v_z(self.geo.x_2)*np.cos(-self.defl), # x2 v_z`


                v_y(self.geo.x_3)*np.cos(-self.defl) + v_z(self.geo.x_3)*np.sin(-self.defl) \
                                                        + self.theta(self.geo.x_3)*self.z_sc, # x3 v_y`
                -v_y(self.geo.x_3)*np.sin(-self.defl) + v_z(self.geo.x_3)*np.cos(-self.defl), # x3 v_z`

                -v_y(self.geo.x_2 - self.geo.x_a/2)*np.sin(-self.defl) + v_z(self.geo.x_2 - self.geo.x_a/2)*np.cos(-self.defl), # actuator I v_z`
                np.array([0,0,0,0,0,0,0,0,0,0,0,0,1])
                ]
        
        return np.vstack(rows)

    @property
    def B(self):
        return np.array([0,0,0,0,0,self.d_1,0,0,0,self.d_3,0,0,1])
    

if __name__ == "__main__":
    pass