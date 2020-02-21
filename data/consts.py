c_a = 0.505 # chord length aileron [m]
l_a = 1.611 # span of the aileron [m]

E = 72.9E9 # elastic modulus [Pa]
G = 27.1E9 # shear modulus [Pa]

x_1 = 0.125 # x-location of hinge 1 [m]
x_2 = 0.498 # x-location of hinge 2 [m]
x_3 = 1.494 # x-location of hinge 3 [m]

x_a = 0.245 # distance between actuator 1 and 2 [m]

h = 0.161 # aileron height [m]
t_sk = 0.0011 # skin thickness [m]
t_sp = 0.0024 # spar thickness [m]

t_st = 0.0012  # stiffener thickness [m]
h_st = 0.013 # height of stiffener [m]
w_st = 0.017 # width of stiffener [m]
n_st = 11 # number of stiffeners [-]

d_1 = 0.00389 # vertical displacement hinge 1 [m]
d_3 = 0.01245 # vertical displacement hinge 3 [m]

max_defl = 30 # maximum upward deflecection [deg]
load = 49200 # Load in actuator 2 [N]

N = 250 #steps for integration spanwise 

parameters_geometry = {"c_a" : c_a,
                       "l_a" : l_a,

                       "x_1" : x_1,
                       "x_2" : x_2,
                       "x_3" : x_3, 
                       "x_a" : x_a,

                       "h"   : h,
                       "t_sk": t_sk,
                       "t_sp": t_sp,

                       "t_st": t_st,
                       "h_st": h_st,
                       "w_st": w_st,
                       "n_st": n_st}


parameters_case = {"d_1" : d_1,
                   "d_3" : d_3,
                   "defl": max_defl,
                   "load": load,
                   "E": E,
                   "G": G}