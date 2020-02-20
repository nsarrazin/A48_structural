c_a = 0.505 # chord length aileron [m]
l_a = 1.611 # span of the aileron [m]

E = 72.9E9 # elastic modulus [Pa]
G = 27.1E9 # shear modulus [Pa]

x_1 = 0.125 # x-location of hinge 1 [m]
x_2 = 0.498 # x-location of hinge 2 [m]
x_3 = 1.494 # x-location of hinge 3 [m]

x_a = 24.5E-2 # distance between actuator 1 and 2 [m]

h = 16.1E-2 # aileron height [m]
t_sk = 1.1E-3 # skin thickness [m]
t_sp = 2.4E-3 # spar thickness [m]

t_st = 1.2E-3 # stiffener thickness [m]
h_st = 1.3E-2 # height of stiffener [m]
w_st = 1.7E-2 # width of stiffener [m]
n_st = 11 # number of stiffeners [-]

d_1 = 0.389E-2 # vertical displacement hinge 1 [m]
d_3 = 1.245E-2 # vertical displacement hinge 3 [m]

max_defl = 30 # maximum upward deflecection [deg]
load = 49.2E3 # Load in actuator 2 [N]

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
                   "E": e_mod,
                   "G": g_mod}