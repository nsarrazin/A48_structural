import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt
from consts import *

# from mpl_toolkits.mplot3d.axes3d import Axes3D
mpl.rcParams['figure.dpi'] = 155
data = np.genfromtxt("aerodynamicloadf100.dat", delimiter=",")
x_ticks = np.linspace(0,1.611,41)
y_ticks = np.linspace(0,0.505,81)

x = range(41)
y = range(81)
xs, ys = np.meshgrid(x, y)
# z = calculate_R(xs, ys)
zs = data

fig = plt.figure()
# ax = Axes3D(fig)
cs = plt.imshow(zs, extent = (0., 1.611, 0., 0.505), cmap='Blues')
cbar= plt.colorbar(cs, shrink =0.5)
cbar.ax.set_ylabel("Aerodynamic load [kN]")
plt.xlabel("Span-wise direction [m]")
plt.ylabel("Chord-wise direction [m]")
plt.title("Aerodynamic loading over the aileron")
plt.tight_layout()
plt.show()

print('test')
