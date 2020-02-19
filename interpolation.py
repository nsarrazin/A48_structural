import numpy as np
import matplotlib.pyplot as plt

## Technique : bilinear interpolation
## NOT FINISHED

class interpolation:
    def __init__(self):
        self.F = np.genfromtxt("data/aerodynamicloadf100.dat", delimiter=",")
        self.x = np.zeros(41)
        self.z = np.zeros(81)

        # Assigning the z coordinates to the z-array
        for i, zi in enumerate(self.z):
            phi_i = ((i+1)-1)/81 * np.pi
            phi_ii = ((i+2)-1)/81 * np.pi
            self.z[i] = -0.5 * (((0.505/2) * (1 - np.cos(phi_i))) + ((0.505/2) * (1 - np.cos(phi_ii))))

        # Assigning the x coordinates to the x-array
        for i, xi in enumerate(self.x):
            phi_i = ((i + 1) - 1) / 41 * np.pi
            phi_ii = ((i + 2) - 1) / 41 * np.pi
            self.x[i] = 0.5 * (((1.611 / 2) * (1 - np.cos(phi_i))) + ((1.611 / 2) * (1 - np.cos(phi_ii))))

    def index_search(self, x, z):
        # searches the beginning index of the interval
        row = 0
        column = 0
        for i, zi in enumerate(self.z):
            if abs(z) >= abs(zi):
                row = i
            else:
                break
        for j, xj in enumerate(self.x):
            if x >= xj:
                column = j
            else:
                break
        return row, column

    def bilinear_interpolation(self, x, z):
        row, column = self.index_search(x, z)
        A = np.matrix([[1, self.x[column]  , self.z[row]  , (self.x[column]*self.z[row])],
                       [1, self.x[column]  , self.z[row+1], (self.x[column]*self.z[row+1])],
                       [1, self.x[column+1], self.z[row]  , (self.x[column+1]*self.z[row])],
                       [1, self.x[column+1], self.z[row+1], (self.x[column+1]*self.z[row+1])]])
        f = np.array([self.F[row, column], self.F[(row+1), column], self.F[row, (column+1)], self.F[(row+1), (column+1)]])
        a = np.linalg.solve(A, f)
        return a[0] + a[1]*x + a[2]*z + a[3]*x*z



## Testing implementation
test = interpolation()
x_test = np.linspace(0, 1.6, 82)
z_test = np.linspace(0, 0.504, 162)
u = np.zeros((len(z_test), len(x_test)))


for row, zi in enumerate(u):
    for column, xi in enumerate(zi):
        u[row,column] = test.bilinear_interpolation(x_test[column], z_test[row])


plt.imshow(u)
plt.show()