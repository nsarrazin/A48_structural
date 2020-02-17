import numpy as np
import matplotlib.pyplot as plt

## NOT FINISHED YET


class intepolation2d:
    def __init__(self):
        self.z = np.zeros(41)
        self.x = np.zeros(81)
        self.f = np.genfromtxt("data/aerodynamicloadf100.dat", delimiter=",")
        self.sigma = 1.0

        # Assigning the z coordinates to the z-array
        for i, zi in enumerate(self.z):
            phi_i = ((i+1)-1)/41 * np.pi
            phi_ii = ((i+2)-1)/41 * np.pi
            self.z[i] = -0.5 * (((0.505/2) * (1 - np.cos(phi_i))) + ((0.505/2) * (1 - np.cos(phi_ii))))

        # Assigning the x coordinates to the x-array
        for i, xi in enumerate(self.x):
            phi_i = ((i + 1) - 1) / 81 * np.pi
            phi_ii = ((i + 2) - 1) / 81 * np.pi
            self.x[i] = 0.5 * (((1.611 / 2) * (1 - np.cos(phi_i))) + ((1.611 / 2) * (1 - np.cos(phi_ii))))

        self.node_array = np.zeros(((81*41), 2))
        for i, xi in enumerate(self.x):
            for j, zi in enumerate(self.z):
                self.node_array[j + (i*41)][0] = xi
                self.node_array[j + (i*41)][1] = zi

    def basisfn(self, x_1, x_2):
        """
        x_1 : vector with coordinates of first point
        x_2 : vector with coordinates of second point
        """
        temp = np.linalg.norm((x_1 - x_2))
        return np.exp(-1* (temp/self.sigma**2))

    def get_w(self):
        """
         nx : the number of points in x direction desired (span-wise)
         nz : the number of points in z direction desired (chord-wise)
        """
        num_points = 81*41
        A = np.zeros((num_points, num_points))
        for i, row in enumerate(A):
            for j, column in enumerate(row):
                print(round((i / num_points * 100), 2), "%")
                A[i, j] = self.basisfn(self.node_array[j], self.node_array[i])

        b = self.f.reshape(num_points)
        w = np.linalg.solve(A, b)
        return w

    def approximatedfn(self, x, weigths):
        """
        :param x: vector with the coordinates
        :return: the value of the interpolated function at x
        """
        f = 0
        for i in range(81*41):
            f += weights[i]*self.basisfn(x, self.node_array[i])
        return f

    def solve(self, nx, nz):
        solve_grid_x = np.linspace(0, 1.611, nx)
        solve_grid_z = np.linspace(0, 0.505, nz)
        solve_node_array = np.zeros(((nx*nz), 2))
        for i, xi in enumerate(solve_grid_x):
            for j, zi in enumerate(solve_grid_z):
                solve_node_array[j + (i*nz)][0] = xi
                solve_node_array[j + (i*nz)][1] = zi

        u = np.zeros(nx*nz)
        weigths = self.get_w()
        for i, ui in enumerate(u):
            u[i] = self.approximatedfn(solve_node_array[i], weigths)

        return u



## Testing implementation
inte = intepolation2d()
u = inte.solve(81, 41).reshape((81, 41))
plt.imshow(u)
plt.show()








