import numpy as np

## Technique : bilinear interpolation
## NOT FINISHED

class interpolation:
    def __init__(self):
        self.F = np.genfromtxt("aerodynamicloadf100.dat", delimiter=",")
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
        for i, xi in enumerate(self.x):
            if abs(x) >= abs(xi):
                column = i
            else:
                break
        return row, column

    ## TO IMPLEMENT:
    ## The actual method :)




## Testing implementation
test = interpolation()
print(test.x[0:5])
print(test.z[0:5])

row, column = test.index_search(0.03, -1.3E-4)
print(row, column)