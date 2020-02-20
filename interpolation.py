import numpy as np
import matplotlib.pyplot as plt
from data.consts import N

"""
Technique   : Bi-linear interpolation
Author      : Group A48, AE year 2020
Date        : February 2020
"""

class Interpolation:
    def __init__(self):
        # Generating the x, z, F arrays
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
        """
        Searches the first index of the square 4 point domain in which point (x, z) lies.
        :param x: scalar; x-position of the point
        :param z: scalar; z-position of the point
        :return: 2 scalars; the row and column index (z, x index) of the first point in the square 4 point domain.
        """
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
        """
        :param x: scalar; x-position of the point
        :param z: scalar; z-position of the point
        :return: scalar, array; the value of q at (x, z) , the weight coefficients of the bi-linear interpolation
        """
        row, column = self.index_search(x, z)
        # If point (x,z) is lies on the right boundary
        if row == 80 or column == 40:
            A = np.matrix([[1, self.x[column-1], self.z[row-1], (self.x[column-1] * self.z[row-1])],
                           [1, self.x[column-1], self.z[row], (self.x[column-1] * self.z[row])],
                           [1, self.x[column], self.z[row-1], (self.x[column] * self.z[row-1])],
                           [1, self.x[column], self.z[row], (self.x[column] * self.z[row])]])
            f = np.array([self.F[row-1, column-1], self.F[row, column-1], self.F[row-1, column],
                          self.F[row, column]])
            a = np.linalg.solve(A, f)
        # Anywhere else in the domain
        else:
            A = np.matrix([[1, self.x[column]  , self.z[row]  , (self.x[column]*self.z[row])],
                           [1, self.x[column]  , self.z[row+1], (self.x[column]*self.z[row+1])],
                           [1, self.x[column+1], self.z[row]  , (self.x[column+1]*self.z[row])],
                           [1, self.x[column+1], self.z[row+1], (self.x[column+1]*self.z[row+1])]])
            f = np.array([self.F[row, column], self.F[(row+1), column], self.F[row, (column+1)], self.F[(row+1), (column+1)]])
            a = np.linalg.solve(A, f)
        return (a[0] + a[1]*x + a[2]*z + a[3]*x*z), a

    def q_intergration_fixed_x(self, x_fixed):
        """
        Integrates the curve at a fixed x location along the chord (z-direction) using an analytical solution to the bilinear integration
        :param x_fixed: scalar; fixed span-wise location (x-position)
        :return: scalar; integrated value along the z-direction
        """
        z_integration = np.zeros(80)
        for i, zi in enumerate(self.z):
            if i == 80:
                break
            z_integration[i] = (self.z[i] + self.z[i+1])/2

        a_list = []
        for i, zi in enumerate(z_integration):
            func, a = self.bilinear_interpolation(x_fixed, zi)
            a_list.append(a)

        result = 0
        for i, zi in enumerate(self.z):
            if i == 80:
                break
            z_begin = self.z[i]
            z_end = self.z[i+1]
            result = result + (((a_list[i][0]*z_end) + a_list[i][1]*x_fixed*z_end + (a_list[i][2]/2)*z_end**2 + (a_list[i][3]/2)*x_fixed*z_end**2) - ((a_list[i][0]*z_begin) + a_list[i][1]*x_fixed*z_begin + (a_list[i][2]/2)*z_begin**2 + (a_list[i][3]/2)*x_fixed*z_begin**2 ))

        return result

    def q_intergrate_double(self, x_begin, x_end, dx):
        """
        :param x_begin: scalar; starting value of the integration (span-wise location)
        :param x_end: scalar; ending value of the integration (span-wise location)
        :param dx: scalar; integration step-size in x-direction
        :return: scalar; integrated value of the integrated values along z, along x (double integral)
        """
        x_intergrate = np.arange(x_begin, (x_end+dx), dx)
        result = 0
        for i, xi in enumerate(x_intergrate):
            if i == (len(x_intergrate) - 1):
                break
            result = result + ((x_intergrate[i+1] - x_intergrate[i])/2) * (self.q_intergration_fixed_x(x_intergrate[i]) + self.q_intergration_fixed_x(x_intergrate[i+1]))
        return result

    def trapezoidalrule(self, y, x):
        """
        :param y: array; array to integrate
        :param x: array; array with the x-values associated with the y-values
        :return: array; value of the integral on each x location + the value of the previous x location
        """
        if len(y) == len(x):
            result = np.zeros(len(x))
            result[0] = 0
            for i, xi in enumerate(x):
                if i == (len(y) - 1):
                    break
                temp = (x[i+1] - x[i]) * ((y[i] + y[i+1])/2)
                result[i+1] = result[i] + temp
        else:
            print("Trapezoidalrule: x and y do not have the same size", ", len(x) =", len(x), ", len(y) =", len(y))
            print("Returning None")
            result = None
        return result

    def integrate_q(self, x, ord=2):
        xs = np.linspace(0, x, N)
        ys = np.array([self.q_intergration_fixed_x(x) for x in xs])
        
        for i in range(ord):
            ys = self.trapezoidalrule(ys, xs)

        return ys

    def tau(self, x_fixed, z_sc):
        """
        blabla this is tau #TODO: add documentation
        """
        z_integration = np.zeros(80)
        for i, zi in enumerate(self.z):
            if i == 80:
                break
            z_integration[i] = (self.z[i] + self.z[i+1])/2

        a_list = []
        for i, zi in enumerate(z_integration):
            func, a = self.bilinear_interpolation(x_fixed, zi)
            a_list.append(a)

        result = 0
        for i, zi in enumerate(self.z):
            if i == 80:
                break
            z_begin = self.z[i]
            z_end = self.z[i+1]

            z_midpoint = z_integration[i]
            val = (((a_list[i][0]*z_end) + a_list[i][1]*x_fixed*z_end + (a_list[i][2]/2)*z_end**2 + (a_list[i][3]/2)*x_fixed*z_end**2) - \
                 ((a_list[i][0]*z_begin) + a_list[i][1]*x_fixed*z_begin + (a_list[i][2]/2)*z_begin**2 + (a_list[i][3]/2)*x_fixed*z_begin**2 )) 

            result += val * (z_midpoint - z_sc)

        return result

    def integrate_tau(self, x, z_sc, ord=2):
        xs = np.linspace(0, x, N)
        ys = np.array([self.tau(x, z_sc) for x in xs])
        
        for i in range(ord):
            ys = self.trapezoidalrule(ys, xs)

        return ys


if __name__ == "__main__":
        
    ## Testing implementations
    test = Interpolation()
    # x_test = np.linspace(0, 1.611, 82)
    # z_test = np.linspace(0, -0.505, 162)

    # y = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    # x = np.linspace(0, 10, 11)
    # plt.plot(x, y)
    # plt.show()
    # plt.plot(x, test.trapezoidalrule(y, x))
    # plt.show()
    # print(test.trapezoidalrule(y, x))



    for i in range(0, 4):
        # plt.plot(test.integrate_q(1.611, ord=i), label=i)

        plt.plot(test.integrate_tau(1.611, -0.225-16.1E-2/2, ord=i), label=i)
    plt.legend()
    plt.show()

    ## Uncomment for full 2D plot
    # u = np.zeros((len(z_test), len(x_test)))
    # for row, zi in enumerate(u):
    #     for column, xi in enumerate(zi):
    #         func, a = test.bilinear_interpolation(x_test[column], z_test[row])
    #         u[row,column] = func
    # plt.imshow(u)
    # plt.show()

    ## Uncomment for x-direction crossection
    # u_cross = np.zeros(len(x_test))
    # for i, ui in enumerate(u_cross):
    #     func, a = test.bilinear_interpolation(x_test[i], test.z[40])
    #     u_cross[i] = func
    # plt.plot(test.x, test.F[40, :], "-")
    # plt.plot(x_test, u_cross, "-" )
    # plt.show()

    ## Uncomment for z-direction crossection
    # u_cross = np.zeros(len(z_test))
    # for i, ui in enumerate(u_cross):
    #     func, a = test.bilinear_interpolation(test.x[21], z_test[i])
    #     u_cross[i] = func
    # plt.plot(test.z, test.F[:, 21], "-")
    # plt.plot(z_test, u_cross, "-" )
    # plt.show()

    ## Uncomment to test integration
    # u_cross = np.zeros(len(z_test))
    # for i, ui in enumerate(u_cross):
    #     func, a = test.bilinear_interpolation(test.x[20], z_test[i])
    #     u_cross[i] = func
    # trapz_result = np.trapz(u_cross, z_test)
    # result = test.q_intergration_fixed_x(test.x[20])
    # print("Numpy trapezoidal implementation:", trapz_result)
    # print("Own analytical implementation   :", result)
    # print("Difference                      :", abs(trapz_result-result))

    ## Uncomment to test double integral
    # result = test.q_intergrate_double(1, 1.5, 0.01)
    # print(result)

