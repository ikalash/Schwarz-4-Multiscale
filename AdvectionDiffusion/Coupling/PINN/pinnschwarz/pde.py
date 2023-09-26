import tensorflow as tf
import numpy as np


tf.keras.backend.set_floatx('float64')


class PDE_1D_Steady_AdvecDiff:
    def __init__(self, nu=0.1, beta=1, n_FD=2**10, order=1):
        # Store parameters of PDE
        self.nu = nu
        self.beta = beta
        self.order = order

        xl = 0
        xr = 1
        x = np.linspace(xl, xr, num=n_FD)
        h = x[1] - x[0]

        u = np.array([[0], [0]])

        a = - nu/(h**2)
        b = (2*nu)/(h**2)
        c = -(nu/(h**2))

        if order == 1:
            b += beta / h
            c += -beta / h
            d = 0.0
        elif order == 2:
            b += 3 / 2 * beta / h
            c += -2 * beta / h
            d = 1 / 2 * beta / h

        else:
            raise ValueError(f"Invalid order: {order}")

        self.coeff = (a, b, c, d)

        self.A = np.diagflat([b]*(n_FD-2)) + np.diagflat([c]*(n_FD - 3), -1) + np.diagflat([a]*(n_FD - 3), 1) + np.diagflat([d] * (n_FD - 4), -2)

        f = np.ones((n_FD-2, 1))

        u_int = np.linalg.solve( self.A, f )

        u = np.hstack((np.hstack((u[0], u_int.flatten())), u[-1]))

        self.x_FD = x
        self.u_FD = u

    # Define boundary condition
    def f_b(self, x):
        return tf.zeros((x.shape[0],1), dtype = 'float64')

    # Define residual of the PDE
    def f_r(self, u_x, u_xx):
        return -self.nu*u_xx + self.beta*u_x - 1

    # Define analytical solution to PDE
    def f(self, x):
        return np.interp(x, self.x_FD, self.u_FD)

