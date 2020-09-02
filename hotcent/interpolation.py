""" Definition of spline functions used by the
AtomicDFT class.

The code below draws heavily from the Hotbit code 
written by Pekka Koskinen (https://github.com/pekkosk/
hotbit/blob/master/box/interpolation.py).
"""
import numpy as np
from scipy.linalg import norm
from scipy.optimize import fminbound, brentq
from scipy.interpolate import splprep, splrep, splev, splint, CubicSpline
try:
    import matplotlib.pyplot as plt
except:
    plt = None


class CubicSplineFunction(CubicSpline):
    def __init__(self, x, y):
        CubicSpline.__init__(self, x, y, bc_type='natural', extrapolate=False)

    def __call__(self, x, der=0):
        if isinstance(x, np.ndarray):
            y = CubicSpline.__call__(self, x, nu=der)
            return np.nan_to_num(y, copy=False)  # NaN -> 0
        else:
            if x < self.x[0] or x > self.x[-1]:
                return 0.
            else:
                return CubicSpline.__call__(self, x, nu=der)


class SplineFunction:
    def __init__(self, x, y, k=3, s=0, name=None):
        """ Simple B-spline function; order k is cubic by default.

        Parameters:
        -----------
        x:  x-grid
        y:  y-values for given x-grid points
        k:  order of spline (cubic by default)
        s:  smoothness parameters (means exact reproduction of x,y values)
        name: name of the function
        """
        if s == -1:
            s = len(x) - np.sqrt(2. * len(x))
        self.tck = splrep(x, y, s=s, k=k)
        self.x = x
        self.y = y
        self.a = x[0]
        self.b = x[-1]
        self.M = len(y)
        self.name=name

    def __call__(self, x, der=0):
        """ Return der'th derivative of f(x)

        Return zero if x beyond the original grid range.
        """
        if isinstance(x, np.ndarray):
            return np.piecewise(x, [x < self.x[0], x > self.x[-1]],
                                [0, 0, lambda x: splev(x, self.tck, der=der)])
        else:
            if x < self.x[0] or x > self.x[-1]:
                return 0.
            else:
                return splev(x, self.tck, der=der)

    def get_name(self):
        """ Return the name of the function. """
        return self.name

    def get_range(self):
        return (self.x[0], self.x[-1])

    def limits(self):
        return self.get_range()

    def solve(self, y, a=None, b=None):
        """ Solve x for f(x)=y, where x in [a, b]. """
        if a is None:
            a = self.x[0]
        if b is None:
            b = self.x[-1]
        assert a < b
        return brentq(lambda x: self(x) - y, a=a, b=b)

    def integrate(self, a, b):
        """ Integrate given function within [a, b]. """
        return splint(a, b, self.tck)

    def plot(self, der=0, filename=None):
        """ Plot the function with matplolib """
        X = np.linspace(self.a, self.b, self.M * 10)
        Y = [self(x, der=der) for x in X]
        plt.plot(X, Y)
        if der == 0:
            plt.scatter(self.x, self.y)

        f1 = SplineFunction(self.x, self.y, k=1)
        plt.plot(X, [f1(x, der=der) for x in X])

        if filename is not None:
            plt.savefig(filename)
        else:
            plt.show()
        plt.clf()

    def max_deviation_from_linear(self):
        """ For given spline (default cubic), return maximum difference
        wrt linear (k=0) interpolation.
        """
        min_dx = np.min(np.array(self.x[1:]) - np.array(self.x[0:-1]))
        M = 10 * (self.b - self.a) / min_dx
        X = np.linspace(self.a, self.b, M)
        f1 = SplineFunction(self.x, self.y, k=1)
        return max([f1(x) - self(x) for x in X])

    def smoothness(self):
        """ Return a measure for the ''smoothness'' of interpolation.

        It is measured by max_deviation_from_linear / average|dy|.
        Smooth interpolations should have <<1.
        If ~1, then interpolation deviates as much as is the variations
        in y and interpolation is not smooth.
        """
        avg = np.average(np.abs(np.array(self.y[1:]) - np.array(self.y[0:-1])))
        return self.max_deviation_from_linear() / avg


class Function:
    def __init__(self, mode, *args, **kwargs):
        assert mode in ['spline', 'cubic_spline']
        if mode == 'spline':
            self.f = SplineFunction(*args, **kwargs)
        elif mode == 'cspline':
            self.f = CubicSplineFunction(*args, **kwargs)

    def __call__(self, x, der=0):
        return self.f(x, der=der)

    def plot(self, der=0, a=None, b=None, npoints=1000, filename=None):
        """ Plot the function with matplolib. """
        a0, b0 = self.f.get_range()
        lower = [a0, a][a is not None]
        upper = [b0, b][b is not None]
        X = np.linspace(lower, upper, npoints)
        Y = [self(x, der=der) for x in X]
        plt.plot(X, Y)
        if filename is not None:
            plt.savefig(filename)
        else:
            plt.show()
        plt.clf()
