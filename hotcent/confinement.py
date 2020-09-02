""" Definition of confinement potentials """
import numpy as np


class Confinement():
    def __init__(self, adjustable=[]):
        """ Keyword argument:

        adjustable: list of confinement potential parameters that
                    are allowed to be adjusted (e.g. during
                    bandstructure or energy fitting procedures).
        """
        # Check validity of 'adjustable' keyword argument;
        # it is only supposed to contain names of parameters
        # that are actually used by the confinement potential
        msg = self.__str__() + ' does not have %s parameter'
        for parameter in adjustable:
            assert parameter in self.__dir__(), msg % parameter
        self.adjustable = adjustable

    def __call__(self, r):
        raise NotImplementedError
    def __str__(self):
        raise NotImplementedError


class ZeroConfinement(Confinement):
    def __call__(self, r):
        return np.zeros_like(r)

    def __str__(self):
        return 'ZeroConfinement'


class PowerConfinement(Confinement):
    def __init__(self, r0=1., s=2, adjustable=[]):
        self.r0 = r0
        self.s = s
        Confinement.__init__(self, adjustable=adjustable)

    def __call__(self, r):
        return (r / self.r0) ** self.s

    def __str__(self):
        return 'PowerConfinement(r0=%.6f, s=%.6f)' % (self.r0, self.s)


class WoodsSaxonConfinement(Confinement):
    def __init__(self, w=1., r0=1., a=1., adjustable=[]):
        self.w = w
        self.r0 = r0
        self.a = a
        Confinement.__init__(self, adjustable=adjustable)

    def __call__(self, r):
        return self.w / (1 + np.exp(self.a * (self.r0 - r)))

    def __str__(self):
        return 'WoodsSaxonConfinement(w=%.6f, r0=%.6f, a=%.6f)' % \
               (self.w, self.r0, self.a)


class SoftConfinement(Confinement):
    """ As in Junquera et al. PRB 64, 23511 (2001). """
    def __init__(self, amp=12., rc=0., x_ri=0.6, adjustable=[]):
        self.amp = amp
        self.rc = rc
        self.x_ri = x_ri
        Confinement.__init__(self, adjustable=adjustable)

    def __call__(self, r):
        ri = self.x_ri * self.rc
        condlist = [r < ri, np.logical_and(ri <= r, r < self.rc), r >= self.rc]
        funclist = [0., lambda x: self.amp * np.exp(-(self.rc - ri) / (x - ri)) \
                                  / (self.rc - x), lambda x: 0., np.inf]
        return np.piecewise(r, condlist, funclist)

    def __str__(self):
        return 'SoftConfinement(amp=%.6f, rc=%.6f, x_ri=%.6f)' % \
               (self.amp, self.rc, self.x_ri)
