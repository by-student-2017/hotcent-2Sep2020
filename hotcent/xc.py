import numpy as np
try:
    from pylibxc import LibXCFunctional
    has_pylibxc = True
except ImportError:
    print('Warning -- could not load LibXC')
    has_pylibxc = False


class LibXC:
    def __init__(self, xcname):
        """ Interface to PyLibXC

        xcname:  str with combination of LibXC functional names,
                 e.g. xcname='GGA_X_PBE+GGA_C_PBE'
        """
        assert has_pylibxc, 'Using XC other than LDA requires PyLibXC!'

        self.xcname = xcname
        self.names = self.xcname.split('+')
        self.functionals = []
        self.types = []
        self.add_gradient_corrections = False

        for name in self.names:
            try:
                self.functionals.append(LibXCFunctional(name, 'unpolarized'))
            except KeyError as err:
                print('KeyError:', err)
                print('>>> Bad XC name. For valid LibXC functional names, see')
                print('>>> https://www.tddft.org/programs/libxc/functionals/')
                raise

            if 'mgga' in name.lower():
                raise ValueError('Meta-GGA functionals not allowed:', name)
            if 'lda' in name.lower():
                self.types.append('LDA')
            elif 'gga' in name.lower():
                self.types.append('GGA')
                self.add_gradient_corrections = True
            else:
                raise ValueError('XC func %s is not LDA or GGA' % name)

    def compute_all(self, rho, sigma):
        inp = {'rho': rho, 'sigma': sigma}
        zk = np.zeros_like(rho)
        vrho = np.zeros_like(rho)
        vsigma = np.zeros_like(rho) if self.add_gradient_corrections else None

        for i, func in enumerate(self.functionals):
            out = func.compute(inp)
            zk += out['zk'][0]
            vrho += out['vrho'][0]
            if self.types[i] == 'GGA':
                vsigma += out['vsigma'][0]

        return {'zk': zk, 'vrho': vrho, 'vsigma': vsigma}

    def evaluate(self, rho, gd):
        """ Returns XC energy and potential given:

        rho:  array-like, the electron density
        gd:   an object that can carry out gradient and divergence
              operations on a grid-based array
        """
        grad = gd.gradient(rho)
        sigma = grad ** 2
        out = self.compute_all(rho, sigma)
        exc = out['zk']
        vxc = out['vrho']
        if self.add_gradient_corrections:
            assert out['vsigma'] is not None
            vxc += -2 * gd.divergence(out['vsigma'] * grad)
        return exc, vxc


class XC_PW92:
    def __init__(self):
        """ The Perdew-Wang 1992 LDA exchange-correlation functional. """
        self.small = 1e-90
        self.a1 = 0.21370
        self.c0 = 0.031091
        self.c1 = 0.046644
        self.b1 = 1.0 / 2.0 / self.c0 * np.exp(-self.c1 / 2.0 / self.c0)
        self.b2 = 2 * self.c0 * self.b1 ** 2
        self.b3 = 1.6382
        self.b4 = 0.49294

    def exc(self, n, der=0):
        """ Exchange-correlation with electron density n. """
        is_scalar = type(n) == np.float64
        n = np.array([n]) if is_scalar else n
        indices = n < self.small
        n[indices] = self.small
        e = self.e_x(n, der=der) + self.e_corr(n, der=der)
        e[indices] = 0.
        if is_scalar:
            e = e[0]
        return e

    def e_x(self, n, der=0):
        """ Exchange. """
        if der == 0:
            return -3. / 4 * (3 * n / np.pi) ** (1. / 3)
        elif der == 1:
            return -3. / (4 * np.pi) * (3 * n / np.pi) ** (-2. / 3)

    def e_corr(self, n, der=0):
        """ Correlation energy. """
        rs = (3. / (4 * np.pi * n)) ** (1. / 3)
        aux = 2 * self.c0
        aux *= self.b1 * np.sqrt(rs) + self.b2 * rs + self.b3 * rs ** (3. / 2) + self.b4 * rs ** 2
        if der == 0:
            return -2 * self.c0 * (1 + self.a1 * rs) * np.log(1 + aux ** -1)
        elif der == 1:
            return (-2 * self.c0 * self.a1 * np.log(1 + aux ** -1) \
                    -2 * self.c0 * (1 + self.a1 * rs) * (1 + aux ** -1) ** -1 * (-aux ** -2) \
                   * 2 * self.c0 * (self.b1 / (2 * np.sqrt(rs)) + self.b2 + 3 * self.b3 * np.sqrt(rs) / 2 \
                   + 2 * self.b4 * rs)) * (-(4 * np.pi * n ** 2 * rs ** 2) ** -1)

    def vxc(self, n):
        """ Exchange-correlation potential (functional derivative of exc). """
        indices = n < self.small
        n[indices] = self.small
        v = self.exc(n) + n * self.exc(n, der=1)
        v[indices] = 0.
        return v

    def evaluate(self, n, *args, **kwargs):
        """ Return the XC energy and potential

            n: array-like, the electron density
        """
        return self.exc(n), self.vxc(n)
