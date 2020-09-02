""" Definition of the AtomicDFT class for atomic
DFT calculations.

The code below draws heavily from the Hotbit code
written by Pekka Koskinen (https://github.com/pekkosk/
hotbit/blob/master/hotbit/parametrization/atom.py).
"""
import pickle
import numpy as np
from ase.data import atomic_numbers, covalent_radii
from ase.units import Bohr
from hotcent.interpolation import CubicSplineFunction
from hotcent.atomic_base import AtomicBase, nl2tuple
from hotcent.confinement import ZeroConfinement
from hotcent.xc import XC_PW92, LibXC
try:
    import matplotlib.pyplot as plt
except:
    plt = None
try:
    import _hotcent
except ModuleNotFoundError:
    print('Warning: C-extensions not available')
    from hotcent.shoot import shoot
    _hotcent = None


class AtomicDFT(AtomicBase):
    def __init__(self,
                 symbol,
                 xc='LDA',
                 convergence={'density':1e-7, 'energies':1e-7},
                 perturbative_confinement=False,
                 **kwargs):
        """ Run Kohn-Sham all-electron calculations for a given atom.

        Example:
        ---------
        from hotcent.atomic_dft import AtomicDFT
        from hotcent.confinement import PowerConfinement
        atom = AtomicDFT('C',
                         xc='GGA_C_PBE+GGA_X_PBE',
                         configuration='[He] 2s2 2p2',
                         valence=['2s', '2p'],
                         confinement=PowerConfinement(r0=3.0, s=2))
        atom.run()

        Parameters:
        -----------
        xc: Name of the XC functional. If 'LDA' or 'PW92' are provided,
            then Hotcent's native LDA implementation will be used.
            For all other functionals, the PyLibXC module is required,
            which is bundled with LibXC.
            The names of the implemented functionals can be found
            on https://www.tddft.org/programs/libxc/functionals/
            Often one needs to combine different LibXC functionals, e.g.
                xc='GGA_X_PBE+GGA_C_PBE'  # for PBE XC

        convergence: convergence criterion dictionary
                        * density: max change for integrated |n_old-n_new|
                        * energies: max change in single-particle energy (Ha)

        perturbative_confinement: determines which type of self-
            consistent calculation is performed when applying each
            of the orbital- or density-confinement potentials:

            False: apply the confinement potential in a conventional
                  calculation with self-consistency between
                  the density and the effective potential,

            True: add the confinement potential to the effective
                  potential of the free (nonconfined) atom and
                  solve for the eigenstate(s)* while keeping this
                  potential fixed.

            * i.e. all valence orbitals when confining the density and
            only the orbital in question in wave function confinement

            The perturbative scheme is e.g. how basis sets are
            generated in GPAW. This option is also faster than the
            self-consistent one, in particular for heavier atoms.
        """
        AtomicBase.__init__(self, symbol, **kwargs)

        print('*******************************************', file=self.txt)
        print('Kohn-Sham all-electron calculator for %s' % self.symbol,
              file=self.txt)
        print('*******************************************', file=self.txt)

        if xc in ['PW92', 'LDA']:
            self.xc = XC_PW92()
        else:
            self.xc = LibXC(xc)

        self.convergence = convergence
        self.perturbative_confinement = perturbative_confinement

        if self.scalarrel:
            print('Using scalar relativistic corrections.', file=self.txt)

        maxnodes = max([n - l - 1 for n, l, nl in self.list_states()])
        self.rmin = 1e-2 / self.Z
        self.N = (maxnodes + 1) * self.nodegpts
        print('max %i nodes, %i grid points' % (maxnodes, self.N),
              file=self.txt)

        self.xgrid = np.linspace(0, np.log(self.rmax / self.rmin), self.N)
        self.rgrid = self.rmin * np.exp(self.xgrid)
        self.grid = RadialGrid(self.rgrid)
        self.timer.stop('init')

    def calculate_energies(self, enl, dens, echo='valence'):
        """ Calculate energy contributions. """
        self.timer.start('energies')
        assert echo in [None, 'valence', 'all']

        self.bs_energy = 0.0
        for n, l, nl in self.list_states():
            self.bs_energy += self.configuration[nl] * enl[nl]

        vhar = self.calculate_hartree_potential(dens)
        self.vhar_energy = 0.5 * self.grid.integrate(vhar * dens, use_dV=True)

        exc, vxc = self.xc.evaluate(dens, self.grid)
        self.vxc_energy = self.grid.integrate(vxc * dens, use_dV=True)
        self.exc_energy = self.grid.integrate(exc * dens, use_dV=True)

        vconf = self.confinement(self.rgrid)
        self.confinement_energy = self.grid.integrate(vconf * dens,
                                                      use_dV=True)

        self.total_energy = self.bs_energy - self.vhar_energy
        self.total_energy += -self.vxc_energy + self.exc_energy

        if echo is not None:
            line = '%s orbital eigenvalues:' % echo
            print('\n'+line, file=self.txt)
            print('-' * len(line), file=self.txt)
            for n, l, nl in self.list_states():
                if echo == 'all' or nl in self.valence:
                    print('  %s:   %.12f' % (nl, enl[nl]), file=self.txt)

            print('\nenergy contributions:', file=self.txt)
            print('----------------------------------------', file=self.txt)
            print('sum of eigenvalues:     %.12f' % self.bs_energy,
                  file=self.txt)
            print('Hartree energy:         %.12f' % self.vhar_energy,
                  file=self.txt)
            print('vxc correction:         %.12f' % self.vxc_energy,
                  file=self.txt)
            print('exchange + corr energy: %.12f' % self.exc_energy,
                  file=self.txt)
            print('----------------------------------------', file=self.txt)
            print('total energy:           %.12f\n' % self.total_energy,
                  file=self.txt)

        self.timer.stop('energies')

    def calculate_density(self, unlg):
        """ Calculate the radial electron density:
        sum_nl occ_nl |Rnl(r)|**2 / (4*pi)
        """
        self.timer.start('density')
        dens = np.zeros_like(self.rgrid)
        for n, l, nl in self.list_states():
            dens += self.configuration[nl] * (unlg[nl] ** 2)

        nel1 = self.grid.integrate(dens)
        nel2 = self.get_number_of_electrons()

        if abs(nel1 - nel2) > 1e-10:
            err = 'Integrated density %.3g' % nel1
            err += ', number of electrons %.3g' % nel2
            raise RuntimeError(err)

        dens = dens / (4 * np.pi * self.rgrid ** 2)

        self.timer.stop('density')
        return dens

    def calculate_hartree_potential(self, dens):
        """ Calculate the Hartree potential. """
        self.timer.start('Hartree')
        dV = self.grid.get_dvolumes()
        r, r0 = self.rgrid, self.grid.get_r0grid()
        N = self.N
        n0 = 0.5 * (dens[1:] + dens[:-1])
        nel = self.get_number_of_electrons()
        n0 *= nel / np.sum(n0 * dV)

        if _hotcent is not None:
            vhar = _hotcent.hartree(n0, dV, r, r0, N)
        else:
            lo, hi, vhar = np.zeros(N), np.zeros(N), np.zeros(N)
            lo[0] = 0.0
            for i in range(1, N):
                lo[i] = lo[i-1] + dV[i-1] * n0[i-1]

            hi[-1] = 0.0
            for i in range(N - 2, -1, -1):
                hi[i] = hi[i + 1] + n0[i] * dV[i] / r0[i]

            for i in range(N):
                vhar[i] = lo[i] / r[i] + hi[i]

        self.timer.stop('Hartree')
        return vhar

    def calculate_veff(self, dens):
        """ Calculate effective potential. """
        self.timer.start('veff')
        vnuc = self.nuclear_potential(self.rgrid)
        vhar = self.calculate_hartree_potential(dens)
        exc, vxc = self.xc.evaluate(dens, self.grid)
        vconf = self.confinement(self.rgrid)
        self.timer.stop('veff')
        return vnuc + vhar + vxc + vconf

    def guess_density(self):
        """ Guess initial density. """
        r2 = 0.02 * self.Z  # radius at which density has dropped to half;
                            # can this be improved?
        dens = np.exp(-self.rgrid / (r2 / np.log(2)))
        nel = self.get_number_of_electrons()
        dens = dens / self.grid.integrate(dens, use_dV=True) * nel
        return dens

    def run(self, write=None):
        """ Execute the required atomic DFT calculations

        Parameters:

        write: None or a filename for saving the rgrid, effective
               potential and electron density.
        """
        def header(*args):
            print('=' * 50, file=self.txt)
            print('\n'.join(args), file=self.txt)
            print('=' * 50, file=self.txt)

        val = self.get_valence_orbitals()
        confinement = self.confinement

        assert all([nl in val for nl in self.wf_confinement])
        nl_x = [nl for nl in val if nl not in self.wf_confinement]
        assert len(nl_x) == 0 or len(nl_x) == len(val), nl_x

        self.enl = {}
        self.unlg = {}
        self.Rnlg = {}
        self.unl_fct = {nl: None for nl in self.configuration}
        self.Rnl_fct = {nl: None for nl in self.configuration}

        if self.perturbative_confinement:
            self.confinement = ZeroConfinement()
            header('Initial run without any confinement',
                   'for pre-converging orbitals and eigenvalues')
            dens_free, veff_free, enl_free, unlg_free, Rnlg_free = \
                                                               self.outer_scf()

        for nl, wf_confinement in self.wf_confinement.items():
            self.confinement = wf_confinement
            if self.confinement is None:
                self.confinement = ZeroConfinement()
            header('Applying %s' % self.confinement,
                   'to get a confined %s orbital' % nl)

            if self.perturbative_confinement:
                veff = veff_free + self.confinement(self.rgrid)
                enl = {nl: enl_free[nl]}
                itmax, enl, d_enl, unlg, Rnlg = self.inner_scf(0, veff, enl,
                                                               {}, solve=[nl])
                print('Confined %s eigenvalue: %.6f' % (nl, enl[nl]),
                      file=self.txt)
            else:
                dens, veff, enl, unlg, Rnlg = self.outer_scf()

            self.enl[nl] = enl[nl]
            self.unlg[nl] = unlg[nl]
            self.Rnlg[nl] = Rnlg[nl]

        self.confinement = confinement
        if self.confinement is None:
            self.confinement = ZeroConfinement()
        extra = '' if len(nl_x) == 0 else '\nand the confined %s orbital(s)' \
                                           % ' and '.join(nl_x)
        header('Applying %s' % self.confinement,
               'to get the confined electron density%s' % extra)

        if self.perturbative_confinement:
            veff = veff_free + self.confinement(self.rgrid)
            enl = {nl_: enl_free[nl_] for nl_ in val}
            itmax, enl, d_enl, unlg, Rnlg = self.inner_scf(0, veff, enl, {},
                                                           solve=val)
        else:
            self.dens, veff, enl, unlg, Rnlg = self.outer_scf()

        for n, l, nl in self.list_states():
            if nl not in self.wf_confinement:
                assert nl in enl or self.perturbative_confinement
                self.enl[nl] = enl[nl] if nl in enl else enl_free[nl]
                self.unlg[nl] = unlg[nl] if nl in enl else unlg_free[nl]
                self.Rnlg[nl] = Rnlg[nl] if nl in enl else Rnlg_free[nl]

        if self.perturbative_confinement:
            self.dens = self.calculate_density(self.unlg)

        self.veff = self.calculate_veff(self.dens)
        self.vhar = self.calculate_hartree_potential(self.dens)
        exc, self.vxc = self.xc.evaluate(self.dens, self.grid)

        if write is not None:
            with open(write, 'w') as f:
                pickle.dump(self.rgrid, f)
                pickle.dump(self.veff, f)
                pickle.dump(self.dens, f)

        self.solved = True
        self.timer.summary()
        self.txt.flush()

    def outer_scf(self):
        """ Solve the self-consistent potential. """
        self.timer.start('outer_scf')
        print('\nStart iteration...', file=self.txt)
        enl = {nl: 0. for n, l, nl in self.list_states()}
        d_enl = {nl: 0. for n, l, nl in self.list_states()}

        dens = self.guess_density()
        veff = self.nuclear_potential(self.rgrid)
        veff += self.confinement(self.rgrid)

        for it in range(self.maxiter):
            veff *= 1. - self.mix
            veff += self.mix * self.calculate_veff(dens)

            dveff = None
            if self.scalarrel:
                spl = CubicSplineFunction(self.rgrid, veff)
                dveff = spl(self.rgrid, der=1)

            itmax, enl, d_enl, unlg, Rnlg = self.inner_scf(it, veff, enl, d_enl,
                                                           dveff=dveff)
            dens0 = dens.copy()
            dens = self.calculate_density(unlg)
            diff = self.grid.integrate(np.abs(dens - dens0), use_dV=True)

            if diff < self.convergence['density'] and it > 5:
                d_enl_max = max(d_enl.values())
                if d_enl_max < self.convergence['energies']:
                    break

            if np.mod(it, 10) == 0:
                line = 'iter %3i, dn=%.1e>%.1e, max %i sp-iter' % \
                       (it, diff, self.convergence['density'], itmax)
                print(line, file=self.txt, flush=True)

            if it == self.maxiter - 1:
                if self.timing:
                    self.timer.summary()
                err = 'Density not converged in %i iterations' % (it + 1)
                raise RuntimeError(err)

        self.calculate_energies(enl, dens, echo='valence')
        print('converged in %i iterations' % it, file=self.txt)
        nel = self.get_number_of_electrons()
        line = '%9.4f electrons, should be %9.4f' % \
               (self.grid.integrate(dens, use_dV=True), nel)
        print(line, file=self.txt, flush=True)

        self.timer.stop('outer_scf')
        return dens, veff, enl, unlg, Rnlg

    def inner_scf(self, iteration, veff, enl, d_enl, dveff=None, itmax=100,
                  solve='all'):
        """ Solve the eigenstates for given effective potential.

        u''(r) - 2*(v_eff(r)+l*(l+1)/(2r**2)-e)*u(r)=0
        ( u''(r) + c0(r)*u(r) = 0 )

        r=r0*exp(x) --> (to get equally spaced integration mesh)

        u''(x) - u'(x) + c0(x(r))*u(r) = 0

        Parameters:

        iteration: iteration number in the SCF cycle
        itmax: maximum number of optimization steps per eigenstate
        solve: which eigenstates to solve: solve='all' -> all states;
               solve = [nl1, nl2, ...] -> only the given subset
        """
        self.timer.start('inner_scf')
        if self.scalarrel and dveff is None:
            spl = CubicSplineFunction(self.rgrid, veff)
            dveff = spl(self.rgrid, der=1)
        elif not self.scalarrel:
            dveff = np.array([])

        rgrid = self.rgrid
        xgrid = self.xgrid
        dx = xgrid[1] - xgrid[0]
        N = self.N
        unlg, Rnlg = {}, {}

        for n, l, nl in self.list_states():
            if solve != 'all' and nl not in solve:
                continue

            nodes_nl = n - l - 1

            if iteration == 0:
                eps = -1.0 * self.Z ** 2 / n ** 2
            else:
                eps = enl[nl]

            if iteration <= 3:
                delta = 0.5 * self.Z ** 2 / n ** 2  # previous!!!!!!!!!!
            else:
                delta = d_enl[nl]

            direction = 'none'
            epsmax = veff[-1] - l * (l + 1) / (2 * self.rgrid[-1] ** 2)
            it = 0
            u = np.zeros(N)
            hist = []

            while True:
                eps0 = eps
                self.timer.start('coeff')
                if _hotcent is not None:
                    c0, c1, c2 = _hotcent.construct_coefficients(l, eps, veff,
                                                             dveff, self.rgrid)
                else:
                    c0, c1, c2 = self.construct_coefficients(l, eps, veff,
                                                             dveff=dveff)
                assert c0[-2] < 0 and c0[-1] < 0
                self.timer.stop('coeff')

                # boundary conditions for integration from analytic behaviour
                # (unscaled)
                # u(r)~r**(l+1)   r->0
                # u(r)~exp( -sqrt(c0(r)) ) (set u[-1]=1
                # and use expansion to avoid overflows)
                u[0:2] = rgrid[0:2] ** (l + 1)
                self.timer.start('shoot')
                if _hotcent is not None:
                    u, nodes, A, ctp = _hotcent.shoot(u, dx, c2, c1, c0, N)
                else:
                    u, nodes, A, ctp = shoot(u, dx, c2, c1, c0, N)
                self.timer.stop('shoot')

                self.timer.start('norm')
                norm = self.grid.integrate(u ** 2)
                u /= np.sqrt(norm)
                self.timer.stop('norm')

                if nodes > nodes_nl:
                    # decrease energy
                    if direction == 'up':
                        delta /= 2
                    eps -= delta
                    direction = 'down'
                elif nodes < nodes_nl:
                    # increase energy
                    if direction == 'down':
                        delta /= 2
                    eps += delta
                    direction = 'up'
                elif nodes == nodes_nl:
                    shift = -0.5 * A / (rgrid[ctp] * norm)
                    if abs(shift) < 1e-8:  # convergence
                        break
                    if shift > 0:
                        direction = 'up'
                    elif shift < 0:
                        direction = 'down'
                    eps += shift

                if eps > epsmax:
                    eps = 0.5 * (epsmax + eps0)
                hist.append(eps)

                it += 1
                if it > 100:
                    print('Epsilon history for %s' % nl, file=self.txt)
                    for h in hist:
                        print(h)
                    print('nl=%s, eps=%f' % (nl,eps), file=self.txt)
                    print('max epsilon', epsmax, file=self.txt)
                    err = 'Eigensolver out of iterations. Atom not stable?'
                    raise RuntimeError(err)

            itmax = max(it, itmax)
            unlg[nl] = u
            Rnlg[nl] = unlg[nl] / self.rgrid
            d_enl[nl] = abs(eps - enl[nl])
            enl[nl] = eps

            if self.verbose:
                line = '-- state %s, %i eigensolver iterations' % (nl, it)
                line += ', e=%9.5f, de=%9.5f' % (enl[nl], d_enl[nl])
                print(line, file=self.txt)

            assert nodes == nodes_nl
            assert u[1] > 0.0

        self.timer.stop('inner_scf')
        return itmax, enl, d_enl, unlg, Rnlg

    def construct_coefficients(self, l, eps, veff, dveff=None):
        """ Construct the coefficients for Numerov's method; see shoot.py """
        c = 137.036
        ll = l * (l + 1)
        c2 = np.ones(self.N)
        if not self.scalarrel:
            c0 = -ll - 2 * self.rgrid ** 2 * (veff - eps)
            c1 = -1. * np.ones(self.N)
        else:
            assert dveff is not None
            # from Paolo Giannozzi: Notes on pseudopotential generation
            ScR_mass = 1 - 0.5 * (veff - eps) / c ** 2
            c0 = -ll - 2 * ScR_mass * self.rgrid ** 2 * (veff - eps)
            c0 -= dveff * self.rgrid / (2 * ScR_mass * c ** 2)
            c1 = dveff * self.rgrid / (2 * ScR_mass * c ** 2) - 1
        return c0, c1, c2


class RadialGrid:
    def __init__(self,grid):
        """
        mode
        ----

        rmin                                                        rmax
        r[0]     r[1]      r[2]            ...                     r[N-1] grid
        I----'----I----'----I----'----I----'----I----'----I----'----I
           r0[0]     r0[1]     r0[2]       ...              r0[N-2]     r0grid
           dV[0]     dV[1]     dV[2]       ...              dV[N-2]         dV

           dV[i] is volume element of shell between r[i] and r[i+1]
        """

        rmin, rmax = grid[0], grid[-1]
        N = len(grid)
        self.N = N
        self.grid = grid
        self.dr = self.grid[1:N] - self.grid[0:N - 1]
        self.r0 = self.grid[0:N - 1] + self.dr / 2
        # first dV is sphere (treat separately), others are shells
        self.dV = 4 * np.pi * self.r0 ** 2 * self.dr
        self.dV *= (4 * np.pi * rmax ** 3 / 3) / sum(self.dV)

    def get_grid(self):
        """ Return the whole radial grid. """
        return self.grid

    def get_N(self):
        """ Return the number of grid points. """
        return self.N

    def get_drgrid(self):
        """ Return the grid spacings (array of length N-1). """
        return self.dr

    def get_r0grid(self):
        """ Return the mid-points between grid spacings
        (array of length N-1).
        """
        return self.r0

    def get_dvolumes(self):
        """ Return dV(r) = 4 * pi * r ** 2 * dr. """
        return self.dV

    def plot(self, screen=True):
        rgrid = self.get_grid()
        plt.scatter(list(range(len(rgrid))), rgrid)
        if screen: 
            plt.show()

    def integrate(self, f, use_dV=False):
        """ Integrate function f (given with N grid points).
        int_rmin^rmax f*dr (use_dv=False) or int_rmin^rmax*f dV (use_dV=True)
        """
        if use_dV:
            return ((f[0:self.N - 1] + f[1:self.N]) * self.dV).sum() * 0.5
        else:
            return ((f[0:self.N - 1] + f[1:self.N]) * self.dr).sum() * 0.5

    def gradient(self, f):
        return np.gradient(f, self.grid)

    def divergence(self, f):
        return (1. / self.grid ** 2) * self.gradient(self.grid ** 2 * f)
