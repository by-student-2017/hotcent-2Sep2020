""" Definition of the base class for atomic DFT
calculations.

The code below draws heavily from the Hotbit code 
written by Pekka Koskinen (https://github.com/pekkosk/
hotbit/blob/master/hotbit/parametrization/atom.py).
"""
import os
import sys
import collections
import numpy as np
from scipy.optimize import minimize
from ase.data import atomic_numbers, covalent_radii
from ase.units import Bohr, Ha
from hotcent.interpolation import CubicSplineFunction
from hotcent.timing import Timer
from hotcent.confinement import Confinement, ZeroConfinement
try:
    import matplotlib.pyplot as plt
except:
    plt = None


not_solved_message = 'A required attribute is missing. ' \
                     'Please call the run() method first.'

class AtomicBase:
    def __init__(self,
                 symbol,
                 configuration='',
                 valence=[],
                 confinement=None,
                 wf_confinement=None,
                 scalarrel=False,
                 mix=0.2,
                 maxiter=200,
                 rmax=100.0,
                 nodegpts=500,
                 timing=False,
                 verbose=False,
                 txt='-'):
        """ Base class for atomic DFT calculators

        symbol:         chemical symbol
        configuration:  e.g. '[He] 2s2 2p2'
        valence:        valence orbitals, e.g. ['2s','2p'].
        confinement:    confinement potential for the electron density 
                        (see hotcent.confinement). The default None
                        means no confinement will be applied.
        wf_confinement: dictionary with confinement potentials for the
                        valence orbitals. If None, the same confinement will
                        be used as for the electron density. If a certain
                        hotcent.confinement.Confinement instance is provided,
                        this will be applied to all valence states. If a
                        dictionary is provided, it is supposed to look like
                        this: {nl: <a certain Confinement instance, or None>
                         for each nl in your set of valence states}.
                        For missing entries, no confinement will be applied.
        scalarrel:      Use scalar relativistic corrections
        mix:            effective potential mixing constant
        maxiter:          maximum number of iterations for self-consistency.
        rmax:           radial cutoff in Bohr
        nodegpts:       total number of grid points is nodegpts times the max number
                        of antinodes for all orbitals
        timing:         output of timing summary
        verbose:        increase verbosity during iterations
        txt:            where output should be printed
                        use '-' for stdout (default), None for /dev/null,
                        any other string for a text file, or a file handle
        """
        self.symbol = symbol
        self.valence = valence
        self.scalarrel = scalarrel
        self.mix = mix
        self.maxiter = maxiter
        self.rmax = rmax
        self.nodegpts = nodegpts
        self.timing = timing
        self.verbose = verbose

        if txt is None:
            self.txt = open(os.devnull, 'w')
        elif isinstance(txt, str):
            if txt == '-':
                self.txt = sys.stdout
            else:
                self.txt = open(txt, 'a')
        else:
            self.txt = txt

        self.set_confinement(confinement)
        self.set_wf_confinement(wf_confinement)

        self.timer = Timer('Atomic', txt=self.txt, enabled=self.timing)
        self.timer.start('init')

        self.Z = atomic_numbers[self.symbol]
        assert len(self.valence) > 0

        assert len(configuration) > 0, "Specify the electronic configuration!"
        self.set_configuration(configuration)

        self.maxl = 9
        self.maxn = 9
        self.unlg = {}
        self.Rnlg = {}
        self.unl_fct = {nl: None for nl in self.configuration}
        self.Rnl_fct = {nl: None for nl in self.configuration}
        self.veff_fct = None
        self.dens_fct = None
        self.vhar_fct = None
        self.solved = False

    def set_configuration(self, configuration):
        """ Set the electron configuration

        configuration: e.g. '[He] 2s2 2p2'
        """
        self.configuration = {}
        noble_conf = {'He':{'1s':2}}
        noble_conf['Ne'] = dict({'2s':2, '2p':6}, **noble_conf['He'])
        noble_conf['Ar'] = dict({'3s':2, '3p':6}, **noble_conf['Ne'])
        noble_conf['Kr'] = dict({'3d':10, '4s':2, '4p':6}, **noble_conf['Ar'])
        noble_conf['Xe'] = dict({'4d':10, '5s':2, '5p':6}, **noble_conf['Kr'])
        noble_conf['Rn'] = dict({'4f':14, '5d':10, '6s':2, '6p':6},
                                **noble_conf['Xe'])
        noble_conf['Og'] = dict({'5f':14, '6d':10, '7s':2, '7p':6},
                                **noble_conf['Rn'])

        for term in configuration.split():
            if term[0] == '[' and term[-1] == ']':
                core = term[1:-1]
                assert core in noble_conf, "[Core] config is not a noble gas!"
                conf = noble_conf[core]
            else:
                conf = {term[:2]: float(term[2:])}
            self.configuration.update(conf)

    def set_confinement(self, confinement):
        if confinement is None:
            self.confinement = ZeroConfinement()
        else:
            self.confinement = confinement

    def set_wf_confinement(self, wf_confinement):
        if wf_confinement is None:
            self.wf_confinement = {}
        elif isinstance(wf_confinement, Confinement):
            self.wf_confinement = {nl: wf_confinement for nl in self.valence}
        elif isinstance(wf_confinement, dict):
            self.wf_confinement = {}
            for nl in self.valence:
                if nl not in wf_confinement or wf_confinement[nl] is None:
                    self.wf_confinement[nl] = ZeroConfinement()
                else:
                    self.wf_confinement[nl] = wf_confinement[nl]
        else:
            msg = "Don't know what to do with the provided wf_confinement:\n"
            msg += str(wf_confinement)
            raise ValueError(msg)

    def run(self, **kwargs):
        """ Child classes must implement a run() method which,
        in turn, is supposed to set the following attributes:

        self.solved: whether the calculations are considered to be done
        self.total_energy: the total energy
        self.rgrid: an array with the radial grid points g
        self.dens: an array with the electron density on the radial grid
        self.vhar: an array with the Hartree potential on the radial grid
        self.veff: an array with the effective potential on the radial grid
                   (note: veff = vnuc + vhar + vxc + vconf)
        self.enl: a {'nl': eigenvalue} dictionary
        self.Rnlg: a {'nl': R_nl(g) array} dictionary
        self.unlg: a {'nl': u_nl(g) array} dictionary (u_nl = R_nl / r)
        """
        raise NotImplementedError('Child class must implement run() method!')

    def __getstate__(self):
        """ Return dictionary of all pickable items. """
        d = self.__dict__.copy()
        for key in self.__dict__:
            if isinstance(d[key], collections.Callable):
                d.pop(key)
        d.pop('out')
        return d

    def get_symbol(self):
        """ Return atom's chemical symbol. """
        return self.symbol

    def get_number_of_electrons(self):
        return sum(self.configuration.values())

    def list_states(self):
        """ List all potential states {(n,l,'nl')}. """
        states = []
        for l in range(self.maxl + 1):
            for n in range(1, self.maxn + 1):
                nl = tuple2nl(n, l)
                if nl in self.configuration:
                    states.append((n, l, nl))
        return states

    def get_valence_orbitals(self):
        """ Get list of valence orbitals, e.g. ['2s','2p'] """
        return self.valence

    def get_energy(self):
        assert self.solved, not_solved_message
        return self.total_energy

    def get_epsilon(self, nl):
        """ E.g. get_eigenvalue('2p') """
        assert self.solved, not_solved_message
        return self.enl[nl]

    def get_valence_energies(self):
        """ Return list of valence eigenenergies. """
        assert self.solved, not_solved_message
        return [(nl, self.enl[nl]) for nl in self.valence]

    def get_eigenvalue(self, nl):
        return self.get_epsilon(nl)

    def get_wf_range(self, nl, fractional_limit=1e-7):
        """ Return the maximum r for which |R(r)| is
        less than fractional_limit * max(|R(r)|) """
        assert self.solved, not_solved_message
        wfmax = np.nanmax(np.abs(self.Rnlg[nl]))
        for r, wf in zip(self.rgrid[-1::-1], self.Rnlg[nl][-1::-1]):
            if abs(wf) > fractional_limit * wfmax:
                return r

    def Rnl(self, r, nl, der=0):
        """ Rnl(r, '2p') """
        assert self.solved, not_solved_message
        if self.Rnl_fct[nl] is None:
            self.Rnl_fct[nl] = CubicSplineFunction(self.rgrid, self.Rnlg[nl])
        return self.Rnl_fct[nl](r, der=der)

    def unl(self, r, nl, der=0):
        """ unl(r, '2p') = Rnl(r, '2p') / r """
        assert self.solved, not_solved_message
        if self.unl_fct[nl] is None:
            self.unl_fct[nl] = CubicSplineFunction(self.rgrid, self.unlg[nl])
        return self.unl_fct[nl](r, der=der)

    def electron_density(self, r, der=0):
        """ Return the all-electron density at r. """
        assert self.solved, not_solved_message
        if self.dens_fct is None:
            self.dens_fct = CubicSplineFunction(self.rgrid, self.dens)
        return self.dens_fct(r, der=der)

    def nuclear_potential(self,r):
        return -self.Z / r

    def effective_potential(self, r, der=0):
        """ Return effective potential at r or its derivatives. """
        assert self.solved, not_solved_message
        if self.veff_fct is None:
            self.veff_fct = CubicSplineFunction(self.rgrid, self.veff)
        return self.veff_fct(r, der=der)

    def hartree_potential(self, r):
        """ Return the Hartree potential at r. """
        assert self.solved, not_solved_message
        if self.vhar_fct is None:
            self.vhar_fct = CubicSplineFunction(self.rgrid, self.vhar)
        return self.vhar_fct(r)

    def plot_Rnl(self, filename=None, only_valence=True):
        """ Plot radial wave functions with matplotlib.
        
        filename:  output file name + extension (extension used in matplotlib)
                   default = <Element>_Rnl.pdf
        only_valence: whether to only plot the valence states or all of them
        """
        assert plt is not None, 'Matplotlib could not be imported!'
        assert self.solved, not_solved_message

        rmax = 3 * covalent_radii[self.Z] / Bohr
        ri = np.where(self.rgrid < rmax)[0][-1]

        if only_valence:
            states = self.valence
        else:
            states = [x[2] for x in self.list_states()]

        p = np.ceil(np.sqrt(len(states)))
        q = 2 * p - 1 if len(states) % 2 == 0 else 2 * p

        fig = plt.figure(dpi=400)
        i = 1
        # as a function of grid points
        for nl in states:
            ax = plt.subplot(q, p, i)
            plt.plot(self.Rnlg[nl])
            plt.xticks(size=5)
            plt.grid(ls='--')

            # annotate
            c = 'k'
            if nl in self.valence:
                c = 'r'
            plt.text(0.5, 0.4, r'$R_{%s}(i)$' % nl, transform=ax.transAxes,
                     size=15, color=c)
            if ax.is_first_col():
                plt.ylabel(r'$R_{nl}(i)$', size=8)
            i += 1

        # as a function of radius
        i = p ** 2 + 1
        for nl in states:
            ax = plt.subplot(q, p, i)
            plt.plot(self.rgrid[:ri], self.Rnlg[nl][:ri])
            plt.xticks(size=5)
            plt.grid(ls='--')

            # annotate
            c = 'k'
            if nl in self.valence: 
                c = 'r'
            plt.text(0.5, 0.4, r'$R_{%s}(r)$' % nl, transform=ax.transAxes,
                     size=15, color=c)
            if ax.is_first_col():
                plt.ylabel(r'$R_{nl}(r)$', size=8)
            if ax.is_last_row():
                plt.xlabel('r (Bohr)', size=8)
            i += 1

        fig.subplots_adjust(hspace=0.2, wspace=0.2)
        plt.figtext(0.4, 0.95, r'$R_{nl}(r)$ for %s' % self.symbol)

        if filename is None:
            filename = '%s_Rnl.pdf' % self.symbol
        plt.savefig(filename, bbox_inches='tight')
        plt.clf()

    def plot_density(self, filename=None):
        """ Plot the electron density and valence orbital densities.

        Note that the plotted electron density (rho_0) generally does
        not correspond to the sum of the valence orbital densities in
        the valence region. For this to be the case, the orbital densities
        would need to be multiplied by their occupation numbers, and
        the same confinement potential would need to be applied throughout.

        filename:  output file name + extension (extension used in matplotlib)
                   default = <Element>_rho.pdf
        """
        assert plt is not None, 'Matplotlib could not be imported!'
        assert self.solved, not_solved_message

        rmax = 3 * covalent_radii[self.Z] / Bohr
        ri = np.where(self.rgrid > rmax)[0][0]

        fig = plt.figure(figsize=(6.4, 4.8), dpi=400)

        core_dens = 0
        colors = ['red', 'green', 'blue']  # s, p , d
        for n, l, nl in self.list_states():
            if nl not in self.valence:
                continue

            dens = self.Rnlg[nl] ** 2 / (4 * np.pi)
            occupied = self.configuration[nl] > 0
            suffix = '' if occupied else '*'
            ls = '-' if occupied else '--'
            label = r'$|R_\mathrm{%s%s}(r) / \sqrt{4\pi}|^2$' % (nl, suffix)

            plt.semilogy(self.rgrid[:ri], dens[:ri], ls=ls, color=colors[l],
                        label=label)

        dens = self.dens[:ri]
        plt.semilogy(self.rgrid[:ri], dens, 'k-', label=r'$\rho_0(r)$')

        ymax = np.exp(np.ceil(np.log(np.max(dens))))
        plt.ylim([1e-7, ymax])
        plt.xlim([-0.05 * rmax, rmax])
        plt.xlabel('r (Bohr)')
        plt.grid(ls='--')
        plt.legend(loc='upper right', ncol=2)
        plt.title('Electron and orbital densities for %s' % self.symbol)

        if filename is None:
            filename = '%s_rho.pdf' % self.symbol
        plt.savefig(filename, bbox_inches='tight')
        plt.clf()

    def plot_rho(self, *args, **kwargs):
        self.plot_density(*args, **kwargs)

    def write_unl(self, filename, only_valence=True, step=20):
        """ Append functions unl=Rnl*r into file.
            Only valence functions by default.

        Parameters:
        -----------
        filename:         output file name (e.g. XX.elm)
        only_valence:     output of only valence orbitals
        step:             step size for output grid
        """
        assert self.solved, not_solved_message
        if only_valence:
            orbitals = self.valence
        else:
            orbitals = [nl for n,l,nl in self.list_states()]

        with open(filename, 'a') as f:
            for nl in orbitals:
                f.write('\n\nu_%s=' % nl)
                for r, u in zip(self.rgrid[::step], self.unlg[nl][::step]):
                    f.write(r, u)
            f.write('\n\n')

    def fit_sto(self, nl, num_exp, num_pow, regularization=1e-6,
                filename=None):
        """ Fit Slater-type orbitals to the one on the grid.
            See self.write_hsd() for more information.

        Parameters:
        -----------
        nl:              the (valence) orbital of interest (e.g. '2p')
        num_exp:         number of exponents to use
        num_pow:         number of r-powers for each exponents
        regularization:  penalty to be used in the L2-regularization
        filename:        filename for a figure with the grid-based
                         and STO-fitted orbitals (to verify that the
                         fit is decent)
        """
        assert self.solved, not_solved_message
        print('Fitting Slater-type orbitals to eigenstate %s' % nl,
              file=self.txt)
        r = self.rgrid
        y = self.Rnlg[nl]
        n, l = nl2tuple(nl)
        num_coeff = num_exp * num_pow
        num_r = len(r)

        def regression(param):
            A = np.zeros((num_r, num_coeff))

            for i in range(num_exp):
                rexp = np.exp(-param[i] * r)[:, None]
                A[:, i * num_pow : (i + 1) * num_pow] = rexp

            for j in range(num_pow):
                rpow = np.power(r, l + j)
                for i in range(num_exp):
                    A[:, num_pow * i + j] *= rpow

            A_loss = np.identity(num_coeff) * regularization
            AA = np.vstack((A, A_loss))

            y_loss = np.zeros(num_coeff)
            yy = np.hstack((y, y_loss))

            coeff, residual, rank, s = np.linalg.lstsq(AA, yy, rcond=None)
            values = np.dot(A, coeff)
            return coeff, values, residual[0]

        def residual(param):
            coeff, values, residual = regression(param)
            return residual

        if num_exp > 1:
            x0, x1 = self.Z, 0.5
            ratio = (x0 / x1) ** (1. / (num_exp - 1.))
            guess = [x0 / (ratio ** i) for i in range(num_exp)]
        else:
            guess = [1.]

        result = minimize(residual, guess, method='COBYLA',
                          options={'rhobeg': 0.1, 'tol': 1e-8})
        exponents = result.x
        coeff, values, residual = regression(exponents)

        integral = np.trapz((r * y) ** 2, x=r)
        if abs(integral - 1) > 1e-1:
            print('Warning -- significant deviation from unity for integral'
                  ' of grid-based %s orbital: %.5f' % (nl, integral),
                  file=self.txt)

        integral = np.trapz((r * values) ** 2, x=r)
        if abs(integral - 1) > 1e-1:
            print('Warning -- significant deviation from unity for integral'
                  ' of STO-based %s orbital: %.5f' % (nl, integral),
                  file=self.txt)

        if filename is not None:
            rmax = 3 * covalent_radii[self.Z] / Bohr
            imax = np.where(r < rmax)[0][-1]
            rmin = 1e-3 * self.Z
            imin = np.where(r < rmin)[0][-1]
            plt.plot(r[imin:imax], y[imin:imax], '-', label='On the grid')
            plt.plot(r[imin:imax], values[imin:imax], '--', label='With STOs')
            plt.xlim([0., rmax])
            plt.grid(ls='--')
            plt.legend(loc='upper right')
            plt.xlabel('r (Bohr radii)')
            plt.ylabel('Psi_%s (a.u.)' % nl)
            plt.savefig(filename)
            plt.clf()

        return exponents, coeff, values, residual

    def write_hsd(self, filename=None, num_exp=None, num_pow=4, wfthr=1e-2):
        """ Writes a HSD-format file with information on the valence
        orbitals. This includes a projection of these orbitals
        on a set of Slater-type orbitals, for post-processing
        purposes (e.g. using the Waveplot tool part of DFTB+).

        The expansion is the same as in DFTB+.
        For an atomic orbital with angular momentum l:

          R_l(r) = \sum_{i=0}^{num_exp-1} \sum_{j=0}^{num_pow-1}
                     coeff_{i,j} * r ^ (l + j) * \exp(-exponent_i * r)

        Note that also the same normalization is used as in DFTB+.
        This means that \int_{r=0}^{\infty} r^2 * |R_l(r)|^2 dr = 1.

        Starting from a reasonable initial guess, the exponents
        are optimized to reproduce the grid-based orbitals,
        with the coefficient matrix being determined by
        (L2-regularized) linear regression at each iteration.

        Parameters:
        -----------
        filename:   output file name. If None, the name
                    defaults to wcf.<Element>.hsd
        num_exp:    number of exponents to use
                    default = highest principal quantum number
        num_pow:    number of powers for each exponent
                    default = 4
        wfthr:      parameter determining the 'Cutoff' radius,
                    which will be where the orbital tail goes
                    below wfthr in absolute value
        """
        assert self.solved, not_solved_message
        if filename is None:
            filename = 'wfc.%s.hsd' % self.symbol

        if num_exp is None:
            num_exp = max([nl2tuple(nl)[0] for nl in self.valence])

        with open(filename, 'a') as f:
            f.write('%s = {\n' % self.symbol)
            f.write('  AtomicNumber = %d\n' % self.Z)

            for nl in self.valence:
                n, l = nl2tuple(nl)
                exp, coeff, values, resid = self.fit_sto(nl, num_exp, num_pow)
                icut = len(values) - 1
                while abs(values[icut]) < wfthr:
                    icut -= 1
                rcut = np.round(self.rgrid[icut + 1], 1)

                f.write('  Orbital = {\n')
                f.write('    AngularMomentum = %d\n' % l)
                f.write('    Occupation = %.6f\n' % self.configuration[nl])
                f.write('    Cutoff = %.1f\n' % rcut)
                f.write('    Exponents = {\n    ')
                for e in exp:
                    f.write('  %.8f' % e)
                f.write('\n    }\n    Coefficients = {\n')
                for c in coeff:
                    f.write('      {: .8E}\n'.format(c))
                f.write('    }\n  }\n')
            f.write('}\n')
        return

    def get_hubbard_value(self, nl, maxstep=1., scheme=None):
        """ Calculates the Hubbard value of an orbital using
        (second order) finite differences.

        nl:      e.g. '2p'
        maxstep: the maximal step size in the orbital occupancy;
                 the default value of 1 means not going further
                 than the monovalent ions
        scheme:  the finite difference scheme, either 'central',
                 'forward' or 'backward' or None. In the last
                 case the appropriate scheme will be chosen
                 based on the orbital in question being empty
                 of partly or entirely filled
        """
        assert scheme in [None, 'central', 'forward', 'backward']

        if scheme is None:
            n, l = nl2tuple(nl)
            max_occup = 2 * (2 * l + 1)
            occup = self.configuration[nl]
            if occup == 0:
                scheme = 'forward'
            elif occup == max_occup:
                scheme = 'backward'
            else:
                scheme = 'central'

        directions = {'forward': [0, 1, 2],
                      'central': [-1, 0, 1],
                      'backward': [-2, -1, 0]}
        delta = maxstep if scheme == 'central' else 0.5 * maxstep

        configuration = self.configuration.copy()

        energies = {}
        bar = '+' * 12
        for direction in directions[scheme]:
            self.configuration = configuration.copy()
            self.configuration[nl] += direction * delta
            s = ' '.join([nl + '%.1f' % self.configuration[nl]
                          for nl in self.valence])
            print('\n%s Configuration %s %s' % (bar, s, bar), file=self.txt)
            self.run()
            energies[direction] = self.total_energy

        if scheme in ['forward', 'central']:
            EA = (energies[0] - energies[1]) / delta
            print('\nElectron affinity = %.5f Ha (%.5f eV)' % (EA, EA * Ha),
                  file=self.txt)

        if scheme in ['backward', 'central']:
            IE = (energies[-1] - energies[0]) / delta
            print('\nIonization energy = %.5f Ha (%.5f eV)' % (IE, IE * Ha),
                  file=self.txt)
        U = 0.
        for i, d in enumerate(directions[scheme]):
            factor = 1 if i % 2 == 0 else -2
            U += energies[d] * factor / (delta ** 2)

        self.configuration = configuration.copy()
        self.solved = False

        return U


angular_momenta = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l']

def nl2tuple(nl):
    """ Transforms e.g. '2p' into (2, 1) """
    return (int(nl[0]), angular_momenta.index(nl[1]))

def tuple2nl(n, l):
    """ Transforms e.g. (2, 1) into '2p' """
    return '%i%s' % (n, angular_momenta[l])
