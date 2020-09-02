""" Tools for optimizing confinement parameters. """
import copy
import numpy as np
from scipy.optimize import minimize
from ase.io import read
from ase.dft.band_structure import BandStructure
from ase.calculators.dftb import Dftb
from hotcent.atomic_base import AtomicBase
from hotcent.slako import SlaterKosterTable
try:
    import matplotlib
    matplotlib.use('agg')
except ImportError:
    print('Warning: could not import matplotlib')


class ConfinementOptimizer:
    def __init__(self, *atoms, verbose=True, **sk_kwargs):
        """ Class which allows to optimize density and wave function
        confinement parameters to minimize a given residual function.

        arguments: one (or several) AtomicBase instance(s).
                   Note: the eigenvalues, hubbardvalues, occupations,
                   and spe kwargs passed to the SlaterKosterTable.write()
                   function (to be included in the SKF files) need to
                   be included in the atom's info (dict) attribute, i.e.
                       atom.info = {'eigenvalues': {'2s': ...}, ...}
        verbose:   whether to print output or not
        sk_kwargs: additional keyword arguments to be passed to the
                   SlaterKosterTable.run() method. If different arguments
                   need to be passed to different element pairs, these
                   can be specified using dictionaries. To e.g. set a
                   different number of grid points N for the C-C pair:
                   N={'default': 600, 'C-C': 450}
        """
        for atom in atoms:
            assert isinstance(atom, AtomicBase)
            assert 'info' in atom.__dir__() and len(atom.info) > 0, \
                   'Provide eigenvalues/Hubbard-values/... in atom.info dict'
        self.atoms = atoms
        self.verbose = verbose

        # Check integrity of sk_kwargs:
        for key, val in sk_kwargs.items():
            msg = '%s or %s pair not specified in %s dict'
            if isinstance(val, dict) and 'default' not in val:
                for i in range(self.atoms):
                    s1 = self.atoms[i].symbol
                    for j in range(i + 1):
                        s2 = self.atoms[j].symbol
                        prefix1, prefix2 = s1 + '-' + s2, s2 + '-' + s1
                        assert prefix1 in val or prefix2 in val, \
                               msg % (pair1, pair2, key)
        self.sk_kwargs = sk_kwargs

    def run(self, func, initial_guess={}, args=(), **opt_kwargs):
        """ Perform the actual optimization starting from an initial guess

        Returns a dict just like initial_guess but with adjusted parameters

        func: function which, when simply called without arguments,
              returns a residual (making use of SKF files generated
              for set of confinement parameters)

        initial_guess: a dict with the starting confinement potentials, e.g.
            {'O_2s': PowerConfinement(r0=3., s=2, adjustable=['r0', 's']),
             'O_2p': PowerConfinement(r0=3., s=2, adjustable=['r0', 's'])}
            for adjusting both the 'r0' and 's' parameters of a powerlaw
            wavefunction confinements for the oxygen 2s and 2p states.
            In case the same confinement should be applied to both, set e.g.
            {'O_2s,O_2p': PowerConfinement(r0=3., s=2, adjustable=['r0'])}

        args: arguments that need to be passed to func 

        opt_kwargs: extra options for scipy.optimize.minimize, see:
            https://docs.scipy.org/doc/scipy/reference/generated/
                    scipy.optimize.minimize.html
            E.g.: method='COBYLA', tol=1e-2, options={'maxiter': 1000}
        """
        assert len(initial_guess) > 0, 'Provide initial guesses for confinement'

        self.opt_param = []  # [('O_2s', 'r0'), ...]-like list
        opt_val = []   # list of parameter values
        self.vconf = initial_guess.copy()

        for key in sorted(self.vconf):
            for param in sorted(self.vconf[key].adjustable):
                value = self.vconf[key].__getattribute__(param)
                self.opt_param.append((key, param))
                opt_val.append(value)

        if self.verbose:
            print('*******************************************')
            print('Starting confinement parameter optimization')
            print('Adjustable parameters and initial values:')
            for (key, param), val in zip(self.opt_param, opt_val):
                print('  ' + key + '.' + param + ' : %.6f' % val)
            print('*******************************************', flush=True)

        def residual(opt_val, *args):
            if self.verbose:
                print('\nOPT_VAL:', '   '.join(map(str, opt_val)), flush=True)

            try:
                self.generate_skf(opt_val)
            except Exception as err:
                print('ERR:', err, flush=True)
                return 1e23

            residual = func(*args)
            if self.verbose:
                print('RESIDUAL:', residual, flush=True)
            return residual

        result = minimize(residual, opt_val, args=args, **opt_kwargs)

        if self.verbose:
            print('\n*******************************************')
            print('Optimization finished after %d iterations' % result.nfev)
            print('Adjustable parameters and final values:')
            for (key, param), val in zip(self.opt_param, result.x):
                print('  ' + key + '.' + param + ' : %.6f' % val)
            print('*******************************************\n', flush=True)

        # Make sure the SKF files and confinement potentials are right:
        self.generate_skf(result.x)
        return self.vconf

    def generate_skf(self, opt_val):
        """ Generates all *-*.skf files for a given set of
        confinement parameters

        opt_val: the (internally defined) list of parameter values used
                 in setting up the confinement potentials (and hence
                 changing the Slater-Koster integrals).
        """
        # Update the confinement parameters:
        for i, (key, param) in enumerate(self.opt_param):
            self.vconf[key].__setattr__(param, opt_val[i])

        if self.verbose:
            print('VCONF:')
            for key in sorted(self.vconf):
                print('  ' + key + ' : ' + str(self.vconf[key]))

        # Run the atomic DFT calculations
        for atom in self.atoms:
            conf, wf_conf = None, {}
            for key, vconf in self.vconf.items():
                if '%s_n' % atom.symbol in key:
                    conf = vconf
                else:
                    for nl in atom.valence:
                        if '%s_%s' % (atom.symbol, nl) in key:
                            wf_conf[nl] = vconf
            if conf is not None:
                atom.set_confinement(conf)
            if wf_conf:
                atom.set_wf_confinement(wf_conf)
            atom.run()

        # Perform the Slater-Koster integrations
        for i, atom1 in enumerate(self.atoms):
            s1 = atom1.symbol
            for j in range(i + 1):
                atom2 = self.atoms[j]
                s2 = atom2.symbol

                prefix1 = s1 + '-' + s2
                prefix2 = s2 + '-' + s1
                sk_kwargs = {}
                for key, val in self.sk_kwargs.items():
                    if isinstance(val, dict):
                        if prefix1 in val:
                            sk_kwargs[key] = val[prefix1]
                        elif prefix2 in val:
                            sk_kwargs[key] = val[prefix2]
                        else:
                            sk_kwargs[key] = val['default']
                    else:
                        sk_kwargs[key] = val

                sk = SlaterKosterTable(atom1, atom2, timing=False, txt=None)
                sk.run(**sk_kwargs)

                filename = '%s-%s.skf' % (s1, s2)
                if s1 == s2:
                    sk.write(filename=filename, pair=(s1, s2), **atom1.info)
                else:
                    sk.write(filename=filename, pair=(s1, s2))
                    filename = '%s-%s.skf' % (s2, s1)
                    sk.write(filename=filename, pair=(s2, s1))
        return


class DftbPlusBandStructure:
    def __init__(self, **dftbplus_kwargs):
        """ Class for calculating band structures with DFTB+
        and computing residuals with respect to reference
        band structures.

        dftbplus_kwargs: extra keyword arguments needed for initializing
                         the ASE DFTB+ calculator (other than k-points
                         and Atoms objects).
        """
        self.dftbplus_kwargs = dftbplus_kwargs
        self.bs_ref = []

    def add_reference_bandstructure(self, bs_ref, atoms=None, kpts_scf=None,
                                reference_level=None, nsemicore=0, weight=1.,
                                distribution={'type': 'Boltzmann', 'kBT': 1e3}):
        """ Adds a bandstructure to be used as a reference.

        bs_ref: an ase.dft.band_structure.BandStructure instance
                containing a path through the Brillouin zone and
                the corresponding eigenvalues calculated with the
                reference method of choice (typically DFT).

                Note: the energy value used as a reference (i.e. the Fermi
                energy or the valence band maximum) should be consistent
                with the reference using in the DFTB calculations.

        atoms: the corresponding Atoms object.
        kpts_scf: the k-points to be used in the DFTB self-consistency
                  cycle before performing a band structure calculation.
        reference_level: the type of energy level to which the eigenvalues
                   are referenced (as stored in the bs_ref.reference
                   attribute). Needs to be one of the following:
                   'vbm': the valence band maximum, generally the best
                          choice, in particular when a band gap is present
                   'fermi': the Fermi level
        nsemicore: if the reference calculation included e.g. semicore
                   states that are not included in the DFTB model, the
                   number of semicore bands need to be specified here.
        weight: the weight relative to other reference band structures
                in the calculation of the total residual.
        distribution: how to weigh the deviations of the DFTB eigenvalues
                with respect to the reference. Currently only a Boltzmann
                distribution is implemented where the contribution to
                the residual of one point in the band path is weighted by:
                    w[k, n] = exp( abs(eigenvalue[k, n] - e_ref) / kBT )
                The default thermal energy kBT amounts to 1000 eV, meaning
                all bands are weighted equally. Reducing kBT increases the
                relative weights attached to the bands closer to the
                reference energy.
        """
        assert isinstance(bs_ref, BandStructure)
        if isinstance(atoms, str):
            bs_ref.atoms = read(atoms)
        else:
            bs_ref.atoms = atoms
        bs_ref.kpts_scf = kpts_scf
        assert reference_level is not None
        assert reference_level.lower() in ['vbm', 'fermi']
        bs_ref.reference_level = reference_level.lower()
        bs_ref.nsemicore = nsemicore
        bs_ref.weight = weight
        assert distribution['type'].lower() == 'boltzmann'
        bs_ref.distribution = distribution
        self.bs_ref.append(bs_ref)

    def get_residual(self):
        """ Returns the total residual for the DFTB band structures
        compared to the reference band structures.
        """
        assert len(self.bs_ref) > 0, 'Add reference band structure(s)!'
        residual = 0.

        for bs_ref in self.bs_ref:
            bs_dftb = self.calculate_bandstructure(bs_ref)
            shape_ref = np.shape(bs_ref.energies)
            shape_dftb = np.shape(bs_dftb.energies)

            assert shape_ref[0] == shape_dftb[0], [shape_ref, shape_dftb]
            assert shape_ref[1] == shape_dftb[1], [shape_ref, shape_dftb]

            nskip = bs_ref.nsemicore // 2
            imin = min(shape_ref[2] - nskip, shape_dftb[2])

            diffs = bs_ref.energies[:, :, nskip:nskip + imin] - \
                    bs_dftb.energies[:, :, :imin]
            diffs -= bs_ref.reference - bs_dftb.reference

            if bs_ref.distribution['type'].lower() == 'boltzmann':
                logw = -1. * np.abs(bs_ref.energies[:, :, nskip:nskip + imin])
                logw /= bs_ref.distribution['kBT']
                w = np.exp(logw)

            residual += ((bs_ref.weight * w * diffs) ** 2).sum()
        return residual

    def calculate_bandstructure(self, bs):
        """ Returns a BandStructure object with a DFTB band structure.
        Assumes that the necessary *-*.skf files have been generated.

        bs: ase.dft.band_structure.BandStructure instance with an
            additional 'atoms' and 'kpts_scf' attributes. The latter
            represents the k-points to be used in the self-consistent
            calculation prior to the bandstructure calculation in DFTB.
        """
        # Perform a regular DFTB calculation first
        # to converge the charges and get the
        # reference energy
        atoms = bs.atoms.copy()
        calc = Dftb(atoms=atoms, kpts=bs.kpts_scf, **self.dftbplus_kwargs)
        atoms.set_calculator(calc)
        etot = atoms.get_potential_energy()

        efermi = calc.get_fermi_level()
        if bs.reference_level == 'vbm':
            eig = np.array(calc.results['eigenvalues'])
            eref = np.max(eig[eig < efermi])
        elif bs.reference_level == 'fermi':
            eref = efermi

        # Now a band structure calculation
        kwargs = self.dftbplus_kwargs.copy()
        kwargs.update({'Hamiltonian_MaxSCCIterations': 1,
                       'Hamiltonian_ReadInitialCharges': 'Yes',
                       'Hamiltonian_SCCTolerance': 1e6})
        calc = Dftb(atoms=atoms, kpts=bs.path.kpts, **kwargs)
        atoms.set_calculator(calc)
        etot = atoms.get_potential_energy()

        # Update the reference level, k-points,
        # and eigenenergies
        bs_new = copy.deepcopy(bs)
        bs_new.reference = eref
        bs_new.kpts = calc.get_ibz_k_points()
        bs_new.energies = []

        for s in range(calc.get_number_of_spins()):
            bs_new.energies.append([calc.get_eigenvalues(kpt=k, spin=s)
                                    for k in range(len(bs_new.kpts))])
        bs_new.energies = np.array(bs_new.energies)

        return bs_new
