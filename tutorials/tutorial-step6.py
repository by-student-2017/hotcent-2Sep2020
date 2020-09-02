#!/usr/bin/python3

from ase.io.jsonio import read_json
from ase.units import Bohr
from ase.build import bulk
from ase.data import atomic_numbers, covalent_radii
from hotcent.atomic_dft import AtomicDFT
from hotcent.confinement import PowerConfinement
from hotcent.tools import ConfinementOptimizer, DftbPlusBandStructure


# Setting up the atomic DFT instance(s) and
# calculating the eigenvalues and Hubbard values
atom = AtomicDFT('Si',
                 configuration='[Ne] 3s2 3p2 3d0',
                 valence=['3s', '3p', '3d'],
                 xc='LDA',
                 scalarrel=False,
                 mix=0.05,
                 maxiter=500,
                 confinement=PowerConfinement(r0=40., s=4),
                 txt=None)
atom.run()
atom.info = {}
atom.info['eigenvalues'] = {nl: atom.get_eigenvalue(nl) for nl in atom.valence}

U_p = atom.get_hubbard_value('3p', scheme='central', maxstep=1.)
atom.info['hubbardvalues'] = {'s': U_p}
atom.info['occupations'] = {'3s': 2, '3p': 2, '3d': 0}

# Creating a DFTB+ band structure evaluator and
# supplying it with a reference (DFT) band structure
dpbs = DftbPlusBandStructure(Hamiltonian_SCC='Yes',
                             Hamiltonian_OrbitalResolvedSCC='No',
                             Hamiltonian_MaxAngularMomentum_='',
                             Hamiltonian_MaxAngularMomentum_Si='d',
                             Hamiltonian_PolynomialRepulsive='SetForAll {Yes}')

bs_gpaw = read_json('bs_gpaw.json')  # the reference band structure (DFT)
atoms = bulk('Si', 'diamond')
# see hotcent.tools.DftbPlusBandStructure for more information
# on the various keyword arguments used below
dpbs.add_reference_bandstructure(bs_gpaw, atoms=atoms, kpts_scf=(5, 5, 5),
                                 reference_level='vbm', nsemicore=0, weight=1.,
                                 distribution={'type': 'Boltzmann', 'kBT': 1.5})

# Setting up and running the actual optimizer
# (the keyword arguments are known from hotcent.slako.SlaterKosterTable)
confopt = ConfinementOptimizer(atom, N=500, rmin=0.4, dr=0.02, stride=4,
                               superposition='density', xc='LDA')

# The initial confinement parameters are the same as in Tutorial #1.
# The additional 'adjustable' keyword argument serves to indicate
# which parameters are allowed to vary. Here we keep the quadratic
# form (s=2) and treat the confinement radii r0 as variable.
rcov = covalent_radii[atomic_numbers['Si']] / Bohr
initial_guess = {'Si_3s,Si_3p': PowerConfinement(r0=2 * rcov, s=2, adjustable=['r0']),
                 'Si_3d': PowerConfinement(r0=2 * rcov, s=2, adjustable=['r0']),
                 'Si_n': PowerConfinement(r0=3 * rcov, s=2, adjustable=['r0'])}
# Note that the 'Si_3s,Si_3p' combination indicates that the same
# confinement potential is used for both states. The Si_3d confinement
# (as well as the density confinement) is defined separately and
# their r0 parameters are hence allowed to vary independently.

# Only two iterations are performed here to limit the computational time.
# The confinement therefore hardly changes from the initial guess.
# In real scenarios, typically on the order of 100 iterations are required.
vconf = confopt.run(dpbs.get_residual, initial_guess=initial_guess, tol=1e-2,
                    method='COBYLA', options={'maxiter': 100, 'rhobeg': 0.2})
# See scipy.optimize.minimize for more information on available optimization methods:
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html

# Make the DFTB band structure plot on the basis of the latest Si-Si.skf file
# (generated with the confinement parameters yielding the lowest residual)
bs_dftb = dpbs.calculate_bandstructure(bs_gpaw)
bs_dftb.plot(filename='bs_dftb.png', emax=bs_dftb.reference + 5, emin=bs_dftb.reference - 15)