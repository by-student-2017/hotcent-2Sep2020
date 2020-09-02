""" Note: currently only works with DFTB+ 17.1
because newer versions do not include eigenvalues,
Fermi levels, ... in the 'results.tag' file,
which the ASE interface relies upon.
"""
import os
import numpy as np
from ase.build import bulk
from ase.io import write
from ase.io.jsonio import write_json, read_json
from ase.units import Bohr
from ase.data import atomic_numbers, covalent_radii
from hotcent.atomic_dft import AtomicDFT
from hotcent.confinement import PowerConfinement
from hotcent.tools import ConfinementOptimizer, DftbPlusBandStructure


atoms = bulk('NO', 'zincblende', a=3.6)
atoms.set_initial_magnetic_moments([1., 0.])
write('NO.traj', atoms)

if not os.path.exists('bs_dft.json'):
    from ase.dft.kpoints import bandpath
    from ase.dft.band_structure import get_band_structure
    from gpaw import GPAW, PW, MixerSum, FermiDirac
    from gpaw.eigensolvers import CG

    calc = GPAW(mode=PW(400),
                maxiter=250,
                spinpol=True,
                kpts=(3, 3, 3),
                xc='LDA',
                txt='-',
                occupations=FermiDirac(0.02, fixmagmom=True),
                mixer=MixerSum(0.02, 5, 100),
                )
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    efermi = calc.get_fermi_level().max()
    Nk = len(calc.get_ibz_k_points())
    Ns = calc.get_number_of_spins()
    eigval = np.array([[calc.get_eigenvalues(kpt=k, spin=s)
                        for k in range(Nk)] for s in range(Ns)])
    evbm = np.max(eigval[eigval < efermi])

    calc.set(nbands=8,  # 4 occupied and 4 unoccupied bands
             fixdensity=True,
             eigensolver=CG(niter=5),
             symmetry='off',
             kpts={'path': 'WGKL', 'npoints': 20},
             convergence={'bands': 'all'},
             )
    calc.get_potential_energy()

    path = bandpath('WGKL', atoms.get_cell(), npoints=20)
    bs_dft = get_band_structure(atoms=atoms, calc=calc, path=path,
                                reference=evbm)
    write_json('bs_dft.json', bs_dft)


bs_dft = read_json('bs_dft.json')
bs_dft.plot(filename='bs_dft.png', show=False, emax=bs_dft.reference + 10,
            emin=bs_dft.reference - 20)

dpbs = DftbPlusBandStructure(Hamiltonian_SCC='Yes',
                             Hamiltonian_OrbitalResolvedSCC='No',
                             Hamiltonian_MaxAngularMomentum_='',
                             Hamiltonian_MaxAngularMomentum_N='p',
                             Hamiltonian_MaxAngularMomentum_O='p',
                             Hamiltonian_PolynomialRepulsive='SetForAll {Yes}',
                             Hamiltonian_SpinPolarisation='Colinear {',
                             Hamiltonian_SpinPolarisation_UnpairedElectrons=1,
                             Hamiltonian_SpinConstants_='',
                             Hamiltonian_SpinConstants_O=-0.028,
                             Hamiltonian_SpinConstants_N=-0.026)

dpbs.add_reference_bandstructure(bs_dft, atoms='NO.traj', kpts_scf=(3, 3, 3),
                                 reference_level='vbm', nsemicore=0, weight=1.,
                                 distribution={'type': 'Boltzmann', 'kBT': 1.5})

elements = ['N', 'O']
atoms = []
for element in elements:
    occ_2p = 3 if element == 'N' else 4
    atom = AtomicDFT(element,
                    configuration='[He] 2s2 2p%d' % occ_2p,
                    valence=['2s', '2p'],
                    xc='LDA',
                    scalarrel=False,
                    confinement=PowerConfinement(r0=40., s=4),
                    txt='atomic.out')
    atom.run()
    atom.info = {}
    atom.info['eigenvalues'] = {nl: atom.get_eigenvalue(nl) for nl in atom.valence}
    U_p = atom.get_hubbard_value('2p', scheme='central', maxstep=1.)
    atom.info['hubbardvalues'] = {'s': U_p}
    atom.info['occupations'] = {'2s': 2, '2p': occ_2p}
    atoms.append(atom)

sk_kwargs = {'N': {'default': 380, 'O-O': 420}, 'rmin': 0.4, 'dr': 0.02,
             'stride': 5, 'superposition': 'density', 'xc': 'LDA'}
confopt = ConfinementOptimizer(*atoms, **sk_kwargs)

initial_guess = {}
for element in elements:
    rcov = covalent_radii[atomic_numbers[element]] / Bohr
    initial_guess['%s_2s,%s_2p' % (element, element)] = \
                  PowerConfinement(r0=2 * rcov, s=2, adjustable=['r0'])
    initial_guess['%s_n' % element] = PowerConfinement(r0=3 * rcov, s=2,
                                                       adjustable=['r0'])

vconf = confopt.run(dpbs.get_residual, initial_guess=initial_guess, tol=1e-2,
                    method='COBYLA', options={'maxiter':1, 'rhobeg': 0.2})

bs_dftb = dpbs.calculate_bandstructure(bs_dft)
bs_dftb.plot(filename='bs_dftb.png', emax=bs_dftb.reference + 10,
             emin=bs_dftb.reference - 20)
