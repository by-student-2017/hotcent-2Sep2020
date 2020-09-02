#!/usr/bin/python3

from ase.units import Bohr
from ase.data import atomic_numbers, covalent_radii
from hotcent.atomic_dft import AtomicDFT
from hotcent.slako import SlaterKosterTable
from hotcent.confinement import PowerConfinement

# Define standard, rule-of-thumb confinement potentials
rcov = covalent_radii[atomic_numbers['Si']] / Bohr
conf = PowerConfinement(r0=3 * rcov, s=2)
wf_conf = {'3s': PowerConfinement(r0=2 * rcov, s=2),
           '3p': PowerConfinement(r0=2 * rcov, s=2)}

atom = AtomicDFT('Si',
                 xc='LDA',
                 configuration='[Ne] 3s2 3p2',
                 valence=['3s', '3p'],
                 scalarrel=False,
                 confinement=conf,
                 wf_confinement=wf_conf,
                 )
atom.run()
atom.plot_Rnl('Si_Rnl_conf.png')
atom.plot_rho('Si_rho_conf.png')

rmin, dr, N = 0.4, 0.02, 600
sk = SlaterKosterTable(atom, atom)
sk.run(rmin, dr, N, superposition='density', xc='LDA')

# Write the Slater-Koster tables to file (without two-body repulsion at this point).
# This file also stores the eigenvalues, Hubbardvalues, occupations, as well as the
# so-called spin-polarization error (the magnetization energy of the atom, which we
# don't need to consider here).
sk.write('Si-Si_no_repulsion.skf', eigenvalues={'3s': -0.397837, '3p': -0.153163},
         hubbardvalues={'s': 0.244242}, occupations={'3s': 2, '3p': 2}, spe=0.)
sk.plot('Si-Si_slako.png')
print('Done!')