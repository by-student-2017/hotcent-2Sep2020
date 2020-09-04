#!/usr/bin/python3

from ase.units import Ha
from hotcent.atomic_dft import AtomicDFT
from hotcent.confinement import PowerConfinement

atom = AtomicDFT('Si',
                 xc='LDA',
                 configuration='[Ne] 3s2 3p2',
                 valence=['3s', '3p'],
                 scalarrel=False,
                 # Add a very weak confinement potential to aid anion convergence:
                 confinement=PowerConfinement(r0=40., s=4),
                 )

# Like above, we use a central difference scheme
# with changes in the occupation of 1 |e|
U = atom.get_hubbard_value('3p', scheme='central', maxstep=1)
print('=======================================')
print('U [Ha]:', U, '[eV]:', U * Ha)