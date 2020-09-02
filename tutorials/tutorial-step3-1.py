#!/usr/bin/python3

from ase.units import Ha
from hotcent.atomic_dft import AtomicDFT
from hotcent.confinement import PowerConfinement

energies = {}
for occupation, kind in zip([2, 1, 3], ['neutral', 'cation', 'anion']):
    atom = AtomicDFT('Si',
                     xc='LDA',
                     configuration='[Ne] 3s2 3p%d' % occupation,
                     valence=['3s', '3p'],
                     scalarrel=False,
                     # Add a very weak confinement potential to aid anion convergence:
                     confinement=PowerConfinement(r0=40., s=4),
                     )
    atom.run()
    energies[kind] = atom.get_energy()

EA = energies['neutral'] - energies['anion']  # the electron affinity
IE = energies['cation'] - energies['neutral']  # the ionization energy
U = IE - EA  # the Hubbard value

print('=======================================')
for value, label in zip([EA, IE, U], ['EA', 'IE', 'U']):
    print(label, '[Ha]:', value, '[eV]:', value * Ha)