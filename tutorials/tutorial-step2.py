#!/usr/bin/python3

from ase.units import Ha
from hotcent.atomic_dft import AtomicDFT

atom = AtomicDFT('Si',
                 xc='LDA',
                 configuration='[Ne] 3s2 3p2',  # the electron configuration we want to use
                 valence=['3s', '3p'],  # these will be our valence states
                 scalarrel=False,  # for Si we don't need (scalar) relativistic effects
                 )
atom.run()
atom.plot_Rnl('Si_Rnl_free.png')  # plot the radial parts of the valence orbitals
atom.plot_rho('Si_rho_free.png')  # plot the valence orbital densities and total electron density

print('=======================================')
for nl in ['3s', '3p']:
    e = atom.get_eigenvalue(nl)
    print(nl, '[Ha]:', e, '[eV]:', e * Ha)