#!/usr/bin/python3

import numpy as np
from ase.build import bulk
from ase.io.jsonio import write_json
from gpaw import GPAW, PW, Mixer, FermiDirac
from gpaw.eigensolvers import CG


# First perform a regular SCF run
calc = GPAW(mode=PW(400),
            maxiter=200,
            spinpol=False,
            kpts=(5, 5, 5),
            xc='LDA',
            txt='-',
            occupations=FermiDirac(0.02),
            mixer=Mixer(0.05, 8, 100),
            )

atoms = bulk('Si', 'diamond')
atoms.set_calculator(calc)
atoms.get_potential_energy()

# Get the valence band maximum
efermi = calc.get_fermi_level()
Nk = len(calc.get_ibz_k_points())
Ns = calc.get_number_of_spins()
eigval = np.array([[calc.get_eigenvalues(kpt=k, spin=s)
                    for k in range(Nk)] for s in range(Ns)])
evbm = np.max(eigval[eigval < efermi])

# Next, a band structure calculation
calc.set(nbands=8,  # 4 occupied and 4 unoccupied bands
         fixdensity=True,
         eigensolver=CG(niter=5),
         symmetry='off',
         kpts={'path': 'LGXUG', 'npoints': 50},
         convergence={'bands': 'all'},
         )
calc.get_potential_energy()

bs_gpaw = calc.band_structure()
#bs_gpaw.reference = evbm
bs_gpaw_reference = evbm
bs_gpaw.plot(filename='bs_gpaw.png', show=False, emax=evbm + 5., emin=evbm - 15.)
write_json('bs_gpaw.json', bs_gpaw)