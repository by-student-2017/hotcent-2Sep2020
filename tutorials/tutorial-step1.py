#!/usr/bin/python3

import os
from ase.build import bulk
from ase.calculators.dftb import Dftb
from ase.optimize import BFGS
from ase.constraints import ExpCellFilter

atoms = bulk('Si', 'diamond')
calc = Dftb(atoms=atoms,
            kpts=(5, 5, 5),
            Hamiltonian_SCC='Yes',  # SCC = self-consistency charges
            Hamiltonian_ShellResolvedSCC='No',  # Use l-dependent Hubbard values?
            Hamiltonian_SCCTolerance=1e-5,  # SCC convergence criterion
            Hamiltonian_MaxSCCIterations=50,
            Hamiltonian_MaxAngularMomentum_Si='p',  # s- and p-states for Si
            Hamiltonian_Charge=0.0,
            Hamiltonian_ReadInitialCharges='No',  # DFTB-equivalent of restarting from saved electron density
            Hamiltonian_Filling='Fermi {',
            Hamiltonian_Filling_empty='Temperature [Kelvin] = 500.0',
            Hamiltonian_PolynomialRepulsive='SetForAll {No}',  # Use polynomial or spline repulsive?
            Hamiltonian_Solver='RelativelyRobust {}',
            )

print('Lattice constant (Angstrom): %.4f (ASE reference, taken from Ashcroft and Mermin)' % (atoms.get_cell()[0,1]*2))
atoms.set_calculator(calc)
f = ExpCellFilter(atoms)
dyn = BFGS(f, logfile='-')
dyn.run(fmax=0.01)
print('Energy (eV):', atoms.get_potential_energy())
print('Lattice constant (Angstrom): %.4f (DFTB with pbc-0-3 parameters)' % (atoms.get_cell()[0,1]*2))