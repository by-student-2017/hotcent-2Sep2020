import os
import numpy as np
from ase.calculators.calculator import kptdensity2monkhorstpack
from ase.calculators.dftb import Dftb

class DftbPlusCalculator(Dftb):
    def __init__(self, atoms, kpts=(1, 1, 1), use_spline=False,
                 maximum_angular_momenta={}, read_chg=False,
                 label='dftb_run', **extra_dftbplus_kwargs):

        if type(kpts) == float or type(kpts) == int:
            mp = kptdensity2monkhorstpack(atoms, kptdensity=kpts, even=False)
            kpts = tuple(mp)

        s = 'No' if use_spline else 'Yes'
        polyrep = 'SetForAll { %s }' % s

        self.read_chg = read_chg

        parameters = {'Hamiltonian_SCC': 'Yes',
                      'Hamiltonian_OrbitalResolvedSCC': 'Yes',
                      'Hamiltonian_SCCTolerance': '1e-5',
                      'Hamiltonian_MaxSCCIterations': 250,
                      'Hamiltonian_MaxAngularMomentum_': '',
                      'Hamiltonian_Charge': '0.000000',
                      'Hamiltonian_ReadInitialCharges': 'No',
                      'Hamiltonian_Filling': 'Fermi {',
                      'Hamiltonian_Filling_empty': 'Temperature [Kelvin] = 500',
                      'Hamiltonian_PolynomialRepulsive': polyrep,
                      'Hamiltonian_Eigensolver': 'RelativelyRobust {}',
                      }

        symbols = atoms.get_chemical_symbols()
        unique_symbols = list(set(symbols))
        for s in unique_symbols:
              key = 'Hamiltonian_MaxAngularMomentum_%s' % s
              maxmom = maximum_angular_momenta[s]
              parameters[key] = 'spd'[maxmom].__repr__()

        parameters.update(extra_dftbplus_kwargs)

        Dftb.__init__(self, label=label, kpts=kpts, atoms=atoms, **parameters)

    def exit(self):
        pass

    def write_dftb_in(self, filename):
        if self.read_chg and os.path.exists('charges.bin'):
            read_charges = 'Yes'
        else:
            read_charges = 'No'
        self.parameters['Hamiltonian_ReadInitialCharges'] = read_charges
        Dftb.write_dftb_in(self, filename)

    def get_potential_energy(self, atoms=None, force_consistent=True):
        return Dftb.get_potential_energy(self, atoms=atoms)
