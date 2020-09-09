''' Wrapper around CP2K Calculator in ASE '''
from ase.calculators.calculator import all_changes, kptdensity2monkhorstpack
from ase.calculators.cp2k import CP2K
from ase.build import niggli_reduce
import numpy as np

template = '''&FORCE_EVAL
  METHOD Quickstep
  STRESS_TENSOR ANALYTICAL
  &DFT
    &MGRID
      CUTOFF 600
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-10
      EXTRAPOLATION USE_GUESS
    &END QS
    &KPOINTS
      FULL_GRID .TRUE.
      SCHEME MONKHORST-PACK __MPMESH__
    &END KPOINTS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      MAX_SCF 300
      ADDED_MOS 20
      &SMEAR 
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 500
      &END SMEAR
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.1
        BETA 0.05
        NBROYDEN 12
      &END MIXING
      &PRINT
        &RESTART OFF
        &END RESTART
      &END PRINT
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
&END FORCE_EVAL
'''

class CP2KCalculator(CP2K):
    def __init__(self, atoms, kpts=(1, 1, 1), label='cp2k_run',
                 run_type=None, **extra_cp2k_kwargs):
        try:
            niggli_reduce(atoms)
        except RuntimeError:
            pass

        kwargs = {k:None for k in CP2K.default_parameters}
        kwargs['debug'] = False
        kwargs['potential_file'] = 'GTH_POTENTIALS'
        kwargs['pseudo_potential'] = 'GTH-PBE'
        kwargs['basis_set_file'] = 'BASIS_MOLOPT'
        kwargs['basis_set'] = 'DZVP-MOLOPT-SR-GTH'

        if type(kpts) == float or type(kpts) == int:
            mp = kptdensity2monkhorstpack(atoms, kptdensity=kpts, even=False)
            kpts = tuple(mp)
        inp = template.replace('__MPMESH__', ' '.join(map(str, kpts)))
        kwargs['inp'] = inp

        kwargs.update(extra_cp2k_kwargs)

        CP2K.__init__(self, atoms=atoms, label=label, **kwargs)

    def calculate(self, atoms=None, properties=None,
                  system_changes=all_changes):
        try:
            niggli_reduce(atoms)
        except RuntimeError:
            pass 
        CP2K.calculate(self, atoms, properties, system_changes)

    def exit(self):
        CP2K.__del__(self)
