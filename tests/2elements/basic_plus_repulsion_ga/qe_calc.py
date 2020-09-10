''' Wrapper around the Espresso Calculator in ASE.
Uses GBRV pseudopotentials, which should be
located in your $ESPRESSO_PSEUDO directory.
'''
import os
from ase.calculators.espresso import Espresso


class QECalculator(Espresso):
    def __init__(self, atoms, kpts=(1,1,1), run_type=None, input_data={},
                 xc='PBE', command=None, tprnfor=True, tstress=True, **kwargs):
        if command is None:
            command = 'mpirun -np $NCORES pw.x -in PREFIX.pwi > PREFIX.pwo'

        gbrv_pp = {}
        ppdir = os.environ['ESPRESSO_PSEUDO']
        sym = list(set(atoms.get_chemical_symbols()))

        for s in sym:
            for f in os.listdir(ppdir):
                keys = f.split('_')
                if keys[0] == s.lower() and keys[1] == xc.lower():
                    gbrv_pp[s] = f

        indat = {'control':{'disk_io': 'none',
                            'restart_mode': 'from_scratch'},
                 'system':{'ecutwfc': 40.,
                           'ecutrho': 200.,
                           'input_dft': xc,
                           'occupations': 'smearing',
                           'smearing': 'gaussian',
                           'degauss': 0.001},
                 'electrons':{'electron_maxstep': 250,
                              'scf_must_converge': False,
                              'mixing_beta': 0.1,
                              'conv_thr': 1e-7},
                 }

        indat.update(input_data)
        Espresso.__init__(self, kpts=kpts, pseudopotentials=gbrv_pp,
                          input_data=indat, tprnfor=tprnfor, tstress=tstress,
                          **kwargs)
        self.command = command

    def get_potential_energy(self, atoms=None, force_consistent=True):
        return Espresso.get_potential_energy(self, atoms=atoms)

    def exit(self):
        pass
