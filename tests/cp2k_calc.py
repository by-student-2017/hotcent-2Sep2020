from ase.calculators.cp2k import CP2K

inp = '''
&GLOBAL
  PRINT_LEVEL LOW 
&END GLOBAL
&FORCE_EVAL
  METHOD FIST
  &MM 
   &FORCEFIELD
     IGNORE_MISSING_CRITICAL_PARAMS .TRUE.
     &NONBONDED
       &TERSOFF
         ATOMS Si Si
       &END TERSOFF
     &END NONBONDED
   &END FORCEFIELD
   &POISSON
     PERIODIC XYZ
     &EWALD
       GMAX 11
     &END EWALD
   &END POISSON
  &END MM
&END FORCE_EVAL
'''

class CP2KCalculator(CP2K):
    def __init__(self, atoms, kpts=(1,1,1), run_type=None):
        '''
        Initializes the CP2K calculator, accepting the 'atoms',
        'kpts' and 'run_type' arguments as required by TANGO
        (see tango.main.py). 
        '''
        kwargs = {k:None for k in CP2K.default_parameters}
        kwargs['inp'] = inp
        CP2K.__init__(self, atoms=atoms, label='cp2k_run',
                      command='mpirun -np $NCORES cp2k_shell.popt', **kwargs)

    def exit(self):
        ''' 
        An exit() method is required by TANGO, specifying 
        what to do when a calculation is over.
        In this CP2K interface, we need to explicitly
        terminate the interactive shell.
        '''
        CP2K.__del__(self)
