import os
import sys
import traceback
from tango import TANGO
from ase.data import atomic_numbers, covalent_radii
#from tango.calculators import DftbPlusCalculator, CP2KCalculator
#from cp2k_calc import CP2KCalculator
from qe_calc import QECalculator
#from dftbplus_calc import DftbPlusCalculator
#from ga import run_ga, prepare_ga, comparator
from ga import run_ga, prepare_ga, comparator, DftbPlusCalc

def run(args): 
    rundir, maxiter = args

    if not os.path.exists(rundir):
        os.mkdir(rundir)
    os.chdir(rundir)

    sys.stdout = open('out.txt', 'a')
    sys.stderr = open('err.txt', 'w')

    try:
        if not os.path.exists('godb.db'): 
            prepare_ga(splits={(2, 2): 1, (2,): 2, (1,): 1}, N=N)
        run_ga(maxiter, kptdensity=kptdensity)
    except:
        traceback.print_exc()
        sys.stderr.flush()
        raise

    os.chdir('..')
    return

def generator(dbfile='godb.db'):
    prepare_ga(dbfile=dbfile, splits={(2, 2): 1, (2,): 2, (1,): 1}, N=20)
    return


if __name__=='__main__':
    elements = ['Xx', 'Yy']
    rcuts = None
    powers = {'Yy-Xx': range(2, 7), 'Yy-Yy': range(2, 7), 'Xx-Xx': range(2, 7)}
    kptdensity = 3.0

    calc = TANGO(elements,
                 DftCalc=QECalculator,
                 #DftCalc=CP2KCalculator,
                 DftbPlusCalc=DftbPlusCalc,
                 #DftbPlusCalc=DftbPlusCalculator,
                 kptdensity=kptdensity,
                 initial_training='random_vc_relax',
                 generator=generator,
                 maximum_angular_momenta={'Xx': MAM1, 'Yy': MAM2},
                 rcuts=rcuts,
                 powers=powers,
                 update_rcuts=False,
                 kBT=1.,
                 comparator=comparator,
                 max_select=20,
                 )

    #os.environ['ASE_ESPRESSO_COMMAND'] = 'mpirun -np NCORES pw.x -in PREFIX.pwi > PREFIX.pwo'
    #os.environ['ASE_CP2K_COMMAND'] = 'mpirun -np NCORES cp2k_shell.popt'

    N = 1
    for i in range(5):
        calc.run(steps=1, go_steps=0, recalculate_dftb=True,
                 run_go=run, number_of_go_runs=NCORES, restart_go=True)

    N = 20
    calc.max_select = 100
    calc.score_limit = 15.
    calc.run(steps=1, go_steps=50, recalculate_dftb=True,
             run_go=run, number_of_go_runs=NCORES, restart_go=True)
    calc.run(steps=1, go_steps=50, recalculate_dftb=True,
             run_go=run, number_of_go_runs=NCORES, restart_go=False)
    calc.run(steps=1, number_of_go_runs=0, recalculate_dftb=True)
