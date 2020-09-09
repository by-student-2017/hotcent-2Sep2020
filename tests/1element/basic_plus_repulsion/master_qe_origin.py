''' Quick test of a Tango run for bulk Xx,
with a Tersoff potential as reference method
for training the DFTB repulsion. First run
the generate_skf.py script to create the
iter000/Xx-Xx_no_repulsion.skf file.
'''
import os
import sys
import traceback
from ase.data import atomic_numbers, covalent_radii
from tango import TANGO
#from cp2k_calc import CP2KCalculator
from qe_calc import QECalculator
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
    prepare_ga(dbfile=dbfile, splits={(2, 2): 1, (2,): 2, (1,): 1}, N=10)
    return


if __name__=='__main__':
    elements = ['Xx']
    rcov = covalent_radii[atomic_numbers['Xx']]
    rcuts = {'Xx-Xx': 1.5 * 2 * rcov}
    rmins = {'Xx-Xx': 0.6 * 2 * rcov}
    powers = {'Xx-Xx': range(2, 4)}
    kptdensity = 1.5

    calc = TANGO(elements,
                 DftCalc=QECalculator,
                 DftbPlusCalc=DftbPlusCalc,
                 kptdensity=kptdensity,
                 initial_training='random_vc_relax',
                 generator=generator,
                 maximum_angular_momenta={'Xx': MAM},
                 mode='exp_poly',
                 rcuts=rcuts,
                 rmins=rmins,
                 powers=powers,
                 kBT=10.,
                 comparator=comparator,
                 max_select=20,
                 )

    N = 1
    calc.run(steps=1, go_steps=0, recalculate_dftb=True,
             run_go=run, number_of_go_runs=2, restart_go=True)
    N = 4
    calc.run(steps=1, go_steps=10, recalculate_dftb=True,
             run_go=run, number_of_go_runs=2, restart_go=True)
