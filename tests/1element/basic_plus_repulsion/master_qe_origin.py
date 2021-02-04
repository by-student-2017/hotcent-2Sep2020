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
    ''' This method specifies how to conduct a global optimization
    search, accepting a (run_directory, maximum_iterations)
    tuple as required by TANGO.
    '''
    rundir, maxiter = args

    # Create the run directory if it doesn't already exist:
    if not os.path.exists(rundir):
        os.mkdir(rundir)
    os.chdir(rundir)

    # Redirect the standard output and error to local files:
    sys.stdout = open('out.txt', 'a')
    sys.stderr = open('err.txt', 'w')

    try:
        # Create a database with N initial random structures
        # (N is here a global variable set in the __main__
        # function):
        if not os.path.exists('godb.db'):
            prepare_ga(N=N)
            #prepare_ga(splits={(2, 2): 1, (2,): 2, (1,): 1}, N=N)
        # Run the genetic algorithm (kptdensity is another
        # global variable set in the __main__ function):
        run_ga(maxiter, kptdensity=kptdensity)
    except:
        # When using the multiprocessing module
        # we need to explicitly print the error
        # traceback:
        traceback.print_exc()
        sys.stderr.flush()
        raise

    os.chdir('..')
    return


def generator(dbfile='godb.db'):
    ''' This function will be called to generate the initial batch
    of random structures, which will then be partially relaxed
    using the reference method (usually DFT). Here, we simply
    use the random structure generator which will also be
    used in the genetic algorithm searches. We wish to have
    twenty of such random structures (N=20):
    '''
    prepare_ga(dbfile=dbfile, N=20)
    #prepare_ga(dbfile=dbfile, splits={(2, 2): 1, (2,): 2, (1,): 1}, N=10)
    return


if __name__=='__main__':
    # The different elements which are present:
    elements = ['Xx']
    rcov = covalent_radii[atomic_numbers['Xx']]
    # The cutoffs for the repulsive interactions for each element
    # (if None, the default cutoffs will be used):
    #rcuts = None
    rcuts = {'Xx-Xx': 1.5 * 2 * rcov}
    rmins = {'Xx-Xx': 0.6 * 2 * rcov}
    # The powers to include in the polynomial fits of the
    # repulsive interactions for each element (if None, the
    # default ranges will be used, which is powers 2 till 6):
    powers = None
    #powers = {'Xx-Xx': range(2, 4)}
    # The k-point density in reciprocal Angstrom (if None,
    # only the Gamma point will be used):
    #kptdensity = None
    kptdensity = 1.5

    # Setting up the TANGO search (see tango.main for more
    # information on the remaining keyword arguments):
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
                 #kBT=10.,
                 kBT=0.1,
                 comparator=comparator,
                 #max_select=20,
                 )
    # Note: an "iter000" directory must be present in the
    # current working directory and must contain a
    # "Xx-Xx_no_repulsion.skf" file containing the
    # electronic DFTB parameters. Such a file can be
    # generated in different ways, one being the Hotcent
    # code (https://gitlab.com/mvdb/hotcent).

    # First, we run a series of (re)parametrizations based on
    # randomly generated structures (see tutorial/si7_example.md
    # for more info):
    N = 1
    calc.max_select = 2
    for i in range(2):
        calc.run(steps=1, go_steps=0, recalculate_dftb=True,
                run_go=run, number_of_go_runs=2, restart_go=True)
    #N = 1
    #calc.run(steps=1, go_steps=0, recalculate_dftb=True,
    #         run_go=run, number_of_go_runs=2, restart_go=True)

    # Now is the time for the actual global optimization runs,
    # for which we increase N and max_select:
    N = 5
    calc.max_select = 10
    calc.score_limit = 15.
    calc.run(steps=1, go_steps=5, recalculate_dftb=True,
            run_go=run, number_of_go_runs=2, restart_go=True)
    calc.run(steps=1, go_steps=5, recalculate_dftb=True,
            run_go=run, number_of_go_runs=2, restart_go=False)
    calc.run(steps=1, number_of_go_runs=0, recalculate_dftb=True)
    #N = 4
    #calc.run(steps=1, go_steps=10, recalculate_dftb=True,
    #         run_go=run, number_of_go_runs=2, restart_go=True)
