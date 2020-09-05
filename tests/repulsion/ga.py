import os
from random import random
import numpy as np
from time import time
from ase import Atoms
from ase.io import write, read
from ase.io.formats import UnknownFileTypeError
from ase.data import atomic_numbers
from ase.build import niggli_reduce
from ase.ga import get_raw_score, set_raw_score
from ase.ga.data import DataConnection, PrepareDB
from ase.ga.population import Population
from ase.ga.utilities import (get_all_atom_types, closest_distances_generator,
                              atoms_too_close)
from ase.ga.offspring_creator import OperationSelector
from ase.ga.ofp_comparator import OFPComparator
from ase.ga.bulk_utilities import CellBounds
from ase.ga.bulk_startgenerator import StartGenerator
from ase.ga.bulk_crossovers import CutAndSplicePairing
from ase.ga.bulk_mutations import *
from ase.ga.standardmutations import RattleMutation
from tango.relax_utils import push_apart, relax_precon, finalize
from tango.calculators import DftbPlusCalculator

comparator = OFPComparator(n_top=None, dE=1.0, cos_dist_max=5e-2,
                           rcut=10., binwidth=0.05, pbc=[True]*3,
                           sigma=0.1, nsigma=4, recalculate=False)

class DftbPlusCalc(DftbPlusCalculator):
    def __init__(self, *args, **kwargs):
        kwargs['Hamiltonian_SCC'] = 'No'
        kwargs['Hamiltonian_ShellResolvedSCC'] = 'No'
        kwargs['Hamiltonian_OrbitalResolvedSCC'] = 'No'
        kwargs['maximum_angular_momenta'] = {'Si': 1}
        DftbPlusCalculator.__init__(self, *args, **kwargs)

def penalize(t):
    # penalize explosion:
    raw_score = get_raw_score(t)
    max_volume_per_atom = 50.
    if t.get_volume()/len(t) >= max_volume_per_atom:
        raw_score -= 1e9
    set_raw_score(t, raw_score)

def singlepoint(t, kptdensity=1.5):
    if get_raw_score(t) < -1e5:
        return t
    try:
        calc = DftbPlusCalc(t, kpts=kptdensity, use_spline=True, read_chg=True)
        t.set_calculator(calc)
        E = t.get_potential_energy()
        F = t.get_forces()
        S = t.get_stress()
        finalize(t, energy=E, forces=F, stress=S)
        penalize(t)
    except (RuntimeError, IOError):
        print('Warning: problems with singlepoint recalculation')
        finalize(t, energy=1e9, forces=None, stress=None)
    return t

def relax_one(t, kptdensity=1.5):
    cellbounds = CellBounds(bounds={'phi':[0.1*180., 0.9*180.],
                                    'chi':[0.1*180., 0.9*180.],
                                    'psi':[0.1*180., 0.9*180.],
                                    'a':[1.5,20], 'b':[1.5,20], 'c':[1.5,20]})

    if not cellbounds.is_within_bounds(t.get_cell()):
        print('Candidate outside cellbounds -- skipping')
        finalize(t, energy=1e9, forces=None, stress=None)
        return t

    pos = t.get_positions()
    numbers = list(set(t.get_atomic_numbers())) 
    blmin = closest_distances_generator(numbers, 0.5)
    t = push_apart(t, blmin, variable_cell=True)

    print('Starting relaxation', flush=True)
    clock = time()
    t.wrap()
    calc = DftbPlusCalc(t, kpts=0.66*kptdensity, use_spline=True, read_chg=True)

    try:
        t = relax_precon(t, calc, fmax=5e-1, smax=5e-2, variable_cell=True,
                         optimizer='LBFGS', a_min=1e-4, cellbounds=cellbounds,
                         logfile='opt_first.log', trajfile='opt_first.traj')
    except IOError as err:
        # probably convergence problem
        print(err.message)
    except TypeError as err:
        if 'Cannot cast array data' in err.message:
            # Rare precon failure
            print(err.message)
        else:
            raise
    except RuntimeError as err:
        print(err.message)

    del t.constraints
    t.wrap()
    calc = DftbPlusCalc(t, kpts=kptdensity, use_spline=True, read_chg=True)
    try:
        t = relax_precon(t, calc, fmax=5e-1, smax=5e-2, variable_cell=True,
                         optimizer='LBFGS', a_min=1e-4, cellbounds=cellbounds,
                         logfile='opt.log', trajfile='opt.traj')
    except (IOError, TypeError, RuntimeError) as err:
        # probably convergence problem
        print(err.message)
        if isinstance(err, TypeError):
            if 'Cannot cast array data' not in err.message:
                raise
        elif isinstance(err, RuntimeError):
            if 'dftb_my in . returned an error:' not in err.message:
                raise
        try:
            t = read('opt.traj@-1')
            energy = t.get_potential_energy()
            forces = t.get_forces()
            stress = t.get_stress()
        except UnknownFileTypeError as err:
            print(err.message)
            energy, forces, stress = (1e9, None, None)
        finalize(t, energy=energy, forces=forces, stress=stress)

    print('Relaxing took %.3f seconds.' % (time() - clock), flush=True)
    os.system('mv opt_first.traj prev_first.traj')
    os.system('mv opt_first.log prev_first.log')
    os.system('mv opt.log prev.log')
    os.system('mv opt.traj prev.traj')
    penalize(t)
    return t

def run_ga(n_to_test, kptdensity=1.5):
    population_size = 10
    da = DataConnection('godb.db')
    atom_numbers_to_optimize = da.get_atom_numbers_to_optimize()
    n_to_optimize = len(atom_numbers_to_optimize)
    slab = da.get_slab()
    all_atom_types = get_all_atom_types(slab, atom_numbers_to_optimize)
    blmin = closest_distances_generator(all_atom_types, 0.05)  # 0.5

    # defining genetic operators:
    mutation_probability = 0.75
    pairing = CutAndSplicePairing(blmin, p1=1., p2=0., minfrac=0.15,
                                  use_tags=False)
    cellbounds = CellBounds(bounds={'phi':[0.2*180., 0.8*180.],
                                    'chi':[0.2*180., 0.8*180.],
                                    'psi':[0.2*180., 0.8*180.]})
    strainmut = StrainMutation(blmin, stddev=0.7, cellbounds=cellbounds,
                               use_tags=False)
    blmin_soft = closest_distances_generator(all_atom_types, 0.1)
    softmut = SoftMutation(blmin_soft, bounds=[2., 5.], use_tags=False)
    rattlemut = RattleMutation(blmin, n_to_optimize, rattle_prop=0.8,
                               rattle_strength=2.5, use_tags=False)
    mutations = OperationSelector([4., 4., 2],
                                  [softmut, strainmut, rattlemut])

    if True:
        # recalculate raw scores
        structures = da.get_all_relaxed_candidates()
        for atoms in structures:
            atoms = singlepoint(atoms, kptdensity=kptdensity)
            da.c.delete([atoms.info['relax_id']])
            if 'data' not in atoms.info:
                atoms.info['data'] = {}
            da.add_relaxed_step(atoms)
        print('Finished recalculating raw scores')

    # relaxing the initial candidates:
    while da.get_number_of_unrelaxed_candidates() > 0:
        a = da.get_an_unrelaxed_candidate()
        a.wrap()
        a = relax_one(a, kptdensity=kptdensity)
        da.add_relaxed_step(a)

    # create the population
    population = Population(data_connection=da,
                            population_size=population_size,
                            comparator=comparator,
                            logfile='log.txt')

    current_pop = population.get_current_population()
    strainmut.update_scaling_volume(current_pop, w_adapt=0.5, n_adapt=4)
    pairing.update_scaling_volume(current_pop, w_adapt=0.5, n_adapt=4)

    # Test n_to_test new candidates
    ga_raw_scores = []
    step = 0
    for step in range(n_to_test):
        print('Starting configuration number %d' % step, flush=True)

        clock = time()
        a3 = None
        r = random()
        if r > mutation_probability:
            while a3 is None:
                a1, a2 = population.get_two_candidates()
                a3, desc = pairing.get_new_individual([a1, a2])
        else:
            while a3 is None:
                a1 = population.get_one_candidate()
                a3, desc = mutations.get_new_individual([a1])

        dt = time() - clock
        op = 'pairing' if r > mutation_probability else 'mutating'
        print('Time for %s candidate(s): %.3f' % (op, dt), flush=True)

        a3.wrap()
        da.add_unrelaxed_candidate(a3, description=desc)

        a3 = relax_one(a3, kptdensity=kptdensity)
        da.add_relaxed_step(a3)

        # Various updates:
        population.update()
        current_pop = population.get_current_population()

        if step % 10 == 0:
            strainmut.update_scaling_volume(current_pop, w_adapt=0.5,
                                            n_adapt=4)
            pairing.update_scaling_volume(current_pop, w_adapt=0.5, n_adapt=4)
            write('current_population.traj', current_pop)

        # Print out information for easy analysis/plotting afterwards:
        if r > mutation_probability:
            print('Step %d %s %.3f %.3f %.3f' % (step, desc,\
                   get_raw_score(a1), get_raw_score(a2), get_raw_score(a3)))
        else:
            print('Step %d %s %.3f %.3f' % (step, desc,\
                   get_raw_score(a1), get_raw_score(a3)))

        print('Step %d highest raw score in pop: %.3f' % \
              (step, get_raw_score(current_pop[0])))
        ga_raw_scores.append(get_raw_score(a3))
        print('Step %d highest raw score generated by GA: %.3f' % \
              (step, max(ga_raw_scores)))

    emin = population.pop[0].get_potential_energy()
    print('GA finished after step %d' % step)
    print('Lowest energy = %8.3f eV' % emin, flush=True)
    write('all_candidates.traj', da.get_all_relaxed_candidates())
    write('current_population.traj', population.get_current_population())


def prepare_ga(dbfile='godb.db', splits={(2,):1}, N=20):

    blocks = [('Si', 4)]  # the building blocks
    volume = 8. * 4 # volume in angstrom^3

    l = [list(Atoms(block).numbers)*count for block, count in blocks]
    stoichiometry = [int(item) for sublist in l for item in sublist]
    atom_numbers = list(set(stoichiometry))

    blmin = closest_distances_generator(atom_numbers=atom_numbers,
                                        ratio_of_covalent_radii=0.6)

    cellbounds = CellBounds(bounds={'phi':[0.2*180., 0.8*180.],
                                    'chi':[0.2*180., 0.8*180.],
                                    'psi':[0.2*180., 0.8*180.],
                                    'a':[2,8], 'b':[2,8], 'c':[2,8]})

    # create the starting population
    sg = StartGenerator(blocks, blmin, volume, cellbounds=cellbounds,
                        splits=splits)

    # create the database to store information in
    da = PrepareDB(db_file_name=dbfile, stoichiometry=stoichiometry)

    for i in range(N):
        a = sg.get_new_candidate()
        niggli_reduce(a)
        da.add_unrelaxed_candidate(a)

    return
