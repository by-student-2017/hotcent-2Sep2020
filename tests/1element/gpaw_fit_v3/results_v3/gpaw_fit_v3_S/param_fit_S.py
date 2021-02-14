#!/usr/bin/python3

from ase.io.jsonio import read_json
from ase.units import Ha
from ase.units import Bohr
from ase.build import bulk
from ase.data import atomic_numbers, covalent_radii
from hotcent.atomic_dft import AtomicDFT
from hotcent.confinement import PowerConfinement
from hotcent.tools import ConfinementOptimizer, DftbPlusBandStructure

#import command #python2
import subprocess #python3
import sys
import ast

#read cif file
import ase.io

element = 'S'
#elem_data = commands.getoutput("awk '{if($1==\""+str(element)+"\"{print $0}}' gpaw_table") #python2
elem_data = subprocess.getoutput("awk '{if($1==\""+str(element)+"\"){print $0}}' gpaw_table") #python3
elem_list = elem_data.split(" | ")

element = elem_list[0]
print("element = ",element)
struct = elem_list[1]
print("struct = ",struct)
nkpts = elem_list[4]
print("kpts = ",nkpts)
kpts = ast.literal_eval(nkpts)

#elem_data = commands.getoutput("awk '{if($1==\""+str(element)+"\"{print $0}}' table_fit") #python2
elem_data = subprocess.getoutput("awk '{if($1==\""+str(element)+"\"){print $0}}' table_fit") #python3
elem_list = elem_data.split(" | ")

#element = 'Si'
element = elem_list[0]
print("element = ",element)
#xc = 'LDA'
xc = 'GGA_X_PBE+GGA_C_PBE'
#configuration = '[Ne] 3s2 3p2 3d0'
configuration = elem_list[1].replace("'","")
print("configuration = ",configuration)
#valence = ['3s', '3p', '3d']
valence = ast.literal_eval(elem_list[2])
print("valence = ",valence)
#occupations = {'3s': 2, '3p': 2, '3d': 0}
occupations = ast.literal_eval(elem_list[3])
print("occupations = ",occupations)
#nls = '3p'
nls = elem_list[4].replace("'","")
print("nls = ",nls)
#rNum = 600
rNum = int(elem_list[5])
print("rNum = ",rNum)

#lmax = valence[-1][1]
lmax = str(elem_list[6])
print("lmax = ",lmax)

# Setting up the atomic DFT instance(s) and
# calculating the eigenvalues and Hubbard values
atom = AtomicDFT(element,
                 configuration=configuration,
                 valence=valence,
                 xc=xc,
                 scalarrel=True,
                 mix=0.005,
                 maxiter=50000,
                 confinement=PowerConfinement(r0=40., s=4),
                 txt=None)
atom.run()
atom.info = {}
atom.info['eigenvalues'] = {nl: atom.get_eigenvalue(nl) for nl in atom.valence}

if (element == "H"):
  # https://webbook.nist.gov/cgi/cbook.cgi?ID=C12385136&Mask=20
  # hubbard value = U = IE - EA = 13.59844 - 0.75497 = 12.84347 [eV]
  U_p = 12.84347/Ha
else:
  U_p = atom.get_hubbard_value(nls, scheme='central', maxstep=1.)
atom.info['hubbardvalues'] = {'s': U_p}
atom.info['occupations'] = occupations

# Creating a DFTB+ band structure evaluator and
# supplying it with a reference (DFT) band structure
dpbs = DftbPlusBandStructure(Hamiltonian_SCC='Yes',
                             Hamiltonian_OrbitalResolvedSCC='No',
                             Hamiltonian_MaxAngularMomentum_='',
                             Hamiltonian_MaxAngularMomentum_S=lmax,
                             Hamiltonian_PolynomialRepulsive='SetForAll {Yes}')

bs_gpaw = read_json('bs_gpaw.json')  # the reference band structure (DFT)
atoms = ase.io.read("./primCIFs/S.cif")
#atoms = bulk(element,struct)
# see hotcent.tools.DftbPlusBandStructure for more information
# on the various keyword arguments used below
dpbs.add_reference_bandstructure(bs_gpaw, atoms=atoms, kpts_scf=kpts,
                                 reference_level='vbm', nsemicore=0, weight=1.,
                                 distribution={'type': 'Boltzmann', 'kBT': 1.5}
                                 )

# Setting up and running the actual optimizer
# (the keyword arguments are known from hotcent.slako.SlaterKosterTable)
confopt = ConfinementOptimizer(atom, N=500, rmin=0.4, dr=0.02, stride=4,
                               superposition='density', xc=xc)

# The initial confinement parameters are the same as in Tutorial #1.
# The additional 'adjustable' keyword argument serves to indicate
# which parameters are allowed to vary. Here we keep the quadratic
# form (s=2) and treat the confinement radii r0 as variable.
rcov = covalent_radii[atomic_numbers[element]] / Bohr

if len(valence) == 1 :
  initial_guess = {element+'_'+valence[0]: PowerConfinement(r0=2 * rcov, s=2, adjustable=['r0']),
                   element+'_n': PowerConfinement(r0=3 * rcov, s=2, adjustable=['r0'])}
elif len(valence) == 2 :
  initial_guess = {element+'_'+valence[0]: PowerConfinement(r0=2 * rcov, s=2, adjustable=['r0']),
                   element+'_'+valence[1]: PowerConfinement(r0=2 * rcov, s=2, adjustable=['r0']),
                   element+'_n': PowerConfinement(r0=3 * rcov, s=2, adjustable=['r0'])}
elif len(valence) == 3 :
  initial_guess = {element+'_'+valence[0]: PowerConfinement(r0=2 * rcov, s=2, adjustable=['r0']),
                   element+'_'+valence[1]: PowerConfinement(r0=2 * rcov, s=2, adjustable=['r0']),
                   element+'_'+valence[2]: PowerConfinement(r0=2 * rcov, s=2, adjustable=['r0']),
                   element+'_n': PowerConfinement(r0=3 * rcov, s=2, adjustable=['r0'])}
# Note that the 'Si_3s,Si_3p' combination indicates that the same
# confinement potential is used for both states. The Si_3d confinement
# (as well as the density confinement) is defined separately and
# their r0 parameters are hence allowed to vary independently.

# Only two iterations are performed here to limit the computational time.
# The confinement therefore hardly changes from the initial guess.
# In real scenarios, typically on the order of 100 iterations are required.
vconf = confopt.run(dpbs.get_residual, initial_guess=initial_guess, tol=1e-2,
                    method='COBYLA', options={'maxiter': 100, 'rhobeg': 0.2})
# See scipy.optimize.minimize for more information on available optimization methods:
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html

# Make the DFTB band structure plot on the basis of the latest Si-Si.skf file
# (generated with the confinement parameters yielding the lowest residual)
bs_dftb = dpbs.calculate_bandstructure(bs_gpaw)
bs_dftb.plot(filename='bs_dftb.png', emax=bs_dftb.reference + 5, emin=bs_dftb.reference - 15)
bs_gpaw.plot(filename='bs_gpaw.png', emax=bs_dftb.reference + 5, emin=bs_dftb.reference - 15)
