#!/usr/bin/python3

import numpy as np
from ase.build import bulk
from ase.io.jsonio import write_json
#from gpaw import GPAW, PW, Mixer, FermiDirac
#from gpaw.eigensolvers import CG
from ase.calculators.espresso import Espresso

#import command #python2
import subprocess #python3
import sys
import ast
import os

element = 'F'
#xc = 'LDA'
xc = 'PBE'
#elem_data = commands.getoutput("awk '{if($1==\""+str(element)+"\"{print $0}}' qe_table") #python2
elem_data = subprocess.getoutput("awk '{if($1==\""+str(element)+"\"){print $0}}' qe_table") #python3
elem_list = elem_data.split(" | ")

element = elem_list[0]
print("element = ",element)
struct = elem_list[1]
print("struct = ",struct)
nbnd = int(elem_list[2])
print("nbands = ",nbnd)
spath = elem_list[3]
print("path = ",spath)
nkpts = elem_list[4]
print("kpts = ",nkpts)
kpts = ast.literal_eval(nkpts)
La = float(elem_list[5])
print("lattice constant a = ",La)
Lb = float(elem_list[6])
print("lattice constant b = ",Lb)
Lc = float(elem_list[7])
print("lattice constant c = ",Lc)
Lalpha = float(elem_list[8])
print("alpha angle = ",Lalpha)

atoms = bulk(element,struct,a=La,b=Lb,c=Lc,alpha=Lalpha)
#atoms = bulk(element,struct)
# sc,fcc,bcc,tetragonal,bct,hcp,rhmbohedral,orthorhombic
# mlc, diamond,zincblende,rocksalt,cesiumchloride, fluorite, wurtzite

# First perform a regular SCF run
gbrv_pp = {}
ppdir = os.environ['ESPRESSO_PSEUDO']
sym = list(set(atoms.get_chemical_symbols()))

for s in sym:
    for f in os.listdir(ppdir):
        keys = f.split('_')
        if keys[0] == s.lower() and keys[1] == xc.lower():
            gbrv_pp[s] = f

input_data = {'control':{'disk_io': 'low',
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

calc = Espresso(kpts=kpts, pseudopotentials=gbrv_pp,
        tstress=True, tprnfor=True,  # kwargs added to parameters
        input_data=input_data)

atoms.calc = calc
atoms.get_potential_energy()

# Get the valence band maximum
efermi = calc.get_fermi_level()
Nk = len(calc.get_ibz_k_points())
Ns = calc.get_number_of_spins()
eigval = np.array([[calc.get_eigenvalues(kpt=k, spin=s)
                    for k in range(Nk)] for s in range(Ns)])
evbm = np.max(eigval[eigval < efermi])

# Next, a band structure calculation
input_data['control'].update({'calculation':'bands',
    'restart_mode':'restart',
    'verbosity':'high'})
calc.set(kpts={'path': spath, 'npoints': 50},
         input_data=input_data
         )
calc.calculate(atoms)

bs_qe = calc.band_structure()
bs_qe.reference = evbm
bs_qe.plot(filename='bs_qe_data.png', show=False, emax=evbm + 5., emin=evbm - 15.)
write_json('bs_qe.json', bs_qe)
