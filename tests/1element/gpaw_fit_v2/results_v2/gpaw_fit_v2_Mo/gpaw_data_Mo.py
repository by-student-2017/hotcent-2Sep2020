#!/usr/bin/python3

import numpy as np
from ase.build import bulk
from ase.io.jsonio import write_json
from gpaw import GPAW, PW, Mixer, FermiDirac
from gpaw.eigensolvers import CG

#import command #python2
import subprocess #python3
import sys
import ast

element = 'Mo'
#xc = 'LDA'
xc = 'PBE'
#elem_data = commands.getoutput("awk '{if($1==\""+str(element)+"\"{print $0}}' gpaw_table") #python2
elem_data = subprocess.getoutput("awk '{if($1==\""+str(element)+"\"){print $0}}' gpaw_table") #python3
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

# First perform a regular SCF run
calc = GPAW(mode=PW(400),
            maxiter=1500,
            spinpol=True,
            kpts=kpts,
            xc=xc,
            txt='-',
            occupations=FermiDirac(0.02),
            mixer=Mixer(0.01, 11, 100.0),
            )

atoms = bulk(element,struct,a=La,b=Lb,c=Lc,alpha=Lalpha)
# sc,fcc,bcc,tetragonal,bct,hcp,rhmbohedral,orthorhombic
# mlc, diamond,zincblende,rocksalt,cesiumchloride, fluorite, wurtzite
atoms.set_calculator(calc)
atoms.get_potential_energy()

# Get the valence band maximum
efermi = calc.get_fermi_level()
Nk = len(calc.get_ibz_k_points())
Ns = calc.get_number_of_spins()
eigval = np.array([[calc.get_eigenvalues(kpt=k, spin=s)
                    for k in range(Nk)] for s in range(Ns)])
evbm = np.max(eigval[eigval < efermi])

# Next, a band structure calculation
calc.set(nbands=nbnd,  # 4 occupied and 4 unoccupied bands
         fixdensity=True,
         eigensolver=CG(niter=5),
         symmetry='off',
         kpts={'path': spath, 'npoints': 50},
         convergence={'bands': 'all'},
         )
calc.get_potential_energy()

bs_gpaw = calc.band_structure()
#bs_gpaw.reference = evbm
bs_gpaw_reference = evbm
bs_gpaw.plot(filename='bs_gpaw_data.png', show=False, emax=evbm + 5., emin=evbm - 15.)
write_json('bs_gpaw.json', bs_gpaw)
