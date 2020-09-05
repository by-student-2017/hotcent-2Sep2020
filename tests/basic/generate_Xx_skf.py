import os
from ase.units import Bohr
from ase.data import covalent_radii, atomic_numbers
from hotcent.atomic_dft import AtomicDFT
from hotcent.confinement import PowerConfinement
from hotcent.slako import SlaterKosterTable

#import command #python2
import subprocess #python3
import sys
import ast

if not os.path.isdir('iter000'):
    os.mkdir('iter000')

element = 'Xx'
#elem_data = commands.getoutput("awk '{if($1==\""+str(element)+"\"{print $0}}' table") #python2
elem_data = subprocess.getoutput("awk '{if($1==\""+str(element)+"\"){print $0}}' table") #python3
elem_list = elem_data.split(" | ")

#element = 'Si'
element = elem_list[0]
print("element = ",element)
xc = 'LDA'
#xc = 'GGA_X_PBEᡠ_C_PBE'
print("xc = "+xc)
#configuration = '[Ne] 3s2 3p2'
configuration = elem_list[1].replace("'","")
print("configuration = ",configuration)
#valence = ['3s', '3p']
valence = ast.literal_eval(elem_list[2])
print("valence = ",valence)
#occupations = {'3s': 2, '3p': 2}
occupations = ast.literal_eval(elem_list[3])
print("occupations = ",occupations)
#nls = '3p'
nls = elem_list[4].replace("'","")
print("nls = ",nls)
#rNum = 600
rNum = int(elem_list[5])
print("rNum = ",rNum)

# ------------------------------------
# Compute eigenvalues of the free atom
# ------------------------------------

atom = AtomicDFT(element,
                 xc=xc,
                 configuration=configuration,
                 valence=valence,
                 scalarrel=False,
                 )
atom.run()
eigenvalues = {nl: atom.get_eigenvalue(nl) for nl in valence}

# ---------------------------------------
# Compute Hubbard values of the free atom
# ---------------------------------------

atom = AtomicDFT(element,
                 xc=xc,
                 configuration=configuration,
                 valence=valence,
                 scalarrel=False,
                 confinement=PowerConfinement(r0=40., s=4),
                 )
U = atom.get_hubbard_value(nls, scheme='central', maxstep=1.)
hubbardvalues = {'s': U}

# -------------------------------
# Compute Slater-Koster integrals
# -------------------------------

r_cov = covalent_radii[atomic_numbers[element]] / Bohr
r_wfc = 2 * r_cov
r_rho = 3 * r_cov
atom = AtomicDFT(element,
                 xc=xc,
                 confinement=PowerConfinement(r0=r_wfc, s=2),
                 wf_confinement=PowerConfinement(r0=r_rho, s=2),
                 configuration=configuration,
                 valence=valence,
                 scalarrel=False,
                 )
atom.run()

# Compute Slater-Koster integrals:
rmin, dr, N = 0.4, 0.02, rNum
sk = SlaterKosterTable(atom, atom)
sk.run(rmin, dr, N, superposition='density', xc=xc, stride=2)
sk.write('iter000/%s-%s_no_repulsion.skf' % (element, element),
         eigenvalues=eigenvalues, spe=0., hubbardvalues=hubbardvalues,
         occupations=occupations)
