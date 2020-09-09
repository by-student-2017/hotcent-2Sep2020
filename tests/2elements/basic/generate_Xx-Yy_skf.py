from ase.units import Bohr
from ase.data import covalent_radii, atomic_numbers
from hotcent.atomic_dft import AtomicDFT
from hotcent.confinement import PowerConfinement
from hotcent.slako import SlaterKosterTable

#import command #python2
import subprocess #python3
import sys
import ast

element1 = 'Xx'
elem1_data = subprocess.getoutput("awk '{if($1==\""+str(element1)+"\"){print $0}}' table") #python3
elem1_list = elem1_data.split(" | ")

element2 = 'Yy'
elem2_data = subprocess.getoutput("awk '{if($1==\""+str(element2)+"\"){print $0}}' table") #python3
elem2_list = elem2_data.split(" | ")

#elements = ['Mg', 'O']
element1 = elem1_list[0]
element2 = elem2_list[0]
elements = [element1,element2]
#xc = 'LDA'
xc = 'GGA_X_PBE+GGA_C_PBE'
#configurations = {'Mg': '[Ne] 3s2', 'O': '[He] 2s2 2p4'}
configuration1 = elem1_list[1].replace("'","")
configuration2 = elem2_list[1].replace("'","")
configurations = {element1 : configuration1, element2 : configuration2}
#valences = {'Mg': ['3s'], 'O': ['2s', '2p']}
valence1 = ast.literal_eval(elem1_list[2])
valence2 = ast.literal_eval(elem2_list[2])
valences = {element1 : valence1, element2 : valence2}
#occupations = {'Mg':{'2s': 2}, 'O': {'2s': 2, '2p': 4}}
occupations1 = ast.literal_eval(elem1_list[3])
occupations2 = ast.literal_eval(elem2_list[3])
occupations = {element1 : occupations1, element2 : occupations2}
eigenvalues, hubbardvalues, atoms = {}, {}, {}

nls1 = elem1_list[4].replace("'","")
nls2 = elem2_list[4].replace("'","")

rNum1 = int(elem1_list[5])
rNum2 = int(elem2_list[5])

# --------------------------------------------------------
# Compute eigenvalues and Hubbard values of the free atoms
# --------------------------------------------------------

for el in elements:
    valence = valences[el]
    atom = AtomicDFT(el,
                     xc=xc,
                     configuration=configurations[el],
                     valence=valence,
                     scalarrel=True,
                     confinement=PowerConfinement(r0=40., s=4),
                     txt='-',
                     )
    nl = nls2 if el == element2 else nls1
    scheme = 'central' if el == element2 else 'backward'
    u = atom.get_hubbard_value(nl, scheme=scheme)
    hubbardvalues[el] = {'s': u}
    atom.run()
    eigenvalues[el] = {nl: atom.get_eigenvalue(nl) for nl in valence}

# --------------------------------------------------
# Get KS all-electron ground state of confined atoms
# --------------------------------------------------

atoms = {}
for el in elements:
    r_cov = covalent_radii[atomic_numbers[el]] / Bohr
    r_wfc = 2 * r_cov
    r_rho = 3 * r_cov
    atom = AtomicDFT(el,
                     xc=xc,
                     configuration=configurations[el],
                     valence=valence,
                     scalarrel=True,
                     confinement=PowerConfinement(r0=r_rho, s=2),
                     wf_confinement=PowerConfinement(r0=r_wfc, s=2),
                     txt='-',
                     )
    atom.run()
    atoms[el] = atom

# -------------------------------
# Compute Slater-Koster integrals
# -------------------------------

for i in range(len(elements)):
    el_a = elements[i]
    rNum = rNum2/2.0 if el == element2 else rNum1/2.0
    for j in range(i + 1):
        rNum = rNum + rNum2/2.0 if el == element2 else rNum + rNum1/2.0
        el_b = elements[j]
        rmin, dr, N = 0.4, 0.02, int(rNum)
        sk = SlaterKosterTable(atoms[el_a], atoms[el_b], txt='-')
        sk.run(rmin, dr, N, superposition='density', xc=xc, stride=4)
        if el_a == el_b:
            sk.write('%s-%s_no_repulsion.skf' % (el_a, el_a), spe=0.,
                     eigenvalues=eigenvalues[el_a],
                     hubbardvalues=hubbardvalues[el_a],
                     occupations=occupations[el_a])
        else:
            for pair in [(el_a, el_b), (el_b, el_a)]:
                sk.write('%s-%s_no_repulsion.skf' % pair, pair=pair)
