from optparse import OptionParser
from ase.units import Bohr, Hartree
from ase.data import covalent_radii, atomic_numbers
from hotcent.slako import SlaterKosterTable
from hotcent.confinement import PowerConfinement
from hotcent.atomic_dft import AtomicDFT

# Run script with --help to see the options
p = OptionParser(usage='%prog')
p.add_option('-f', '--functional', default='LDA',
             help='Which density functional to apply? '
                  'E.g. LDA (default), GGA_X_PBE+GGA_C_PBE, ...')
p.add_option('-s', '--superposition',
             default='potential',
             help='Which superposition scheme? '
                  'Choose "potential" (default) or "density"')
p.add_option('-t', '--stride', default=1, type=int,
             help='Which SK-table stride length? Default = 1. '
                  'See hotcent.slako.run for more information.')
opt, args = p.parse_args()

# Get KS all-electron ground state of confined atom:
element = 'C'
r0 = 1.85 * covalent_radii[atomic_numbers[element]] / Bohr
atom = AtomicDFT(element,
                 xc=opt.functional,
                 confinement=PowerConfinement(r0=r0, s=2),
                 configuration='[He] 2s2 2p2',
                 valence=['2s', '2p'],
                 timing=True,
                 )
atom.run()
atom.plot_Rnl(only_valence=False)
atom.plot_density()

# Compute Slater-Koster integrals:
rmin, dr, N = 0.5, 0.05, 250
sk = SlaterKosterTable(atom, atom, timing=True)
sk.run(rmin, dr, N, superposition=opt.superposition,
       xc=opt.functional, stride=opt.stride)
sk.write('%s-%s_no_repulsion.par' % (element, element))
sk.write('%s-%s_no_repulsion.skf' % (element, element))
sk.plot()
