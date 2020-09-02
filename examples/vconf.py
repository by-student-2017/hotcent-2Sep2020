from hotcent.atomic_base import AtomicBase
from hotcent.confinement import ZeroConfinement, PowerConfinement
from ase.data import covalent_radii, atomic_numbers
from ase.units import Bohr

element = 'C'
r0 = 1.85 * covalent_radii[atomic_numbers[element]] / Bohr
vc = PowerConfinement(r0=r0, s=2)
v0 = ZeroConfinement()

for conf in [None, vc]:
    for wf_conf in [None, vc, {}, {'2s': vc}, {'2s': vc, '2p': None},
                    {'2s': vc, '2p': vc}, {'2s': None, '2p':None},
                    {'1s': v0, '2s': None, '2p': vc}]:
        atom = AtomicBase(element,
                          configuration='[He] 2s2 2p2',
                          valence=['2s', '2p'],
                          timing=False,
                          confinement=conf,
                          wf_confinement=wf_conf,
                          )
        print('-------------')
        print(conf, atom.confinement)
        print(wf_conf, atom.wf_confinement)
print('-------------')
