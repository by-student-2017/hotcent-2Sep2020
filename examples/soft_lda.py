from ase.units import Ha
from hotcent.confinement import SoftConfinement
from hotcent.atomic_dft import AtomicDFT

kwargs = {'xc': 'LDA',
          'configuration': '[Ne] 3s2 3p2',
          'valence': ['3s', '3p'],
          'confinement': None,
          'scalarrel': True,
          'timing': True,
          'txt': '-'}

atom = AtomicDFT('Si',
                 wf_confinement=None,
                 **kwargs)
atom.run()
eps_free = {nl: atom.get_eigenvalue(nl) for nl in atom.valence}

wf_confinement = {'3s': SoftConfinement(amp=12., rc=6.74, x_ri=0.6),
                  '3p': SoftConfinement(amp=12., rc=8.70, x_ri=0.6)}
atom = AtomicDFT('Si',
                 wf_confinement=wf_confinement,
                 perturbative_confinement=True,
                 **kwargs)
atom.run()
eps_conf = {nl: atom.get_eigenvalue(nl) for nl in atom.valence}

print('\nChecking eigenvalues and their shifts upon confinement:')
# gpaw-setup Si -f LDA -a
eps_ref = {'3s': -0.399754, '3p': -0.152954}
# gpaw-basis Si -f LDA -t sz
shift_ref = {'3s': 0.104 / Ha, '3p': 0.103 / Ha}

for nl in atom.valence:
    ok_abs = abs(eps_free[nl] - eps_ref[nl]) < 1e-4
    items = (nl, eps_ref[nl], eps_free[nl], 'OK' if ok_abs else 'FAIL')
    print('  %s  %.6f  %.6f\t| %s' % items)

    shift = eps_conf[nl] - eps_free[nl]
    ok_shift = abs(shift - shift_ref[nl]) < 1e-4
    items = (nl, shift_ref[nl], shift, 'OK' if ok_shift else 'FAIL')
    print('  %s  %.6f  %.6f\t| %s' % items)
