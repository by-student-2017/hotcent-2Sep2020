from ase.units import Ha
from hotcent.confinement import SoftConfinement
from hotcent.atomic_dft import AtomicDFT
try:
    import pylibxc
except ImportError:
    print('Need PyLibXC to run this test!')
    raise

kwargs = {'xc': 'GGA_X_PBE+GGA_C_PBE',
          'configuration': '[Ar] 3d6 4s2 4p0',
          'valence': ['3d', '4s'],
          'confinement': None,
          'scalarrel': True,
          'timing': True,
          'txt': '-'}

atom = AtomicDFT('Fe',
                 wf_confinement=None,
                 **kwargs)
atom.run()
eps_free = {nl: atom.get_eigenvalue(nl) for nl in atom.valence}

wf_confinement = {'3d': SoftConfinement(amp=12., rc=5.11, x_ri=0.6),
                  '4s': SoftConfinement(amp=12., rc=8.85, x_ri=0.6)}
atom = AtomicDFT('Fe',
                 wf_confinement=wf_confinement,
                 perturbative_confinement=True,
                 **kwargs)
atom.run()
eps_conf = {nl: atom.get_eigenvalue(nl) for nl in atom.valence}

print('\nChecking eigenvalues and their shifts upon confinement:')
# gpaw-setup Fe -f PBE -a
eps_ref = {'4s': -0.194442, '3d': -0.275800}
# gpaw-basis Fe -f PBE -t sz
shift_ref = {'4s': 0.098 / Ha, '3d': 0.100 / Ha}

for nl in atom.valence:
    ok_abs = abs(eps_free[nl] - eps_ref[nl]) < 5e-4
    items = (nl, eps_ref[nl], eps_free[nl], 'OK' if ok_abs else 'FAIL')
    print('  %s  %.6f  %.6f\t| %s' % items)

    shift = eps_conf[nl] - eps_free[nl]
    ok_shift = abs(shift - shift_ref[nl]) < 1e-4
    items = (nl, shift_ref[nl], shift, 'OK' if ok_shift else 'FAIL')
    print('  %s  %.6f  %.6f\t| %s' % items)
