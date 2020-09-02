""" Defintion of the SlaterKosterTable class (and
supporting methods) for calculating Slater-Koster
integrals.

The code below draws heavily from the Hotbit code
written by Pekka Koskinen (https://github.com/pekkosk/
hotbit/blob/master/hotbit/parametrization/slako.py).
"""
import os
import sys
import numpy as np
from scipy.interpolate import SmoothBivariateSpline
from ase.units import Bohr
from ase.data import atomic_numbers, atomic_masses, covalent_radii
from hotcent.timing import Timer
from hotcent.xc import XC_PW92, LibXC
from hotcent.interpolation import CubicSplineFunction
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None
try:
    import _hotcent
except ModuleNotFoundError:
    print('Warning: C-extensions not available')
    _hotcent = None


class SlaterKosterTable:
    def __init__(self, ela, elb, txt='-', timing=False):
        """ Construct Slater-Koster table for given elements.

        Parameters:
        -----------
        ela:    AtomicDFT object
        elb:    AtomicDFT object
        txt:    where output should be printed
                use '-' for stdout (default), None for /dev/null,
                any other string for a text file, or a file handle
        timing: output of timing summary after calculation
        """
        self.ela = ela
        self.elb = elb

        if ela.get_symbol() != elb.get_symbol():
            self.nel = 2
            self.pairs = [(ela, elb), (elb, ela)]
            self.elements = [ela, elb]
        else:
            self.nel = 1
            self.pairs = [(ela, elb)]
            self.elements = [ela]

        if txt is None:
            self.txt = open(os.devnull, 'w')
        elif isinstance(txt, str):
            if txt == '-':
                self.txt = sys.stdout
            else:
                self.txt = open(txt, 'a')
        else:
            self.txt = txt

        self.timer = Timer('SlaterKosterTable', txt=self.txt, enabled=timing)

    def __del__(self):
        self.timer.summary()

    def write(self, filename=None, pair=None, eigenvalues={},
              hubbardvalues={}, occupations={}, spe=0.):
        """ Write SK tables to a file

        filename: str with name of file to write to.
                  The file format is selected by the extension
                  (.par or .skf).
                  Defaults to self.ela-self.elb_no_repulsion.skf

        pair: either (symbol_a, symbol_b) or (symbol_b, symbol_a)
              to select which of the two SK tables to write

        other kwargs: {nl: value}-dictionaries with eigenvalues,
              hubbardvalues and valence orbital occupations, as well
              as the spin-polarization error (all typically calculated
              on the basis of atomic DFT calculations). These will be
              written to the second line of a homo-nuclear .skf file.
              Examples: hubbardvalues={'2s': 0.5}, spe=0.2,
                        occupations={'3d':10, '4s': 1}, etc.
        """
        if pair is None:
            pair = (self.ela.get_symbol(), self.elb.get_symbol())

        fn = '%s-%s_no_repulsion.skf' % pair if filename is None else filename

        ext = fn[-4:]

        assert ext in ['.par', '.skf'], \
               "Unknown format: %s (-> choose .par or .skf)" % ext

        with open(fn, 'w') as handle:
            if ext == '.par':
                self._write_par(handle)
            elif ext == '.skf':
                self._write_skf(handle, pair, eigenvalues, hubbardvalues,
                                occupations, spe)

    def _write_skf(self, handle, pair, eigval, hubval, occup, spe):
        """ Write to SKF file format; this function
        is an adaptation of hotbit.io.hbskf 
        """
        symbols = (self.ela.get_symbol(), self.elb.get_symbol())
        if pair == symbols:
             index = 0
        elif pair == symbols[::-1]:
             index = 1
        else:
             msg = 'Requested ' + str(pair) + ' pair, but this calculator '
             msg += 'is restricted to the %s-%s pair.' % symbols
             raise ValueError(msg)

        grid_dist = self.Rgrid[1] - self.Rgrid[0]
        grid_npts = len(self.tables[index])
        grid_npts += int(self.Rgrid[0] / (self.Rgrid[1] - self.Rgrid[0]))
        print("%.12f, %d" % (grid_dist, grid_npts), file=handle)

        el1, el2 = self.ela.get_symbol(), self.elb.get_symbol()
        if el1 == el2:
            line = 'E_d   E_p   E_s   SPE   U_d   U_p   U_s   f_d   f_p   f_s  '
            labels = {x: 0 for x in line.split()}
            labels['SPE'] = spe
            for prefix, d in zip(['E', 'U', 'f'], [eigval, hubval, occup]):
                keys = list(d.keys())
                for l in ['s', 'p', 'd']:
                    check = [l in key for key in keys]
                    if any(check):
                        key = keys[check.index(True)]
                        labels['%s_%s' % (prefix, l)] = d[key]

            for label, value in labels.items():
                s = '%d' % value if isinstance(value, int) else '%.6f' % value
                line = line.replace(label + '  ', s)
            print(line, file=handle)

        m = atomic_masses[atomic_numbers[symbols[index]]]
        print("%.3f, 19*0.0" % m, file=handle)

        # Integral table containing the DFTB Hamiltonian

        if self.Rgrid[0] != 0:
            n = int(self.Rgrid[0] / (self.Rgrid[1] - self.Rgrid[0]))
            for i in range(n):
                print('%d*0.0,' % len(self.tables[0][0]), file=handle)

        ct, theader = 0, ''
        for i in range(len(self.tables[index])):
            line = ''

            for j in range(len(self.tables[index][i])):
                if self.tables[index][i, j] == 0:
                    ct += 1
                    theader = str(ct) + '*0.0 '
                else:
                    ct = 0
                    line += theader
                    theader = ''
                    line += '{0: 1.12e}  '.format(self.tables[index][i, j])

            if theader != '':
                ct = 0
                line += theader

            print(line, file=handle)

    def _write_par(self, handle):
        for p, (e1, e2) in enumerate(self.pairs):
            line = '%s-%s_table=' % (e1.get_symbol(), e2.get_symbol())
            print(line, file=handle)

            for i, R in enumerate(self.Rgrid):
                print('%.6e' % R, end=' ', file=handle)

                for t in range(20):
                    x = self.tables[p][i, t]
                    if abs(x) < 1e-90:
                        print('0.', end=' ', file=handle)
                    else:
                        print('%.6e' % x, end=' ', file=handle)
                print(file=handle)

            print('\n\n', file=handle)

    def plot(self, filename=None):
        """ Plot the Slater-Koster table with matplotlib.

        parameters:
        ===========
        filename:     name for the figure
        """
        self.timer.start('plotting')
        assert plt is not None, 'Matplotlib could not be imported!'

        fig = plt.figure()
        fig.subplots_adjust(hspace=1e-4, wspace=1e-4)

        el1 = self.ela.get_symbol()
        rmax = 6 * covalent_radii[atomic_numbers[el1]] / Bohr
        ymax = max(1, self.tables[0].max())
        if self.nel == 2:
            el2 = self.elb.get_symbol()
            rmax = max(rmax, 6 * covalent_radii[atomic_numbers[el2]] / Bohr)
            ymax = max(ymax, self.tables[1].max())

        for i in range(10):
            name = integrals[i]
            ax = plt.subplot(5, 2, i + 1)

            for p, (e1, e2) in enumerate(self.pairs):
                s1, s2 = e1.get_symbol(), e2.get_symbol()

                if p == 0: 
                    s = '-'
                    lw = 1
                    alpha = 1.0
                else: 
                    s = '--'
                    lw = 4
                    alpha = 0.2

                if np.all(abs(self.tables[p][:, i]) < 1e-10):
                    ax.text(0.03, 0.5 + p * 0.15,
                            'No %s integrals for <%s|%s>' % (name, s1, s2),
                            transform=ax.transAxes, size=10, va='center')

                    if not ax.is_last_row():
                        plt.xticks([], [])
                    if not ax.is_first_col():
                        plt.yticks([], [])
                else:
                    plt.plot(self.Rgrid, self.tables[p][:, i] , c='r',
                             ls=s, lw=lw, alpha=alpha)
                    plt.plot(self.Rgrid, self.tables[p][:, i + 10], c='b',
                             ls=s, lw=lw, alpha=alpha)
                    plt.axhline(0, c='k', ls='--')
                    ax.text(0.8, 0.1 + p * 0.15, name, size=10,
                            transform=ax.transAxes)

                    if ax.is_last_row():
                        plt.xlabel('r (Bohr)')
                    else:
                        plt.xticks([], [])
                    if not ax.is_first_col():
                        plt.yticks([],[])

                plt.xlim([0, rmax])
                plt.ylim(-ymax, ymax)

        plt.figtext(0.3, 0.95, 'H', color='r', size=20)
        plt.figtext(0.34, 0.95, 'S', color='b', size=20)
        plt.figtext(0.38, 0.95, ' Slater-Koster tables', size=20)
        e1, e2 = self.ela.get_symbol(), self.elb.get_symbol()
        plt.figtext(0.3, 0.92, '(thin solid: <%s|%s>, wide dashed: <%s|%s>)' \
                    % (e1, e2, e2, e1), size=10)
        
        if filename is None:
            filename = '%s-%s_slako.pdf' % (e1, e2)
        plt.savefig(filename, bbox_inches='tight')
        plt.clf()
        self.timer.stop('plotting')

    def get_range(self, fractional_limit):
        """ Define ranges for the atoms: largest r such that Rnl(r)<limit. """
        wf_range = 0.

        for el in self.elements:
            r = max([el.get_wf_range(nl, fractional_limit) 
                     for nl in el.get_valence_orbitals()])
            print('Wave function range for %s = %.5f a0' % (el.get_symbol(), r),
                  file=self.txt)
            wf_range = max(r, wf_range)

        assert wf_range < 20, 'Wave function range exceeds 20 Bohr radii. ' \
                              'Decrease wflimit?'
        return wf_range
        
    def run(self, rmin=0.4, dr=0.02, N=None, ntheta=150, nr=50, wflimit=1e-7,
            superposition='potential', xc='LDA', stride=1):
        """ Calculate the Slater-Koster table.

        parameters:
        ------------
        rmin, dr, N: parameters defining the equidistant grid of interatomic
                separations: the shortest distance rmin and grid spacing dr
                (both in Bohr radii) and the number of grid points N.
        ntheta: number of angular divisions in polar grid
                (more dense towards bonding region).
        nr:     number of radial divisions in polar grid
                (more dense towards origins).
                with p=q=2 (powers in polar grid) ntheta~3*nr is 
                optimal (with fixed grid size)
                with ntheta=150, nr=50 you get~1E-4 accuracy for H-elements
                (beyond that, gain is slow with increasing grid size)
        wflimit: value below which the radial wave functions are considered
                to be negligible. This determines how far the polar grids
                around the atomic centers extend in space.
        superposition: 'density' or 'potential': whether to use the density
                superposition or potential superposition approach for the
                Hamiltonian integrals.
        xc:     name of the exchange-correlation functional to be used
                in calculating the effective potential in the density
                superposition scheme. If the PyLibXC module is available,
                any LDA or GGA (but not hybrid or MGGA) functional available
                via LibXC can be specified. E.g. for using the N12
                functional, set xc='XC_GGA_X_N12+XC_GGA_C_N12'.
                If PyLibXC is not available, only the local density
                approximation xc='PW92' (alias: 'LDA') can be chosen.
        stride: the desired SK-table typically has quite a large number
                of points (N=500-1000), even though the integrals
                themselves are comparatively smooth. To speed up the
                construction of the SK-table, one can therefore restrict
                the expensive integrations to a subset N' = N // stride,
                and map the resulting curves on the N-grid afterwards.
                The default stride = 1 means that N' = N (no shortcut).
        """
        print('\n\n', file=self.txt)
        print('***********************************************', file=self.txt)
        print('Slater-Koster table construction for %s and %s' % \
              (self.ela.get_symbol(), self.elb.get_symbol()), file=self.txt)
        print('***********************************************', file=self.txt)
        self.txt.flush()

        assert N is not None, 'Need to set number of grid points N!'
        assert rmin >= 1e-3, 'For stability, please set rmin >= 1e-3'
        assert superposition in ['density', 'potential']

        self.timer.start('calculate_tables')
        self.wf_range = self.get_range(wflimit)
        Nsub = N // stride
        Rgrid = rmin + stride * dr * np.arange(Nsub)
        tables = [np.zeros((Nsub, 20)) for i in range(self.nel)]
        dH = 0.
        Hmax = 0.

        for p, (e1, e2) in enumerate(self.pairs):
            print('Integrals:', end=' ', file=self.txt)
            selected = select_integrals(e1, e2)
            for s in selected:
                print(s[0], end=' ', file=self.txt)
            print(file=self.txt, flush=True)

        for i, R in enumerate(Rgrid):
            if R > 2 * self.wf_range:
                break

            grid, area = self.make_grid(R, nt=ntheta, nr=nr)

            if  i == Nsub - 1 or Nsub // 10 == 0 or i % (Nsub // 10) == 0:
                print('R=%8.2f, %i grid points ...' % (R, len(grid)),
                      file=self.txt, flush=True)

            for p, (e1, e2) in enumerate(self.pairs):
                selected = select_integrals(e1, e2)
                S, H, H2 = 0., 0., 0.
                if len(grid) > 0:
                    S, H, H2 = self.calculate_mels(selected, e1, e2, R, grid,
                                                   area, xc=xc,
                                                   superposition=superposition)
                    Hmax = max(Hmax, max(abs(H)))
                    dH = max(dH, max(abs(H - H2)))
                tables[p][i, :10] = H
                tables[p][i, 10:] = S

        if superposition == 'potential':
            print('Maximum value for H=%.2g' % Hmax, file=self.txt)
            print('Maximum error for H=%.2g' % dH, file=self.txt)
            print('     Relative error=%.2g %%' % (dH / Hmax * 100),
                  file=self.txt)

        self.Rgrid = rmin + dr * np.arange(N)

        if stride > 1:
            self.tables = [np.zeros((N, 20)) for i in range(self.nel)]
            for p in range(self.nel):
                for i in range(20):
                    spl = CubicSplineFunction(Rgrid, tables[p][:, i])
                    self.tables[p][:, i] = spl(self.Rgrid)
        else:
            self.tables = tables

        # Smooth the curves near the cutoff
        for p in range(self.nel):
            for i in range(20):
                self.tables[p][:, i] = tail_smoothening(self.Rgrid,
                                                        self.tables[p][:, i])

        self.timer.stop('calculate_tables')

    def calculate_mels(self, selected, e1, e2, R, grid, area,
                       superposition='potential', xc='LDA'):
        """ Perform integration for selected H and S integrals.

        parameters:
        -----------
        selected: list of [('dds', '3d', '4d'), (...)]
        e1: <bra| element
        e2: |ket> element
        R: e1 is at origin, e2 at z=R
        grid: list of grid points on (d, z)-plane
        area: d-z areas of the grid points.
        superposition: 'density' or 'potential' superposition scheme
        xc: exchange-correlation functional (see description in self.run())

        return:
        -------
        List of H,S and H2 for selected integrals. In the potential
        superposition scheme, H2 is calculated using a different technique
        and can be used for error estimation. This is not available
        for the density superposition scheme, where simply H2=0 is returned.

        S: simply R1 * R2 * angle_part

        H: operate (derivate) R2 <R1 | t + Veff - Conf1 - Conf2 | R2>.
           With potential superposition: Veff = Veff1 + Veff2
           With density superposition: Veff = Vxc(n1 + n2)

        H2: operate with full h2 and hence use eigenvalue of | R2>
            with full Veff2:
              <R1 | (t1 + Veff1) + Veff2 - Conf1 - Conf2 | R2>
            = <R1 | h1 + Veff2 - Conf1 - Conf2 | R2> (operate with h1 on left)
            = <R1 | e1 + Veff2 - Conf1 - Conf2 | R2>
            = e1 * S + <R1 | Veff2 - Conf1 - Conf2 | R2>
            -> H and H2 can be compared and error estimated
        """
        self.timer.start('calculate_mels')
        Sl, Hl, H2l = np.zeros(10), np.zeros(10), np.zeros(10)

        # common for all integrals (not wf-dependent parts)
        self.timer.start('prelude')
        N = len(grid)
        x = grid[:N, 0]
        y = grid[:N, 1]
        r1 = np.sqrt(x ** 2 + y ** 2)
        r2 = np.sqrt(x ** 2 + (R - y) ** 2)
        t1 = np.arccos(y / r1)
        t2 = np.arccos((y - R) / r2)
        radii = np.array([r1, r2]).T
        gphi = g(t1, t2).T

        if superposition == 'potential':
            self.timer.start('vrho')
            v1 = e1.effective_potential(r1) - e1.confinement(r1)
            v2 = e2.effective_potential(r2) - e2.confinement(r2)
            veff = v1 + v2
            self.timer.stop('vrho')
        elif superposition == 'density':
            self.timer.start('vrho')
            rho = e1.electron_density(r1) + e2.electron_density(r2)
            veff = e1.nuclear_potential(r1) + e1.hartree_potential(r1)
            veff += e2.nuclear_potential(r2) + e2.hartree_potential(r2)
            if xc in ['LDA', 'PW92']:
                xc = XC_PW92()
                veff += xc.vxc(rho)
                self.timer.stop('vrho')
            else:
                xc = LibXC(xc)
                drho1 = e1.electron_density(r1, der=1)
                drho2 = e2.electron_density(r2, der=1)
                grad_x = drho1 * np.sin(t1)
                grad_x += drho2 * np.sin(t2)
                grad_y = drho1 * np.cos(t1)
                grad_y += drho2 * np.cos(t2)
                sigma = np.sqrt(grad_x ** 2 + grad_y ** 2) ** 2
                out = xc.compute_all(rho, sigma)
                veff += out['vrho']
                self.timer.stop('vrho')
                self.timer.start('vsigma')
                # add gradient corrections to vxc
                # provided that we have enough points
                # (otherwise we get "dfitpack.error:
                # (m>=(kx+1)*(ky+1)) failed for hidden m")
                if out['vsigma'] is not None and len(x) > 16:
                    splx = SmoothBivariateSpline(x, y, out['vsigma'] * grad_x)
                    sply = SmoothBivariateSpline(x, y, out['vsigma'] * grad_y)
                    veff += -2. * splx(x, y, dx=1, dy=0, grid=False)
                    veff += -2. * sply(x, y, dx=0, dy=1, grid=False)
                self.timer.stop('vsigma')

        assert np.shape(gphi) == (N, 10)
        assert np.shape(radii) == (N, 2)
        assert np.shape(veff) == (N,)
        self.timer.stop('prelude')

        # calculate all selected integrals
        for integral, nl1, nl2 in selected:
            index = integrals.index(integral)
            S, H, H2 = 0., 0., 0.
            l2 = angular_momentum[nl2[1]]

            nA = len(area)
            r1 = radii[:nA, 0]
            r2 = radii[:nA, 1]
            d, z = grid[:nA, 0], grid[:nA, 1]
            aux = gphi[:nA, index] * area * d
            Rnl1 = e1.Rnl(r1, nl1)
            Rnl2 = e2.Rnl(r2, nl2)
            ddunl2 = e2.unl(r2, nl2, der=2)

            S = np.sum(Rnl1 * Rnl2 * aux)
            H = np.sum(Rnl1 * (-0.5 * ddunl2 / r2 + (veff + \
                       l2 * (l2 + 1) / (2 * r2 ** 2)) * Rnl2) * aux)

            if superposition == 'potential':
                H2 = np.sum(Rnl1 * Rnl2 * aux * (v2 - e1.confinement(r1)))
                H2 += e1.get_epsilon(nl1) * S
            elif superposition == 'density':
                H2 = 0

            Sl[index] = S
            Hl[index] = H
            H2l[index] = H2

        self.timer.stop('calculate_mels')
        return Sl, Hl, H2l

    def make_grid(self, Rz, nt, nr, p=2, q=2, view=False):
        """ Construct a double-polar grid.

        Parameters:
        -----------
        Rz: element 1 is at origin, element 2 at z=Rz
        nt: number of theta grid points
        nr: number of radial grid points
        p: power describing the angular distribution of grid points
           (larger puts more weight towards theta=0)
        q: power describing the radial disribution of grid points
           (larger puts more weight towards centers)
        view: view the distribution of grid points with matplotlib

        Plane at R/2 divides two polar grids.


         ^ (z-axis)
         |--------_____               phi_j
         |       /     ----__         *
         |      /            \       /  *
         |     /               \    /  X *                X=coordinates of the center of area element(z,d),
         |    /                  \  \-----* phi_(j+1)     area=(r_(i+1)**2-r_i**2)*(phi_(j+1)-phi_j)/2
         |   /                    \  r_i   r_(i+1)
         |  /                      \
         | /                       |
         *2------------------------|           polar centered on atom 2
         | \                       |
         |  \                     /                                                     1
         |   \                   /                                                     /  \
         |-------------------------- z=h -line         ordering of sector slice       /     \
         |   /                   \                                      points:      /        \
         |  /                     \                                                 /          \
         | /                       |                                               /     0       4
         *1------------------------|--->      polar centered on atom 1            2            /
         | \                       |    (r_perpendicular (xy-plane) = 'd-axis')    \        /
         |  \                      /                                                 \   /
         |   \                    /                                                    3
         |    \                  /
         |     \               /
         |      \           /
         |       \ ___ ---
         |---------
        """
        self.timer.start('make_grid')
        rmin, rmax = 1e-7, self.wf_range
        h = Rz / 2

        if _hotcent is not None:
            grid, area = _hotcent.make_grid(Rz, rmin, rmax, nt, nr, p, q)
        else:
            max_range = self.wf_range
            T = np.linspace(0, 1, nt) ** p * np.pi
            R = rmin + np.linspace(0, 1, nr) ** q * (rmax - rmin)

            area = np.array([])
            d = np.array([])
            z = np.array([])

            # first calculate grid for polar centered on atom 1:
            # the z=h-like starts cutting full elements starting from point (1)
            Tj0 = T[:nt - 1]
            Tj1 = T[1: nt]

            for i in range(nr - 1):
                # corners of area element
                d1 = R[i + 1] * np.sin(Tj0)
                z1 = R[i + 1] * np.cos(Tj0)
                d2 = R[i] * np.sin(Tj0)
                z2 = R[i] * np.cos(Tj0)
                d3 = R[i] * np.sin(Tj1)
                z3 = R[i] * np.cos(Tj1)
                d4 = R[i + 1] * np.sin(Tj1)
                z4 = R[i + 1] * np.cos(Tj1)

                cond_list = [z1 <= h,  # area fully inside region
                     (z1 > h) * (z2 <= h) * (z4 <= h),  # corner 1 outside region
                     (z1 > h) * (z2 > h) * (z3 <= h) * (z4 <= h),  # 1 & 2 outside
                     (z1 > h) * (z2 > h) * (z3 <= h) * (z4 > h),  # only 3 inside
                     (z1 > h) * (z2 <= h) * (z3 <= h) * (z4 > h),  # 1 & 4 outside
                     (z1 > h) * (z3 > h) * ~((z2 <= h) * (z4 > h))]

                r0_list = [0.5 * (R[i] + R[i + 1]),
                           0.5 * (R[i] + R[i + 1]),
                           0.5 * (R[i] + R[i + 1]),
                           lambda x: 0.5 * (R[i] + h / np.cos(x)),
                           lambda x: 0.5 * (R[i] + h / np.cos(x)),
                           0,
                           np.nan]
                r0 = np.piecewise(Tj1, cond_list, r0_list)

                Th0 = np.piecewise(h / R[i], [np.abs(h / R[i]) > 1],
                                   [np.nan, lambda x: np.arccos(x)])
                Th1 = np.piecewise(h / R[i + 1], [np.abs(h / R[i + 1]) > 1],
                                   [np.nan, lambda x: np.arccos(x)])

                t0_list = [lambda x: 0.5 * x,
                           0.5 * Th1,
                           0.5 * Th1,
                           0.5 * Th0,
                           lambda x: 0.5 * x,
                           0,
                           np.nan]
                t0 = 0.5 * Tj1
                t0 += np.piecewise(Tj0, cond_list, t0_list)

                rr = 0.5 * (R[i + 1] ** 2 - R[i] ** 2)
                A_list0 = [lambda x: rr * -x,
                           lambda x: rr * -x - 0.5 * R[i + 1] ** 2 * (Th1 - x) \
                                     + 0.5 * h ** 2 * (np.tan(Th1) - np.tan(x)),
                           lambda x: rr * -x - (rr * -x + 0.5 * R[i + 1] ** 2 \
                                     * (Th1 - Th0)),
                           0.,
                           lambda x: 0.5 * h ** 2 * -np.tan(x) \
                                     - 0.5 * R[i] ** 2 * -x,
                           -1,
                           np.nan]
                A = np.piecewise(Tj0, cond_list, A_list0)

                A_list1 = [lambda x: rr * x,
                           lambda x: rr * x,
                           lambda x: rr * x - (rr * Th0 - 0.5 * h ** 2 \
                                     * (np.tan(Th1) - np.tan(Th0))),
                           lambda x: 0.5 * h ** 2 * (np.tan(x) - np.tan(Th0)) \
                                     - 0.5 * R[i] ** 2 * (x - Th0),
                           lambda x: 0.5 * h ** 2 * np.tan(x) \
                                     - 0.5 * R[i] ** 2 * x,
                           0,
                           np.nan]
                A += np.piecewise(Tj1, cond_list, A_list1)

                dd = r0 * np.sin(t0)
                zz = r0 * np.cos(t0)
                select = np.sqrt(dd ** 2 + zz ** 2) < max_range
                select *= np.sqrt(dd ** 2 + (Rz - zz) ** 2) < max_range
                select *= A > 0
                area = np.hstack((area, A[select]))
                d = np.hstack((d, dd[select]))
                z = np.hstack((z, zz[select]))
            grid = np.array([d, z]).T

        self.timer.start('symmetrize')
        # calculate the polar centered on atom 2 by mirroring the other grid
        grid2 = grid.copy()
        grid2[:, 1] = -grid[:, 1]
        shift = np.zeros_like(grid)
        shift[:, 1] = 2 * h
        grid = np.concatenate((grid, grid2 + shift))
        area = np.concatenate((area, area))
        self.timer.stop('symmetrize')

        if view:
            assert plt is not None, 'Matplotlib could not be imported!'
            plt.plot([h, h ,h])
            plt.scatter(grid[:, 0], grid[:, 1], s=10 * area / np.max(area))
            plt.show()
            plt.clf()

        self.timer.stop('make_grid')
        return grid, area


integrals = ['dds', 'ddp', 'ddd', 'pds', 'pdp', 'pps', 'ppp', 
             'sds', 'sps', 'sss']
angular_momentum = {'s': 0, 'p': 1, 'd': 2}


def select_orbitals(val1, val2, integral):
    """ Select orbitals from given valences to calculate given integral.
    e.g. ['2s', '2p'], ['4s', '3d'], 'sds' --> ('2s', '3d')
    """
    nl1 = None
    for nl in val1:
        if nl[1] == integral[0]:
            nl1 = nl

    nl2 = None
    for nl in val2:
        if nl[1] == integral[1]: 
            nl2 = nl

    return nl1, nl2


def select_integrals(e1, e2):
    """ Return list of integrals (integral,nl1,nl2) 
    to be done for element pair e1, e2. """
    selected = []
    val1, val2 = e1.get_valence_orbitals(), e2.get_valence_orbitals()

    for integral in integrals:
        nl1, nl2 = select_orbitals(val1 , val2 , integral)
        if nl1 is None or nl2 is None:
            continue
        else:
            selected.append((integral, nl1, nl2))

    return selected


def g(t1, t2):
    """ Return the angle-dependent part of the two-center
    integral (it) with t1=theta_1 (atom at origin)
    and t2=theta2 (atom at z=Rz). These dependencies
    come after integrating analytically over phi.
    """
    c1, c2, s1, s2 = np.cos(t1), np.cos(t2), np.sin(t1), np.sin(t2)
    return np.array([5. / 8 * (3 * c1 ** 2 - 1) * (3 * c2 ** 2 - 1),\
                     15. / 4 * s1 * c1 * s2 * c2,\
                     15. / 16 * s1 ** 2 * s2 ** 2,\
                     np.sqrt(15.) / 4 * c1 * (3 * c2 ** 2 - 1),\
                     np.sqrt(45.) / 4 * s1 * s2 * c2,\
                     3. / 2 * c1 * c2,\
                     3. / 4 * s1 * s2,\
                     np.sqrt(5.) / 4 * (3 * c2 ** 2 - 1),\
                     np.sqrt(3.) / 2 * c2,\
                     0.5*np.ones_like(t1)])


def tail_smoothening(x, y):
    """ For given grid-function y(x), make smooth tail.

    Aim is to get (e.g. for Slater-Koster tables and repulsions) smoothly
    behaving energies and forces near cutoff region.
    
    Make is such that y and y' go smoothly exactly to zero at last point.
    Method: take largest neighboring points y_k and y_(k+1) (k<N-3) such
    that line through them passes zero below x_(N-1). Then fit
    third-order polynomial through points y_k, y_k+1 and y_N-1.

    Return:
    smoothed y-function on same grid.
    """
    if np.all(abs(y) < 1e-10):
        return y

    Nzero = 0
    for i in range(len(y) - 1, 1, -1):
        if abs(y[i]) < 1e-60:
            Nzero += 1
        else:
            break

    N = len(y) - Nzero
    y = y[:N]
    xmax = x[:N][-1]

    for i in range(N - 3, 1, -1):
        x0i = x[i] - y[i] / ((y[i + 1] - y[i]) /(x[i + 1] - x[i]))
        if x0i < xmax:
            k = i
            break
    else:
        print('N:', N, 'len(y):', len(y))
        for i in range(len(y)):
            print(x[i], y[i])
        raise RuntimeError('Problem with tail smoothening')

    if k < N / 4:
        for i in range(N):
            print(x[i], y[i])
        msg = 'Problem with tail smoothening: requires too large tail.'
        raise RuntimeError(msg)

    if k == N - 3:
        y[-1] = 0.
        y = np.append(y, np.zeros(Nzero))
        return y
    else:
        # g(x)=c2*(xmax-x)**m + c3*(xmax-x)**(m+1) goes through 
        # (xk,yk),(xk+1,yk+1) and (xmax,0)
        # Try different m if g(x) should change sign (this we do not want)
        sgn = np.sign(y[k])
        for m in range(2, 10):
            a1, a2 = (xmax - x[k]) ** m, (xmax - x[k]) ** (m + 1)
            b1, b2 = (xmax-  x[k + 1]) ** m, (xmax - x[k + 1]) ** (m + 1)
            c3 = (y[k] - a1 * y[k + 1] / b1) / (a2 - a1 * b2 / b1)
            c2 = (y[k] - a2 * c3) / a1

            for i in range(k + 2,N):
                y[i] = c2 * (xmax - x[i]) ** 2 + c3 * (xmax - x[i]) ** 3

            y[-1] = 0.  # once more explicitly

            if np.all(y[k:] * sgn >= 0):
                y = np.append(y, np.zeros(Nzero))
                break

            if m == 9:
                msg = 'Problems with smoothening; need for new algorithm?'
                raise RuntimeError(msg)

    return y
