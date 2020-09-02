from hotcent.atomic_dft import AtomicDFT
 
atom = AtomicDFT('Sn',
                 xc='LDA',
                 configuration='[Kr] 4d10 5s2 5p2',
                 valence=['5s', '5p', '4d'],
                 scalarrel=True,
                 nodegpts=150,
                 mix=0.2,
                 txt='-',
                 timing=True,
                 )

atom.run()
atom.plot_density()
atom.plot_Rnl()
atom.fit_sto('5s', 5, 4, filename='Sn_5s_STO.png')
atom.fit_sto('5p', 5, 4, filename='Sn_5p_STO.png')
atom.fit_sto('4d', 5, 4, filename='Sn_4d_STO.png')
atom.write_hsd(filename='Sn_wfc.hsd')
