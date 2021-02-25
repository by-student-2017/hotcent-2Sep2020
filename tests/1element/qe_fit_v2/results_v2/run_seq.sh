#!/bin/bash

mkdir results

# GBRV potentials do not include He, Ne, Ar, Kr and Xe.
elements="H Li Be B C N O F Na Mg Al Si P S Cl K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Cs Ba Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi"

#elements="H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi"

for elm in ${elements} ;do
  cp -r qe_fit_v2 qe_fit_v2_${elm}
  cd qe_fit_v2_${elm}
  ./run.sh ${elm} | tee output_${elm}.txt
  #rm -f -r primCIFs CIFs
  cd ..
  mv qe_fit_v2_${elm} ./results/qe_fit_v2_${elm}
done

ls -ltr ./results
