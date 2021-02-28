#!/bin/bash

mkdir results

# GBRV potentials do not include He, Ne, Ar, Kr and Xe.
elements="H Be C N Mg Sc Ti Co Se Y Zr Cd In Te Hf W Tl"

#elements="H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi"

for elm in ${elements} ;do
  cp -r qe_fit_v3 qe_fit_v3_${elm}
  cd qe_fit_v3_${elm}
  ./run.sh ${elm} | tee output_${elm}.txt
  rm -f -r primCIFs CIFs
  cd ..
  mv qe_fit_v3_${elm} ./results/qe_fit_v3_${elm}
done

ls -ltr ./results
