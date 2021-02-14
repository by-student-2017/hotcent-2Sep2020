#!/bin/bash

mkdir results

elements="H Be C N O F Na Mg P S Cl Sc Ti Co Zn Ga Se Br Y Zr Ru Cd In Te I Hf Re Os Tl "

#elements="H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi"

for elm in ${elements} ;do
  cp -r gpaw_fit_v3 gpaw_fit_v3_${elm}
  cd gpaw_fit_v3_${elm}
  ./run.sh ${elm} | tee output_${elm}.txt
  rm -f -r primCIFs CIFs
  cd ..
  mv gpaw_fit_v3_${elm} ./results/gpaw_fit_v3_${elm}
done

ls -ltr ./results
