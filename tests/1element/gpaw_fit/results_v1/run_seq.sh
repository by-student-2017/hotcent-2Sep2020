#!/bin/bash

elements="Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi"

for elm in ${elements} ;do
  cp -r gpaw_fit gpaw_fit_${elm}
  cd gpaw_fit_${elm}
  ./run.sh ${elm} | tee output_${elm}.txt
  cd ..
  cp -r gpaw_fit_${elm} ./results/gpaw_fit_${elm}
done

ls -ltr ./results
