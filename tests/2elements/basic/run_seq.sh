#!/bin/bash

elements1="H He"

elements2="H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi"

for elm1 in ${elements1} ;do
  for elm2 in ${elements2} ;do
    ./run.sh ${elm1} ${elm2}
    cp ./*_no_repulsion.skf ./results/
  done
done

ls -ltr ./results
