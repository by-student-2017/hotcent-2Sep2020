#!/bin/bash

if [ -e iter000 ]; then
  echo "read iter000 file"
else
  echo "import *.skf file to iter000"
  cp generate_Xx_skf.py generate_$1_skf.py
  sed -i "s/Xx/$1/g" generate_$1_skf.py
  python3 generate_$1_skf.py
  mkdir iter000
  cp $1-$1_no_repulsion.skf ./iter000/$1-$1_no_repulsion.skf
fi

if [ -z "$1" ]; then
  echo "please, input name of element after command"
else
  cp master_cp2k_origin.py master_cp2k.py
  cp master_qe_origin.py   master_qe.py
  cp ga_origin.py   ga.py

  sed -i "s/Xx/$1/g" master_cp2k.py
  sed -i "s/Xx/$1/g" master_qe.py
  sed -i "s/Xx/$1/g" ga.py

  LMAM0=("H" "He" "Li" "Be")
  for i in "${LMAM0[@]}";
  do
    if [ $1 == $i ]; then
      LMAM=0
    fi
  done

  LMAM1=("B" "C" "N" "O" "F" "Ne" "Na" "Mg" "Al" "Si" "P" "S" "Cl" "Ar" "K" "Ca")
  for i in "${LMAM1[@]}";
  do
    if [ $1 == $i ]; then
      LMAM=1
    fi
  done
  
  LMAM2=("Sc" "Ti" "V" "Cr" "Mn" "Fe" "Co" "Ni" "Cu" "Zn" "Ga" "Ge" "As" "Se" "Br" "Kr" "Rb" "Sr" "Y" "Zr" "Nb" "Mo" "Tc" "Ru" "Rh" "Pd" "Ag" "Cd" "In" "Sn" "Sb" "Te" "I" "Xe" "Cs" "Ba" "La")
  for i in "${LMAM2[@]}";
  do
    if [ $1 == $i ]; then
      LMAM=2
    fi
  done

  LMAM3=("Hf" "Ta" "W" "Re" "Os" "Ir" "Pt" "Au" "Hg" "Tl" "Pb" "Bi")
  for i in "${LMAM3[@]}";
  do
    if [ $1 == $i ]; then
      LMAM=2
    fi
  done

  sed -i "s/MAM/$LMAM/g" master_cp2k.py
  sed -i "s/MAM/$LMAM/g" master_qe.py  
  sed -i "s/MAM/$LMAM/g" ga.py

  sed -i "s/Xx/$1/g" cp2k_calc.py

fi

if [ -z "$2" ]; then
  echo "END"
else
  echo ""
  echo "***********************************************"
  echo "Repulsion on tango"
  echo "***********************************************"
  ESPRESSO_PSEUDO=$(pwd)/pseudo
  export "ESPRESSO_PSEUDO=${ESPRESSO_PSEUDO}"
  ncpu=`grep physical.id /proc/cpuinfo | sort -u | wc -l`
  export "NCORES=${ncpu}"
  python3 master_$2.py
fi
