#!/bin/bash

if [ -e iter000 ]; then
  echo "delete previous file (iter000)"	
  rm -f -r iter000
fi

if [ -z "$1" ]; then
  echo "please, input name of element after command"
else
  cp generate_Xx_skf.py generate_$1_skf.py
  cp master_cp2k_origin.py master_cp2k.py
  cp master_qe_origin.py   master_qe.py
  cp ga_origin.py   ga.py

  sed -i "s/Xx/$1/g" generate_$1_skf.py
  sed -i "s/Xx/$1/g" master_cp2k.py
  sed -i "s/Xx/$1/g" master_qe.py

  python3 generate_$1_skf.py

  sed -i "s/Xx/$1/g" master_cp2k.py
  sed -i "s/Xx/$1/g" master_qe.py

  LMAM00=("H" "He" "Li" "Be")
  for i in "${LMAM00[@]}";
  do
    if [ $1 == $i ]; then
      LMAM=0
    fi
  done

  LMAM01=("B" "C" "N" "O" "F" "Ne" "Na" "Mg" "Al" "Si" "P" "S" "Cl" "Ar" "K" "Ca")
  for i in "${LMAM01[@]}";
  do
    if [ $1 == $i ]; then
      LMAM=1
      echo "Si",$i
    fi
  done
  
  LMAM02=("Sc" "Ti" "V" "Cr" "Mn" "Fe" "Co" "Ni" "Cu" "Zn" "Ga" "Ge" "As" "Se" "Br" "Kr" "Rb" "Sr" "Y" "Zr" "Nb" "Mo" "Tc" "Ru" "Rh" "Pd" "Ag" "Cd" "In" "Sn" "Sb" "Te" "I" "Xe" "Cs" "Ba" "La")
  for i in "${LMAM02[@]}";
  do
    if [ $1 == $i ]; then
      LMAM=2
    fi
  done

  LMAM03=("Hf" "Ta" "W" "Re" "Os" "Ir" "Pt" "Au" "Hg" "Tl" "Pb" "Bi")
  for i in "${LMAM03[@]}";
  do
    if [ $1 == $i ]; then
      LMAM=3
    fi
  done

  sed -i "s/MAM/$LMAM/g" master_cp2k.py
  sed -i "s/MAM/$LMAM/g" master_qe.py
  sed -i "s/MAM/$LMAM/g" ga.py

  #sed -i "s/Xx/$1/g" cp2k_calc.py

  sed -i "s/Xx/$1/g" ga.py
fi

if [ -z "$2" ]; then
  echo "END"
else
  echo ""
  echo "***********************************************"
  echo "Repulsion on tango"
  echo "***********************************************"
  ESPRESSO_PSEUDO=$(pwd)/../../pseudo
  export "ESPRESSO_PSEUDO=${ESPRESSO_PSEUDO}"
  ncpu=`grep physical.id /proc/cpuinfo | sort -u | wc -l`
  export "NCORES=${ncpu}"
  #export ASE_ESPRESSO_COMMAND="mpirun -np $NCORES pw.x -in PREFIX.pwi > PREFIX.pwo"
  #export ASE_CP2K_COMMAND="mpirun -np $NCORES cp2k_shell.popt"
  #export ASE_NWCHEM_COMMAND="mpirun -np $NCORES nwchem PREFIX.nwi > PREFIX.nwo"
  #export DFTB_PREFIX=$(pwd)
  python3 master_$2.py
fi
