#!/bin/bash

if [ -e iter000 ]; then
  echo "read iter000 file"
else
  echo "import *.skf file to iter000"
  mkdir iter000
  cp $1-$1_no_repulsion.skf ./iter000/$1-$1_no_repulsion.skf
fi

if [ -z "$1" ]; then
  echo "please, input name of element after command"
else
  cp master_cp2k_origin.py master_cp2k.py
  cp master_qe_origin.py   master_qe.py

  sed -i "s/Xx/$1/g" master_cp2k.py
  sed -i "s/Xx/$1/g" master_qe.py

  sed -i "s/Xx/$1/g" cp2k_calc.py

  sed -i "s/Xx/$1/g" ga.py
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
