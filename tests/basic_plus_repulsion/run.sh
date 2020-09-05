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

  sed -i "s/Xx/$1/g" generate_$1_skf.py
  sed -i "s/Xx/$1/g" master_cp2k.py
  sed -i "s/Xx/$1/g" master_qe.py

  python3 generate_$1_skf.py
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
