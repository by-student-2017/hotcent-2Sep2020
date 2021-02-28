#!/bin/bash

if [ -z "$1" ]; then
  echo "please, input name of element after command"
else
  ESPRESSO_PSEUDO=$(pwd)/../pseudo
  export "ESPRESSO_PSEUDO=${ESPRESSO_PSEUDO}"
  ncpu=`grep physical.id /proc/cpuinfo | sort -u | wc -l`
  export "NCORES=${ncpu}"
  export ASE_ESPRESSO_COMMAND="mpirun -np $NCORES pw.x -in PREFIX.pwi > PREFIX.pwo"

  export DFTB_PREFIX=$(pwd)

  cp qe_data_origin.py qe_data_$1.py
  cp param_fit_origin.py param_fit_$1.py

  sed -i "s/Xx/$1/g" qe_data_$1.py
  sed -i "s/Xx/$1/g" param_fit_$1.py

  python3 qe_data_$1.py
  python3 param_fit_$1.py

  #mkdir iter000
  #mv $1-$1.skf ./iter000/$1-$1_no_repulsion.skf
  mv $1-$1.skf $1-$1_no_repulsion.skf

  python3 compare.py
fi

rm -f -r pwscf.save

#if [ -z "$2" ]; then
#  echo "END"
#else
#  echo ""
#  echo "***********************************************"
#  echo "Repulsion on tango"
#  echo "***********************************************"
#  ESPRESSO_PSEUDO=$(pwd)/pseudo
#  export "ESPRESSO_PSEUDO=${ESPRESSO_PSEUDO}"
#  ncpu=`grep physical.id /proc/cpuinfo | sort -u | wc -l`
#  export "NCORES=${ncpu}"
#  python3 master_$2.py
#fi
