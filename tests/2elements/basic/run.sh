#!/bin/bash

if [ -e iter000 ]; then
  echo "delete previous file (iter000)"	
  rm -f -r iter000
fi

if [ -z "$1" ]; then
  echo "please, input name of element 1 after command"
elif [ -z "$2" ]; then
  echo "please, input name of element 2 after command"
else
  cp generate_Xx-Yy_skf.py generate_$1-$2_skf.py

  sed -i "s/Xx/$1/g" generate_$1-$2_skf.py
  sed -i "s/Yy/$2/g" generate_$1-$2_skf.py

  python3 generate_$1-$2_skf.py
fi
