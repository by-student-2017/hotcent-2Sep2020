#!/bin/bash

if [ -e iter000 ]; then
  echo "delete previous file (iter000)"	
  rm -f -r iter000
fi

if [ -z "$1" ]; then
  echo "please, input name of element after command"
else
  cp generate_Xx_skf.py generate_$1_skf.py

  sed -i "s/Xx/$1/g" generate_$1_skf.py

  python3 generate_$1_skf.py
fi
