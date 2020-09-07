#!/bin/bash

echo "Install libraries"
sudo apt update
sudo apt install -y gfortran
sudo apt install -y gcc
sudo apt install -y g++
sudo apt install -y build-essential
sudo apt install -y libopenmpi-dev
sudo apt install -y libscalapack-openmpi-dev
#sudo apt install -y libscalapack-openmpi2.0
sudo apt install -y libblas-dev 
sudo apt install -y liblapack-dev
sudo apt install -y libopenblas-dev
sudo apt install -y libarpack2-dev
sudo apt install -y libfftw3-dev
#sudo apt install -y libxc-dev
#sudo apt install -y git
sudo apt install -y wget
sudo apt install -y make
sudo apt install -y cmake
sudo apt install -y ase
sudo apt install -y python3-ase
sudo apt install -y python3-dev
sudo apt install -y python3-distutils
sudo apt install -y python3-setuptools
sudo apt install -y python3-numpy
sudo apt install -y python3-scipy
sudo apt install -y python3-f2py
sudo apt install -y python3-mpmath
sudo apt install -y python3-matplotlib
sudo apt install -y python3-sympy
#sudo apt install -y grace
#sudo apt install -y jmol
#sudo apt install -y gnuplot
sudo apt install -y gpaw cp2k quantum-espresso

if [ -d $HOME/tango/libxc-4.3.4 ]; then
  echo " "
  echo "skip libxc-4.3.4 installation"
else
  echo " "
  echo "libxc-4.3.4 install"
  tar zxvf libxc-4.3.4.tar.gz
  cd libxc-4.3.4
  sudo python3 setup.py install
fi

cd ~/hotcent

if [ -d $HOME/tango/dftbplus-19.1.x86_64-linux ]; then
  echo " "
  echo "skip DFTB+19.1 installation"
else
  echo " "
  echo "DFTB+ 17.1 install"
  tar xf dftbplus-17.1.x86_64-linux.tar.xz
  tar xf pbc-0-3.tar.xz
  echo "dftb+17.1 environment settings"
  echo ' ' >> ~/.bashrc
  echo '# dftb+17.1' >> ~/.bashrc
  echo 'export DFTB_COMMAND=$HOME/hotcent/dftbplus-17.1.x86_64-linux/bin/dftb+' >> ~/.bashrc
  echo 'export DFTB_PREFIX=$HOME/hotcent/pbc-0-3/' >> ~/.bashrc
fi

if [ -d $HOME/tango/ase-3.19.3 ]; then
  echo " "
  echo "skip ASE v.3.19.3 installation"
else
  echo " "
  echo "ASE v.3.19.3 install"
  tar zxvf ase-3.19.3.tar.gz
  echo "ASE v.3.19.3 environment settings"
  echo ' ' >> ~/.bashrc
  echo '# ASE v.3.19.3' >> ~/.bashrc
  echo 'export PYTHONPATH=$HOME/hotcent/ase-3.19.3:$PYTHONPATH' >> ~/.bashrc
  echo 'export PATH=$HOME/hotcent/ase-3.19.3/bin:$PATH' >> ~/.bashrc
fi

echo " "
echo "hotcent install"
export PYTHONPATH=~/hotcent:$PYTHONPATH
python3 setup.py build_ext --inplace
echo "hotcent environment settings"
echo ' ' >> ~/.bashrc
echo '# hotcent' >> ~/.bashrc
echo 'export PYTHONPATH=$HOME/hotcent:$PYTHONPATH' >> ~/.bashrc

echo "Installation, END"
