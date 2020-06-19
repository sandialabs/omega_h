#!/bin/bash -e
if [ -f $1.cpp ]
then
  echo "$1.cpp"
  git mv $1.cpp Omega_h_$1.cpp
  sed -E "s/ $1.cpp/ Omega_h_$1.cpp/g" -i CMakeLists.txt
fi
if [ -f $1.hpp ]
then
  echo "$1.hpp"
  git mv $1.hpp Omega_h_$1.hpp
  for f in *pp
  do
    sed -E "s/\"$1.hpp/\"Omega_h_$1.hpp/g" -i $f
    sed -E "s/<$1.hpp/<Omega_h_$1.hpp/g" -i $f
  done
fi
