#!/bin/bash

mkdir results/PhysiCell

rm -r packages/PhysiCell/projects/*
cp -r simulations_code/PhysiCell/projects packages/PhysiCell

cd packages/PhysiCell

make reset
cp -f projects/growth/config/* config
cp -f projects/growth/custom_modules/* custom_modules
cp -f projects/growth/main.cpp main.cpp
cp -f projects/growth/Makefile Makefile
make
./growth > ../../results/PhysiCell/time_growth_1.txt

make reset
cp -f projects/growth-death/config/* config
cp -f projects/growth-death/custom_modules/* custom_modules
cp -f projects/growth-death/main.cpp main.cpp
cp -f projects/growth-death/Makefile Makefile
make
./growth_death > ../../results/PhysiCell/time_growth_death_1.txt
./growth_death > ../../results/PhysiCell/time_growth_death_2.txt
./growth_death > ../../results/PhysiCell/time_growth_death_3.txt

cd ../..