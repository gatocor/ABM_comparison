#!/bin/bash

echo "Downloading PhysiCell"
rm PhysiCell_V.1.10.4.zip
wget https://github.com/MathCancer/PhysiCell/releases/download/1.10.4/PhysiCell_V.1.10.4.zip
unzip PhysiCell_V.1.10.4.zip
rm PhysiCell_V.1.10.4.zip
mv PhysiCell packages 
cp -r simulations_code/PhysiCell/projects packages/PhysiCell
rm -r PhysiCell

echo "Downloading Yalla"
rm v1.0.zip
wget https://github.com/germannp/yalla/archive/refs/tags/v1.0.zip
unzip v1.0.zip
rm v1.0.zip
mv yalla-1.0/ yalla 
mv yalla packages 
cp -r simulations_code/yalla/projects packages/yalla
rm -r yalla
