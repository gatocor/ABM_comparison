#!/bin/bash

nvcc -std=c++11 random_walk.cu

./a.out

sudo chown gabriel output/
sudo chown gabriel output/*
sudo chmod +r output/*
sudo chmod +w output/*
