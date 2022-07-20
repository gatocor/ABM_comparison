#!/bin/bash

mkdir results/Agents.jl

julia simulations_code/Agents.jl/growth/growth.jl > results/Agents.jl/time_growth_1.txt
