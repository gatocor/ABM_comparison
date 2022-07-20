#!/bin/bash

mkdir results/AgentBasedModels.jl

julia simulations_code/AgentBasedModels.jl/growth/growth.jl > results/AgentBasedModels.jl/time_growth_1.txt
