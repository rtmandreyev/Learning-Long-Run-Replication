#!/bin/bash 
#SBATCH --cpus-per-task 16 

matlab -nodesktop -nodisplay  <ucSimulation2.m> ucSimulation2.out
