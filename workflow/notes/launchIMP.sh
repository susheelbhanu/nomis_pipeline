#!/bin/bash -l

sbatch --time=120:00:0 -N1 -n8 -p bigmem --qos qos-bigmem  -J "$sample" ./"$sample".runIMP.sh
