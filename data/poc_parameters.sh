#!/bin/bash
# Minimal POC parameters - designed for ~6 hour completion

# Valid system sizes (L² must be divisible by 16)
declare -A SYSTEM_PARAMS
SYSTEM_PARAMS[8]="2048"    # L=8,  S=64  ✓
SYSTEM_PARAMS[12]="2048"   # L=12, S=144 ✓  
SYSTEM_PARAMS[16]="2048"   # L=16, S=256 ✓
SYSTEM_PARAMS[20]="2048"   # L=20, S=400  ✓

# Minimal but sufficient parameters
K_MIN=0.28              # Start of critical region
K_MAX=0.30              # End of critical region
NK=24                   # 24 temperatures (reduced from 24)

MC_PER_BLOCK[8]=66      # For L=8, each block = 66 Metropolis+PT passes
MC_PER_BLOCK[12]=84     # For L=12, each block = 84 passes
MC_PER_BLOCK[16]=102    # For L=16, each block = 102 passes
MC_PER_BLOCK[20]=120    # For L=20, each block = 120 passes

NCONF=100                # 40 configurations (enough for statistics)
NMEAS=50              # 25 blocks between configs
NREP=6                  # 10 replicas per disorder
NDIS=10                  # 10 disorder samples per L

# Thread counts for each L (must be power of 2)
declare -A THREAD_COUNTS
THREAD_COUNTS[8]=1024    # L=8: 8×8×2048=131072 spins → 1024 threads (128 spins/thread)
THREAD_COUNTS[12]=256    # L=12: 12×12×2048=294912 → 256 threads (1152 spins/thread)
THREAD_COUNTS[16]=512    # L=16: 16×16×2048=524288 → 512 threads (1024 spins/thread)
THREAD_COUNTS[20]=512    # L=20: 20×20×2048=819200 → 512 threads (1600 spins/thread)
