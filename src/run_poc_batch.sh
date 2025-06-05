#!/bin/bash
#
#SBATCH --job-name=mcqsg_full
#SBATCH --output=/ada/ptmp/pks/munbe/Documents/QISG_PROJECT/logs/mcqsg_full_%A_%a.out
#SBATCH --error=/ada/ptmp/pks/munbe/Documents/QISG_PROJECT/logs/mcqsg_full_%A_%a.err
#SBATCH --time=10-00:00:00
#SBATCH --partition=p.ada
#SBATCH --gres=gpu:4
#SBATCH --array=0-299
#SBATCH --nodes=1
#SBATCH --cpus-per-task=72
#SBATCH --qos=longrun
#SBATCH --mem=128G

module purge
module load cuda/12.6
module load gcc/11

# === 1. Source global simulation parameters ===
source /ada/ptmp/pks/munbe/Documents/QISG_PROJECT/scripts/poc_parameters.sh

# === 2. Decode SLURM_ARRAY_TASK_ID → (L, ibit, rep) ===
TASK_ID=${SLURM_ARRAY_TASK_ID}

REPLICA=$(( TASK_ID % NREP ))                 # 0–5
IBIT=$(( (TASK_ID / NREP) % NDIS ))           # 0–5
LIDX=$(( TASK_ID / (NREP * NDIS) ))           # 0–2

# === 3. Set L, Lz, threads ===
LVALS=(8 12 16 20)
L=${LVALS[$LIDX]}
LZ=${SYSTEM_PARAMS[$L]}
THREAD_COUNT=${THREAD_COUNTS[$L]}

# === 4. File paths ===
ISAMPLE=0
COFF=0
BIT=$(printf "%03d" $IBIT)
RR=$(printf "%02d" $REPLICA)

K_FILE=/ada/ptmp/pks/munbe/Documents/MCQSG/input/k.minimal
INPUT_FILE=/ada/ptmp/pks/munbe/Documents/MCQSG/input/input.${L}x${LZ}_minimal
LUT_FILE=/ada/ptmp/pks/munbe/Documents/MCQSG/input/LUT_for_PRNG_nbits10_NK${NK}.bin

# === 5. Output directory ===
OUTDIR=/ada/ptmp/pks/munbe/Documents/MCQSG/output/minimal/${L}x${LZ}/I000/BIT${BIT}/R${RR}
mkdir -p "$OUTDIR"

# === 6. Assign GPU ===
GPU_ID=$(( SLURM_ARRAY_TASK_ID % 4 ))
export CUDA_VISIBLE_DEVICES=$GPU_ID

echo "▶ Starting task $TASK_ID: L=$L, ibit=$IBIT, rep=$REPLICA, GPU=$GPU_ID"
echo "▶ Output: $OUTDIR"
echo "▶ Input: $INPUT_FILE"

# === 7. Run ===
srun /ada/ptmp/pks/munbe/Documents/MCQSG/MCQSG $ISAMPLE $IBIT $REPLICA $THREAD_COUNT $COFF $K_FILE $INPUT_FILE $LUT_FILE
