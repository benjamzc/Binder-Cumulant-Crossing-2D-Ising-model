#!/bin/bash
#
#SBATCH --job-name=overlap
#SBATCH --output=/ada/ptmp/pks/munbe/Documents/QISG_PROJECT/logs/overlap_%A_%a.out
#SBATCH --error=/ada/ptmp/pks/munbe/Documents/QISG_PROJECT/logs/overlap_%A_%a.err
#SBATCH --time=10-00:00:00
#SBATCH --partition=p.ada
#SBATCH --gres=gpu:4
#SBATCH --array=0
#SBATCH --nodes=1
#SBATCH --cpus-per-task=72
#SBATCH --qos=longrun
#SBATCH --mem=128G

module purge
module load cuda/12.6
module load gcc/11

cd /ada/ptmp/pks/munbe/Documents/MCQSG/analysis/OVERLAP

for ibit in 0 1 2 3 4 5; do
  /ada/ptmp/pks/munbe/Documents/MCQSG/analysis/OVERLAP/coverlap_pt_L8_Lz2048_NK24_NTHXBLO1024 \
    /ada/ptmp/pks/munbe/Documents/MCQSG/output/minimal/ \
    0 ${ibit} 100 0.5 0 0 0
  echo "→ Finished L=8, IBIT=$ibit at $(date)"

  /ada/ptmp/pks/munbe/Documents/MCQSG/analysis/OVERLAP/coverlap_pt_L12_Lz2048_NK24_NTHXBLO1024 \
    /ada/ptmp/pks/munbe/Documents/MCQSG/output/minimal/ \
    0 ${ibit} 100 0.5 0 0 0
  echo "→ Finished L=12, IBIT=$ibit at $(date)"

  /ada/ptmp/pks/munbe/Documents/MCQSG/analysis/OVERLAP/coverlap_pt_L16_Lz2048_NK24_NTHXBLO1024 \
    /ada/ptmp/pks/munbe/Documents/MCQSG/output/minimal/ \
    0 ${ibit} 100 0.5 0 0 0
  echo "→ Finished L=16, IBIT=$ibit at $(date)"

  /ada/ptmp/pks/munbe/Documents/MCQSG/analysis/OVERLAP/coverlap_pt_L20_Lz2048_NK24_NTHXBLO1024 \
    /ada/ptmp/pks/munbe/Documents/MCQSG/output/minimal/ \
    0 ${ibit} 100 0.5 0 0 0
  echo "→ Finished L=20, IBIT=$ibit at $(date)"

done
